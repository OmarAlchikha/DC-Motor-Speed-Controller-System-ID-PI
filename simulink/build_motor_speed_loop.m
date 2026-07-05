function build_motor_speed_loop()
%BUILD_MOTOR_SPEED_LOOP Script-build the closed-loop motor speed model.
%
%   Creates and saves 'motor_speed_pi.slx': the PI speed loop on the
%   IDENTIFIED plant, block-for-block equivalent to the discrete
%   simulation in control-design/pi_design.py:
%
%     wref --> [Sum] --> [Discrete PI, Ts=10ms, limits 0..Vbus,     ]
%               ^-       [back-calculation anti-windup, Kb = Ki/Kp ]--+
%               |                                                     |
%               |     +--------------- v_dist (hill step) ---(+)<-----+
%               |     v
%               |   [Transport Delay theta] --> [K/(tau.s+1)] --> omega
%               |                                                  |
%               +-- [15 Hz discrete LPF] <-- [Quantizer] <-- [ZOH] -+
%                   (encoder: speed quantized to 2*pi/CPR/Ts steps)
%
%   Parameters are read from ../control-design/results/controller.json
%   (written by pi_design.py) so Python remains the single source of
%   truth; if the JSON is missing, the transcribed defaults below are
%   used. Run this once, then run_and_compare.m to overlay the Python
%   waveforms.
%
%   Requires: Simulink only (no Simscape). Written against R2021b+.

% ---- design values (fallback = transcribed from design_report.txt) ----
p.K     = 83.68;    p.tau   = 0.254;   p.theta = 0.009;   % identified plant
p.Kp    = 0.0446;   p.Ki    = 0.3493;                     % PI design
p.Ts    = 0.010;    p.fmeas = 15;      p.Vbus  = 5.0;
p.wref  = 300;      p.tstep = 0.25;                       % scenario
p.vdist = -0.75;    p.tdist = 3.0;     p.Tstop = 5.0;
p.cpr   = 48;

jsonfile = fullfile(fileparts(mfilename('fullpath')), '..', ...
    'control-design', 'results', 'controller.json');
if exist(jsonfile, 'file')
    j = jsondecode(fileread(jsonfile));
    p.K = j.plant.K;  p.tau = j.plant.tau;  p.theta = j.plant.theta;
    p.Kp = j.Kp;      p.Ki = j.Ki;          p.Ts = j.ts_s;
    p.fmeas = j.f_meas_Hz;                  p.Vbus = j.vbus_V;
    p.wref = j.w_ref; p.tstep = j.t_step;
    p.vdist = j.v_dist;                     p.tdist = j.t_dist;
    fprintf('Loaded design values from controller.json\n');
end
alpha = 1 - exp(-2*pi*p.fmeas*p.Ts);        % discrete LPF coefficient
qstep = 2*pi/p.cpr/p.Ts;                    % encoder speed quantum [rad/s]

mdl = 'motor_speed_pi';
close_system(mdl, 0);
if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end
new_system(mdl);
open_system(mdl);

    function h = blk(libpath, name, pos, varargin)
        h = add_block(libpath, [mdl '/' name], 'Position', pos, varargin{:});
    end

% ---------------- reference, error, controller ----------------
blk('simulink/Sources/Step', 'Wref', [40 100 80 140], ...
    'Time', num2str(p.tstep), 'Before', '0', 'After', num2str(p.wref));
blk('simulink/Math Operations/Sum', 'Err', [130 105 160 135], ...
    'Inputs', '+-');
% Discrete PI, Backward Euler (same difference equation as the firmware),
% output clamped to the buck rail, back-calculation anti-windup Kb = 1/Tt
% with Tt = Kp/Ki -- identical to pi_design.py / the Arduino sketch.
blk('simulink/Discrete/Discrete PID Controller', 'PI', [200 90 270 150], ...
    'Controller', 'PI', ...
    'P', num2str(p.Kp), 'I', num2str(p.Ki), ...
    'IntegratorMethod', 'Backward Euler', ...
    'SampleTime', num2str(p.Ts), ...
    'LimitOutput', 'on', ...
    'UpperSaturationLimit', num2str(p.Vbus), ...
    'LowerSaturationLimit', '0', ...
    'AntiWindupMode', 'back-calculation', ...
    'Kb', num2str(p.Ki / p.Kp));

% ---------------- disturbance + plant ----------------
blk('simulink/Sources/Step', 'Vdist', [200 190 240 230], ...
    'Time', num2str(p.tdist), 'Before', '0', 'After', num2str(p.vdist));
blk('simulink/Math Operations/Sum', 'PlantIn', [310 105 340 135], ...
    'Inputs', '++');
blk('simulink/Continuous/Transport Delay', 'Theta', [380 100 420 140], ...
    'DelayTime', num2str(p.theta));
blk('simulink/Continuous/Transfer Fcn', 'Plant', [460 95 540 145], ...
    'Numerator', sprintf('[%g]', p.K), ...
    'Denominator', sprintf('[%g 1]', p.tau));

% ---------------- encoder measurement path ----------------
% Count-difference speed over one Ts window ~ ZOH sample + quantize to
% one-count resolution, then the firmware's single-pole IIR at 15 Hz.
blk('simulink/Discrete/Zero-Order Hold', 'EncZOH', [460 220 500 260], ...
    'SampleTime', num2str(p.Ts));
blk('simulink/Discontinuities/Quantizer', 'EncQuant', [390 220 430 260], ...
    'QuantizationInterval', num2str(qstep));
blk('simulink/Discrete/Discrete Transfer Fcn', 'MeasLPF', ...
    [290 215 350 265], ...
    'Numerator', sprintf('[%g]', alpha), ...
    'Denominator', sprintf('[1 %g]', -(1 - alpha)), ...
    'SampleTime', num2str(p.Ts));

% ---------------- logging ----------------
blk('simulink/Sinks/To Workspace', 'SpeedLog', [620 95 690 135], ...
    'VariableName', 'omega_log', 'SaveFormat', 'Timeseries');
blk('simulink/Sinks/To Workspace', 'VcmdLog', [310 20 380 60], ...
    'VariableName', 'vcmd_log', 'SaveFormat', 'Timeseries');
blk('simulink/Sinks/Scope', 'Scope', [620 170 680 220], ...
    'NumInputPorts', '2');

wires = {
 'Wref/1',    'Err/1'
 'Err/1',     'PI/1'
 'PI/1',      'PlantIn/1'
 'PI/1',      'VcmdLog/1'
 'Vdist/1',   'PlantIn/2'
 'PlantIn/1', 'Theta/1'
 'Theta/1',   'Plant/1'
 'Plant/1',   'SpeedLog/1'
 'Plant/1',   'Scope/1'
 'Plant/1',   'EncZOH/1'
 'EncZOH/1',  'EncQuant/1'
 'EncQuant/1','MeasLPF/1'
 'MeasLPF/1', 'Err/2'
 'PI/1',      'Scope/2'
};
for k = 1:size(wires, 1)
    add_line(mdl, wires{k,1}, wires{k,2}, 'autorouting', 'on');
end

set_param(mdl, 'StopTime', num2str(p.Tstop), ...
    'Solver', 'ode23tb', 'RelTol', '1e-5', ...
    'MaxStep', num2str(p.Ts / 5), ...
    'SaveTime', 'on', 'SaveOutput', 'on');

save_system(mdl);
fprintf('Built and saved %s.slx\n', mdl);
fprintf('Simulate it directly or run run_and_compare.m\n');
end
