function run_and_compare()
%RUN_AND_COMPARE Simulate motor_speed_pi.slx and overlay the Python sim.
%
%   Loads control-design/results/closedloop_sim.csv (written by
%   pi_design.py) and overlays speed and command voltage. The two
%   implementations share nothing but the design numbers, so agreement
%   validates both: the Python difference equations and the Simulink
%   block diagram. Expected agreement:
%     - rise time / overshoot within a few percent (the Python sim
%       quantizes the encoder inside the loop *count-exactly*, Simulink
%       uses a Quantizer approximation, so small limit-cycle detail may
%       differ -- the envelope must match)
%     - disturbance dip and recovery within a few percent

here = fileparts(mfilename('fullpath'));
csvf = fullfile(here, '..', 'control-design', 'results', ...
                'closedloop_sim.csv');
assert(exist(csvf, 'file') == 2, ...
    'Run control-design/pi_design.py first (writes closedloop_sim.csv)');
d = readmatrix(csvf);                  % t, ref, omega, vcmd
tp = d(:,1); refp = d(:,2); wp = d(:,3); up = d(:,4);

if ~bdIsLoaded('motor_speed_pi')
    if exist(fullfile(here, 'motor_speed_pi.slx'), 'file') ~= 2
        run(fullfile(here, 'build_motor_speed_loop.m'));
    end
    load_system(fullfile(here, 'motor_speed_pi'));
end
out = sim('motor_speed_pi', 'ReturnWorkspaceOutputs', 'on');
ws = out.get('omega_log');  vs = out.get('vcmd_log');

% Project 1 palette, fixed assignment
cP = [0.180 0.435 0.718];   % Python
cS = [0.851 0.471 0.176];   % Simulink

figure('Name', 'Python vs Simulink closed loop');
subplot(2,1,1); hold on; grid on;
plot(tp, wp, '-',  'Color', cP, 'LineWidth', 1.4);
plot(ws.Time, ws.Data, '--', 'Color', cS, 'LineWidth', 1.4);
plot(tp, refp, ':', 'Color', [0.43 0.43 0.43]);
ylabel('speed [rad/s]');
legend({'Python firmware-exact sim', 'Simulink', 'setpoint'}, ...
       'Location', 'southeast');
title('Closed-loop step + disturbance: Python vs Simulink');

subplot(2,1,2); hold on; grid on;
plot(tp, up, '-',  'Color', cP, 'LineWidth', 1.2);
plot(vs.Time, vs.Data, '--', 'Color', cS, 'LineWidth', 1.2);
ylabel('v\_cmd [V]'); xlabel('time [s]');

% quantitative cross-checks
w_sl = interp1(ws.Time, ws.Data, tp, 'previous', 'extrap');
wref = max(refp);
fprintf('\n--- cross-checks (Python | Simulink) ---\n');
seg = tp >= 0.25 & tp < 3.0;
fprintf('step overshoot : %5.1f %% | %5.1f %%\n', ...
    100*(max(wp(seg))-wref)/wref, 100*(max(w_sl(seg))-wref)/wref);
seg = tp >= 3.0;
fprintf('dist. dip      : %5.1f rad/s | %5.1f rad/s\n', ...
    wref - min(wp(seg)), wref - min(w_sl(seg)));
fprintf('final speed    : %6.1f | %6.1f rad/s (setpoint %g)\n', ...
    mean(wp(tp > 4.5)), mean(w_sl(tp > 4.5)), wref);
end
