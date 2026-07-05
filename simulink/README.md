# simulink — independent re-implementation of the closed loop

```matlab
build_motor_speed_loop   % constructs and saves motor_speed_pi.slx
run_and_compare          % simulates + overlays the Python waveforms
```

`build_motor_speed_loop.m` script-builds the model (nothing hand-drawn,
same convention as Project 1): identified plant `K/(τs+1)` behind a
`Transport Delay θ`, a Discrete PID block (PI, Backward Euler, Ts = 10 ms,
output clamped 0..Vbus, **back-calculation anti-windup with
Kb = Ki/Kp = 1/Tt**) and an encoder measurement path (ZOH → quantizer at
13.1 rad/s → 15 Hz discrete IIR) — block-for-block the difference
equations in `control-design/pi_design.py` and the Arduino sketch.

Design numbers are read from `../control-design/results/controller.json`
at build time so Python stays the single source of truth; transcribed
fallbacks are inlined for running the folder standalone.

`run_and_compare.m` overlays the Simulink run on
`../control-design/results/closedloop_sim.csv` and prints overshoot /
disturbance-dip / final-value cross-checks. The two implementations share
only the design numbers, so agreement (envelope within a few percent;
the encoder limit-cycle fine structure may differ, since Python counts
encoder edges exactly while Simulink uses a Quantizer block) validates
both. Requires Simulink only — no Simscape, no toolboxes.
