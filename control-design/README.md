# control-design — PI from the identified plant

```
python3 pi_design.py     # reads ../system-id/results/identified_plant.json
```

Same design discipline as Project 1's buck compensator: pick a crossover
and phase margin, evaluate the (identified) plant there, solve for the PI
zero and gain in closed form, then verify on a dense grid and in a
nonlinear time-domain simulation.

**Targets: fc = 2.5 Hz, PM ≥ 60°.** That is ~4× the identified open-loop
pole (0.63 Hz) — a real "throttle response" improvement — while spending
only ~8° on the 9 ms identified dead time and ~9.5° on the 15 Hz
measurement filter. The filter is not optional hygiene: the encoder gives
13.1 rad/s per count at 100 Hz, and its phase lag is designed for
*inside* the loop, not discovered later on the bench.

Result (synthetic-data run): **Kp = 0.0446 V/(rad/s),
Ki = 0.349 V/(rad/s)/s** (zero at 1.25 Hz) → verified fc = 2.50 Hz,
PM = 60.0°, GM = 18.1 dB. Robustness sweep ±30 % on K and τ (hot
winding, added payload) keeps PM ≥ 55.9° with no retuning.

Verification is firmware-exact, not just linear: `simulate_firmware()`
runs the same difference equations the Arduino executes — 0..5 V
saturation against the buck rail, encoder quantization, 15 Hz IIR,
back-calculation anti-windup — against the identified plant:

| Metric | Value |
|---|---|
| 0 → 300 rad/s rise (10–90 %) | 250 ms |
| Overshoot, anti-windup ON | 5.6 % |
| Overshoot, anti-windup OFF | 30.6 % (command rails at 5 V for ~0.3 s — this is why the firmware has AW) |
| −0.75 V input disturbance ("the hill") | 11.8 rad/s dip (3.9 %), 2 %-band recovery 0.19 s |

Outputs in `results/`: `design_report.txt` (full derivation),
`bode_loop.png`, `closedloop_step.png`, `disturbance_rejection.png`,
`controller.json` (consumed by `simulink/build_motor_speed_loop.m`),
`closedloop_sim.csv` (consumed by `simulink/run_and_compare.m`), and the
firmware `#define` block to transcribe into
`firmware/motor_speed_pi/motor_speed_pi.ino`.
