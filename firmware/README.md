# firmware — Arduino Mega, 100 Hz PI with anti-windup

`motor_speed_pi/motor_speed_pi.ino`. Gains are **outputs** of the design
pipeline (`control-design/results/design_report.txt`) — never hand-tune
them here; re-run the pipeline on new data instead.

## Wiring

```
buck 5 V (Project 1 output) ──┬── motor (+)
                              │
                        SB540 K (freewheel,
                              │  cathode to 5 V)
                              │
        motor (−) ────────────┴── FET drain
Mega pin 11 ── 150 Ω ── gate  (IRLB8721 / IRLZ44N, logic-level N-FET)
        FET source ── GND (star point, common with Mega GND and buck GND)

encoder A ── Mega pin 2 (INT4)     encoder VCC ── Mega 5 V (USB rail,
encoder B ── Mega pin 3 (INT5)     NOT the motor rail — keep PWM ripple
encoder GND ── Mega GND            out of the sensor supply)
```

| Item | Choice | Why |
|---|---|---|
| PWM | pin 11, Timer1, **31.4 kHz** | above audio; ≫ motor electrical pole (~1 ms) so the winding sees clean duty-average voltage; Timer0/millis untouched |
| Drive | low-side N-FET, single quadrant | same perfboard-survival logic as Project 1's P-FET choice: no bootstrap, no shoot-through possible |
| Freewheel | SB540 Schottky | same part as Project 1's BOM — one order |
| Loop rate | 100 Hz, drift-free micros() scheduler | 40× the 2.5 Hz crossover; matches the Ts the plant was identified at |
| Encoder | ×4 quadrature, two CHANGE ISRs, table decode | 48 counts/rev → 13.1 rad/s per count per 10 ms window, hence the 15 Hz IIR whose lag is in the loop design |
| Anti-windup | back-calculation, Tt = Kp/Ki | the 0→300 rad/s step rails the command at 5 V for ~0.3 s; without AW overshoot goes 6 % → 31 % |
| Stall guard | >0.9·Vbus commanded, <20 rad/s, 1.5 s → latch off | locked rotor draws 1.85 A continuously — legal for the buck, pointless and hot |

## Serial interface (115200 baud)

| Command | Action |
|---|---|
| `s300` | closed-loop RUN at 300 rad/s |
| `i` | run the system-ID staircase (open loop), stream CSV |
| `x` | stop / clear fault |
| `t` | toggle 10 Hz telemetry (`ref  speed  volts`) |
| `?` | help + active gains |

## Collecting real step-response data

The `i` command runs the same staircase as the synthetic generator and
streams a CSV in the exact format `system-id/fit_plant.py --csv` consumes
(capture it with `pio device monitor -b 115200 | tee bench.csv` or any
serial logger — keep the `#`/header lines, they're part of the format).
Re-running `fit_plant.py` and `pi_design.py` on that file regenerates the
identified plant and gains directly from bench data, and transcribing the
resulting `#define` block back into this sketch closes the loop. Because
the staircase is byte-identical to the synthetic one, the fitted synthetic
vs. bench parameters are directly comparable — the difference *is* the
model error.
