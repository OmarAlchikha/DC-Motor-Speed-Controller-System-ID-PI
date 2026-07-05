# system-id — fit the plant you actually have

```
python3 motor_parameters.py          # print the ground-truth "answer key"
python3 generate_synthetic_data.py   # simulate the bench -> data/*.csv
python3 fit_plant.py                 # fit FO & SO models -> results/
python3 fit_plant.py --csv bench.csv # same pipeline on your real data
```

`generate_synthetic_data.py` plays the role of the real bench: it
integrates the full nonlinear motor (Coulomb friction, stiction, brush
drop), applies the firmware's one-sample command latency, and measures
speed the way the Arduino ISR does — encoder count differences over each
10 ms window, quantized to 13.1 rad/s per count. `fit_plant.py` never
sees the true parameters; it sees only the CSV, exactly as it will see
your bench log.

`fit_plant.py` fits both candidates with `scipy.optimize.curve_fit`
(each candidate is ZOH-discretized and simulated over the recorded input
inside the objective):

- first order + dead time `K·e^(−θs)/(τs+1)`
- second order + dead time `K·e^(−θs)/((τ₁s+1)(τ₂s+1))`

and selects on parsimony: the second pole must buy ≥ 1 fit-point or it is
rejected. On the synthetic set it buys 0.00 — the electrical pole
(~1.1 ms) is invisible at Ts = 10 ms, and its lag lands in θ where the
controller design correctly treats it as phase loss at crossover.

Outputs in `results/`:

| File | Contents |
|---|---|
| `sysid_report.txt` | parameters ± std errors, fit %, model-selection verdict |
| `fit_overlay.png` | staircase input, measured vs both fits, residuals |
| `fit_zoom.png` | largest step zoomed: dead time + encoder quantization visible |
| `identified_plant.json` | the hand-off consumed by `control-design/pi_design.py` |

Synthetic-run scorecard (fit vs hidden truth): K 83.7 vs 88.1 rad/s/V
(the ~5 % droop is Coulomb friction — a linear model averages it into the
gain), τ 254 vs 237 ms, θ 9 ms ≈ 1 sample latency + ZOH + electrical
pole. The residual trace shows structure only at the low-voltage levels —
that *is* the friction nonlinearity, and the PI integrator absorbs it in
closed loop.
