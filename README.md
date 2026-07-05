# DC motor speed controller — system ID + PI, the miniature EV powertrain

Project 2 in the series. Project 1
([buck-converter](../../../buck-converter)) built the "charging / DC-link"
stage: a 12 V → 5 V, 2 A closed-loop buck. This project bolts the
"traction" stage onto it: a brushed DC motor speed loop powered from that
buck output, with the plant model **measured, not assumed** — and every
number below produced by the scripts in this repo, nothing hand-waved.

## Power budget: Project 1's actual output

| Buck output (measured design values from `buck-converter/python-model/results/design_report.txt`) | Value |
|---|---|
| Rail voltage | **5.0 V** regulated (voltage-mode TL494 PI loop, fc = 7.77 kHz, PM = 49.2°) |
| Continuous current | **2.0 A** |
| Output ripple | ≈ 33 mV pk-pk at 100 kHz |
| Load-step behavior | 1.5 A step → ~90 mV dip, ~0.3 ms recovery |

The motor stage is designed *into* that budget: a Pololu-25D-class
brushed gearmotor (R ≈ 2.7 Ω) whose **stall current at 5 V is 1.85 A <
2.0 A** — the worst case the motor can do is legal for the supply, no
current loop needed. The buck's 33 mV ripple is 0.7 % of the rail and
two decades above the speed loop's bandwidth: invisible. Its ~0.3 ms
load-step recovery means motor current transients (τ_elec ≈ 1 ms) don't
fight the voltage loop. The stages are separable because their
bandwidths are: 7.8 kHz ≫ 2.5 Hz.

## Why system ID — the actual point of this project

Project 1's plant came from physics: L and C are purchased numbers, the
averaged model follows from Kirchhoff, and the bench matches to a few
percent. A motor + drivetrain is not like that:

- **The datasheet doesn't know your load.** τ_mech = JR/(Rb + KtKe)
  depends on the *total* inertia and friction of whatever you bolted on.
  A real EV's mass, tire, and gearbox losses exist nowhere on the motor
  datasheet — same for this rig's flywheel.
- **Half the parameters aren't on the datasheet at all.** Coulomb
  friction, stiction, brush drop, gearbox drag — and they're not even
  linear. The fit on synthetic data recovers K = 83.7 rad/s/V where the
  linearized truth is 88.1: that 5 % gap *is* the friction, correctly
  averaged into an effective gain.
- **Latency is a measured quantity.** The identified 9 ms dead time is
  firmware pipeline + ZOH + the un-resolvable electrical pole, in one
  number, measured through the exact signal path the closed loop uses —
  8° of phase at crossover that an assumed model would have missed.
- **Parameters drift.** Winding heats (K drops), payload changes (τ
  grows). The design carries a ±30 % robustness sweep so the identified
  model's error bars are budgeted, not ignored.

So the pipeline is: **measure → fit → design → verify → flash**, and it
runs end-to-end on a synthetic bench first, so every script is proven
before real data exists.

## The numbers (synthetic-data run; replace after the bench run)

| Stage | Result |
|---|---|
| Identified plant (FO + dead time, 90.7 % fit) | **K = 83.7 ± 0.1 rad/s/V, τ = 254 ± 3 ms, θ = 9 ± 2 ms** |
| Model selection | 2nd pole rejected: buys 0.00 fit-points at Ts = 10 ms (τ_elec ≈ 1 ms is sub-sample; its lag lives in θ) |
| Design targets | fc = 2.5 Hz (4× open-loop pole), PM ≥ 60° |
| PI gains | **Kp = 0.0446 V/(rad/s), Ki = 0.349 V/(rad/s)/s** (zero at 1.25 Hz) |
| Verified margins | fc = 2.50 Hz, PM = 60.0°, GM = 18.1 dB (dead time + 15 Hz meas. filter included) |
| Robustness (K, τ ± 30 %) | PM ≥ 55.9° everywhere — survives hot windings and payload |
| Step 0 → 300 rad/s | 250 ms rise, 5.6 % overshoot (31 % without anti-windup — command rails at 5 V for 0.3 s) |
| "Hill" disturbance (−0.75 V ≈ +2.8 mN·m) | 3.9 % dip, recovered in 0.19 s |

### Measured plant parameters (fill in after the bench run)

| Parameter | Synthetic run | **Real bench** |
|---|---|---|
| K [rad/s per V] | 83.7 ± 0.1 | *TBD* |
| τ [ms] | 254 ± 3 | *TBD* |
| θ [ms] | 9 ± 2 | *TBD* |
| fit [%] | 90.7 | *TBD* |
| → Kp, Ki | 0.0446, 0.349 | *TBD (re-run pi_design.py)* |

Procedure: `firmware/README.md` §"Collecting real step-response data" —
the firmware's `i` command streams a CSV in exactly the format
`fit_plant.py --csv` eats; the design pipeline then regenerates gains and
this table.

## The EV analogy, made explicit

| This rig | Full-size EV powertrain |
|---|---|
| Buck converter (Project 1) | DC-DC link / on-board charger stage |
| 5 V rail, 2 A budget | DC bus with a current limit the drive must respect by design |
| Brushed motor + flywheel | Traction motor + vehicle mass |
| Encoder, 13.1 rad/s per count at 100 Hz | Resolver/encoder with finite resolution forcing a filtered speed estimate |
| Step-response system ID | Drive commissioning run ("auto-tune") — real drives *measure* the plant on installation for exactly the reasons above |
| PI + anti-windup, command clamped to the rail | Torque/speed loop that saturates against bus voltage on every hard launch |
| −0.75 V disturbance step | The hill |
| Stall guard in firmware | Locked-rotor / traction fault protection |

## Repository map

| Folder | Contents |
|---|---|
| `system-id/` | synthetic bench generator, `curve_fit` plant fitting, fit/residual plots, `identified_plant.json` hand-off |
| `control-design/` | PI design from the identified plant (closed-form fz/Kp/Ki + dense-grid verification + firmware-exact nonlinear sim), all plots and `design_report.txt` |
| `simulink/` | `build_motor_speed_loop.m` (script-built .slx, reads `controller.json`), `run_and_compare.m` (overlays Python waveforms) |
| `firmware/` | Arduino Mega sketch: 100 Hz PI, ×4 quadrature decode, back-calculation anti-windup, stall guard, sysid streaming mode; wiring in its README |

## Run it

```bash
cd system-id
python3 generate_synthetic_data.py   # or skip, once you have bench.csv
python3 fit_plant.py                 # -> results/identified_plant.json
cd ../control-design
python3 pi_design.py                 # -> gains, plots, firmware constants
# optional cross-check in MATLAB:  cd ../simulink; build_motor_speed_loop; run_and_compare
# then flash firmware/motor_speed_pi/motor_speed_pi.ino and type: s300
```

The three implementations (Python firmware-exact sim, Simulink block
model, Arduino sketch) share only the design numbers — the same
cross-validation discipline as Project 1. If they agree and the bench
doesn't, the discrepancy is physical, and the place to fix it is the
model (re-run the ID), never the gains by hand.
