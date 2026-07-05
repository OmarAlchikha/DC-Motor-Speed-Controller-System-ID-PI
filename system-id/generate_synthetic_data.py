"""Synthetic bench: simulate the step-response experiment the firmware runs.

Produces data/step_response_synthetic.csv with exactly the columns the
Arduino sysid mode streams (see firmware/README.md):

    t_s, v_cmd_V, omega_meas_rad_s

so the whole pipeline (fit_plant.py -> control-design/pi_design.py) can be
exercised end-to-end before any real data exists. Swap in your bench CSV
with `python3 fit_plant.py --csv <file>` later -- nothing else changes.

What the generator deliberately includes (and the fit must cope with):
  * full nonlinear motor ODEs at a 10 kHz inner step (RK-free semi-implicit
    Euler; the electrical pole at ~1.1 ms is resolved with 9 sub-steps)
  * Coulomb friction + stiction + brush drop  -> gain droops at low speed
  * one control-period command latency + ZOH  -> shows up as dead time
  * encoder quantization: speed is the *count difference* over each 10 ms
    window, exactly like the ISR-based measurement in the firmware
"""

import csv
import math
import os

import numpy as np

import motor_parameters as mp

DT = 1e-4                      # inner simulation step [s]
SUB = int(round(mp.TS / DT))   # inner steps per control sample


def simulate_staircase(seed=1):
    """Run the staircase experiment; return (t, v_cmd, omega_meas, omega_true)."""
    rng = np.random.default_rng(seed)
    n_per = int(round(mp.STEP_HOLD / mp.TS))
    v_cmd = np.repeat(mp.STEP_LEVELS, n_per)          # commanded voltage
    n = len(v_cmd)
    t = np.arange(n) * mp.TS

    # command reaches the PWM one control period late (ISR pipeline)
    v_applied_seq = np.concatenate([np.zeros(mp.LATENCY_SAMPLES),
                                    v_cmd[:n - mp.LATENCY_SAMPLES]])

    i = 0.0          # armature current [A]
    w = 0.0          # true shaft speed [rad/s]
    theta = 0.0      # true shaft angle [rad]
    counts_prev = 0

    omega_meas = np.empty(n)
    omega_true = np.empty(n)

    for k in range(n):
        v = v_applied_seq[k]                          # ZOH over the period
        for _ in range(SUB):
            # electrical: semi-implicit in i (stable for dt >> L/R too)
            v_eff = v - mp.V_BRUSH * np.sign(i)
            i = (i + DT / mp.L * (v_eff - mp.KE * w)) / (1.0 + DT * mp.R / mp.L)
            # mechanical with stiction: motor stays stuck until torque
            # exceeds breakaway
            tau_m = mp.KT * i
            if w == 0.0 and abs(tau_m) < mp.STICTION * mp.TAU_C:
                dw = 0.0
            else:
                tau_fric = mp.B * w + mp.TAU_C * np.sign(w if w != 0 else tau_m)
                dw = DT / mp.J * (tau_m - tau_fric)
            w = w + dw
            if w < 0.0:        # single-quadrant drive: freewheel, no reverse
                w = 0.0
            theta += w * DT

        # encoder measurement: counts accumulated this window
        counts = math.floor(theta / (2 * math.pi) * mp.CPR)
        # +/-1 count edge jitter (comparator/edge-timing noise)
        jitter = int(rng.integers(-1, 2))
        omega_meas[k] = ((counts - counts_prev) + jitter) \
            * (2 * math.pi / mp.CPR) / mp.TS
        counts_prev = counts
        omega_true[k] = w

    return t, v_cmd, np.maximum(omega_meas, 0.0), omega_true


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    outdir = os.path.join(here, "data")
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "step_response_synthetic.csv")

    t, v, wm, wt = simulate_staircase()

    with open(path, "w", newline="") as f:
        wr = csv.writer(f)
        # metadata header consumed by fit_plant.py (rig facts travel with
        # the data, so control-design never re-declares them)
        wr.writerow(["# vbus_V=%.2f" % mp.VBUS,
                     "ts_s=%.4f" % mp.TS,
                     "cpr=%d" % mp.CPR,
                     "source=synthetic"])
        wr.writerow(["t_s", "v_cmd_V", "omega_meas_rad_s"])
        for k in range(len(t)):
            wr.writerow(["%.3f" % t[k], "%.3f" % v[k], "%.4f" % wm[k]])

    q = 2 * math.pi / mp.CPR / mp.TS
    print("Wrote %s  (%d samples, %.1f s)" % (path, len(t), t[-1] + mp.TS))
    print("  staircase levels: %s V, %.1f s each" % (mp.STEP_LEVELS, mp.STEP_HOLD))
    print("  encoder quantization step: %.1f rad/s "
          "(the fit must average through this)" % q)
    print("  peak measured speed: %.0f rad/s" % wm.max())


if __name__ == "__main__":
    main()
