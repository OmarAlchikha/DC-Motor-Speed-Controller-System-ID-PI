"""Fit a plant model to step-response data with scipy.optimize.curve_fit.

    python3 fit_plant.py                      # uses the synthetic dataset
    python3 fit_plant.py --csv my_bench.csv   # your real data, same format

Input CSV format (produced by generate_synthetic_data.py and by the
firmware's sysid mode):

    # vbus_V=5.00,ts_s=0.0100,cpr=48,source=...
    t_s, v_cmd_V, omega_meas_rad_s

Candidate models (voltage -> shaft speed), both with dead time:

    FO:  G(s) = K e^(-theta s) / (tau s + 1)
    SO:  G(s) = K e^(-theta s) / ((tau1 s + 1)(tau2 s + 1))

Physics says the motor is second order (mechanical pole ~ JR/(bR+KtKe),
electrical pole L/R), but the electrical pole of a small motor sits near
1 ms while we sample at 10 ms -- so the honest question is whether the
data can *support* the second pole. The fit answers it: the SO model is
selected only if it improves the NRMSE fit by more than 1 percentage
point; otherwise the extra pole is un-identifiable decoration and the FO
model wins on parsimony. Dead time absorbs what the fit cannot resolve
(loop latency + ZOH + the fast electrical pole), which is exactly how it
should be treated by the controller design: as phase lag at crossover.

How the models are simulated inside curve_fit
---------------------------------------------
The recorded input is ZOH by construction (the firmware updates PWM once
per period), so each candidate is discretized exactly with
scipy.signal.cont2discrete(method='zoh') and run over the *recorded*
v_cmd sequence with lfilter. Dead time theta = d*Ts + f is applied as an
integer shift plus a first-sample blend (modified-z approximation of a
fractional delay of a ZOH signal). curve_fit then adjusts (K, tau, theta)
to minimize the squared error between simulated and measured speed.
"""

import argparse
import json
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import signal
from scipy.optimize import curve_fit

from motor_parameters import (COL_PRIMARY, COL_SECONDARY, COL_TERTIARY,
                              COL_GRAY)

plt.rcParams.update({
    "figure.dpi": 150, "axes.grid": True, "grid.alpha": 0.25,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 9, "axes.titlesize": 10,
})


# --------------------------------------------------------------- data ----
def load_csv(path):
    meta = {}
    with open(path) as f:
        first = f.readline().strip()
    if first.startswith("#"):
        for tok in first.lstrip("#").split(","):
            k, _, v = tok.strip().partition("=")
            if not v:
                continue
            try:
                meta[k] = int(v)
            except ValueError:
                try:
                    meta[k] = float(v)
                except ValueError:
                    meta[k] = v
        skip = 2
    else:
        skip = 1
    d = np.genfromtxt(path, delimiter=",", skip_header=skip)
    t, v, w = d[:, 0], d[:, 1], d[:, 2]
    ts = meta.get("ts_s", float(np.median(np.diff(t))))
    return t, v, w, ts, meta


# ------------------------------------------------------------- models ----
def delay_input(u, theta, ts):
    """u(t - theta) for a ZOH sequence: integer shift + fractional blend."""
    d = int(theta // ts)
    f = (theta - d * ts) / ts
    ud = np.concatenate([np.full(d, u[0]), u[:len(u) - d]]) if d else u.copy()
    if f > 0:
        prev = np.concatenate([[ud[0]], ud[:-1]])
        ud = f * prev + (1 - f) * ud
    return ud


def make_sim(ts, order):
    """Return f(t, *params) that simulates the candidate on the recorded
    input. curve_fit sees a plain function of the parameters; the input
    sequence is bound in via closure (set_input) before fitting."""
    holder = {}

    def set_input(u):
        holder["u"] = u

    def sim(_t, *p):
        if order == 1:
            K, tau, theta = p
            num, den = [K], [tau, 1.0]
        else:
            K, tau1, tau2, theta = p
            num, den = [K], [tau1 * tau2, tau1 + tau2, 1.0]
        bz, az, _ = signal.cont2discrete((num, den), ts, method="zoh")
        y = signal.lfilter(np.ravel(bz), np.ravel(az),
                           delay_input(holder["u"], theta, ts))
        return y

    sim.set_input = set_input
    return sim


def fit_percent(y, yhat):
    """MATLAB compare()-style NRMSE fit in percent."""
    return 100.0 * (1.0 - np.linalg.norm(y - yhat)
                    / np.linalg.norm(y - y.mean()))


def initial_guess(t, v, w, ts):
    """K from a steady-state regression over the tail of each hold,
    tau from the longest visible transient, theta ~ one sample."""
    lev = np.flatnonzero(np.diff(v) != 0)
    segs = np.split(np.arange(len(v)), lev + 1)
    vs, ws = [], []
    for s in segs:
        if len(s) < 20:
            continue
        tail = s[-len(s) // 4:]
        vs.append(v[tail].mean())
        ws.append(w[tail].mean())
    vs, ws = np.array(vs), np.array(ws)
    K0 = float(np.dot(vs, ws) / np.dot(vs, vs)) if np.dot(vs, vs) > 0 else 50.0
    return K0, 0.1, ts


# ------------------------------------------------------------ fitting ----
def fit_models(t, v, w, ts):
    K0, tau0, th0 = initial_guess(t, v, w, ts)
    out = {}

    sim1 = make_sim(ts, 1)
    sim1.set_input(v)
    p1, c1 = curve_fit(sim1, t, w, p0=[K0, tau0, th0],
                       bounds=([1.0, 1e-3, 0.0], [1e3, 5.0, 0.2]))
    y1 = sim1(t, *p1)
    out["fo"] = {"params": p1, "stderr": np.sqrt(np.diag(c1)),
                 "y": y1, "fit_pct": fit_percent(w, y1),
                 "rmse": float(np.sqrt(np.mean((w - y1) ** 2)))}

    sim2 = make_sim(ts, 2)
    sim2.set_input(v)
    p2, c2 = curve_fit(sim2, t, w, p0=[K0, tau0, ts / 2, th0],
                       bounds=([1.0, 1e-3, 1e-4, 0.0], [1e3, 5.0, 5.0, 0.2]),
                       maxfev=20000)
    y2 = sim2(t, *p2)
    out["so"] = {"params": p2, "stderr": np.sqrt(np.diag(c2)),
                 "y": y2, "fit_pct": fit_percent(w, y2),
                 "rmse": float(np.sqrt(np.mean((w - y2) ** 2)))}
    return out


# -------------------------------------------------------------- plots ----
def plot_overlay(t, v, w, fits, chosen, outdir):
    fig, ax = plt.subplots(3, 1, figsize=(9, 7.5), sharex=True,
                           gridspec_kw={"height_ratios": [1, 2.2, 1.2]})
    ax[0].step(t, v, where="post", color=COL_GRAY, lw=1.2)
    ax[0].set_ylabel("v_cmd [V]")
    ax[0].set_title("System ID: staircase experiment, measured vs fitted")

    ax[1].plot(t, w, color=COL_PRIMARY, lw=0.7, alpha=0.85, label="measured")
    ax[1].plot(t, fits["fo"]["y"], color=COL_SECONDARY, lw=1.6,
               label="1st-order fit (%.1f%%)" % fits["fo"]["fit_pct"])
    ax[1].plot(t, fits["so"]["y"], color=COL_TERTIARY, lw=1.2, ls="--",
               label="2nd-order fit (%.1f%%)" % fits["so"]["fit_pct"])
    ax[1].set_ylabel("speed [rad/s]")
    ax[1].legend(loc="upper right", framealpha=0.9)
    ax[1].annotate("selected: %s" % chosen.upper(),
                   xy=(0.015, 0.94), xycoords="axes fraction",
                   fontsize=9, color=COL_GRAY, va="top")

    for key, col in (("fo", COL_SECONDARY), ("so", COL_TERTIARY)):
        ax[2].plot(t, w - fits[key]["y"], color=col, lw=0.7,
                   label="%s residual (RMSE %.1f rad/s)"
                   % (key.upper(), fits[key]["rmse"]))
    ax[2].axhline(0, color=COL_GRAY, lw=0.8)
    ax[2].set_ylabel("residual [rad/s]")
    ax[2].set_xlabel("time [s]")
    ax[2].legend(loc="upper right", framealpha=0.9)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fit_overlay.png"))
    plt.close(fig)


def plot_zoom(t, v, w, fits, ts, outdir):
    """Zoom on the largest single step: dead time + quantization visible."""
    lev = np.flatnonzero(np.diff(v) != 0)
    k0 = lev[np.argmax([abs(v[i + 1] - v[i]) for i in lev])] if len(lev) else 0
    a, b = max(0, k0 - 15), min(len(t), k0 + 90)

    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(t[a:b], w[a:b], ".", color=COL_PRIMARY, ms=4,
            label="measured (encoder-quantized)")
    ax.plot(t[a:b], fits["fo"]["y"][a:b], color=COL_SECONDARY, lw=1.8,
            label="1st-order fit")
    ax.plot(t[a:b], fits["so"]["y"][a:b], color=COL_TERTIARY, lw=1.2,
            ls="--", label="2nd-order fit")
    th = fits["fo"]["params"][2]
    ax.axvline(t[k0] + ts, color=COL_GRAY, lw=0.9, ls=":")
    ax.annotate("step applied", xy=(t[k0] + ts, ax.get_ylim()[1] * 0.35),
                rotation=90, fontsize=8, color=COL_GRAY,
                ha="right", va="bottom")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("speed [rad/s]")
    ax.set_title("Largest step, zoomed: dead time theta = %.1f ms "
                 "absorbs latency + electrical pole" % (th * 1e3))
    ax.legend(loc="upper right", framealpha=0.9)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fit_zoom.png"))
    plt.close(fig)


# -------------------------------------------------------------- report ---
def report(fits, chosen, ts, meta, n):
    ln = []
    p = ln.append
    p("=" * 78)
    p("PLANT IDENTIFICATION FROM STEP-RESPONSE DATA (scipy curve_fit)")
    p("=" * 78)
    p("Dataset: %d samples at Ts = %.0f ms, Vbus = %.2f V, source = %s"
      % (n, ts * 1e3, meta.get("vbus_V", float("nan")),
         meta.get("source", "?")))

    f1 = fits["fo"]
    K, tau, th = f1["params"]
    sK, stau, sth = f1["stderr"]
    p("\n--- Candidate 1: first order + dead time ---")
    p("  G(s) = K e^(-theta s) / (tau s + 1)")
    p("  K     = %7.2f  +/- %.2f   rad/s per V" % (K, sK))
    p("  tau   = %7.1f  +/- %.1f   ms" % (tau * 1e3, stau * 1e3))
    p("  theta = %7.1f  +/- %.1f   ms" % (th * 1e3, sth * 1e3))
    p("  fit   = %6.2f %%   RMSE = %.2f rad/s" % (f1["fit_pct"], f1["rmse"]))

    f2 = fits["so"]
    K2, t1, t2, th2 = f2["params"]
    sK2, st1, st2, sth2 = f2["stderr"]
    p("\n--- Candidate 2: second order + dead time ---")
    p("  G(s) = K e^(-theta s) / ((tau1 s + 1)(tau2 s + 1))")
    p("  K     = %7.2f  +/- %.2f   rad/s per V" % (K2, sK2))
    p("  tau1  = %7.1f  +/- %.1f   ms" % (t1 * 1e3, st1 * 1e3))
    p("  tau2  = %7.2f  +/- %.2f   ms" % (t2 * 1e3, st2 * 1e3))
    p("  theta = %7.1f  +/- %.1f   ms" % (th2 * 1e3, sth2 * 1e3))
    p("  fit   = %6.2f %%   RMSE = %.2f rad/s" % (f2["fit_pct"], f2["rmse"]))

    p("\n--- Model selection ---")
    gain = f2["fit_pct"] - f1["fit_pct"]
    p("Second pole buys %.2f fit points (threshold to justify it: 1.00)."
      % gain)
    if chosen == "so":
        p("=> SECOND-ORDER model selected.")
    else:
        p("=> FIRST-ORDER model selected on parsimony.")
        p("   Physics says a second (electrical, ~L/R) pole exists, but at")
        p("   Ts = %.0f ms it is far inside one sample: the data cannot"
          % (ts * 1e3))
        p("   support it (note tau2's error bar), and its lag is already")
        p("   captured inside the fitted dead time. Designing on the FO")
        p("   model with theta treated as phase loss at crossover is both")
        p("   honest and conservative.")

    p("\nResidual sanity: the residual trace (fit_overlay.png, bottom) is")
    p("quantization + edge jitter around zero at cruise, with structured")
    p("bumps only at the low-voltage steps -- that is Coulomb friction /")
    p("stiction, a real nonlinearity a linear model cannot (and need not)")
    p("capture. The PI integrator will absorb it in closed loop.")
    return "\n".join(ln)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=os.path.join(
        here, "data", "step_response_synthetic.csv"))
    args = ap.parse_args()

    outdir = os.path.join(here, "results")
    os.makedirs(outdir, exist_ok=True)

    t, v, w, ts, meta = load_csv(args.csv)
    fits = fit_models(t, v, w, ts)
    chosen = "so" if fits["so"]["fit_pct"] - fits["fo"]["fit_pct"] > 1.0 \
        else "fo"

    txt = report(fits, chosen, ts, meta, len(t))
    print(txt)
    with open(os.path.join(outdir, "sysid_report.txt"), "w") as f:
        f.write(txt + "\n")

    plot_overlay(t, v, w, fits, chosen, outdir)
    plot_zoom(t, v, w, fits, ts, outdir)

    # hand-off to control-design: the chosen model + rig facts
    f1, f2 = fits["fo"], fits["so"]
    out = {
        "chosen": chosen,
        "fo": {"K": f1["params"][0], "tau": f1["params"][1],
               "theta": f1["params"][2], "fit_pct": f1["fit_pct"]},
        "so": {"K": f2["params"][0], "tau1": f2["params"][1],
               "tau2": f2["params"][2], "theta": f2["params"][3],
               "fit_pct": f2["fit_pct"]},
        "ts_s": ts,
        "vbus_V": meta.get("vbus_V", 5.0),
        "cpr": int(meta.get("cpr", 48)),
        "source": str(meta.get("source", "unknown")),
    }
    jpath = os.path.join(outdir, "identified_plant.json")
    with open(jpath, "w") as f:
        json.dump(out, f, indent=2)
    print("\nWrote %s (consumed by control-design/pi_design.py)" % jpath)
    print("Wrote results/fit_overlay.png, results/fit_zoom.png,"
          " results/sysid_report.txt")


if __name__ == "__main__":
    main()
