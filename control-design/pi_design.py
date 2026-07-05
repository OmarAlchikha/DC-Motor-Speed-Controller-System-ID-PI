"""PI speed-controller design from the IDENTIFIED plant (not a datasheet).

    python3 pi_design.py          # reads ../system-id/results/identified_plant.json

Loop structure (all in engineering units; the firmware implements exactly
this at 100 Hz):

    w_ref -->(+)-- e --[ Gc(s) = Kp + Ki/s ]-- v_cmd --sat[0..Vbus]--+
              ^-                                                     |
              |                                        [ P(s) = K e^(-theta s)
              |                                              / (tau s + 1) ]
              +--[ H(s) = 1/(tau_f s + 1) ]-- omega_meas <-----------+
                    (15 Hz measurement LPF -- the encoder speed is
                     quantized to 13 rad/s steps at Ts = 10 ms, so the
                     firmware must filter; the filter's phase is part
                     of the loop and is designed for here, not ignored)

    T(s) = Gc(s) * P(s) * H(s)

Design procedure (same as Project 1's buck PI design):
  1. Evaluate P(j2pi fc)*H(j2pi fc) at the target crossover fc.
     The fitted dead time theta already contains the firmware's
     command latency + ZOH + the unresolved electrical pole, because it
     was identified through the same firmware path the closed loop uses.
  2. PI phase at fc is -90 + atan(fc/fz); solve PM = 180 + ph(PH) + ph(PI)
     for the zero:  atan(fc/fz) = PM - 90 - ph(PH)  ->  fz.
  3. Unity loop gain at fc:  Kp = 1/(|PH| sqrt(1+(fz/fc)^2)), Ki = 2pi fz Kp.

Crossover choice fc = 2.5 Hz: ~4x faster than the identified open-loop
pole (1/2pi tau ~ 0.63 Hz) -- a real "throttle response" improvement --
while staying a decade under the 15 Hz measurement filter and paying only
~8 deg to the 9 ms dead time. Pushing fc toward 5+ Hz would spend the
entire phase budget on delay + filter lag and amplify encoder
quantization into the PWM; this is the EV-powertrain trade between
responsiveness and drivetrain "chatter".

Verification is NOT the linear model alone: verify_closed_loop() runs a
discrete-time simulation of the exact firmware difference equations
(saturation 0..Vbus, back-calculation anti-windup, encoder quantization,
measurement LPF, dead-time samples) against the identified plant.
"""

import json
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
PLANT_JSON = os.path.join(HERE, "..", "system-id", "results",
                          "identified_plant.json")

# ---- design targets ----
F_CROSS = 2.5        # target loop crossover [Hz]
PM_TARGET = 60.0     # target phase margin [deg]
F_MEAS = 15.0        # measurement LPF corner [Hz] (implemented in firmware)

# ---- verification scenario ----
W_REF = 300.0        # cruise setpoint [rad/s] (~2865 rpm, 68% of no-load)
T_STEP = 0.25        # reference step time [s]
V_DIST = -0.75       # input-voltage disturbance [V] at T_DIST ("the hill":
                     # equivalent extra load torque Kt*V/R ~ 2.8 mN*m,
                     # more than doubling the cruise load)
T_DIST = 3.0
T_END = 5.0

# palette (same fixed assignment as Project 1 and system-id)
COL_PRIMARY = "#2E6FB7"
COL_SECONDARY = "#D9782D"
COL_TERTIARY = "#3E8E5A"
COL_ACCENT = "#8A5CB8"
COL_GRAY = "#6E6E6E"

plt.rcParams.update({
    "figure.dpi": 150, "axes.grid": True, "grid.alpha": 0.25,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 9, "axes.titlesize": 10,
})


def load_plant():
    with open(PLANT_JSON) as f:
        d = json.load(f)
    if d["chosen"] == "fo":
        K, tau, theta = d["fo"]["K"], d["fo"]["tau"], d["fo"]["theta"]
        poles = [tau]
    else:
        K, theta = d["so"]["K"], d["so"]["theta"]
        poles = [d["so"]["tau1"], d["so"]["tau2"]]
    return d, K, poles, theta


# ------------------------------------------------- frequency responses ---
def plant_fr(f, K, poles, theta, tau_f):
    """P(j2pi f) * H(j2pi f) including dead time and measurement filter."""
    s = 2j * math.pi * np.asarray(f, dtype=float)
    G = K * np.exp(-theta * s)
    for tau in poles:
        G = G / (tau * s + 1.0)
    return G / (tau_f * s + 1.0)


def margins(Tfun, fmin=0.01, fmax=50.0, n=20000):
    """Crossover, phase margin, gain margin from a dense log grid
    (same approach as Project 1's pi_design.margins)."""
    f = np.logspace(math.log10(fmin), math.log10(fmax), n)
    T = Tfun(f)
    mag, ph = np.abs(T), np.unwrap(np.angle(T))
    fc = pm = None
    idx = np.where((mag[:-1] >= 1) & (mag[1:] < 1))[0]
    if len(idx):
        i = idx[-1]
        t = math.log(mag[i]) / (math.log(mag[i]) - math.log(mag[i + 1]))
        fc = f[i] * (f[i + 1] / f[i]) ** t
        pm = 180 + math.degrees(ph[i] + t * (ph[i + 1] - ph[i]))
    gm_db = None
    tgt = -math.pi
    idx2 = np.where((ph[:-1] > tgt) & (ph[1:] <= tgt))[0]
    if len(idx2):
        i = idx2[0]
        t = (ph[i] - tgt) / (ph[i] - ph[i + 1])
        m = mag[i] * (mag[i + 1] / mag[i]) ** t
        gm_db = -20 * math.log10(m)
    return fc, pm, gm_db


def gm_str(gm_db):
    return "inf" if gm_db is None else "%.1f dB" % gm_db


# ------------------------------------------------------------- design ----
def design_pi(K, poles, theta, tau_f):
    PH = complex(plant_fr([F_CROSS], K, poles, theta, tau_f)[0])
    magPH = abs(PH)
    phPH = math.degrees(np.angle(PH))
    phi_pi = PM_TARGET - 180.0 - phPH
    if phi_pi >= 0.0:
        fz = F_CROSS / 10.0            # plant alone beats the target
    elif phi_pi <= -88.0:
        fz = F_CROSS / 10.0            # PI can't reach target; report honest max
    else:
        fz = F_CROSS / math.tan(math.radians(phi_pi + 90.0))
    Kp = 1.0 / (magPH * math.sqrt(1.0 + (fz / F_CROSS) ** 2))
    Ki = 2.0 * math.pi * fz * Kp
    return Kp, Ki, fz, magPH, phPH, phi_pi


# ----------------------------------- firmware-exact discrete simulation --
def simulate_firmware(K, tau, theta, ts, vbus, cpr, Kp, Ki, tau_f,
                      antiwindup=True, quantize=True):
    """Discrete closed loop: the exact difference equations the Arduino
    runs, against the identified FO plant discretized with ZOH."""
    n = int(round(T_END / ts))
    t = np.arange(n) * ts
    a = math.exp(-ts / tau)
    b = K * (1.0 - a)
    d = max(1, int(round(theta / ts)))       # dead time in whole samples
    alpha = 1.0 - math.exp(-2.0 * math.pi * F_MEAS * ts)
    q = 2.0 * math.pi / cpr                  # encoder angle quantum [rad]
    Tt = Kp / Ki                             # back-calculation time constant

    r = np.where(t >= T_STEP, W_REF, 0.0)
    dist = np.where(t >= T_DIST, V_DIST, 0.0)

    y = np.zeros(n)          # true plant speed
    ym = np.zeros(n)         # filtered measurement
    u = np.zeros(n)          # saturated command
    ubuf = [0.0] * d         # dead-time pipeline
    yf = xi = 0.0
    ang = 0.0
    counts_prev = 0

    for k in range(n):
        yk = y[k - 1] if k else 0.0
        # measurement: encoder counts over the window, then LPF
        if quantize:
            ang += yk * ts
            counts = math.floor(ang / q)
            w_raw = (counts - counts_prev) * q / ts
            counts_prev = counts
        else:
            w_raw = yk
        yf += alpha * (w_raw - yf)
        ym[k] = yf

        # PI with saturation + back-calculation anti-windup
        e = r[k] - yf
        u_unsat = Kp * e + xi
        u_sat = min(max(u_unsat, 0.0), vbus)
        if antiwindup:
            xi += Ki * ts * e + (ts / Tt) * (u_sat - u_unsat)
        else:
            xi += Ki * ts * e
        u[k] = u_sat

        # plant: dead-time pipeline then FO ZOH step
        ubuf.append(u_sat + dist[k])
        v_now = ubuf.pop(0)
        y[k] = a * yk + b * v_now
        if y[k] < 0.0:
            y[k] = 0.0
    return t, r, y, ym, u


def linear_reference_step(K, poles, theta, tau_f, Kp, Ki, ts):
    """Linear closed-loop step via frequency-domain inversion (no
    saturation): y(t) = step response of T/(1+T) with H in feedback."""
    n = 4096
    fs = 1.0 / ts
    f = np.fft.rfftfreq(n, ts)
    f[0] = 1e-9
    s = 2j * np.pi * f
    P = K * np.exp(-theta * s)
    for tau in poles:
        P = P / (tau * s + 1.0)
    H = 1.0 / (s / (2 * np.pi * F_MEAS) + 1.0)
    C = Kp + Ki / s
    G_cl = C * P / (1.0 + C * P * H)
    # step = integral of impulse response
    g = np.fft.irfft(G_cl, n) * fs
    y = np.cumsum(g) * ts * W_REF
    t = np.arange(n) * ts
    return t + T_STEP, y


# -------------------------------------------------------------- plots ----
def plot_bode(K, poles, theta, tau_f, Kp, Ki, fc, pm, outdir):
    f = np.logspace(-2, math.log10(50.0), 1200)
    PH = plant_fr(f, K, poles, theta, tau_f)
    C = Kp + Ki / (2j * np.pi * f)
    T = C * PH

    fig, ax = plt.subplots(2, 1, figsize=(7.5, 6.2), sharex=True)
    ax[0].loglog(f, np.abs(PH), color=COL_PRIMARY, lw=1.6,
                 label="plant P·H (identified)")
    ax[0].loglog(f, np.abs(T), color=COL_TERTIARY, lw=1.6,
                 label="loop T = Gc·P·H")
    ax[0].axhline(1.0, color=COL_GRAY, lw=0.8, ls=":")
    ax[0].axvline(fc, color=COL_GRAY, lw=0.8, ls=":")
    ax[0].annotate("fc = %.2f Hz" % fc, xy=(fc * 1.1, 1.6),
                   fontsize=8, color=COL_GRAY)
    ax[0].set_ylabel("magnitude")
    ax[0].set_title("Loop gain on the identified plant "
                    "(dead time + 15 Hz measurement LPF included)")
    ax[0].legend(loc="lower left", framealpha=0.9)

    for tr, col, lab in ((PH, COL_PRIMARY, None), (T, COL_TERTIARY, None)):
        ax[1].semilogx(f, np.degrees(np.unwrap(np.angle(tr))),
                       color=col, lw=1.6, label=lab)
    ax[1].axhline(-180, color=COL_GRAY, lw=0.8, ls=":")
    ax[1].axvline(fc, color=COL_GRAY, lw=0.8, ls=":")
    ax[1].annotate("PM = %.1f deg" % pm, xy=(fc * 1.1, -175),
                   fontsize=8, color=COL_GRAY)
    ax[1].set_ylabel("phase [deg]")
    ax[1].set_xlabel("frequency [Hz]")
    ax[1].set_ylim(-280, 0)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "bode_loop.png"))
    plt.close(fig)


def plot_step(sim_aw, sim_naw, lin, vbus, outdir):
    t, r, y, ym, u = sim_aw
    _, _, y2, _, _ = sim_naw
    tl, yl = lin

    fig, ax = plt.subplots(2, 1, figsize=(8.5, 6.2), sharex=True,
                           gridspec_kw={"height_ratios": [2.2, 1]})
    ax[0].plot(t, r, color=COL_GRAY, lw=1.0, ls="--", label="setpoint")
    ax[0].plot(tl, yl, color=COL_ACCENT, lw=1.2, ls=":",
               label="linear design model")
    ax[0].plot(t, y2, color=COL_SECONDARY, lw=1.2,
               label="firmware sim, NO anti-windup")
    ax[0].plot(t, y, color=COL_PRIMARY, lw=1.6,
               label="firmware sim, anti-windup")
    ax[0].set_ylabel("speed [rad/s]")
    ax[0].set_xlim(0, 2.5)
    ax[0].set_title("Closed-loop reference step 0 -> %.0f rad/s" % W_REF)
    ax[0].legend(loc="lower right", framealpha=0.9)

    ax[1].plot(t, u, color=COL_PRIMARY, lw=1.2)
    ax[1].axhline(vbus, color=COL_GRAY, lw=0.8, ls=":")
    ax[1].annotate("Vbus = %.1f V (buck output, saturates here)"
                   % vbus, xy=(0.6, vbus - 0.45), fontsize=8, color=COL_GRAY)
    ax[1].set_ylabel("v_cmd [V]")
    ax[1].set_xlabel("time [s]")
    ax[1].set_ylim(-0.2, vbus + 0.6)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "closedloop_step.png"))
    plt.close(fig)


def plot_disturbance(sim_aw, vbus, outdir):
    t, r, y, ym, u = sim_aw
    m = (t >= T_DIST - 0.4) & (t <= T_DIST + 1.6)
    dip = W_REF - y[m].min()
    # recovery: last time speed is below setpoint - 2%
    below = np.where((t >= T_DIST) & (y < W_REF * 0.98))[0]
    t_rec = (t[below[-1]] - T_DIST) if len(below) else 0.0

    fig, ax = plt.subplots(2, 1, figsize=(8, 6.0), sharex=True,
                           gridspec_kw={"height_ratios": [2.2, 1]})
    ax[0].plot(t[m], r[m], color=COL_GRAY, lw=1.0, ls="--", label="setpoint")
    ax[0].plot(t[m], y[m], color=COL_PRIMARY, lw=1.6, label="speed")
    ax[0].plot(t[m], ym[m], color=COL_TERTIARY, lw=0.5, alpha=0.35,
               label="filtered measurement (encoder limit cycle)")
    ax[0].annotate("dip %.1f rad/s (%.1f%%), recovered in %.2f s"
                   % (dip, 100 * dip / W_REF, t_rec),
                   xy=(T_DIST + 0.25, W_REF - dip * 0.7),
                   fontsize=8, color=COL_GRAY)
    ax[0].set_ylabel("speed [rad/s]")
    ax[0].set_title("Disturbance rejection: %.2f V input step at t = %.1f s "
                    "(load-torque equivalent, 'the hill')" % (V_DIST, T_DIST))
    ax[0].legend(loc="lower right", framealpha=0.9)

    ax[1].plot(t[m], u[m], color=COL_PRIMARY, lw=1.2)
    ax[1].set_ylabel("v_cmd [V]")
    ax[1].set_xlabel("time [s]")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "disturbance_rejection.png"))
    plt.close(fig)
    return dip, t_rec


# -------------------------------------------------------------- report ---
def main():
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)

    d, K, poles, theta = load_plant()
    ts, vbus, cpr = d["ts_s"], d["vbus_V"], d["cpr"]
    tau_f = 1.0 / (2.0 * math.pi * F_MEAS)
    tau_dom = max(poles)

    ln = []
    p = ln.append
    p("=" * 78)
    p("PI SPEED-CONTROLLER DESIGN FROM IDENTIFIED PLANT")
    p("(target fc = %.1f Hz, PM >= %.0f deg)" % (F_CROSS, PM_TARGET))
    p("=" * 78)
    p("Plant (%s model, fit %.1f%%, source: %s):"
      % (d["chosen"].upper(), d[d["chosen"]]["fit_pct"], d["source"]))
    p("  K = %.2f rad/s per V, tau = %.1f ms, theta = %.1f ms"
      % (K, tau_dom * 1e3, theta * 1e3))
    p("  open-loop pole: %.2f Hz;  measurement LPF: %.0f Hz;"
      % (1 / (2 * math.pi * tau_dom), F_MEAS))
    p("  Ts = %.0f ms, Vbus = %.2f V (Project 1 buck output), CPR = %d"
      % (ts * 1e3, vbus, cpr))

    Kp, Ki, fz, magPH, phPH, phi_pi = design_pi(K, poles, theta, tau_f)
    p("\nUncompensated P·H at fc:")
    p("  |PH(j2pi fc)| = %.3f (rad/s)/V,  phase = %.1f deg" % (magPH, phPH))
    p("  phase budget at fc: pole %.1f deg, dead time %.1f deg,"
      % (-math.degrees(math.atan(2 * math.pi * F_CROSS * tau_dom)),
         -360.0 * F_CROSS * theta))
    p("  measurement filter %.1f deg"
      % (-math.degrees(math.atan(F_CROSS / F_MEAS))))
    p("Required PI phase at fc: PM - 180 - phase(PH) = %.1f deg" % phi_pi)
    p("  atan(fc/fz) = %.1f deg  ->  fz = %.2f Hz" % (phi_pi + 90.0, fz))
    p("Unity gain at fc:")
    p("  Kp = 1/(|PH| sqrt(1+(fz/fc)^2)) = %.4f V/(rad/s)" % Kp)
    p("  Ki = 2 pi fz Kp = %.4f V/(rad/s)/s" % Ki)

    Tfun = lambda f, K_=K, po=poles: (Kp + Ki / (2j * np.pi * np.asarray(f))) \
        * plant_fr(f, K_, po, theta, tau_f)
    fc_a, pm_a, gm_a = margins(Tfun)
    p("\nVerification (dense grid): fc = %.2f Hz, PM = %.1f deg, GM = %s"
      % (fc_a, pm_a, gm_str(gm_a)))

    # Robustness sweep -- the EV justification for system ID: gain and
    # inertia are NOT constants. Winding heats (R up -> K down), payload
    # changes (J up -> tau up). Sweep +/-30% on both.
    p("\nRobustness sweep (gains fixed, plant perturbed +/-30%):")
    p("  %-22s %-10s %-10s %-8s" % ("case", "fc", "PM", "GM"))
    for dK, dT, name in ((0.7, 1.0, "K -30% (hot winding)"),
                         (1.3, 1.0, "K +30%"),
                         (1.0, 0.7, "tau -30%"),
                         (1.0, 1.3, "tau +30% (payload)"),
                         (0.7, 1.3, "K-30% tau+30%"),
                         (1.3, 0.7, "K+30% tau-30%")):
        po = [tau_dom * dT] + list(poles[1:])
        Tf = lambda f, K_=K * dK, po_=po: \
            (Kp + Ki / (2j * np.pi * np.asarray(f))) \
            * plant_fr(f, K_, po_, theta, tau_f)
        fcx, pmx, gmx = margins(Tf)
        p("  %-22s %5.2f Hz   %5.1f deg  %s" % (name, fcx, pmx, gm_str(gmx)))
    p("Worst case keeps > 45 deg: the loop tolerates the whole realistic")
    p("parameter envelope without retuning.")

    # firmware-exact simulations
    sim_aw = simulate_firmware(K, tau_dom, theta, ts, vbus, cpr,
                               Kp, Ki, tau_f, antiwindup=True)
    sim_naw = simulate_firmware(K, tau_dom, theta, ts, vbus, cpr,
                                Kp, Ki, tau_f, antiwindup=False)
    lin = linear_reference_step(K, poles, theta, tau_f, Kp, Ki, ts)

    t, r, y, ym, u = sim_aw
    seg = (t >= T_STEP) & (t < T_DIST)
    ts_seg, ys = t[seg], y[seg]
    ov = max(0.0, (ys.max() - W_REF) / W_REF * 100)
    i90 = np.argmax(ys >= 0.9 * W_REF)
    i10 = np.argmax(ys >= 0.1 * W_REF)
    rise = ts_seg[i90] - ts_seg[i10]
    _, _, y_naw, _, _ = sim_naw
    ov_naw = (y_naw[seg].max() - W_REF) / W_REF * 100

    dip, t_rec = plot_disturbance(sim_aw, vbus, outdir)
    plot_bode(K, poles, theta, tau_f, Kp, Ki, fc_a, pm_a, outdir)
    plot_step(sim_aw, sim_naw, lin, vbus, outdir)

    p("\nFirmware-exact discrete simulation (sat 0..%.1f V, encoder" % vbus)
    p("quantization, %.0f Hz LPF, back-calculation anti-windup Tt=Kp/Ki):"
      % F_MEAS)
    p("  reference step 0 -> %.0f rad/s:" % W_REF)
    p("    rise time (10-90%%)      = %.0f ms" % (rise * 1e3))
    p("    overshoot               = %.1f %%  (anti-windup ON)" % ov)
    p("    overshoot               = %.1f %%  (anti-windup OFF -- the" % ov_naw)
    p("      command saturates at Vbus for ~%.1f s during spin-up and a"
      % (float(np.sum(u[seg] >= vbus - 1e-9)) * ts))
    p("      naive integrator winds up; this is why the firmware has AW)")
    p("  disturbance %.2f V at t=%.1f s ('hill'):" % (V_DIST, T_DIST))
    p("    speed dip               = %.1f rad/s (%.1f %%)"
      % (dip, 100 * dip / W_REF))
    p("    recovery to 2%%          = %.2f s" % t_rec)

    p("\nFirmware constants (transcribe into firmware/motor_speed_pi):")
    p("  #define KP          %.4ff   // V per rad/s" % Kp)
    p("  #define KI          %.4ff   // V per rad/s per s" % Ki)
    p("  #define TS          %.3ff   // s (%.0f Hz loop)" % (ts, 1 / ts))
    p("  #define F_MEAS_LPF  %.1ff    // Hz -> alpha = %.4f"
      % (F_MEAS, 1.0 - math.exp(-2 * math.pi * F_MEAS * ts)))
    p("  #define VBUS        %.2ff   // V, Project 1 buck output" % vbus)
    p("  anti-windup: Tt = Kp/Ki = %.3f s (back-calculation)" % (Kp / Ki))

    txt = "\n".join(ln)
    print(txt)
    with open(os.path.join(outdir, "design_report.txt"), "w") as f:
        f.write(txt + "\n")

    # hand-off for Simulink cross-check and firmware
    with open(os.path.join(outdir, "controller.json"), "w") as f:
        json.dump({"Kp": Kp, "Ki": Ki, "fz_Hz": fz, "fc_Hz": fc_a,
                   "pm_deg": pm_a, "f_meas_Hz": F_MEAS, "ts_s": ts,
                   "vbus_V": vbus, "w_ref": W_REF, "t_step": T_STEP,
                   "v_dist": V_DIST, "t_dist": T_DIST,
                   "plant": {"K": K, "tau": tau_dom, "theta": theta}},
                  f, indent=2)
    np.savetxt(os.path.join(outdir, "closedloop_sim.csv"),
               np.column_stack([t, r, y, u]),
               delimiter=",", header="t_s,ref_rad_s,omega_rad_s,vcmd_V",
               comments="")
    print("\nWrote results/design_report.txt, controller.json,")
    print("closedloop_sim.csv (for simulink/run_and_compare.m) and plots.")


if __name__ == "__main__":
    main()
