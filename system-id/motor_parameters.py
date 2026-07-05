"""Ground-truth motor used ONLY by the synthetic-data generator.

============================  READ THIS FIRST  ============================
The identification pipeline (fit_plant.py) must never import the physical
constants below. They exist so that generate_synthetic_data.py can play the
role of the real bench: a motor with parasitics the datasheet doesn't tell
you about (Coulomb friction, brush drop, encoder quantization, loop
latency). fit_plant.py sees only the CSV, exactly as it will see your real
bench log. If the pipeline recovers a good model *without* peeking, it will
also work on the real motor.
===========================================================================

The "true" motor is modeled on a small brushed DC gearmotor of the Pololu
25D class (6 V winding run from the 5 V buck rail, 48 CPR quadrature
encoder on the motor shaft). Chosen so the worst case fits Project 1's
buck output:

    stall current at 5 V = VBUS / R = 5.0 / 2.7 = 1.85 A  <  2.0 A budget

Electromechanical model integrated by the generator (10 kHz inner step):

    L  di/dt = v_applied - R i - Ke w        v_applied = v_cmd - Vbrush*sgn(i)
    J  dw/dt = Kt i - b w - tau_c*sgn(w) - tau_load

Linearized (no Coulomb/brush terms) this is the classic second-order plant

    w(s)/v(s) = Kt / ( (J s + b)(L s + R) + Kt Ke )

whose DC gain / poles are printed by `python3 motor_parameters.py`. Those
numbers are the *answer key* for judging the fit -- not an input to it.
"""

import math

# ---------------------------------------------------------------- rig ----
VBUS = 5.0        # motor supply [V] = Project 1 buck output (12V->5V, 2A)
TS = 0.010        # controller / logging sample period [s]  (100 Hz)
CPR = 48          # encoder counts per motor-shaft rev (quadrature-decoded)
LATENCY_SAMPLES = 1   # command-to-PWM latency of the firmware loop [samples]

# ------------------------------------------------- true motor (hidden) ---
R = 2.7           # winding resistance [ohm]
L = 3.0e-3        # winding inductance [H]
KT = 0.010        # torque constant [N*m/A]
KE = 0.010        # back-EMF constant [V*s/rad]  (= KT in SI)
J = 1.0e-5        # rotor + flywheel inertia [kg*m^2]
B = 5.0e-6        # viscous friction [N*m*s/rad]
TAU_C = 2.0e-4    # Coulomb (kinetic) friction torque [N*m]
STICTION = 1.2    # breakaway torque = STICTION * TAU_C
V_BRUSH = 0.15    # brush contact drop [V]

# ------------------------------------------------ sysid input sequence ---
# Staircase of voltage steps: multiple levels up AND down so the fit sees
# several operating points (Coulomb friction makes the low end nonlinear).
# Each level is held 2.5 s ~ 10 mechanical time constants.
STEP_LEVELS = [0.0, 1.5, 3.0, 4.5, 2.5, 5.0, 1.0, 3.5, 0.0]   # [V]
STEP_HOLD = 2.5                                               # [s]

# ------------------------------------- plot palette (same as Project 1) --
COL_PRIMARY = "#2E6FB7"    # measured data / primary trace
COL_SECONDARY = "#D9782D"  # first-order fit / comparison trace
COL_TERTIARY = "#3E8E5A"   # second-order fit / third series
COL_ACCENT = "#8A5CB8"     # fourth series
COL_GRAY = "#6E6E6E"       # annotations, limits, references


def linearized_truth():
    """Answer key: DC gain, poles and equivalent (K, tau) of the true motor."""
    den0 = R * B + KT * KE                 # s^0 coefficient
    K_dc = KT / den0                       # [rad/s per V]
    # (J s + b)(L s + R) + Kt Ke = JL s^2 + (JR + bL) s + (bR + KtKe)
    a2, a1, a0 = J * L, J * R + B * L, den0
    disc = math.sqrt(a1 * a1 - 4 * a2 * a0)
    p1, p2 = (-a1 + disc) / (2 * a2), (-a1 - disc) / (2 * a2)
    return {
        "K_dc": K_dc,
        "tau_mech": -1.0 / p1,             # dominant pole
        "tau_elec": -1.0 / p2,             # fast electrical pole
        "no_load_speed_5V": K_dc * VBUS,
    }


if __name__ == "__main__":
    t = linearized_truth()
    print("Linearized ground truth (answer key for the fit):")
    print("  DC gain          K      = %.1f rad/s per V" % t["K_dc"])
    print("  mechanical pole  tau_m  = %.1f ms" % (t["tau_mech"] * 1e3))
    print("  electrical pole  tau_e  = %.2f ms  (%.1fx faster than Ts=%.0f ms"
          % (t["tau_elec"] * 1e3, TS / t["tau_elec"], TS * 1e3)
          + " -- expect it to be invisible to the fit)")
    print("  no-load speed at %.1f V = %.0f rad/s (%.0f rpm)"
          % (VBUS, t["no_load_speed_5V"], t["no_load_speed_5V"] * 60 / (2 * math.pi)))
    print("  stall current at %.1f V = %.2f A (buck budget 2.0 A)"
          % (VBUS, VBUS / R))
    print("  encoder speed resolution at Ts: %.1f rad/s per count"
          % (2 * math.pi / CPR / TS))
