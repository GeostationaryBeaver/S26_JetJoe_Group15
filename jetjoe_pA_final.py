import math
from dataclasses import dataclass


@dataclass
class Inputs:
    # Ambient / cycle states
    Tt2: float = 288          # K
    Patm: float = 101325.0    # Pa
    Tt5: float = 953.15       # K

    # Gas properties
    cp_c: float = 1004.5
    cp_h: float = 1277.0
    gamma_c: float = 1.4
    gamma_h: float = 1.3
    R: float = 287.0

    # Efficiencies
    eta_comp: float = 0.75
    eta_turb: float = 0.89
    eta_mech: float = 0.99
    eta_p_comb: float = 1

    # Fuel / geometry
    Q_LHV: float = 43e6
    A_stator: float = 5.915e-4

    # Mean-line turbine variables
    rpm: float = 160000.0      # rev/min
    r_mean: float = 20.05e-3   # m
    alpha1_deg: float = 57.0   # deg
    beta2_deg: float = -70.1   # deg

    # Target
    F_target: float = 62.3


@dataclass
class Results:
    wturb: float
    Tt4: float
    f: float
    Tt3: float
    pi_c: float
    Pt3: float
    Pt4: float
    Pt5: float
    M6: float
    T6: float
    c6: float
    m_dot_t: float
    thrust: float
    spec_thrust_nd: float
    sfc: float


def solve_cycle(inp: Inputs) -> Results:

    # 0) Mean-line turbine specific work
    omega = 2.0 * math.pi * inp.rpm / 60.0
    alpha1 = math.radians(inp.alpha1_deg)
    beta2 = math.radians(inp.beta2_deg)
    wturb = -(omega * inp.r_mean) ** 2 * math.tan(alpha1) / math.tan(beta2)

    # 1) Turbine inlet temp
    Tt4 = inp.Tt5 + wturb / inp.cp_h

    # 2) Solve f and Tt3
    f = 0.022
    for _ in range(200):
        Tt3 = inp.Tt2 + inp.eta_mech * (1 + f) * inp.cp_h * (Tt4 - inp.Tt5) / inp.cp_c
        f_new = (inp.cp_h * Tt4 - inp.cp_c * Tt3) / (inp.Q_LHV - inp.cp_h * Tt4)
        if abs(f_new - f) < 1e-10:
            f = f_new
            break
        f = f_new

    # 3) Compressor pressure ratio
    Tt3_ideal = inp.Tt2 + inp.eta_comp * (Tt3 - inp.Tt2)
    pi_c = (Tt3_ideal / inp.Tt2) ** (inp.gamma_c / (inp.gamma_c - 1))
    Pt3 = pi_c * inp.Patm
    Pt4 = Pt3  # assume 100% combustor efficiency

    # 4) Turbine pressure drop
    Tt5s = Tt4 - (Tt4 - inp.Tt5) / inp.eta_turb
    Pt5 = Pt4 * (Tt5s / Tt4) ** (inp.gamma_h / (inp.gamma_h - 1))

    # 5) Find c6
    P6 = inp.Patm
    Pt6 = Pt5

    pressure_ratio = Pt6 / P6

    M6_sq = (2 / (inp.gamma_h - 1)) * (
        pressure_ratio ** ((inp.gamma_h - 1) / inp.gamma_h) - 1
    )
    M6 = math.sqrt(max(M6_sq, 0.0))
    T6 = inp.Tt5 / (1 + (inp.gamma_h - 1) / 2 * M6**2)
    c6 = M6 * math.sqrt(inp.gamma_h * inp.R * T6)

    # total mass flow
    critical = (2 / (inp.gamma_h + 1)) ** (
        (inp.gamma_h + 1) / (2 * (inp.gamma_h - 1))
    )

    m_dot_t = (
        (Pt4 * inp.A_stator / math.sqrt(Tt4))
        * math.sqrt(inp.gamma_h / inp.R)
        * critical
    )

    # thrust
    thrust = m_dot_t * c6

    # Non-dimensional specific thrust
    a0 = math.sqrt(inp.gamma_c * inp.R * inp.Tt2)
    m_dot_a = m_dot_t / (1 + f)
    spec_thrust_nd = thrust / (m_dot_a * a0)

    # SFC
    m_dot_f = f * m_dot_a
    sfc = m_dot_f / thrust

    return Results(
        wturb=wturb,
        Tt4=Tt4,
        f=f,
        Tt3=Tt3,
        pi_c=pi_c,
        Pt3=Pt3,
        Pt4=Pt4,
        Pt5=Pt5,
        M6=M6,
        T6=T6,
        c6=c6,
        m_dot_t=m_dot_t,
        thrust=thrust,
        spec_thrust_nd=spec_thrust_nd,
        sfc=sfc
    )


def main():
    inp = Inputs()
    r = solve_cycle(inp)

    print("=== Model Version 500000 ===")
    print(f"w_turb = {r.wturb/1000:.2f} kJ/kg")
    print(f"Tt4 = {r.Tt4:.2f} K")
    print(f"f = {r.f:.5f}")
    print(f"Tt3 = {r.Tt3:.2f} K")
    print(f"pi_c = {r.pi_c:.3f}")
    print(f"Pt5 = {r.Pt5:.0f} Pa")
    print(f"M6 = {r.M6:.3f}")
    print(f"T6 = {r.T6:.2f} K")
    print(f"c6 = {r.c6:.2f} m/s")
    print(f"m_dot_t = {r.m_dot_t:.5f} kg/s")
    print(f"Thrust = {r.thrust:.2f} N")
    print(f"Spec thrust (nd) = {r.spec_thrust_nd:.3f}")
    print(f"SFC = {r.sfc:.6e} kg/(N*s)")
    print(f"Error vs 62.3N = {r.thrust - inp.F_target:.2f} N")


if __name__ == "__main__":
    main()

# Note: this model has been updated many times with different approaches to reach a better c6.
# we have concluded that the lower c6 is likely linked to the initial geometric parameters.
