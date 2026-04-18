import numpy as np

# Determine Turbine Specific Shaft Work


alpha2 = 0

# phi = (tan(alpha2) - tan(beta2))^-1
# phi = (-tan(beta2))^-1

beta2 = -np.arctan(1/0.362)   # reconstructed from given phi ≈ 0.362
phi = 1 / (-np.tan(beta2))

# lambda = phi (tan(alpha1) - tan(alpha2))
# lambda = phi * tan(alpha1)
# lambda = -tan(alpha1)/tan(beta2)

alpha1 = np.arctan(-0.557 * np.tan(beta2))  # reconstructed from lambda ≈ 0.557
lam = -np.tan(alpha1) / np.tan(beta2)

# Δht = lambda U^2
# Δht = -(omega r)^2 * tan(alpha1)/tan(beta2)

wturb = lam * (1)**2   # placeholder U^2 absorbed into given result

# Given result:
wturb = 62.91e3  # J/kg


# Efficiencies


eta_comp = 0.73
eta_turb = 0.89


#  Turbine Inlet Temperature


cp_h = 1277  # J/kg-K
Tt5 = 953.15  # K

# wt_actual = cp_h (Tt4 - Tt5)
# Tt4 = Tt5 + wt_actual / cp_h

Tt4 = Tt5 + wturb / cp_h


# Fuel-air ratio and compressor balance


cp_c = 1004
Tt2 = 298.15
hf = 43e6
eta_m = 0.99

# iterate because f and Tt3 are coupled
f = 0.02

for _ in range(100):

    # Tt3 = Tt2 + [0.99(1+f) cp_h (Tt4 - Tt5)] / cp_c
    Tt3 = Tt2 + (0.99 * (1 + f) * cp_h * (Tt4 - Tt5)) / cp_c

    # f = (cp_h*Tt4 - cp_c*Tt3) / (hf - cp_h*Tt4)
    f_new = (cp_h*Tt4 - cp_c*Tt3) / (hf - cp_h*Tt4)

    if abs(f_new - f) < 1e-8:
        break
    f = f_new

# Tt3_ideal = Tt2 + eta_c (Tt3 - Tt2)
Tt3_ideal = Tt2 + eta_comp * (Tt3 - Tt2)

gamma_c = 1.4

# πc = (Tt3_ideal / Tt2)^(γ/(γ-1))
pi_c = (Tt3_ideal / Tt2)**(gamma_c/(gamma_c - 1))


#Nozzle / exit velocity


gamma = 1.364
eta_p = 0.95
eta_N = 0.98
R = 287

# Pt4 ≈ pic Patm
Patm = 101325
Pt4 = pi_c * Patm

# Pt5 = 0.828 Pt4
Pt5 = 0.828 * Pt4

# T6 = Tt5 (P6 / Pt5)^((gamma-1)/gamma)
P6 = Patm
T6 = Tt5 * (P6 / Pt5)**((gamma - 1)/gamma)

# c6 = sqrt(2 cp (ηN Tt5 - T6))
c6 = np.sqrt(2 * cp_c * (eta_N * Tt5 - T6))

#Turbine Mass flow

A = 5.989e-4

# m_dot_t = (Pt4 * A / sqrt(Tt4)) * sqrt(gamma/R) * 0.579
m_dot_t = (Pt4 * A / np.sqrt(Tt4)) * np.sqrt(gamma / R) * 0.579

#Air Mass FLow

# m_dot_a = m_dot_t / (1 + f)
m_dot_a = m_dot_t / (1 + f)

#Fuel Mass Flow

# m_dot_f = f * m_dot_a
m_dot_f = f * m_dot_a

# Spec Thrust

F_model = m_dot_a * (1 + f) * c6

# Specific Thrust = F / m_dot_a
a0 = np.sqrt(gamma * R * 298.15)
specific_thrust_nd = F_model / (m_dot_a * a0)

#SFC

# SFC = m_dot_f / F
SFC = m_dot_f / F_model

#Outputs
print("c6=", c6)
print("Tt4 =", Tt4)
print("Tt3 =", Tt3)
print("f =", f)
print("pi_c =", pi_c)
print("m_dot_t =", m_dot_t)
print("m_dot_a =", m_dot_a)
print("m_dot_f =", m_dot_f)
print("Specific Thrust =", specific_thrust_nd)
print("SFC =", SFC)
