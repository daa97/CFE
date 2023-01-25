
import numpy as np
import matplotlib.pyplot as plt
from fluids import H2

import matplotlib as mpl
mpl.rc('font', family='Times New Roman',size="10")
mpl.rc('figure', figsize=(4.8,3.6))
mpl.rc('savefig', dpi=800)
mpl.rc('lines', linewidth=1.2)
mpl.rc('axes', grid=True)
mpl.rc('grid', linewidth=0.25)
mpl.rc('mathtext', fontset="dejavuserif")
mpl.rc('xtick.minor', visible=True, size=1.5, width=0.5)
mpl.rc('ytick.minor', visible=True, size=1.5, width=0.5)
plt.rcParams['figure.constrained_layout.use'] =  True


def lewis_Nu(Re_w):
    R = np.log10(Re_w)
    lewis_eta = 0.15999/0.22085
    lewis_G_factor = 4*np.pi*lewis_eta/((1-lewis_eta)*(1-lewis_eta**2))
    C0 = np.log10(lewis_G_factor)
    print(C0)
    if 2600<=Re_w<=13e3:
        exp = .2005*R**3 - 1.970*R**2 + (7.775-1)*R - 5.516-C0
    elif 13e3<Re_w<=1e6:
        exp = -.006360*R**3 + .1349*R**2 + (.8850-1)*R + 1.610-C0
    else:
        raise ValueError("Rotational Reynolds number out of range!")
    return 10**exp


def TTCPF_torque(R1, R2, L, omega, state, mdot):
    nu = state.mu/state.rho
    Q = mdot/state.rho
    A = np.pi * (R2**2 - R1**2)
    U = Q/A
    Re_w = Re_omega(omega, R1, R2, nu)
    Re_a = Re_axial(U, R1, R2, nu)
    Nu = 1.1* lewis_Nu(Re_w)
    Mlam = 4*np.pi*state.mu*L*omega / (R1**(-2) - R2**(-2))
    M = Nu*Mlam
    F = M/R1
    shear = F / (2*np.pi*R1*L)
    print(f"Rotational Nusselt Number: {Nu:.2f}")
    print(f"Azimuthal wall shear: {shear:.3f} Pa")
    print(f"Laminar azimuthal wall shear: {shear/Nu:.3f} Pa")
    return M

def Re_omega(omega, R1, R2, nu):
    return omega * R1 * (R2 - R1) / nu

def Re_axial(U, R1, R2, nu):
    return U * (R2 - R1) / nu

def bearing_torque(D,F_axial,F_radial):
    mu_f = .0025
    F_min = 450
    F = np.sqrt(F_axial**2+F_radial**2)
    return mu_f*(F+F_min)*D/2
state = H2(P=12.5e6, T=450)

M = TTCPF_torque(.056, .058, .95, 733, state, .108)
print(f"TTCPF Torque: {M:.5f} N*m")


# vals = []
# rmin = .05
# rmax = .066
# tmin = .0005
# tmax = .010
# rs = np.arange(rmin, rmax, (rmax-rmin)/4)
# ts = np.arange(tmin, tmax, (tmax-tmin)/100)

# for ri in rs:
#     vals.append([])
#     for ti in ts:
#         A = TTCPF_torque(ri, ri+ti, .84, 7000*np.pi/30, state, 0.108)
#         vals[-1].append(A)

# for i in range(len(vals)):
#     plt.plot(ts*1000, vals[i], label=f"$r_i$ = {rs[i]*1000:.0f} mm")
# plt.xlabel("Channel Width $r_o - r_i$ (mm)")
# plt.ylabel("Power Dissipated (W)")

# plt.legend()
# plt.savefig("TTCPF Power")
# plt.show()