# optimize_f_param.py
# Minimal F-parameterization optimizer for ring electrodes
import numpy as np
import nevergrad as ng
import matplotlib.pyplot as plt
from weak_e_reversal_opt import AB_from_params  # uses your established A,B definitions

# --------- Constants / simple helpers ---------
K_E_PER_VSTEP = 0.42  # [V/cm per V] => E_i = k * (V_{i+} - V_{i-})
EPS = 1e-12

def build_uniform_times(N_gaps, d_gap, v):
    return np.arange(N_gaps, dtype=float) * (d_gap / v)

def build_S(t_gap, tau):
    """Sech kernel S_ij = sech((t_i - t_j)/tau)."""
    diff = (t_gap[:, None] - t_gap[None, :]) / float(tau)
    return 1.0 / np.cosh(diff)

def center_selector(N_gaps, center_count):
    """Selection matrix R ∈ R^{N_gaps × center_count} that picks the central gaps."""
    mid0 = (N_gaps - center_count) // 2
    mid1 = mid0 + center_count
    R = np.zeros((N_gaps, center_count), dtype=float)
    for k, i in enumerate(range(mid0, mid1)):
        R[i, k] = 1.0
    return R, (mid0, mid1)

def nullspace_basis_1d(a):
    """Orthonormal basis B for {x : a^T x = 0}, with a shape (m,)."""
    a = np.asarray(a, dtype=float).reshape(1, -1)  # (1, m)
    # SVD of a^T; rows 1..m-1 of V^T span the nullspace
    _, _, Vt = np.linalg.svd(a, full_matrices=True)
    return Vt[1:, :].T  # (m, m-1)

def solve_E_from_F(S, F):
    """Solve S E = F (no explicit inverse)."""
    return np.linalg.solve(S, F)

def voltages_from_E(E, k):
    """Reconstruct ring voltages with V0=0; V_j = (1/k) * sum_{m<j} E_m."""
    V = np.zeros(len(E) + 1, dtype=float)
    V[1:] = np.cumsum(E) / float(k)
    return V

def snr_from_E(E, A, B, w, k, eps_noise=0.0):
    """SNR from E. Since dV = E/k, we have C = (w^T E)/k."""
    C = (w @ E) / float(k)
    gamma = A * np.conj(B)
    num = 4.0 * abs(np.real(gamma * np.conj(C)))          # |4 Re(AB* C*)|
    denom = np.sqrt(2.0 * (abs(A)**2 + abs(B)**2 * abs(C)**2) + eps_noise**2)
    return num / (denom + 1e-300)

def build_u_and_projector(S):
    """Return u solving S u = 1 and a projector P(y) = y - u*(u^T y)/(u^T u)
    so that (1^T S^{-1}) F = 0 is enforced via F = P(y).
    """
    N_gaps = S.shape[0]
    ones = np.ones(N_gaps, dtype=float)
    u = np.linalg.solve(S, ones)  # u = S^{-1} 1
    uTu = float(u @ u)
    def project(y):
        y = np.asarray(y, dtype=float)
        return y - u * float(u @ y) / (uTu + EPS)
    return u, project


def make_center_mask(N_gaps, center_count):
    mid0 = (N_gaps - center_count) // 2
    mid1 = mid0 + center_count
    mask = np.zeros(N_gaps, dtype=bool)
    mask[mid0:mid1] = True
    return mask, (mid0, mid1)


def reconstruct_all_from_F(F, S, k, Delta, W, T, d_dip, tau, t_gap):
    """Given F (gap fields, V/cm), reconstruct E, dV, V, and SNR pieces."""
    E = solve_E_from_F(S, F)
    V = voltages_from_E(E, k)
    w = np.exp(-1j * Delta * t_gap)
    A, B = AB_from_params(W, Delta, T, d_dip, tau)
    C = (w @ E) / float(k)
    gamma = A * np.conj(B)
    num = 4.0 * abs(np.real(gamma * np.conj(C)))
    denom = np.sqrt(2.0 * (abs(A)**2 + abs(B)**2 * abs(C)**2))
    return E, V, E / float(k), C, num, denom

def optimize_f_param(
    N_rings=30,
    total_length=183.2e-3,   # [m]
    v=600.0,                 # [m/s]
    Delta=1500.0*2*np.pi,    # [rad/s]
    W=2*np.pi*5,             # [rad/s]
    T=None,                  # [s] total interaction window; if None, infer from geometry
    d_dip=-33.6*2*np.pi,     # (your existing model’s dipole factor)
    tau=None,                # default: d_gap/v
    center_count=8,          # number of active central gaps (looser bounds)
    Vmax=5.0,                # [V] hard bound for each ring
    eps_field=0.02,          # [V/cm] tolerance outside center
    eps_noise=0.2,           # noise floor in SNR denominator
    budget=4000,             # Nevergrad eval budget
    seed=42,                 # RNG seed
    optimizer_name="NGOpt",  # or "CMA", "PSO", ...
    center_F_bound_factor=4.0  # center bound ≈ factor * (k * Vmax)
):
    # geometry/time grid
    N_gaps = N_rings - 1
    d_gap = total_length / N_rings
    if tau is None:
        tau = d_gap / v
    t_gap = build_uniform_times(N_gaps, d_gap, v)
    if T is None:
        # Approx: length spanning all gaps; center of window not critical for SNR constants
        T = (t_gap[-1] - t_gap[0]) + 8*tau

    # kernel and linear objects
    S = build_S(t_gap, tau)         # N_gaps x N_gaps
    k = K_E_PER_VSTEP                # E = k * dV

    # phase weights and physics constants
    w = np.exp(-1j * Delta * t_gap)
    A, B = AB_from_params(W, Delta, T, d_dip, tau)

    # build projector to enforce 1^T S^{-1} F = 0 (i.e., V_N=0 with V_0=0)
    u, project = build_u_and_projector(S)

    # masks and per-gap bounds in F-space (V/cm)
    center_mask, (mid0, mid1) = make_center_mask(N_gaps, center_count)
    F_center_max = center_F_bound_factor * (k * Vmax)  # heuristic magnitude scale

    lower = np.where(center_mask, -F_center_max, -eps_field).astype(float)
    upper = np.where(center_mask,  F_center_max,  eps_field).astype(float)

    # Decision variable: y ∈ R^{N_gaps}; actual F is its projection onto the constraint plane
    init_y = np.zeros(N_gaps, dtype=float)
    param = ng.p.Array(init=init_y).set_bounds(lower, upper)

    # Hard constraints via cheap constraints
    def _field_tolerance_ok(y):
        F = project(y)
        return float(np.max(np.abs(F[~center_mask]))) <= (eps_field + EPS)

    def _voltage_ok(y):
        F = project(y)
        E = solve_E_from_F(S, F)
        V = voltages_from_E(E, k)  # V0=0; VN follows from projection
        return float(np.max(np.abs(V))) <= (Vmax + EPS)

    param.register_cheap_constraint(_field_tolerance_ok)
    param.register_cheap_constraint(_voltage_ok)

    # optimizer
    Optimizer = getattr(ng.optimizers, optimizer_name)
    opt = Optimizer(parametrization=param, budget=budget, num_workers=1)
    opt.random_state = np.random.RandomState(seed)

    # objective: maximize SNR => minimize -SNR
    def loss_fn(y):
        F = project(y)
        E = solve_E_from_F(S, F)
        snr = snr_from_E(E, A, B, w, k, eps_noise)
        return -float(snr)

    rec = opt.minimize(loss_fn)
    y_best = np.asarray(rec.value, dtype=float)

    # reconstruct optimal objects
    F_opt = project(y_best)
    E_opt = solve_E_from_F(S, F_opt)
    V_opt = voltages_from_E(E_opt, k)
    dV_opt = E_opt / float(k)
    C_opt = (w @ E_opt) / float(k)
    snr_opt = snr_from_E(E_opt, A, B, w, k, eps_noise)

    info = dict(
        A=A, B=B, C_opt=C_opt, snr_opt=float(snr_opt),
        t_gap=t_gap, mid_range=(mid0, mid1), S=S, u=u,
        F_opt=F_opt, E_opt=E_opt, center_mask=center_mask,
        eps_field=eps_field, Vmax=Vmax
    )
    return V_opt, dV_opt, info

# --------- Minimal plotting (optional) ---------
def plot_results(V_opt, dV_opt, info, tau):
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=False)

    # voltages
    ax = axes[0]
    ax.bar(np.arange(1, len(V_opt)+1), V_opt, width=0.85)
    ax.set_ylabel("Voltage [V]")
    ax.set_title("Optimized Ring Voltages")
    ax.grid(alpha=0.35)

    # E(t) via sech superposition (for a quick look)
    t_gap = info["t_gap"]
    t = np.linspace(t_gap[0] - 4*tau, t_gap[-1] + 4*tau, 2000)
    diff = (t[:, None] - t_gap[None, :]) / tau
    # E(t) in V/cm: sum dV_i * sech * k
    E_t = (dV_opt[None, :] * (1.0/np.cosh(diff)) * K_E_PER_VSTEP).sum(axis=1)

    ax = axes[1]
    ax.plot(t*1e6, E_t, lw=1.6)
    ax.set_xlabel("Time [µs]")
    ax.set_ylabel("E field [V/cm]")
    ax.set_title("Field profile (from optimized dV)")
    ax.grid(alpha=0.35)

    fig.suptitle(f"|C_opt|={np.abs(info['C_opt']):.4g}, SNR_opt={info['snr_opt']:.4g}", y=0.98)
    fig.tight_layout()
    plt.show()

# --------- Run a simple case ---------
if __name__ == "__main__":
    # Geometry & physics
    N_rings = 30
    total_length = 183.2e-3
    v = 600.0
    d_gap = total_length / N_rings
    Delta = 1500.0 * 2*np.pi
    W = 2*np.pi*5
    d_dip = -33.6*2*np.pi

    # Constraints / optimizer
    Vmax = 5.0         # ±5 V
    center_count = 8   # central active gaps
    eps_noise = 0.2
    budget = 4000
    seed = 42
    optimizer_name = "NGOpt"
    tau = 7.6e-3 / v

    V_opt, dV_opt, info = optimize_f_param(
        N_rings=N_rings, total_length=total_length, v=v,
        Delta=Delta, W=W, d_dip=d_dip, tau=tau,
        center_count=center_count, Vmax=Vmax,
        eps_field=0.02, eps_noise=eps_noise,
        budget=budget, seed=seed,
        optimizer_name=optimizer_name,
        center_F_bound_factor=4.0,
    )

    print(f"|C_opt| = {np.abs(info['C_opt']):.6g}")
    print(f"SNR_opt = {info['snr_opt']:.6g}")
    plot_results(V_opt, dV_opt, info, tau=tau)