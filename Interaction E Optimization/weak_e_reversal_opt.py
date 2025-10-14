import numpy as np
from typing import Callable, Tuple, Sequence

"""
Weak-force E-reversal optimization utilities
============================================
Core model (user-specified physics):
    c_plus(C) = A + B*C
where A and B are complex constants determined by your configuration.

In this cleaned version, A and B are computed **only** from the five physical inputs:
    W (weak matrix element), Delta (detuning), T (interaction time), d (dipole matrix element), tau (sech width).

Definitions used:
    A(W,Δ,T)   = (2 W / Δ) * sin(Δ T / 2) * exp(-i Δ T / 2)   [Δ→0 limit: A ≈ W T]
    B(Δ,τ,d)   = d * (π τ) / cosh(π Δ τ / 2)

Key functions:
  - c_plus(C, A, B)
  - S(C, A, B) : |c_plus|^2
  - asymmetry(C, A, B)
  - snr_objective(C, A, B, epsilon)
  - AB_from_params(W, Delta, T, d, tau) -> (A,B)
  - c_plus_from_params(C, W, Delta, T, d, tau)
  - objective_from_params(C, W, Delta, T, d, tau, epsilon)
  - sample_C_boundary_from_rings(...) -> boundary of feasible C
  - optimize_C_from_rings(A, B, ...)
"""

# ---------- Core algebra ----------

def c_plus(C: complex, A: complex, B: complex) -> complex:
    """Return c_plus(C) = A + B*C."""
    return A + B*C

def S(C: complex, A: complex, B: complex) -> float:
    """Return S(C) = |c_plus(C)|^2."""
    val = c_plus(C, A, B)
    return float(np.abs(val)**2)

def delta_S(C: complex, A: complex, B: complex) -> float:
    """Return ΔS = S(C) - S(-C) = 4 * Re(A * conj(B) * conj(C))."""
    return float(4.0 * np.real(A * np.conj(B) * np.conj(C)))

def sigma_S(C: complex, A: complex, B: complex) -> float:
    """Return ΣS = S(C) + S(-C) = 2*(|A|^2 + |B|^2*|C|^2)."""
    return float(2.0 * (np.abs(A)**2 + (np.abs(B)**2) * (np.abs(C)**2)))

def asymmetry(C: complex, A: complex, B: complex, eps: float = 0.0) -> float:
    """
    Return asymmetry = ΔS / (ΣS + eps).
    By default eps=0; you can include a small additive stabilizer if desired.
    """
    sig = sigma_S(C, A, B) + float(eps)
    if sig == 0.0:
        return 0.0
    return delta_S(C, A, B) / sig

def snr_objective(C: complex, A: complex, B: complex, epsilon: float = 0.0) -> float:
    """
    SNR-like objective (statistically principled for Poisson + additive tech noise):
        |ΔS| / sqrt(ΣS + epsilon^2)
    epsilon is an additive (count-like) technical noise floor (per measurement).
    """
    num = abs(delta_S(C, A, B))
    den = np.sqrt(sigma_S(C, A, B) + float(epsilon)**2)
    # den = abs(sigma_S(C, A, B))
    if den == 0.0:
        return 0.0
    return float(num / den)

# ---------- Optimization over a boundary ----------

def optimize_over_boundary(
    boundary: Sequence[complex],
    objective_fn: Callable[[complex], float]
) -> Tuple[complex, float]:
    """
    Given a sequence of complex boundary points (polygon/curve), evaluate objective_fn(C)
    and return the best point and its objective value.
    """
    best_C = None
    best_val = -np.inf
    for C in boundary:
        val = float(objective_fn(C))
        if val > best_val:
            best_val = val
            best_C = C
    return best_C, best_val

# ---------- Geometry of feasible C from ring voltages ----------

def _extreme_C_in_direction(theta: float, w: np.ndarray, Vmax: float) -> complex:
    """
    For a given direction angle theta (unit phasor u=e^{i theta}),
    compute the extreme point of C = sum E_i w_i with E_i in [-Vmax, Vmax]
    by bang-bang choice:
       E_i = Vmax * sign( Re( e^{-i theta} * w_i ) )
    This realizes the support function / extreme point of the zonotope.
    """
    u_star = np.exp(-1j * theta)
    proj = np.real(u_star * w)   # projection of each w_i onto direction theta
    E = Vmax * np.sign(proj)     # bang-bang
    return complex(np.sum(E * w))

def sample_C_boundary_from_rings(
    N: int,
    T_f1: float,
    Vmax: float,
    d: float,
    v: float,
    Delta: float,
    num_angles: int = 720
) -> np.ndarray:
    """
    Approximate the feasible boundary of C in the complex plane for an array of N rings,
    each with voltage E_i in [-Vmax, Vmax], spaced by distance d along a beam with velocity v.

    Model:
        t_i = (i-1) * d / v
        w_i = exp(-i * Delta * t_i)
        C = sum_i E_i * w_i ,  E_i in [-Vmax, Vmax]

    The feasible set is a zonotope. We sample its boundary by sweeping direction theta
    and taking the extreme point in that direction (bang-bang construction).

    Returns
    -------
    boundary : np.ndarray of shape (num_angles,)
        Complex points approximating the boundary in CCW order.
    """
    t = (np.arange(N)) * (d / v) + T_f1
    w = np.exp(-1j * Delta * t)  # complex weights
    thetas = np.linspace(0.0, 2.0 * np.pi, num_angles, endpoint=False)
    boundary = np.array([_extreme_C_in_direction(th, w, Vmax) for th in thetas], dtype=np.complex128)
    return boundary

# ---------- Convenience wrapper for optimizing over ring-derived boundary ----------

def optimize_C_from_rings(
    T_f1: float,
    A: complex,
    B: complex,
    N: int,
    Vmax: float,
    d: float,
    v: float,
    Delta: float,
    epsilon: float = 0.0,
    num_angles: int = 720,
    objective: str = "snr"
) -> Tuple[complex, float, np.ndarray]:
    """
    Build the C-feasible boundary from ring parameters and optimize the chosen objective.

    objective: "snr" or "asym"
    Returns: (C_opt, obj_val, boundary_points)
    """
    boundary = sample_C_boundary_from_rings(N, T_f1, Vmax, d, v, Delta, num_angles=num_angles)

    if objective == "snr":
        def obj(C):
            return snr_objective(C, A, B, epsilon)
    elif objective == "asym":
        def obj(C):
            return abs(asymmetry(C, A, B))  # maximize magnitude of asymmetry
    else:
        raise ValueError("objective must be 'snr' or 'asym'")

    C_opt, val = optimize_over_boundary(boundary, obj)
    return C_opt, val, boundary

# ---------- Physics helpers: build A and B from (W, Delta, T, d, tau) ----------

def AB_from_params(W: float, Delta: float, T: float, d: float, tau: float) -> Tuple[complex, complex]:
    """
    Build (A,B) from physical parameters using standard perturbative formulas.

    Parameters
    ----------
    W : float
        Weak matrix element magnitude (use sign via complex phase if needed).
    Delta : float
        Detuning (angular frequency units, rad/s).
    T : float
        Total interaction time window for the weak-only evolution (s).
    d : float
        Dipole matrix element (units such that d*E is a Rabi frequency).
    tau : float
        sech pulse width parameter for each gap (s).

    Returns
    -------
    A, B : complex
    """
    # A term: weak-only off-resonant evolution over duration T
    # A = (2 W / Δ) * sin(Δ T / 2) * exp(-i Δ T / 2)
    if Delta == 0.0:
        # Use the Δ→0 limit: A ≈ W T (since sin(x)/x → 1 and the exponential → 1)
        A = W * T
    else:
        A = (2.0 * W / Delta) * np.sin(0.5 * Delta * T) * np.exp(-1j * 0.5 * Delta * T)

    # B term: envelope from sech(t/τ) Fourier transform (real & positive)
    # B = -i * d * π τ / cosh(π Δ τ / 2)
    B = -1j * d * np.pi * tau / np.cosh(0.5 * np.pi * Delta * tau) 
    return A, B

def c_plus_from_params(C: complex, W: float, Delta: float, T: float, d: float, tau: float) -> complex:
    """Convenience: compute c_plus given C and (W,Δ,T,d,τ) using AB_from_params."""
    A, B = AB_from_params(W, Delta, T, d, tau)
    return c_plus(C, A, B)

def objective_from_params(C: complex, W: float, Delta: float, T: float, d: float, tau: float, epsilon: float = 0.0) -> float:
    """SNR-like objective using AB_from_params."""
    A, B = AB_from_params(W, Delta, T, d, tau)
    return snr_objective(C, A, B, epsilon=epsilon)

# ---------- Example usage ----------
if __name__ == "__main__":
    # Physical parameters (EDIT with your experiment)
    W = 2 * np.pi * 5         # rad/s (weak amplitude)
    T = 87.4 * 1e-6        # s (interaction time)
    d = -33.6 * 2 * np.pi        # dipole matrix element (arb. so that d*E is Rabi freq)
    v = 600           # m/s (beam speed, for ring spacing)
    gap = 7.6 * 1e-3    # m sech width sigma in meters
    tau = gap/v         # s (sech width in seconds)

    import matplotlib.pyplot as plt

    # Loop over Delta values and compute asymmetry for each
    Delta_values = np.linspace(-4000, 4000, 200) * 2 * np.pi
    asym_values = []
    for Delta in Delta_values:
        A, B = AB_from_params(W, Delta, T, d, tau)
        N = 31
        voltages = np.array([
                -0.0,   # 1
                0.05,    # 2
                -0.1,    # 3
                0.2,    # 4
                -0.25,    # 5
                0.4,    # 6
                -0.5,    # 7
                0.70,    # 8
                -0.7,    # 9
                0.5,    # 10
                -0.2,    # 11
                0.75,    # 12
                1,    # 13
                -1.75,    # 14
                -1.25,    # 15
                -2.8,    # 16
                -2.8,    # 17
                -1.25,    # 18
                -1.75,    # 19
                1,    # 20
                0.75,    # 21
                -0.2,    # 22
                0.5,    # 23
                -0.7,    # 24
                0.70,    # 25
                -0.50,    # 26
                0.4,    # 27
                -0.25,    # 28
                0.20,    # 29
                -0.10,    # 30
                0.05     # 31
                ])
        dV = np.diff(voltages, prepend=0.0)
        t = np.arange(N) * 183.2e-3 / 30 / v
        C = np.sum(dV * np.exp(-1j * Delta * t))
        asym_values.append(asymmetry(C, A, B))

    # Plot asymmetry vs Delta
    plt.figure(figsize=(8, 5))
    plt.plot(Delta_values / (2*np.pi), asym_values, 'b-')
    plt.xlabel("Delta [Hz]")
    plt.ylabel("Asymmetry")
    plt.title("Asymmetry vs Detuning")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
