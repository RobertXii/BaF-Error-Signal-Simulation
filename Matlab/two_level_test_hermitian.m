clear; clc;

%% Parameters
d12   = -33.6*2*pi;        % dipole moment (arb. units)
E0    = 40;               % field amplitude
omega = 2*pi*11.44e3;      % angular frequency (rad/s)

Delta = 2*pi*2e3;         % real part of detuning
gamma_decay = 2.7e2 * 2*pi * 10;   % decay rate (choose as needed)
Delta_complex = Delta - 1i*gamma_decay;  % non-Hermitian detuning

t_start = -51.1e-6;       % start time (s)
t_end   = 160e-6;         % end time (s)
tspan   = [t_start, t_end];

%% Define time-dependent coupling: active only when -43.7 µs < t < 43.7 µs.
E_stark = @(t) ((t > -43.7e-6 & t < 43.7e-6) .* ( d12*E0 * sin(omega*(t + 43.7e-6)) ));

%% Non-Hermitian Hamiltonian for Two-Level System
H = @(t) [ 0,            E_stark(t);
           E_stark(t),   Delta_complex ];

%% Solve the Right-State Evolution
psi0 = [1; 0];   % initial state: all in state 1
stateODE_right = @(t, psi) -1i * H(t) * psi;
opts = odeset('RelTol',1e-13,'AbsTol',1e-12);
[t_sol, psi_sol] = ode45(stateODE_right, tspan, psi0, opts);

%% Solve the Left-State Evolution (Adjoint Equation)
% The left state |chi> satisfies: d|chi>/dt = + i H^\dagger(t)|chi>
% (we choose chi0 so that <chi(0)|psi(0)> = 1; here we take chi0 = [1; 0])
chi0 = [1; 0];
stateODE_left = @(t, chi) +1i * (H(t))' * chi;  % note: (H)' gives H^\dagger for complex matrices
[t_sol_left, chi_sol] = ode45(stateODE_left, tspan, chi0, opts);

% For simplicity assume t_sol and t_sol_left are (nearly) the same time grid.
N = length(t_sol);

%% Biorthonormalize: Compute normalized right (v) and left (w) states
% Here, N_k = <chi(t)|psi(t)>; then define:
%     v(t) = psi(t)/sqrt(N_k)
%     w(t) = chi(t)/sqrt(N_k)
v_states = zeros(N, 2);   % right normalized states (each row is a 1x2 vector)
w_states = zeros(N, 2);   % left normalized states
N_factors = zeros(N,1);

for k = 1:N
    psi_k = psi_sol(k, :).';   % column vector (2x1)
    chi_k = chi_sol(k, :).';   % column vector (2x1)
    N_k = chi_k' * psi_k;      % biorthogonal normalization factor
    N_factors(k) = N_k;
    
    % Use sqrt with proper branch; here we assume N_k is nonzero.
    v_states(k, :) = (psi_k / sqrt(N_k)).';
    w_states(k, :) = (chi_k / sqrt(N_k)).';
end

%% Compute Dynamic and Berry Phases via Finite Differences
% The dynamic phase (gamma_D) is:
%    gamma_D(t) = ∫ dt ( w(t)' * H(t) * v(t) )
% The Berry (geometric) phase (gamma_B) is:
%    gamma_B(t) = i ∫ dt ( w(t)' * dv/dt )
dynamic_phase = zeros(N,1);
berry_phase   = zeros(N,1);

for k = 2:N
    dt = t_sol(k) - t_sol(k-1);
    
    % Dynamic phase integrand at t(k)
    v_k = v_states(k, :).';
    w_k = w_states(k, :).';
    Hk = H(t_sol(k));
    integrand_dyn = w_k' * Hk * v_k;  
    dynamic_phase(k) = dynamic_phase(k-1) + integrand_dyn * dt;
    
    % Berry phase integrand via finite difference:
    v_prev = v_states(k-1, :).';
    % Approximate derivative: dv/dt
    dv = (v_k - v_prev) / dt;
    integrand_berry = 1i * (w_k' * dv);
    berry_phase(k) = berry_phase(k-1) + integrand_berry * dt;
end

%% Construct the Predicted State Using the Non-Adiabatic Phase Factors
% The predicted state is:
%    psi_pred(t) = exp(i*berry_phase(t)) * exp(-i*dynamic_phase(t)) * v(t)
psi_pred = zeros(N, 2);
for k = 1:N
    phaseFactor = exp(1i * berry_phase(k)) * exp(-1i * dynamic_phase(k));
    v_k = v_states(k, :).';
    psi_pred(k, :) = (phaseFactor * v_k).';
end

%% (Optional) Compare Predicted Evolution with the Direct ODE Solution
pop_full = abs(psi_sol).^2;  % populations from the full (right-state) evolution
pop_pred = abs(psi_pred).^2;   % populations from the predicted state

figure;
plot(t_sol*1e6, pop_full(:,1), 'b', 'LineWidth',2); hold on;
plot(t_sol*1e6, pop_full(:,2), 'r', 'LineWidth',2);
plot(t_sol*1e6, pop_pred(:,1), 'k--', 'LineWidth',2);
plot(t_sol*1e6, pop_pred(:,2), 'k--', 'LineWidth',2);
xlabel('Time (\mus)');
ylabel('Population');
legend('State 1 (full)', 'State 2 (full)', 'State 1 (pred)', 'State 2 (pred)');
title('Two-Level Non-Hermitian System: Full vs. Berry-Phase Reconstruction');
grid on;