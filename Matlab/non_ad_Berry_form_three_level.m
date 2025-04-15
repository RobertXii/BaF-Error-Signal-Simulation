% clear; clc;

%% Parameters
% Coupling for 1-2 transition (via the oscillatory field)
d12   = -33.6 * 2*pi;         % dipole moment for 1-2 (arbitrary units)
E0    = 40;                   % field amplitude for stark
omega = 2*pi*11.44e3;         % angular frequency for the sine field

% Coupling for 1-3 transition (Gaussian pulse)
d13    = -2.15e4 * 2*pi;      % dipole moment for 1-3 (as in previous code)
E_L2   = 8.514e2*0.02;            % electric field amplitude for L2
E0_nr  = 12;                 % additional amplitude for non-resonant part
t0L2   = 52.6e-6;            % center of the Gaussian pulse (s)
t0nr   = 70.1e-6;
sigma  = 0.52e-6;            % pulse width (s)

% Detunings and decay
Delta1 = 2*pi*2e3;          % detuning for state 2 (1-2 transition)
Delta2 = 2*pi*1e6*3;        % L2 detuning
Gamma  = 2.7e6 * 2*pi;      % decay rate for state 3

% Time span (in seconds)
t_start = -51.1e-6;        % start time
t_end   = 160e-6;          % end time
tspan   = [t_start, t_end];

%% Define time-dependent couplings
% 1-2 coupling: active only when -43.7 µs < t < 43.7 µs plus a non-resonant term
E_stark = @(t) ( (t > -43.7e-6 & t < 43.7e-6) .* ( d12 .* E0 .* sin(omega*(t + 43.7e-6)) ) );
Enr = @(t) d12 * E0_nr * sech(616 * (t - t0nr) / 0.76e-2);
Omega = @(t) E_stark(t) + Enr(t);
% 1-3 coupling: Omega'_{13}(t) = 0.5*d13*E_L2*exp(-((t-t0L2)^2)/(2*sigma^2))
Omega13prime = @(t) 0.5 * d13 .* E_L2 .* exp(-((t-t0L2).^2)/(2*sigma^2));

%% Define the 3x3 non-Hermitian Hamiltonian H(t)
H = @(t) [ 0,               Omega(t),         Omega13prime(t);
           Omega(t),        Delta1,           0;
           Omega13prime(t), 0,   Delta2 - 1i*(Gamma/2) ];

%% Solve the Right-State Evolution: dψ/dt = -i H(t) ψ,  ψ(t_start) = [1; 0; 0]
psi0 = [1; 0; 0];  % initial state
stateODE_right = @(t, psi) -1i * H(t) * psi;
opts = odeset('RelTol',1e-14,'AbsTol',1e-15);
[t_sol, psi_sol] = ode45(stateODE_right, tspan, psi0, opts);

% (For visualization) Extract populations from full simulation
pop1 = abs(psi_sol(:,1)).^2;
pop2 = abs(psi_sol(:,2)).^2;
pop3 = abs(psi_sol(:,3)).^2;

%% Solve the Left-State Evolution (Adjoint Equation):
% dχ/dt = + i H†(t) χ, choose chi0 so that <chi(0)|psi(0)> = 1, e.g. chi0 = [1; 0; 0]
chi0 = [1; 0; 0];
stateODE_left = @(t, chi) +1i * (H(t))' * chi;   % H(t)' is H†(t)
[t_sol_left, chi_sol] = ode45(stateODE_left, tspan, chi0, opts);

% Assume t_sol and t_sol_left are essentially the same grid.
N_time = length(t_sol);

%% Biorthonormalization: at each time t, compute:
%   N(t) = <chi(t)|psi(t)>
%   v(t) = psi(t)/sqrt(N(t))
%   w(t) = chi(t)/sqrt(N(t))
v_states = zeros(N_time, 3);   % normalized right states (each row is 1x3)
w_states = zeros(N_time, 3);   % normalized left states
N_factors = zeros(N_time, 1);

for k = 1:N_time
    psi_k = psi_sol(k, :).';  % 3x1 column vector
    chi_k = chi_sol(k, :).';  % 3x1 column vector
    N_k = chi_k' * psi_k;     % scalar (in general complex)
    N_factors(k) = N_k;
    
    % Use the square root of N_k. (Choose branch consistently.)
    v_states(k, :) = (psi_k / sqrt(N_k)).';
    w_states(k, :) = (chi_k / sqrt(N_k)).';
end

%% Compute Dynamic and Berry Phases
% Dynamic phase: gamma_D(t) = ∫ [w(t)'*H(t)*v(t)] dt
% Berry phase:   gamma_B(t) = i ∫ [w(t)'*(dv/dt)] dt
dynamic_phase = zeros(N_time, 1);
berry_phase   = zeros(N_time, 1);

for k = 2:N_time
    dt = t_sol(k) - t_sol(k-1);
    
    % -- Dynamic phase integrand --
    v_k = v_states(k, :).';
    w_k = w_states(k, :).';
    Hk  = H(t_sol(k));
    integrand_dyn = w_k' * Hk * v_k;  % may be complex in non-Hermitian case
    dynamic_phase(k) = dynamic_phase(k-1) + integrand_dyn * dt;
    
    % -- Berry phase integrand --
    v_prev = v_states(k-1, :).';
    % Finite difference for dv/dt:
    dv = (v_k - v_prev) / dt;
    integrand_berry = 1i * (w_k' * dv);
    berry_phase(k) = berry_phase(k-1) + integrand_berry * dt;
end

%% Construct the Predicted Wavefunction from the Phases
% psi_pred(t) = exp(i * berry_phase(t)) * exp(-i * dynamic_phase(t)) * v(t)
psi_pred = zeros(N_time, 3);
for k = 1:N_time
    phaseFactor = exp(1i * berry_phase(k)) * exp(-1i * dynamic_phase(k));
    v_k = v_states(k, :).';
    psi_pred(k, :) = (phaseFactor * v_k).';
end

%% (Optional) Compare populations of the direct integration and predicted state
pop_pred = abs(psi_pred).^2;

figure;
plot(t_sol*1e6, pop1, 'b', 'LineWidth',2); hold on;
plot(t_sol*1e6, pop2, 'r', 'LineWidth',2);
plot(t_sol*1e6, pop3, 'g', 'LineWidth',2);
plot(t_sol*1e6, pop_pred(:,1), 'r--', 'LineWidth',2);
plot(t_sol*1e6, pop_pred(:,2), 'g--', 'LineWidth',2);
plot(t_sol*1e6, pop_pred(:,3), 'b--', 'LineWidth',2);
% plot(t_sol*1e6, Omega(t_sol), 'k', 'LineWidth',2);
xlabel('Time (\mus)');
ylabel('Population');
legend('State 1 (full)', 'State 2 (full)', 'State 3 (full)', 'State 1 (pred)', 'State 2 (pred)', 'State 3 (pred)');
title('Three-Level Non-Hermitian System: Full vs. Berry-Phase Reconstruction');
grid on;


