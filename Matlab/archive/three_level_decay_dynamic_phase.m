%% Three-Level System: Non-Hermitian Hamiltonian & Dynamic Phase Calculation
clear; clc;

%% Parameters
% Coupling for 1-2 transition (via the oscillatory field)
d12   = -33.6 * 2*pi;               % dipole moment for 1-2 (arbitrary units)
E0    = 40;                       % field amplitude for E_{nL2}
omega = 2*pi*11.44e3;             % angular frequency for the sine field

% Coupling for 1-3 transition (Gaussian pulse)
d13    = -2.15e4 * 2*pi;          % dipole moment for 1-3 (as in previous code)
E_L2   = 8.514e2 * 0.1;           % electric field amplitude for L2
E0_nr  = 12;                     % additional amplitude for non-resonant part
t0L2   = 52.6e-6;                % center of the Gaussian pulse (s)
t0nr   = 70.1e-6;
sigma  = 0.52e-6;                % pulse width (s)

% Detunings and decay
Delta1 = 2*pi*2e3;               % detuning for state 2 (1-2 transition)
Delta2 = 2*pi*1e6*3;             % L2 detuning
Gamma  = 2.7e6 * 2*pi;           % decay rate for state 3

% Time span (in seconds)
t_start = -51.1e-6;             % start time
t_end   = 160e-6;               % end time
tspan   = [t_start, t_end];

%% Define time-dependent couplings
% 1-2 coupling: active only when -43.7 µs < t < 43.7 µs plus a non-resonant term
E_stark = @(t) ( (t > -43.7e-6 & t < 43.7e-6) .* ( d12 .* E0 .* sin(omega*(t + 43.7e-6)) ) );
Enr = @(t) d12 * E0_nr * sech(616 * (t - t0nr) / 0.76e-2);
Omega = @(t)E_stark(t) + Enr(t);
% 1-3 coupling: Omega'_{13}(t) = 0.5*d13*E_L2*exp(-((t-t0L2)^2)/(2*sigma^2))
Omega13prime = @(t) 0.5 * d13 .* E_L2 .* exp(-((t-t0L2).^2)/(2*sigma^2));

%% Define the 3x3 non-Hermitian Hamiltonian function H(t)
H = @(t) [ 0,           Omega(t),         Omega13prime(t);
           Omega(t),    Delta1,           0;
           Omega13prime(t), 0,   Delta2 - 1i*(Gamma/2) ];

%% Solve the state-vector dynamics:
% We solve: dpsi/dt = -i*H(t)*psi, with initial state psi(t_start) = [1;0;0]
psi0 = [1; 0; 0];  % initial state

% Define the ODE for the state vector (note: state may be complex)
stateODE = @(t, psi) -1i * H(t) * psi;

% Solve using ode45
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t_sol, psi_sol] = ode45(stateODE, tspan, psi0, opts);

% Extract populations (squared modulus)
pop1 = abs(psi_sol(:,1)).^2;
pop2 = abs(psi_sol(:,2)).^2;
pop3 = abs(psi_sol(:,3)).^2;

%% Compute the Dynamic Phase of the relevant eigenstate branch
N_phase = 100000;
t_phase = linspace(tspan(1), tspan(2), N_phase);
% Preallocate arrays for the eigenvalue, right eigenvector, and left eigenvector.
lambda_branch = zeros(1, length(t_phase));   % eigenvalue to be tracked
v_branch = zeros(3, length(t_phase));          % right eigenvector (|phi_n>)
w_branch = zeros(3, length(t_phase));          % left eigenvector (|chi_n> as column)

% --- At initial time, choose the branch with maximum weight on state 1 ---
H0 = H(t_phase(1));
[V0, D0] = eig(H0);
eigs0 = diag(D0);
[~, idx0] = max(abs(V0(1,:)));
lambda_branch(1) = eigs0(idx0);
v_branch(:,1) = V0(:,idx0);
% Compute left eigenvectors from the inverse of V0 (rows of inv(V0) are left eigenvectors)
W0 = inv(V0);
w_branch(:,1) = (W0(idx0, :)).';  % transpose to a column

prev_lambda = lambda_branch(1);
disp(prev_lambda);

% --- Loop over time to track the continuous eigenstate branch ---
for k = 2:length(t_phase)
    Hk = H(t_phase(k));
    [Vk, Dk] = eig(Hk);
    eigs_k = diag(Dk);
    % Select the eigenvalue closest to the previous one (branch tracking)
    [~, idx] = min(abs(eigs_k - prev_lambda));
    lambda_branch(k) = eigs_k(idx);
    v_branch(:, k) = Vk(:, idx);
    
    % Compute left eigenvectors: rows of inv(Vk)
    Wk = inv(Vk);
    w_branch(:, k) = (Wk(idx, :)).';
    
    % (Optional) Normalize the biorthogonal pair such that <w|v> = 1
    norm_factor = w_branch(:, k)' * v_branch(:, k);
    v_branch(:, k) = v_branch(:, k) / sqrt(norm_factor);
    w_branch(:, k) = w_branch(:, k) / sqrt(norm_factor);
    
    prev_lambda = lambda_branch(k);
end

% --- Dynamic phase: integrate lambda_branch over time ---
dynamic_phase = cumtrapz(t_phase, lambda_branch);

% --- Berry phase: compute the non-Hermitian connection A(t) ---
% A(t) = <chi(t)| (d/dt)|phi(t)>
A_connection = zeros(1, length(t_phase));
dt = t_phase(2) - t_phase(1);
for k = 1:length(t_phase)-1
    dv = (v_branch(:, k+1) - v_branch(:, k)) / dt;
    A_connection(k) = w_branch(:, k)' * dv;
end
A_connection(end) = A_connection(end-1);  % assign last point

% Integrate A(t) over time to obtain the Berry phase: gamma_B = i*integral A(t) dt.
berry_phase = 1i * cumtrapz(t_phase, A_connection);

% --- Total phase factor including both dynamic and Berry phases ---
total_phase_factor = exp(1i.*berry_phase) .* exp(-1i .* dynamic_phase);
pop_predicted = abs(total_phase_factor).^2;

%% Plot the state evolution and phase prediction
figure;
plot(t_sol*1e6, pop1, 'b', 'LineWidth', 2); hold on;
plot(t_sol*1e6, pop2, 'r', 'LineWidth', 2);
plot(t_sol*1e6, pop3, 'g', 'LineWidth', 2);
% Here we plot the predicted evolution (from full phase factor) vs time.
plot(t_phase*1e6, pop_predicted, 'm', 'LineWidth', 2);
xlabel('Time (\mus)');
ylabel('Population');
legend('\psi_1','\psi_2','\psi_3','Full Phase Prediction');
title('Population Dynamics for a 3-Level System');
grid on;

%% Plot the Berry connection A(t)
figure;
subplot(2,1,1);
plot(t_phase*1e6, real(A_connection), 'LineWidth',2);
xlabel('Time (\mus)');
ylabel('Re[A(t)]');
title('Real Part of the Berry Connection');
grid on;
subplot(2,1,2);
plot(t_phase*1e6, imag(A_connection), 'LineWidth',2);
xlabel('Time (\mus)');
ylabel('Im[A(t)]');
title('Imaginary Part of the Berry Connection');
grid on;

%% Plot Eigenvalues vs Time (Real and Imaginary Parts) ---
% Preallocate array to store all three eigenvalues at each time point.
num_eigs = 3;
eigenvalues_all = zeros(num_eigs, length(t_phase));
for k = 1:length(t_phase)
    % Compute and sort the eigenvalues at time t_phase(k)
    eigenvalues_all(:, k) = sort(eig(H(t_phase(k))));
end

% Plot Real parts of eigenvalues
figure;
subplot(2,1,1);
plot(t_phase*1e6, real(eigenvalues_all(1,:)), 'r', 'LineWidth',2); hold on;
plot(t_phase*1e6, real(eigenvalues_all(2,:)), 'g', 'LineWidth',2);
plot(t_phase*1e6, real(eigenvalues_all(3,:)), 'b', 'LineWidth',2);
xlabel('Time (\mus)');
ylabel('Real Part of Eigenvalue');
legend('eig1','eig2','eig3','Location','Best');
title('Real Parts of Eigenvalues vs. Time');
grid on;

% Plot Imaginary parts of eigenvalues
subplot(2,1,2);
plot(t_phase*1e6, imag(eigenvalues_all(1,:)), 'r', 'LineWidth',2); hold on;
plot(t_phase*1e6, imag(eigenvalues_all(2,:)), 'g', 'LineWidth',2);
plot(t_phase*1e6, imag(eigenvalues_all(3,:)), 'b', 'LineWidth',2);
xlabel('Time (\mus)');
ylabel('Imaginary Part of Eigenvalue');
legend('eig1','eig2','eig3','Location','Best');
title('Imaginary Parts of Eigenvalues vs. Time');
grid on;