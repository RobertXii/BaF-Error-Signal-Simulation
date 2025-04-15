% Parameters
% W = 5 * 2 * pi;                      % Weak interaction matrix element [rad/s]
W = 0;                                 % Coherence decay parameter 
Delta = 2 * pi * 2e3;                  % Energy gap for |b> [rad/s]
Delta13 = 2 * pi * 348.69e12;          % Energy gap for |c> [rad/s]
d12 = -3360 * 2 * pi;                  % Dipole moment between ab  [Hz/(V/cm)]
d13 = 3360 * 2 * pi;                   % Dipole moment between ac  [Hz/(V/cm)] 
d23 = 0;                               % Dipole moment between bc  [Hz/(V/cm)]
E0_stark = 0.1;                        % Stark field amplitude   [V/cm]
E0_nr = 0.42e-1;                       % non-reversing field amplitude [V/cm]
E0_L2 = 1;                             % Depletion laser L2 amplitude  [V/cm]         
omega_stark = 2 * pi * 11.4e3;         % stark field frequency [rad/s]
detuning_L2 = 0;                       % L2 laser detuning [rad/s]
gamma_c = 1e3;                         % Decay rate for |c> 
v = 616;                               % molecular beam velocity [m/s]
sigma_u = 0.76e-2;                     % Width of the non-reversing field [m]
t0 = 70.1e-6;                          % Center of the non-reversing field [s]
% tspan = linspace(-43.7e-6, 43.7e-6);    % Time span for the simulation [s]
tspan = linspace(-51.1e-6, 100e-6,100000);

% Initial state
rho0 = zeros(3);
rho0(1,1) = 1; 
% rho0(2,2) = 1; 

% Define the electric field as a function of time
E_stark = @(t) (t > -43.7e-6 & t < 43.7e-6) .* (E0_stark * sin(omega_stark * (t + 43.7e-6)));
E_nr = @(t) E0_nr * sech(v * (t - t0) / sigma_u);
E_L2 = @(t) (t > 43.7e-6 & t < 61.5e-6) .* (E0_L2 * sin((Delta13 + detuning_L2) * (t - 43.7e-6)));

E_t = @(t) E_stark(t) + E_nr(t) + E_L2(t);    % total field
% E_t = @(t) E_stark(t);

% Define the time-dependent Hamiltonian
H_t = @(t) [
    0, 1i * W + d12 * E_t(t), d13 * E_t(t);
    -1i * W + d12 * E_t(t), Delta, d23 * E_t(t);
    d13 * E_t(t), d23 * E_t(t), Delta13
];

% Define the decay (Lindblad) operators
L_a = sqrt(gamma_c) * [0; 0; 1] * [0, 0, 1]; % Decay from |3>

% Solve the ODE
% options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10); % Tighten tolerances
[t, rho_sol] = ode15s(@(t, rho_vec) master_equation(t, rho_vec, H_t, L_a), tspan, rho0(:));

% Extract populations
rho_aa = real(rho_sol(:, 1));     % Population of |a>
rho_bb = real(rho_sol(:, 5));     % Population of |b>
rho_cc = real(rho_sol(:, 9));     % Population of |c>

% Plot the populations
figure;
% plot(t*1e6, rho_aa, '-k', 'LineWidth', 2); 
hold on;
plot(t*1e6, rho_bb, '-r', 'LineWidth', 2);
plot(t*1e6, rho_cc, '-b', 'LineWidth', 2);
xlabel('Time(\mus)');
ylabel('Probability');
legend('|b>', '|c>', '|c>');
title('Probability Evolution');
lgd = legend('show');
set(lgd, 'Box', 'off','Color', 'none','FontName', 'Times New Roman','FontSize', 15); 
set(gca, 'FontName', 'Times New Roman','FontSize', 15);
grid on;
hold off;

% Master equation function
function drho_dt = master_equation(t, rho_vec, H_t, L_a)
    rho = reshape(rho_vec, 3, 3); % Reshape vector to 3x3 matrix
    % Evaluate the time-dependent Hamiltonian
    H = H_t(t);
    % Commutator term
    commutator = -1i * (H * rho - rho * H);
    % Dissipator term
    L_a_dagger_L_a = L_a' * L_a; 
    dissipator = L_a * rho * L_a' - 0.5 * (L_a_dagger_L_a * rho + rho * L_a_dagger_L_a);    
    % Total derivative
    drho_dt = commutator + dissipator;
    drho_dt = drho_dt(:); % Flatten to a vector for the solver
    % disp("a")
end
