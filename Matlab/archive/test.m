% Parameters
W = 5 * 2 * pi;                      % Weak interaction matrix element [rad/s]
% W = 0;                                 
Delta = 2 * pi * 0.2e3;                  % Energy gap for |b> [rad/s]-
Delta13 = 2 * pi * 348.69e12;          % Energy gap for |c> [rad/s]-
d12 = -33.6 * 2 * pi;                  % Dipole moment between ab  [rad/(s·V/m)]-
d13 = 0*100 * 2 * pi;                    % Dipole moment between ac  [rad/(s·V/m)] 
d23 = 0;                               % Dipole moment between bc  [rad/(s·V/m)]
gamma_c = 2 * pi * 2.7e6;              % Decay rate for |c| in [rad/s], translate into 59ns decay time constant-

% Electric field amplitudes (in V/m)
E0_stark = 10;                        % Stark field amplitude-
E0_nr = 42*0;                       % non-reversing field amplitude-
E0_L2 = 8.514e3;                         % Depletion laser L2 amplitude      
omega_stark = 2 * pi * 11.4e3;         % stark field frequency [rad/s]-
detuning_L2 = 0 * pi * 2e6;            % L2 laser detuning [rad/s]-
v = 616;                               % molecular beam velocity [m/s]-
sigma_u = 0.76e-2;                     % Width constant of the non-reversing field [m]-
t0 = 70.1e-6;                          % Center of the non-reversing field [s]-
% tspan = linspace(-43.7e-6, 43.7e-6);    % Time span for the simulation [s]
tspan = linspace(-51.1e-6, 100e-6,100000);

omega_L2 = Delta13 + detuning_L2;      % Laser L2 frequency
diff = omega_L2 - Delta;

% Initial state
rho0 = zeros(3);
rho0(1,1) = 1; 
% rho0(2,2) = 1; 

% Define the electric field as a function of time
E_stark = @(t) (t > -43.7e-6 & t < 43.7e-6) .* (E0_stark * sin(omega_stark * (t + 43.7e-6)));
E_nr = @(t) E0_nr * sech(v * (t - t0) / sigma_u);
E_L2 = @(t) (t > 52.35e-6 & t < 52.85e-6) .* (E0_L2 * sin(omega_L2 * (t - 43.7e-6)));
E_t = @(t) E_stark(t) + E_nr(t) + E_L2(t);    % total field

% With 2,2 unitary transformation as exp(i * Delta * t)
H_prime = @(t) [
    0, (1i * W + d12 * E_t(t)) * exp(-1i * Delta * t), d13 * E_t(t) * exp(-1i * omega_L2 * t);
    (-1i * W + d12 * E_t(t)) * exp(1i * Delta * t), 0, d23 * E_t(t) * exp(-1i * diff * t);
    d13 * E_t(t) * exp(1i * omega_L2 * t), d23 * E_t(t) * exp(1i * diff * t), -detuning_L2
];


E_t_reversed = @(t) -E_stark(t) + E_nr(t) + E_L2(t);    % total reversed field

H_prime_reversed = @(t) [
    0, (1i * W + d12 * E_t_reversed(t)) * exp(-1i * Delta * t), d13 * E_t_reversed(t) * exp(-1i * omega_L2 * t);
    (-1i * W + d12 * E_t_reversed(t)) * exp(1i * Delta * t), 0, d23 * E_t_reversed(t) * exp(-1i * diff * t);
    d13 * E_t_reversed(t) * exp(1i * omega_L2 * t), d23 * E_t_reversed(t) * exp(1i * diff * t), -detuning_L2
];

% Solve the ODE
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10); % Tighten tolerances
[t, rho_sol] = ode45(@(t, rho_vec) master_equation(t, rho_vec, H_prime, gamma_c), tspan, rho0(:));
[t_reversed, rho_sol_reversed] = ode45(@(t, rho_vec) master_equation(t, rho_vec, H_prime_reversed, gamma_c), tspan, rho0(:));

% Extract the final density matrix
format long;
rho_final_vec = rho_sol(end, :); % The last row of rho_sol
rho_final = reshape(rho_final_vec, 3, 3); % Reshape to 3x3 matrix
% disp(rho_final);
rho_final_vec_reversed = rho_sol_reversed(end, :); % The last row of rho_sol
rho_final_reversed = reshape(rho_final_vec_reversed, 3, 3); % Reshape to 3x3 matrix
% disp(rho_final_reversed);


% Extract final values of rho_22
rho_final = reshape(rho_sol(end, :), 3, 3); % Final density matrix for forward field
rho_final_reversed = reshape(rho_sol_reversed(end, :), 3, 3); % Final density matrix for reversed field

S_plus = real(rho_final(2,2)); % Population of |b> for forward field
S_minus = real(rho_final_reversed(2,2)); % Population of |b> for reversed field

% Calculate the asymmetry
asymmetry = (S_plus - S_minus) / (S_plus + S_minus);

% Display the asymmetry
fprintf('S+: %.6f, S-: %.6f\n', S_plus, S_minus);
disp('Asymmetry:');
disp(asymmetry);

disp('W:');
W_extracted = d12 * Delta * E0_stark * asymmetry / (2 * omega_stark);
disp(W_extracted/(2*pi));


% Extract populations
rho_aa = real(rho_sol(:, 1));     % Population of |a>
rho_bb = real(rho_sol(:, 5));     % Population of |b>
rho_cc = real(rho_sol(:, 9));     % Population of |c>

% Plot the populations
figure;
plot(t*1e6, rho_aa, '-g', 'LineWidth', 2); 
hold on;
plot(t*1e6, rho_bb, '-r', 'LineWidth', 2);
plot(t*1e6, rho_cc, '-b', 'LineWidth', 2);
plot(t*1e6, rho_aa + rho_bb + rho_cc, '-k', 'LineWidth', 2);

xlabel('Time(\mus)');
ylabel('Probability');
legend('|a>', '|b>', '|c>');
title('Probability Evolution');
lgd = legend('show');
set(lgd, 'Box', 'off','Color', 'none','FontName', 'Times New Roman','FontSize', 15); 
set(gca, 'FontName', 'Times New Roman','FontSize', 15);
grid on;
hold off;



% Master equation function
function drho_dt = master_equation(t, rho_vec, H_prime, gamma_c)
    rho = reshape(rho_vec, 3, 3); % Reshape vector to 3x3 matrix
    % Evaluate the time-dependent Hamiltonian
    H = H_prime(t);
    % Commutator term
    commutator = -1i * (H * rho - rho * H);
    % Dissipator term
    dissipator = gamma_c / 2  * [0, 0, rho(1,3); 0, 0, rho(2,3); rho(3,1), rho(3,2), 2 * rho(3,3)];    
    % Total derivative
    drho_dt = commutator - dissipator;
    drho_dt = drho_dt(:); % Flatten to a vector for the solver
    % disp("a")
end
