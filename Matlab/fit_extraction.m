clearvars; % Clears all variables from the workspace

%% Parameters and Initialization
% Constants
W = 0 * 5 * 2 * pi;                             % Weak interaction matrix element [rad/s]
Delta = 2 * pi * 2e3;                           % Energy gap for |b> [rad/s]
% Delta13 = 2 * pi * 348.69e12;                 % Energy gap for |c> [rad/s]
d12 = -33.6 * 2 * pi;                           % Dipole moment between ab  [rad/(s·V/m)]
d13 = -2.15e4 * 2 * pi;                         % Dipole moment between ac  [rad/(s·V/m)]
gamma_c = 2.7e6 * 2 * pi;                       % Decay rate for |c| in [rad/s]

% Field parameters
E0_stark = 40;                                  % Stark field amplitude [V/m]
E0_nr = -12;                                    % Non-reversing field amplitude [V/m]
E0_L2 = 8.514e2 * 0.09;                         % Depletion laser L2 amplitude [V/m]
omega_stark = 2 * pi * 11.4e3;                  % Stark field frequency [rad/s]
detuning_L2 = 2 * pi * 1e6 * 1;                 % L2 laser detuning [rad/s]
v = 616;                                        % Molecular beam velocity [m/s]
sigma_u = 0.76e-2;                              % Width constant of the non-reversing field [m]
t0 = 70.1e-6;                                   % Center of the non-reversing field [s]
tspan  = linspace(-51.1e-6, 100e-6, 100000);    % Time span for the simulation [s]
% tspan = linspace(-43.7e-6, 43.7e-6,10000);    
T_e = 87.4e-6;                                  % duration for stark interference [s]
T_f1 = 7.4e-6;                                  % free evolution time after L1[s]
T_f2 = 8.9e-6;                                  % free evolution time before L2[s]
% T_f2 = 56.3e-6;                               % free evolution time before L2[s]

t_center_L2 = 52.6e-6;                            % Center of the laser pulse [s]
sigma_L2 = 0.52e-6;                             % variance of the laser [s]

% Define Variables
detuning_range = linspace(-4000, 4000, 80) * 2 * pi; % Detuning values in rad/s
% detuning_L2_range = linspace(-4e6, 4e6, 16) * 2 * pi; % Detuning values in rad/s
% E0_nr_range = linspace(-12, 12, 23); % [V/m]

% Initialize asymmetry array
asymmetry = zeros(size(detuning_range));
W_val = zeros(size(detuning_range));

%% Initial Condition
% Initial state
rho0 = zeros(3);
rho0(1,1) = 1; % Start in state |a>
rho0(2,2) = 0; 

%% Start the parallel computing
% Enable parallel pool
if isempty(gcp('nocreate'))
    parpool; % Start a parallel pool if one is not running
end

% Parallel loop for detuning range
parfor i = 1:length(detuning_range)
    % detuning_L2 = detuning_L2_range(i);       
    % E0_nr = E0_nr_range(i); 
    Delta = detuning_range(i);

    % omega_L2 = Delta13 + detuning_L2;           % Laser L2 frequency
    Rabi_half = E0_L2 * d13 / 2;                % half Rabi frequency between a, c

    % Define electric fields as functions of time
    E_stark = @(t) (t > -43.7e-6 & t < 43.7e-6) .* (E0_stark * sin(omega_stark * (t + 43.7e-6))); % Stark Interference field
    E_nr = @(t) E0_nr * sech(v * (t - t0) / sigma_u);                                             % non-reversing field
    % E_L2 = @(t) (t > 52.35e-6 & t < 52.85e-6) .* (E0_L2 * cos(omega_L2 * (t - 52.35e-6)));      % L2 Laser
    % Rabi = @(t) (t > 51.85e-6 & t < 53.35e-6) * Rabi_half;
    Rabi = @(t) Rabi_half * exp(-((t - t_center_L2).^2) / (2 * sigma_L2^2));
    
    E_t_nL2 = @(t) E_stark(t) + E_nr(t);                        % Total field without L2 laser
    E_t_nL2_reverse = @(t) - E_stark(t) + E_nr(t);              % Total field without L2 laser reveresed stark

    % Define Hamiltonians for forward and reversed fields
    H_prime = @(t) [
        0,  1i * W + d12 * E_t_nL2(t), Rabi(t);
         -1i * W + d12 * E_t_nL2(t), Delta, 0;
        Rabi(t), 0, -detuning_L2
    ];

    H_prime_reversed = @(t) [
        0,  1i * W + d12 * E_t_nL2_reverse(t), Rabi(t);
         -1i * W + d12 * E_t_nL2_reverse(t), Delta, 0;
        Rabi(t), 0, -detuning_L2
    ];

    % Solve the ODE for forward and reversed fields
    options = odeset('RelTol', 1e-13, 'AbsTol', 1e-12); % Tighten tolerances
    [~, rho_sol] = ode45(@(t, rho_vec) master_equation(t, rho_vec, H_prime, gamma_c), tspan, rho0(:), options);
    [~, rho_sol_reversed] = ode45(@(t, rho_vec) master_equation(t, rho_vec, H_prime_reversed, gamma_c), tspan, rho0(:),options);
    
    % Extract the final density matrices
    rho_final = reshape(rho_sol(end, :), 3, 3); % Final density matrix for forward field
    rho_final_reversed = reshape(rho_sol_reversed(end, :), 3, 3); % Final density matrix for reversed field
    
    % Calculate populations
    S_plus = real(rho_final(2,2)); % Population of |b> for forward field
    S_minus = real(rho_final_reversed(2,2)); % Population of |b> for reversed field
    
    % Calculate the asymmetry
    asymmetry(i) = (S_plus - S_minus) / (S_plus + S_minus); % The minus term is to cancel the a_0 const
    % W_val(i) = d12 * Delta * E0_stark * asymmetry(i) / (2 * omega_stark * 2 * pi);
    % format long;
    % disp(omega_L2);
    % A_Delta = 2 * d12 * E0_stark * omega_stark ./ (omega_stark^2 - Delta.^2); % Driving amplitude
    % sin_term = sin(Delta * T_e / 2); % Sine term
    % denominator = sin(Delta * (T_f2+T_f1+T_e)/2) .* cos(Delta * (T_f1-T_f2)/2); % Denominator term
    % denominator(denominator == 0) = eps;
    % W_val(i) = (asymmetry(i) + 0.0344 * E0_nr) .* (A_Delta .* (Delta .* sin_term) ./ denominator);

    rho_solutions{i} = rho_sol; % Store the time evolution for each detuning
end

%% Visualization of the Result

% Convert detuning to Hz for plotting
detuning_Hz = detuning_range / (2 * pi);

% Plot the asymmetry vs detuning
figure;
% scatter(E0_nr_range*10, asymmetry, 50, 'k', 'filled');
plot(detuning_Hz / 1000, asymmetry, 'k', 'LineWidth', 1.5);
xlabel('\Delta / 2\pi [kHz]'  );
% xlabel('\delta_{L2} / 2\pi [MHz]');
% xlabel('E_{nr0} [mV/cm]');
% ylabel('W/2\pi (Hz)');
ylabel('Asymmetry')
% title('E_{nr0} vs. W');
% title('\Delta vs. Asymmetry');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
grid on;

rho_sol_mid = rho_solutions{1};

% display density matrix before and after the pulse
% index_before_pulse = round(100000 * 101.1 / 151.1);
% row_vector = rho_sol_mid(index_before_pulse, :);
% rho_matrix = reshape(row_vector, [3, 3]);
% disp('rho before pulse:');
% disp(rho_matrix);
% 
% rho_sol_mid = rho_solutions{1};
% index_after_pulse = round(100000 * 106.1 / 151.1);
% row_vector = rho_sol_mid(index_after_pulse, :);
% rho_matrix = reshape(row_vector, [3, 3]);
% disp('rho after pulse:');
% disp(rho_matrix);
% 
% disp(rho_sol_mid(end, 5))


time = tspan;
pop_a = rho_sol_mid(:, 1); % Population of |a> (rho_aa)
pop_b = rho_sol_mid(:, 5); % Population of |b> (rho_bb)
pop_c = rho_sol_mid(:, 9); % Population of |c> (rho_cc)

% Plot population time evolution
figure;
plot(time/1e6, pop_a, 'g', 'LineWidth', 1.5); % Red for |a>
hold on;
plot(time/1e6, pop_b, 'r', 'LineWidth', 1.5); % Green for |b>
plot(time/1e6, pop_c, 'b', 'LineWidth', 1.5); % Blue for |c>

% Add labels and legend
xlabel('Time [\mus]');
ylabel('Population');
title('Time Evolution of Energy State Populations');
legend('|a>', '|b>', '|c>', 'Location', 'Best');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);

%% fit function
% Exclude data in the range -1 kHz to 1 kHz
exclude_indices = abs(detuning_Hz) > 1000; % Logical array for excluding the range
detuning_fit = detuning_Hz(exclude_indices) * 2 * pi; % Exclude detuning in [-1 kHz, 1 kHz]
asymmetry_fit = asymmetry(exclude_indices); % Exclude corresponding asymmetry data

% Define the fitting function A(Delta)
fit_function = @(params, Delta) ...
    (2 * params(1) ./ Delta) .* ...
    ((omega_stark^2 - Delta.^2) ./ (d12 * E0_stark * omega_stark)) .* ...
    (sin((Delta / 2) .* (T_e + T_f1 + T_f2)) ./ sin((Delta / 2) * T_e)) .* ...
    cos((Delta / 2) * (T_f1 - T_f2)) + ...
    params(2) + params(3) * Delta;

% Initial guess for parameters: [W, a_0, a_1]
initial_params = [1, 0.01, 1e-8];

% Define the error function for fitting
error_function = @(params) ...
    sum((asymmetry_fit - fit_function(params, detuning_fit)).^2);

% Use optimization to find best-fit parameters
options = optimset('Display', 'off', 'TolFun', 1e-12, 'TolX', 1e-12);
best_fit_params = fminsearch(error_function, initial_params, options);

% Extract the parameters
W = best_fit_params(1);
a_0 = best_fit_params(2);
a_1 = best_fit_params(3);

% Display the results
fprintf('Fitted parameters (excluding -1 kHz to 1 kHz)[W, a_0, a_1]:\n');
fprintf('%.5e\n', W/(2*pi));
fprintf('%.5e\n', a_0);
fprintf('%.5e\n', a_1);

% Plot the fitted curve against the data
Delta_fit = linspace(min(detuning_fit), max(detuning_fit), 1000);
asymmetry_fit_curve = fit_function(best_fit_params, Delta_fit);
exclude_fit_indices = abs(Delta_fit / (1000 * 2 * pi)) > 1; % Logical mask for excluding the range
Delta_fit_below = Delta_fit(Delta_fit / (1000 * 2 * pi) < -1);
asymmetry_fit_curve_below = asymmetry_fit_curve(Delta_fit / (1000 * 2 * pi) < -1);
Delta_fit_above = Delta_fit(Delta_fit / (1000 * 2 * pi) > 1);
asymmetry_fit_curve_above = asymmetry_fit_curve(Delta_fit / (1000 * 2 * pi) > 1);

figure;
plot(detuning_Hz/1000, asymmetry, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'All Data');
hold on;
plot(detuning_fit /(1000 * 2 * pi), asymmetry_fit, 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Data Used for Fit');
plot(Delta_fit_below / (1000 * 2 * pi), asymmetry_fit_curve_below, 'r-', 'LineWidth', 2, 'DisplayName', 'Fit'); % Fit below -1 kHz
plot(Delta_fit_above / (1000 * 2 * pi), asymmetry_fit_curve_above, 'r-', 'LineWidth', 2, 'HandleVisibility', 'off'); % Fit above 1 kHz
xlabel('\Delta / 2\pi [kHz]');
ylabel('Asymmetry');
title('Asymmetry vs Detuning E_{nr}=12 V/m;E_{0}=40 V/m; \delta_{L2}=3MHz');
legend('show');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);

%% Define the Master equation
% Master equation function
function drho_dt = master_equation(t, rho_vec, H_prime, gamma_c)
    rho = reshape(rho_vec, 3, 3); % Reshape vector to 3x3 matrix
    % Evaluate the time-dependent Hamiltonian
    H = H_prime(t);
    % Commutator term
    commutator = - 1i * (H * rho - rho * H);
    % Dissipator term
    dissipator =  - gamma_c / 2  * [0, 0, rho(1,3); 0, 0, rho(2,3); rho(3,1), rho(3,2), 2 * rho(3,3)];  
    % Total derivative
    drho_dt = commutator + dissipator; % No dissipator since gamma_c = 0
    drho_dt = drho_dt(:); % Flatten to a vector for the solver
end


