%% Parameters and Initialization
% Constants
W = 0 * 5 * 2 * pi;                             % Weak interaction matrix element [rad/s]
Delta = 2 * pi * 2e3;                           % Energy gap for |b> [rad/s]
% Delta13 = 2 * pi * 348.69e12;                 % Energy gap for |c> [rad/s]
d12 = -33.6 * 2 * pi;                           % Dipole moment between ab  [rad/(s·V/m)]
d13 = 2.15e4 * 2 * pi;                          % Dipole moment between ac  [rad/(s·V/m)]
gamma_c = 2 * pi * 2.7e6;                       % Decay rate for |c| in [rad/s]

% Field parameters
E0_stark = 40;                                  % Stark field amplitude [V/m]
E0_nr = -12;                                     % Non-reversing field amplitude [V/m]
E0_L2 = 8.514e2 * 0.11;                         % Depletion laser L2 amplitude [V/m]
omega_stark = 2 * pi * 11.4e3;                  % Stark field frequency [rad/s]
detuning_L2 = 2 * pi * 1e6 * 3;                 % L2 laser detuning [rad/s]
v = 616;                                        % Molecular beam velocity [m/s]
sigma_u = 0.76e-2;                              % Width constant of the non-reversing field [m]
t0 = 70.1e-6;                                   % Center of the non-reversing field [s]
tspan  = linspace(-51.1e-6, 100e-6, 100000);    % Time span for the simulation [s]
% tspan = linspace(-43.7e-6, 43.7e-6,10000);    
T_e = 87.4e-6;                                  % duration for stark interference [s]
T_f1 = 7.4e-6;                                  % free evolution time after L1[s]
T_f2 = 8.9e-6;                                  % free evolution time before L2[s]

t_center_L2 = 52.6e-6;                          % Center of the laser pulse [s]
sigma_L2 = 0.52e-6;                             % variance of the laser [s]

%%% Define Variables
variable = linspace(-4e6, 4e6, 80) * 2 * pi; % variable
% detuning_range = linspace(-4000, 4000, 80) * 2 * pi; % ground state Detuning values in rad/s
% detuning_L2_range = linspace(-4e6, 4e6, 200) * 2 * pi; % Detuning values in rad/s
% E0_nr_range = linspace(-12, 12, 23); % [V/m]
% E0_L2_range = linspace(0, 3, 8) * 8.514e2; % [V/m]

% Initialize asymmetry array
asymmetry = zeros(size(variable));
W_val = zeros(size(variable));

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
parfor i = 1:length(variable)
    detuning_L2 = variable(i);       
    % E0_nr = E0_nr_range(i); 
    % Delta = variable(i);
    % E0_L2 = variable(i);

    % omega_L2 = Delta13 + detuning_L2;           % Laser L2 frequency
    Rabi_half = E0_L2 * d13 / 2;                % half Rabi frequency between a, c

    % Define electric fields as functions of time
    E_stark = @(t) (t > -43.7e-6 & t < 43.7e-6) .* (E0_stark * sin(omega_stark * (t + 43.7e-6))); % Stark Interference field
    E_nr = @(t) E0_nr * sech(v * (t - t0) / sigma_u);                                             % non-reversing field
    % Rabi = @(t) (t > 51.85e-6 & t < 53.35e-6) * Rabi_half;
    Rabi = @(t) Rabi_half * exp(-((t - t_center_L2).^2) / (2 * sigma_L2^2));
    
    E_t_nL2 = @(t) E_stark(t) + E_nr(t);                        % Total field without L2 laser
    E_t_nL2_reverse = @(t) -E_stark(t) + E_nr(t);              % Total field without L2 laser reveresed stark

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
    asymmetry(i) = (S_plus - S_minus) / (S_plus + S_minus); 

    rho_solutions{i} = rho_sol; % Store the time evolution for each detuning
end

%% Visualization of the Result

% Convert detuning to Hz for plotting
% detuning_Hz = detuning_range / (2 * pi);

% Plot the asymmetry vs detuning
figure;

%%% L2 Detuning 
plot(variable / (1e6 * 2 * pi) , asymmetry, 'k', 'LineWidth', 1.5);
xlabel('\delta_{L2} / 2\pi [MHz]');
ylabel('Asymmetry');
title('Asymmetry vs. \delta_{L2} \Delta/2\pi = 1kHz; E_{nr}=12V/m');

%%% energy Detuning
% plot(variable / (1e3 * 2 * pi) , asymmetry, 'k', 'LineWidth', 1.5);
% xlabel('\Delta / 2\pi [kHz]'  );
% ylabel('Asymmetry');

%%% E0_nr
% scatter(E0_nr_range*10, asymmetry, 50, 'k', 'filled');
% xlabel('E_{nr0} [mV/cm]');
% ylabel('Asymmetry');
% title('E_{nr0} vs. asymmetry');

%%% E_L2
% plot(variable/8.514e2, asymmetry, 'k', 'LineWidth', 1.5);
% xlabel({'$ \times 8.514 \times 10^2\quad\mathcal{E}_{L2}$[V/m]'}, 'Interpreter', 'latex');
% ylabel('Asymmetry');
% title('Asymmetry vs. E_{L2} E_{nr0}=120mV/cm; E_{0,stark} =400mV/cm; \delta_{L2}=3MHz;\Delta=2\times2\pi kHz');

set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
grid on;

rho_sol_mid = rho_solutions{1};
time = tspan;
pop_a = rho_sol_mid(:, 1); % Population of |a> (rho_aa)
pop_b = rho_sol_mid(:, 5); % Population of |b> (rho_bb)
pop_c = rho_sol_mid(:, 9); % Population of |c> (rho_cc)

% Plot population time evolution
figure;
% set(gcf, 'Position', [100, 100, 250, 200]); % Set figure size (Width x Height)
plot(time/1e-6, pop_a, 'g', 'LineWidth', 1.5); % Red for |a>
hold on;
plot(time/1e-6, pop_b, 'r', 'LineWidth', 1.5); % Green for |b>
plot(time/1e-6, pop_c, 'b', 'LineWidth', 1.5); % Blue for |c>

% Add labels and legend
% xlabel('Time [\mus]');
% ylabel('Population');
% title('Time Evolution of Energy State Populations');
legend('|a>', '|b>', '|c>', 'Location', 'west');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);

% %% diagnose wiggle behavior by plotting the final b state population
% % Extract the fifth element from the 100000th row for each cell
% finalValues = zeros(1, 200);
% 
% for i = 1:200
%     finalValues(i) = rho_solutions{i}(70000, 1);  % Access row 5, column 100000/ 70000 after L2
% end
% 
% % Plot the extracted values
% figure;
% plot(variable / (1e6 * 2 * pi), finalValues, '-', 'LineWidth', 2, Color='k');
% xlabel('\delta_{L2} / 2\pi [MHz]');
% ylabel('Population');
% title('Odd Parity Ground State Population after L2 at 20% E_{L2} intensity');
% grid on;
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);

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


