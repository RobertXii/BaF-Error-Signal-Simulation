% Parameters
W = 5 * 2 * pi;             % Coupling strength (Hz)
E0 = 0.1;                   % Electric field strength (V/cm)
omega = 2 * pi * 11.4e3;    % Driving frequency of the electric field (Hz)
d = -3360 * 2 * pi;         % Dipolar matrix element (Hz/(V/cm))
tspan = [0, 87.4e-6];       % Time span for the simulation (s)
rho0 = [1, 0, 0, 0];        % Initial conditions of the density matrix: [rho_11, rho_12_real, rho_12_imag, rho_22]

% Range of detuning values
detuning_range = linspace(-4000, 4000, 1000) * 2 * pi; % Convert from Hz to angular frequency (rad/s)
asymmetry = zeros(size(detuning_range)); % Initialize asymmetry array
asymmetry_TDPT = zeros(size(detuning_range)); % Initialize asymmetry TDPT array

% Define the Bloch Optical Equation as ODE
dynamics = @(t, rho, E_field, Delta) [
    % d(rho_11)/dt
    W * (rho(3) + rho(2)) - 1i * d * E_field * sin(omega * t) * (rho(3) - rho(2));
    % d(rho_12)/dt
    (W - 1i * d * E_field * sin(omega * t)) * (rho(4) - rho(1)) + 1i * Delta * rho(2);
    % d(rho_21)/dt
    (- W - 1i * d * E_field * sin(omega * t)) * (rho(1) - rho(4)) - 1i * Delta * rho(3);
    % d(rho_22)/dt
    -W * (rho(3) + rho(2)) - 1i * d * E_field * sin(omega * t) * (rho(2) - rho(3));
];

% Solver options for numerical integration
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-20);

% Loop over detuning values
for i = 1:length(detuning_range)
    Delta = detuning_range(i); % Current detuning
    
    % Compute S(+E0)
    [~, rho_plus] = ode45(@(t, rho) dynamics(t, rho, E0, Delta), tspan, rho0, opts);
    S_plus = rho_plus(end, 4); % Final value of \rho_{22} for +E0
    
    % Compute S(-E0)
    [~, rho_minus] = ode45(@(t, rho) dynamics(t, rho, -E0, Delta), tspan, rho0, opts);
    S_minus = rho_minus(end, 4); % Final value of \rho_{22} for -E0
    
    % Calculate asymmetry
    asymmetry(i) = (S_plus - S_minus) / (S_plus + S_minus);

    % Skip calculation for the middle 20% of the detuning range
    % n = length(detuning_range);
    % start_skip = floor(n * 0.5); % Start index of middle 20%
    % end_skip = ceil(n * 0.5);   % End index of middle 20%
    % 
    % if i >= start_skip && i <= end_skip
    %     asymmetry_TDPT(i) = NaN; % Yield NA for this range
    %     continue;
    % end

    % Calculate |c_+(t)|^2 based on
    % % equation (1.52) in Emine's Thesis, derived using TDPT
    % sin_term_plus = sin(Delta * tspan(2) / 2).^2; % Sinusoidal term
    % % interaction_term_plus = 2 * (W / Delta) * (d * E0 / omega) + (d * E0 / omega)^2; % Interaction term (1.52)
    % interaction_term_plus = 2 * (W / Delta) * (d * E0 /sqrt(omega^2 - Delta^2)) + (d * E0 /sqrt(omega^2 - Delta^2))^2; % Interaction term (1.99)
    % c_plus_squared_plus = 4 * sin_term_plus * interaction_term_plus; % Probability |c_+(t)|^2
    first_term = 2 * W ./ Delta .* sin(Delta * t/2) .* exp(-1i * Delta * t/2);
    prop = -1i * d * E0 .* exp(-1i * Delta * t)./ (omega.^2 - Delta.^2);

    c_plus = first_term + prop .* (omega .* exp(1i * Delta * t) - 1i * Delta.* sin(omega * t) - omega .* cos(omega * t));
    c_minus = first_term - prop .* (omega .* exp(1i * Delta * t) - 1i * Delta.* sin(omega * t) - omega .* cos(omega * t));
    c_plus_squared_plus = abs(c_plus).^2;
    c_plus_squared_minus = abs(c_minus).^2;

    %reverse the field
    % sin_term_minus = sin(Delta * tspan(2) / 2).^2; % Sinusoidal term
    % % interaction_term_minus = 2 * (W / Delta) * (d * -E0 / omega) + (d * -E0 / omega)^2; % Interaction term (1.52)
    % interaction_term_minus = 2 * (W / Delta) * (d * -E0 / sqrt(omega^2 - Delta^2)) + (d * -E0 / sqrt(omega^2 - Delta^2))^2; % Interaction term (1.99)
    % c_plus_squared_minus = 4 * sin_term_minus * interaction_term_minus; % Probability |c_+(t)|^2
    asymmetry_TDPT(i) = (c_plus_squared_plus(end) - c_plus_squared_minus(end)) / (c_plus_squared_plus(end) + c_plus_squared_minus(end));
end

% Convert detuning to Hz for the plot
detuning_Hz = detuning_range / (2 * pi);

% Plot Detuning vs. Asymmetry
figure;
plot(detuning_Hz/1000, asymmetry_TDPT, 'LineWidth', 2, 'Color', 'r');
hold on;
plot(detuning_Hz/1000, asymmetry, 'LineWidth', 2, 'Color','k');
xlabel('\Delta / 2\pi [kHz]');
ylabel('Asymmetry (A)');
title('Detuning vs. Asymmetry');
lgd = legend('TDPT','Analytic','Location', 'best'); 
set(lgd, 'Box', 'off','Color', 'none','FontName', 'Times New Roman','FontSize', 20); 
set(gca, 'FontName', 'Times New Roman','FontSize', 15);
grid on;
hold off;

