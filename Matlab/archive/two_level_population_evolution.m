% Parameters
W = 0 * 5 * 2 * pi;             % Coupling strength (Hz)
E0 = 40;                   % Electric field strength (V/cm)
omega = 2 * pi * 11.4e3;    % Driving frequency of the electric field (Hz)
d = -33.60 * 2 * pi;         % Dipolar matrix element (Hz/(V/cm))
Delta = 2000 * 2 * pi;      % Detuning/energy difference between two levels (Hz)
tspan = [0, 87.4e-6];       % Time span for the simulation (s)
rho0 = [1, 0, 0, 0];        % Initial conditions of the density matrix: [rho_11, rho_12, rho_21, rho_22]

% Bloch Optical Equation as ODE
dynamics = @(t, rho) [
    % d(rho_11)/dt
    W * (rho(3) + rho(2)) - 1i * d * E0 * sin(omega * t) * (rho(3) - rho(2));
    % d(rho_12)/dt
    (W - 1i * d * E0 * sin(omega * t)) * (rho(4) - rho(1)) + 1i * Delta * rho(2);
    % d(rho_21)/dt
    (- W - 1i * d * E0 * sin(omega * t)) * (rho(1) - rho(4)) - 1i * Delta * rho(3);
    % d(rho_22)/dt
    -W * (rho(3) + rho(2)) - 1i * d * E0 * sin(omega * t) * (rho(2) - rho(3));
];

% Solver options for numerical integration
opts = odeset('RelTol', 1e-18, 'AbsTol', 1e-20);

% Solve the ODEs
[t, rho] = ode45(dynamics, tspan, rho0, opts);

% Extract populations and coherences
rho_11 = rho(:, 1);          % Population of state 1
rho_22 = rho(:, 4);          % Population of state 2
rho_12_real = real(rho(:, 2)); % Real part of coherence
rho_12_imag = imag(rho(:, 2)); % Imaginary part of coherence

% Calculate |c_+(t)|^2 based on
% equation (1.48) in Emine's Thesis, derived using TDPT
first_term = 2 * W ./ Delta .* sin(Delta * t/2) .* exp(-1i * Delta * t/2);
prop = -1i * d * E0 .* exp(-1i * Delta * t)./ (omega.^2 - Delta.^2);
c_plus = first_term + prop .* (omega .* exp(1i * Delta * t) - 1i * Delta.* sin(omega * t) - omega .* cos(omega * t));
c_plus_squared = abs(c_plus).^2;

% Plot results
figure;

% Plot |c_+(t)|^2
plot(t/1e-6, c_plus_squared, '--', 'LineWidth', 3, 'DisplayName','TDPT','Color','r');
hold on;

% Plot population of state 2 (\rho_{22})
plot(t/1e-6, rho_22, 'LineWidth', 2, 'DisplayName', 'Analytic','Color','k');
xlabel('Time (\mus)');
ylabel('Probability');
title('Time Evolution of Even Parity Ground State');
lgd = legend('show');
set(lgd, 'Box', 'off','Color', 'none','FontName', 'Times New Roman','FontSize', 20); 
set(gca, 'FontName', 'Times New Roman','FontSize', 15);
grid on;

% Calculation from Equation (A)
% c_+^(1)(T_e) = (2*omega*d*E0/(omega^2-Delta))*exp(-1i*Delta*t/2).*sin(Delta*t/2)
c_plus1 = (2 * omega * d * E0 .* exp(-1i * Delta * t/2) .* sin(Delta * t/2)./ (omega^2 - Delta^2));
c_plus1_squared = abs(c_plus1).^2;

% Plot |c_+^(1)(T_e)|^2 on the same figure (using a dashed blue line)
plot(t/1e-6, c_plus1_squared, '--', 'LineWidth', 2, 'DisplayName', 'TDPT (Eq. A)', 'Color', 'b');

% Plot coherences
% figure;
% plot(t, rho_12_real, 'LineWidth', 2, 'DisplayName', 'Re(\rho_{12})'); hold on;
% plot(t, rho_12_imag, 'LineWidth', 2, 'DisplayName', 'Im(\rho_{12})');
% xlabel('Time (s)');
% ylabel('Coherence');
% title('Coherence (Re(\rho_{12}) and Im(\rho_{12}))');
% legend('show');
% grid on;

