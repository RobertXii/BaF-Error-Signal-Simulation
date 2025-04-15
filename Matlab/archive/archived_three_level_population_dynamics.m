% Parameters
Delta = 2 * pi * 2e3;          % Detuning between levels 1 and 2 (Hz)
% Delta13 = 2 * pi * 348.69e12;        % Detuning between levels 1 and 3 (Hz)
Delta13 = 0;        % Detuning between levels 1 and 3 (Hz)
d12 = -3360 * 2 * pi;                      % Dipole matrix element between levels 1 and 2 (arb. units)
d13 = 0;                    % Dipole matrix element between levels 1 and 3 (arb. units)
d23 = 0;                    % Dipole matrix element between levels 2 and 3 (arb. units)
tspan = linspace(-51.1e-6, 100e-6); % Time span for the simulation (s)
rho0 = [1, 0, 0, 0, 0, 0, 0, 0, 0]; % Initial density matrix: [rho_11, rho_12, rho_13, rho_21, rho_22, rho_23, rho_31, rho_32, rho_33]

% Stark field parameters
E0_stark = 0.1;                 % Amplitude of the Stark field (arb. units)
omega = 2 * pi * 11.4e3;                      % frequency of the stark field 

% Non-resonant field parameters
% E0_nr = 0.42;                  % Amplitude of non-resonant field (arb. units)
% sigma_u = 0.76e-2;               % Width of the non-resonant field (s)
% t0 = 70.1e-6;                       % Center of the non-resonant field (s)

% Define the Stark field \mathcal{E}_\text{stark}(t)
E_stark = @(t) (t > -43.7e-6 & t < 43.7e-6) .* (E0_stark * sin(omega * (t + 43.7e-6)));

% Define the non-resonant field \mathcal{E}_{nr}(t)
% E_nr = @(t) E0_nr * sech(v * (t - t0) / sigma_u);

% Define the total electric field \mathcal{E}(t)
% E = @(t) E_stark(t) + E_nr(t);
E = @(t) E_stark(t);

% Define the system of ODEs
% dynamics = @(t, rho) reshape(...
%     1i * [
%         % Row 1
%         E(t) * d13 * (rho(3) - rho(7)) + E(t) * d12 * (rho(2) - rho(4)), ...
%         Delta * rho(2) + E(t) * (d23 * rho(3) - d13 * rho(8) + d12 * (rho(1) - rho(5))), ...
%         Delta13 * rho(3) + E(t) * (d13 * (rho(1) - rho(9)) + d23 * rho(2) - d12 * rho(6)); ...
%         % Row 2
%         -Delta * rho(4) + E(t) * (d13 * rho(6) - d23 * rho(7) + d12 * (rho(5) - rho(1))), ...
%         E(t) * d23 * (rho(6) - rho(8)) + E(t) * d12 * (rho(4) - rho(2)), ...
%         (Delta13 - Delta) * rho(6) + E(t) * (d13 * rho(4) + d23 * (rho(5) - rho(9)) - d12 * rho(3)); ...
%         % Row 3
%         -Delta13 * rho(7) + E(t) * (d13 * (rho(9) - rho(1)) - d23 * rho(4) + d12 * rho(8)), ...
%         (Delta - Delta13) * rho(8) + E(t) * (-d13 * rho(2) + d23 * (rho(9) - rho(5)) + d12 * rho(7)), ...
%         E(t) * (d13 * (rho(7) - rho(3)) + d23 * (rho(8) - rho(6)))
%     ], [], 1);
dynamics = @(t, rho) reshape(...
    1i * [
        % Row 1
        E(t) * d12 * (rho(2) - rho(4)), ...
        Delta * rho(2) + E(t) * (d12 * (rho(1) - rho(5))), ...
        0; ...
        % Row 2
        -Delta * rho(4) + E(t) * ( d12 * (rho(5) - rho(1))), ...
        E(t) * d12 * (rho(4) - rho(2)), ...
        0; ...
        % Row 3
        0, ...
        0, ...
        0
    ], [], 1);

% Solve the ODEs
opts = odeset('RelTol', 1e-18, 'AbsTol', 1e-20);
[t, rho] = ode45(dynamics, tspan, rho0, opts);

% Extract populations (diagonal elements of the density matrix)
rho_11 = rho(:, 1); % Population of level 1
rho_12_real = real(rho(:, 2)); % Population of level 1
rho_22 = rho(:, 5); % Population of level 2
rho_33 = rho(:, 9); % Population of level 3

% Plot results
figure;

% Plot populations
% plot(t * 1e6, rho_11, 'LineWidth', 2, 'DisplayName', '\rho_{11} (Level 1)');
hold on;
plot(t * 1e6, rho_22, 'LineWidth', 2, 'DisplayName', '\rho_{22} (Level 2)');
% plot(t * 1e6, rho_33, 'LineWidth', 2, 'DisplayName', '\rho_{33} (Level 3)');
xlabel('Time (\mus)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Population', 'FontName', 'Times New Roman', 'FontSize', 14);
title('Population Dynamics of a Three-Level System', 'FontName', 'Times New Roman', 'FontSize', 16);
legend('show', 'FontName', 'Times New Roman', 'FontSize', 12);
grid on;
