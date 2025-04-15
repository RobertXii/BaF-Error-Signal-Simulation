% Parameters
Omega1 = 2*pi*8e6; % Rabi frequency (5 MHz)
omega0 = 2*pi*6.8e9; % Transition frequency (6.8 GHz)
omega = 2*pi*6.8e9; % Resonance driving frequency
tspan = [0, 2e-7]; % Time span (2 microseconds)

% Initial conditions [sigma_bb, sigma_aa, real(sigma_ab), imag(sigma_ab)]
sigma0 = [1, 0, 0, 0];

% Define the system of ODEs
dynamics = @(t, sigma) [
    1i * Omega1 * cos(omega * t) * (sigma(4) - sigma(3)); % dsigma_bb/dt
    -1i * Omega1 * cos(omega * t) * (sigma(4) - sigma(3)); % dsigma_aa/dt
    1i * omega0 * sigma(3) - 1i * Omega1 * cos(omega * t) * (sigma(1) - sigma(2)); % Re(dsigma_ab/dt)
    -1i * omega0 * sigma(4) + 1i * Omega1 * cos(omega * t) * (sigma(1) - sigma(2)); % Im(dsigma_ab/dt)
];

% Solve the ODEs
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[t, sigma] = ode45(dynamics, tspan, sigma0,opts);

% Extract results
sigma_bb = sigma(:, 1);
sigma_aa = sigma(:, 2);
sigma_ab_real = sigma(:, 3);
sigma_ab_imag = sigma(:, 4);

% Plot the results
figure;
plot(t, sigma(:, 1), '-', 'DisplayName', '\sigma_{bb}'); % sigma_bb
hold on;
plot(t, sigma(:, 2), '-', 'DisplayName', '\sigma_{aa}'); % sigma_aa
hold off;
xlabel('Time');
ylabel('Population');
title('Time Evolution of \sigma_{bb} and \sigma_{aa}');
legend('show'); 

% plot(t, sigma_ab_real, 'LineWidth', 2);
% xlabel('Time');
% ylabel('Re(\sigma_{aa})');
% title('\sigma_{aa}');
