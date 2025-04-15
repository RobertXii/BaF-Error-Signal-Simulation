%% Two-Level System with Decay: OBE Simulation and Dynamic Phase

% Parameters
d13    = -2.15e4 * 2*pi;
E_L2   = 8.514e2 * 0.15;
t0L2   = 0.0;
sigma  = 0.52e-6;
Delta2 = 2*pi*1e6*2;
Gamma  = 2.7e6 * 2*pi;
tspan  = [-1e-5 1e-5];
y0     = [1; 0; 0; 0];

% Define and solve the ODE
blochODE = @(t,y) [
    -1i * (0.5*d13*E_L2*exp(-((t-t0L2)^2)/(2*sigma^2)))*(y(3)-y(2));
    -1i * ((0.5*d13*E_L2*exp(-((t-t0L2)^2)/(2*sigma^2)))*(y(4)-y(1)) - Delta2*y(2)) - (Gamma/2)*y(2);
    -1i * ((0.5*d13*E_L2*exp(-((t-t0L2)^2)/(2*sigma^2)))*(y(1)-y(4)) + Delta2*y(3)) - (Gamma/2)*y(3);
    -1i * (0.5*d13*E_L2*exp(-((t-t0L2)^2)/(2*sigma^2)))*(y(2)-y(3)) - Gamma*y(4)
];
options = odeset('RelTol',1e-13,'AbsTol',1e-12);
[t, y] = ode45(blochODE, tspan, y0, options);
rho11 = real(y(:,1));
rho33 = real(y(:,4));

% Compute dynamic phase (numerical)
z = Delta2 - 1i*(Gamma/2);
t_phase = linspace(tspan(1), tspan(2), 10000);
Omega_13 = @(t) 0.5*d13*E_L2*exp(-((t-t0L2).^2)/(2*sigma^2));
lambda_plus = 0.5 * (z - sqrt(z.^2 + 4*(Omega_13(t_phase)).^2));
phase_integral = cumtrapz(t_phase, lambda_plus);
psi_dyn = exp(-1i*phase_integral);
pop_dyn = abs(psi_dyn).^2;

% Compute analytic dynamic phase (perturbative expansion up to second order)
Phi_dyn_analytic = -((d13*E_L2)^2/(4*(Delta2-1i*(Gamma/2)))) .* ((sqrt(pi)/2)*sigma .* (erf(t_phase/sigma)-erf(tspan(1)/sigma))) ...
    + ((d13*E_L2)^4/(16*(Delta2-1i*(Gamma/2))^3)) .* ((sqrt(pi)/(2*sqrt(2)))*sigma .* (erf(sqrt(2)*t_phase/sigma)-erf(sqrt(2)*tspan(1)/sigma)));
pop_analytic = abs(exp(-1i*Phi_dyn_analytic)).^2;

% Compute dynamic phase expansion (up to n = 3)
phi_dyn = 0;
for n = 1:5
    binom = gamma(1/2+1) / (gamma(n+1)*gamma(1/2-n+1));
    term = -binom * ((d13*E_L2)^(2*n) * sigma) / (2*z^(2*n-1)) * sqrt(pi/n);
    phi_dyn = phi_dyn + term;
end
pop_pred = abs(exp(-1i*phi_dyn))^2;

% Plot results
figure; hold on;
plot(t, rho11, 'b', 'LineWidth', 2);
plot(t, rho33, 'r', 'LineWidth', 2);
plot(t_phase, pop_dyn, 'm', 'LineWidth', 2);
yline(pop_pred, 'k--', 'LineWidth', 2);
% plot(t_phase, pop_analytic, 'k', 'LineWidth', 2);

xlabel('Time (s)'); ylabel('Population');
grid on; legend('\rho_{11}', '\rho_{33}', 'Dynamic Phase (num)', 'Expansion (n=1:3)');