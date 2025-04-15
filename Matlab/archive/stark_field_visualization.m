E0_stark = 0.1; % Amplitude of the field
omega = 2 * pi * 11.4e3; % Frequency (rad/s)
t = linspace(-50e-6, 50e-6, 1000); % Time vector

E_stark = (t > -43.7e-6 & t < 43.7e-6) .* (E0_stark * sin(omega * (t + 43.7e-6)));

plot(t * 1e6, E_stark, 'LineWidth', 2); % Time in microseconds
xlabel('Time (\mus)');
ylabel('Electric Field (arb. units)');
title('Time-Limited Sinusoidal Electric Field');
grid on;