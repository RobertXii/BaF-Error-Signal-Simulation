%% Parameters
% Stark field (NR) parameters
d12    = -33.6 * 2 * pi;               % dipole moment for Stark transition
E0     = 0;                         % Stark field amplitude, 𝓔_0
Delta  = 2*pi*2e3;                  % Stark field detuning
T_e    = 87.4e-6;                    % Stark field duration
omega  = 2*pi*11.44e3;                % Stark field angular frequency

% NR field parameters
E_nr0  = -10;                        % NR field amplitude, 𝓔_{nr0}
t0     = 70.1e-6;                  % NR field center time
t1     = 52.6e-6;                   % L2 laser applied time (t_diff = t1-t0 < 0)
t_diff = t1 - t0;                  

% L2 (dynamic phase) parameters
d13    = -2.15e4 * 2 * pi;           % dipole moment for L2 field
E_L2   = 8.514e2 * 0.0;              % L2 field amplitude
sigma  = 0.52e-6;                   % pulse width
Delta2 = 2*pi*1e6*-3;                % detuning for non-Hermitian Hamiltonian
Gamma  = 2.7e6 * 2*pi;              % decay rate
z = Delta2 - 1i*(Gamma/2);

% NR field envelope parameter (from the Gaussian form)
v = 616;                           % velocity parameter
sigma_u = 0.76e-2;                 % width parameter for NR field envelope
a = v^2/(2*sqrt(2)*sigma_u^2);

% Order of dynamic phase expansion
N = 5;

%% Compute the Stark field (first term)
term1 = 2 * d12 * omega * E0 * exp(-1i * Delta*T_e/2) / (omega^2 - Delta^2) * ...
    sin(Delta*T_e/2);

%% Compute the dynamic phase expansion (L2 field) up to N=3
phi_dyn = 0;
for n = 1:N
    binom = gamma(1/2+1) / ( gamma(n+1) * gamma(1/2 - n + 1) );
    term = -binom * ((d13 * E_L2)^(2*n) * sigma) / (2 * z^(2*n-1)) * sqrt(pi/n);
    phi_dyn = phi_dyn + term;
end
exp_factor = exp(1i * phi_dyn);

%% Compute the NR field integration factor with error function (updated)
NR_factor = (exp(-1i*Delta*t0)/2) * sqrt(pi/a) * exp(-1i*Delta^2/(4*a)) * ...
            (1 - erf_complex((2*a*t_diff+1i*Delta)/(2*sqrt(a))));

%% Compute the NR field contribution (second term) for c_+^(1,(0+2),N)(t2)
% The extra term (in square brackets) is the second-order correction.
term2 = -1i * ( 1 + (d12^2 * E0^2/(omega^2-Delta^2)) * ( (1i*Delta*T_e/2) - (omega^2*(1 - exp(1i*Delta*T_e))/(omega^2-Delta^2) ) ) ) ...
         * exp_factor * d12 * E_nr0 * NR_factor;

%% Total c_+^(1,(0+2),N)(t2)
c_plus = term1 + term2;
disp('c_+^(1,(0+2),N)(t2) = ');
disp(c_plus);

prob = abs(c_plus)^2;
disp('Probability = ');
disp(prob);

%% Compute second order correction c_-^(2)(T_e) ---
c_minus2 = 1 + (d12^2 * E0^2/(omega^2 - Delta^2)) * ( (1i*Delta*T_e/2) - (omega^2*(1 - exp(1i*Delta*T_e))/(omega^2 - Delta^2) ) );
disp('c_-^(2)(T_e) = ');
disp(c_minus2);

prob_odd = abs(c_minus2)^2;
disp('odd Probability = ');
disp(prob_odd);

%% Custom function: Erf for complex arguments using series expansion
function y = erf_complex(z)
    % Series expansion for erf(z):
    %   erf(z) = (2/sqrt(pi)) * sum_{k=0}^∞ [ (-1)^k * z^(2k+1) / (k!(2k+1)) ]
    maxIter = 1000;
    tol = 1e-12;
    y = zeros(size(z));
    term = (2/sqrt(pi)) * z;  % first term (k=0)
    y = y + term;
    k = 1;
    while k <= maxIter
        term_new = (2/sqrt(pi)) * ((-1)^k * z.^(2*k+1)) / (factorial(k)*(2*k+1));
        y = y + term_new;
        if all(abs(term_new(:)) < tol)
            break;
        end
        k = k + 1;
    end
end

