classdef QuantumSimulation
    % QuantumSimulation Simulates a three-level quantum system.
    %
    %   This class simulates a three-level quantum system using two different
    %   formalisms:
    %     - The Optical Bloch Equation (OBE) formalism.
    %     - The non-Hermitian Schrödinger equation approach.
    %
    %   The simulation includes:
    %       - Definition of simulation constants (interaction parameters, field
    %         parameters, decay rates, etc.)
    %       - Time evolution using both the OBE and Schrödinger approaches.
    %       - A detuning scan to study asymmetry in the system.
    %       - Fitting routines for analyzing asymmetry data.
    %       - Plotting routines for visualizing the time evolution and asymmetry.
    %
    %   Usage example:
    %       sim = QuantumSimulation();
    %       sim = sim.runVary_Detuning('OBE');
    %       sim.plotAsymmetry();
    %       fitParams = sim.fitAsymmetry(true);
    %       sim.plotPopulationEvolution(3, 'OBE', 'yes');
    
    %% Properties
    properties
        %% Simulation Constants
        % Interaction and decay parameters
        W = 0 * 5 * 2*pi;            % Weak interaction matrix element [rad/s]
        % Delta = 2*pi*2e3;            % Energy gap for |b> [rad/s]
        d12 = -33.6 * 2*pi;          % Dipole moment between |a> and |b> [rad/(s·V/m)]
        d13 = -2.15e4 * 2*pi;        % Dipole moment between |a> and |c> [rad/(s·V/m)]
        gamma_c = 2.7e6 * 2*pi;      % Decay rate for |c> 2.7e6 * 2*pi [rad/s]
        
        % Field Parameters
        E0_stark = 10;               % Stark field amplitude [V/m]
        E0_nr    = 0;                % Non-reversing field amplitude [V/m]
        E0_L2    = 1204.1 * 0.1;    % 8.514e2 * 0.11;   % Depletion laser (L2) amplitude 8.514e2 [V/m]
        omega_stark = 2*pi*11.44e3;  % Stark field frequency [rad/s]
        detuning_L2 = 2*pi*1e6*2;    % L2 laser detuning [rad/s]
        v = 616;                     % Molecular beam velocity 616[m/s]
        sigma_u = 0.76e-2;           % Width constant of the non-reversing field 0.76e-2[m]
        t0 = 70.1e-6;                % Center time of the non-reversing field 70.1e-6 [s]
          
        % Time Parameters
        timesteps = 10000;          % Number of time steps for integration
        T_e = 87.4e-6;               % Duration for Stark interference [s]
        T_f1 = 7.4e-6;               % Free evolution time after exiting stark field [s]
        T_f2 = 8.9e-6;               % Free evolution time before enter stark field [s]
        t_center_L2 = 52.6e-6;       % Center of the laser pulse 52.6e-6[s]
        sigma_L2 = 0.373e-6;          % Variance of the laser pulse 0.52e-6 [s]
        tspan = linspace(-51.1e-6, 250e-6, 10000);  % Simulation time span (-51.1e-6, 400e-6, 10000)[s]
        
        % Detuning Scan Parameters
        detuning_range = linspace(-4000, 4000, 8) * 2 * pi; % Detuning values for |b> [rad/s]
        
        %% Simulation Results Storage
        initialState      % Vectorized initial density matrix (state |a>)
        asymmetry         % Asymmetry values vs. detuning
        analytic_asymmetry
        rhoSolutions      % Cell array storing time evolution for each detuning
    end
    
    %% Methods
    methods
        %% Constructor
        function obj = QuantumSimulation()
            % QuantumSimulation Construct an instance of this class.
            % Initializes the density matrix for state |a> and vectorizes it.
            
            rho0 = zeros(3);
            rho0(1,1) = 1;  % Initial state |a>
            obj.initialState = rho0(:);  % Vectorize the density matrix
        end
        
        %% OBE Formalism: Master Equation
        function drho_dt = masterEquation(obj, t, rho_vec, H_prime)
            % masterEquation Compute the time derivative of the density matrix.
            %
            %   drho_dt = masterEquation(obj, t, rho_vec, H_prime) computes the 
            %   derivative of the density matrix (vectorized) using the Optical 
            %   Bloch Equation formalism.
            %
            %   Inputs:
            %       t       - Current time [s]
            %       rho_vec - Vectorized density matrix (9x1 vector)
            %       H_prime - Function handle that returns the Hamiltonian matrix at time t
            %
            %   Outputs:
            %       drho_dt - Vectorized derivative of the density matrix
            
            % Reshape vectorized density matrix into 3x3 matrix
            rho = reshape(rho_vec, 3, 3);
            % Get time-dependent Hamiltonian
            H = H_prime(t);
            % Compute commutator term: -i[H, rho]
            commutator = -1i * (H * rho - rho * H);
            % Compute dissipator term modeling decay from |c>
            dissipator = -obj.gamma_c/2 * [ 0,         0,         rho(1,3);
                                             0,         0,         rho(2,3);
                                             rho(3,1),  rho(3,2),  2*rho(3,3) ];
            % Sum contributions and vectorize the result
            drho_dt = commutator + dissipator;
            drho_dt = drho_dt(:);
        end
        
        %% Asymmetry Calculation (OBE Approach)
        function [asym, sol] = computeOBEAsymmetry(obj, H_prime, H_prime_reversed)
            % computeOBEAsymmetry Compute asymmetry between field configurations using OBE.
            %
            %   [asym, sol] = computeOBEAsymmetry(obj, H_prime, H_prime_reversed)
            %   calculates the asymmetry between the forward and reversed field 
            %   configurations using the Optical Bloch Equation approach.
            %
            %   Inputs:
            %       H_prime          - Function handle for the forward configuration Hamiltonian.
            %       H_prime_reversed - Function handle for the reversed configuration Hamiltonian.
            %
            %   Outputs:
            %       asym - Asymmetry computed from the populations of state |b>
            %       sol  - The full simulation solution (density matrix evolution)
            
            options = odeset('RelTol',1e-13, 'AbsTol',1e-12);
            [~, rho_sol] = ode45(@(t, rho_vec) obj.masterEquation(t, rho_vec, H_prime), ...
                                 obj.tspan, obj.initialState, options);
            [~, rho_sol_rev] = ode45(@(t, rho_vec) obj.masterEquation(t, rho_vec, H_prime_reversed), ...
                                     obj.tspan, obj.initialState, options);
            % Reshape the final density matrices for forward and reversed cases
            rho_final     = reshape(rho_sol(end, :), 3, 3);
            rho_final_rev = reshape(rho_sol_rev(end, :), 3, 3);
            % Extract populations from state |b> (element (2,2))
            S_plus  = real(rho_final(1,1));
            S_minus = real(rho_final_rev(1,1));
            asym = (S_plus - S_minus) / (S_plus + S_minus);
            
            % Return the full simulation result for further analysis
            sol = rho_sol;
        end

        %% Non-Hermitian Schrödinger Equation Evolution
        function psi_sol = schrodingerEvolution(obj, H_prime)
            % schrodingerEvolution Evolve the state vector via the Schrödinger equation.
            %
            %   psi_sol = schrodingerEvolution(obj, H_prime) computes the time
            %   evolution of the state vector using a non-Hermitian effective
            %   Hamiltonian that incorporates decay.
            %
            %   Inputs:
            %       H_prime - Function handle returning the Hamiltonian at time t.
            %
            %   Outputs:
            %       psi_sol - Matrix of state vectors over time (each row corresponds
            %                 to the state at a given time)
            
            psi0 = [0; 1; 0];      % Initial state |a>
            proj3 = [0, 0, 0; 0, 0, 0; 0, 0, 1];
            % Define the effective Hamiltonian including decay of |c>
            H_eff = @(t) H_prime(t) - 1i*(obj.gamma_c/2)*proj3;
            % Define the ODE function for the Schrödinger evolution
            odeFun = @(t, psi) -1i * H_eff(t) * psi;
            options = odeset('RelTol',1e-13, 'AbsTol',1e-12);
            % Solve the time-dependent Schrödinger equation
            [~, psi_sol] = ode45(odeFun, obj.tspan, psi0, options);
        end
        
        %% Asymmetry Calculation (Schrödinger Approach)
        function [asym, sol] = computeSchrodingerAsymmetry(obj, H_prime, H_prime_reversed)
            % computeSchrodingerAsymmetry Compute asymmetry using the Schrödinger approach.
            %
            %   [asym, sol] = computeSchrodingerAsymmetry(obj, H_prime, H_prime_reversed)
            %   calculates the asymmetry between forward and reversed field 
            %   configurations using the non-Hermitian Schrödinger evolution.
            %
             %   Inputs:
            %       H_prime          - Function handle for the forward configuration Hamiltonian.
            %       H_prime_reversed - Function handle for the reversed configuration Hamiltonian.
            %
            %   Outputs:
            %       asym - Asymmetry computed from the populations of state |b>
            %       sol  - The simulation solution (state vector evolution)
            
            psi_forward  = obj.schrodingerEvolution(H_prime);
            psi_reversed = obj.schrodingerEvolution(H_prime_reversed);
            % Compute populations for state |b> (second component)
            S_plus  = abs(psi_forward(end,1))^2;
            S_minus = abs(psi_reversed(end,1))^2;
            asym = (S_plus - S_minus) / (S_plus + S_minus);
            
            % Return the forward simulation result for further analysis
            sol = psi_forward;
        end

        %% Run Detuning Scan
        function obj = runVary_Detuning(obj, method)
            warning('off','MATLAB:parfor:UninitializedTemporaries')

            % runVary_Detuning Perform a detuning scan using a chosen evolution method.
            %
            %   obj = runVary_Detuning(obj, method) runs the simulation over a
            %   range of detuning values using either the 'OBE' or 'Schrodinger'
            %   method. If no method is specified, 'OBE' is used by default.
            %
            %   Inputs:
            %       method - 'OBE' or 'Schrodinger'
            %
            %   Outputs:
            %       obj - Updated object with computed asymmetry and state evolutions
            
            if nargin < 2
                method = 'OBE';
            end
            
            nDetuning = numel(obj.detuning_range);
            asymmetry   = zeros(size(obj.detuning_range));
            rhoSolutions = cell(size(obj.detuning_range));
            
            % Start parallel pool if not already active
            if isempty(gcp('nocreate'))
                parpool;
            end
            
            parfor i = 1:nDetuning
                % Set the local detuning for state |b>
                Delta_local = obj.detuning_range(i);
                Rabi_half = obj.E0_L2 * obj.d13 / 2;
                
                % Define time-dependent field functions
                T = 43.7e-6;
                E_stark = @(t) (abs(t) <= T) .* (obj.E0_stark * sin(obj.omega_stark * (t + T)));
                E_nr    = @(t) obj.E0_nr * sech(obj.v * (t - obj.t0) / obj.sigma_u);
                Rabi    = @(t) Rabi_half * exp(-((t - obj.t_center_L2).^2) / (2*obj.sigma_L2^2));

                % E_nr field before/after L2 removed
                % E_nr = @(t) obj.E0_nr * sech(obj.v * (t - obj.t0) / obj.sigma_u) .* (0.5*(1 + tanh((t - obj.t_center_L2)/1e-6))); %before
                % E_nr = @(t) obj.E0_nr * sech(obj.v * (t - obj.t0) / obj.sigma_u) .* (0.5*(1 + tanh((obj.t_center_L2 - t)/1e-6))); %after

                % % Define the flat value for E_nr at 70 µs
                % flatVal = obj.E0_nr * sech(obj.v * (obj.t_center_L2 - obj.t0) / obj.sigma_u);
                % delta = 1e-6;  % Adjust delta to control transition smoothness
                % 
                % % Create smooth transition masks: S1 goes from 0 to 1 around 65 µs, 
                % % and S2 goes from 1 to 0 around 75 µs.
                % S1 = @(t) 0.5 * (1 + tanh((t - obj.t_center_L2 + 5e-6) / delta));
                % S2 = @(t) 0.5 * (1 + tanh((obj.t_center_L2 + 5e-6 - t) / delta));
                % mask = @(t) S1(t) .* S2(t);  % mask is ~1 between 65 and 75 µs, 0 outside, with smooth transitions
                % 
                % % Blend the original E_nr with the flat value using the mask
                % E_nr_modified = @(t) mask(t) * flatVal + (1 - mask(t)) .* E_nr(t);


                % Combined field configurations: forward and reversed
                E_t_nL2          = @(t) E_stark(t) + E_nr(t);
                E_t_nL2_reversed = @(t) -E_stark(t) + E_nr(t);
                
                % Define Hamiltonian functions for forward and reversed configurations
                H_prime_handle = @(t) [  0,                     1i*obj.W + obj.d12*E_t_nL2(t),          0;
                                         -1i*obj.W + obj.d12*E_t_nL2(t),   Delta_local,                    Rabi(t);
                                          0,                          Rabi(t),                     -obj.detuning_L2 ];
                H_prime_rev_handle = @(t) [  0,                     1i*obj.W + obj.d12*E_t_nL2_reversed(t), 0;
                                            -1i*obj.W + obj.d12*E_t_nL2_reversed(t), Delta_local,                    Rabi(t);
                                             0,                           Rabi(t),                     -obj.detuning_L2 ];

                % Evaluate evolution using the chosen method and obtain both asymmetry and simulation solution
                if strcmpi(method, 'Schrodinger')
                    [asymmetry(i), rhoSolutions{i}] = obj.computeSchrodingerAsymmetry(H_prime_handle, H_prime_rev_handle);
                elseif strcmpi(method, 'OBE')
                    [asymmetry(i), rhoSolutions{i}] = obj.computeOBEAsymmetry(H_prime_handle, H_prime_rev_handle);
                else
                    error('Unknown method. Use "Schrodinger" or "OBE".');
                end
            end
            
            % Store simulation results in object properties
            obj.asymmetry = asymmetry;
            obj.rhoSolutions = rhoSolutions;
        end
        
        %% Fit Asymmetry Data
        function fitResults = fitAsymmetry(obj, plotFlag)
            % fitAsymmetry Fit the asymmetry data to a theoretical model.
            %
            %   fitResults = fitAsymmetry(obj, plotFlag) fits the computed
            %   asymmetry data (for detunings |Delta| > 1000 Hz) to a user-defined
            %   model. Optionally, the fit curve is plotted.
            %
            %   Inputs:
            %       plotFlag - (optional) Boolean flag to display the plot. Default is false.
            %
            %   Outputs:
            %       fitResults - Best-fit parameters [W/2pi, a0, a1]
            
            if nargin < 2
                plotFlag = false;
            end
            
            % Convert detuning to Hz for the fitting procedure
            detuning_Hz = obj.detuning_range / (2*pi);
            validIdx = abs(detuning_Hz) > 1000;  % Use only data with |detuning| > 1000 Hz
            detuning_fit = detuning_Hz(validIdx) * 2*pi;
            asymmetry_fit = obj.asymmetry(validIdx);
            % asymmetry_fit = obj.analytic_asymmetry(validIdx);
            
            % Define the theoretical fit function (as from thesis)
            fitFunc = @(params, Delta) (2*params(1)./Delta) .* ...
                ((obj.omega_stark^2 - Delta.^2) / (obj.d12*obj.E0_stark*obj.omega_stark)) .* ...
                (sin((Delta/2)*(obj.T_e + obj.T_f1 + obj.T_f2)) ./ sin((Delta/2)*obj.T_e)) .* ...
                cos((Delta/2)*(obj.T_f1 - obj.T_f2)) + params(2) + params(3)*Delta;
            
            % Initial parameter guess: [W, a0, a1]
            initial_params = [1, 0.01, 1e-8];
            errFun = @(params) sum((asymmetry_fit - fitFunc(params, detuning_fit)).^2);
            options = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);
            best_fit_params = fminsearch(errFun, initial_params, options);
            
            % Display the fitted parameters
            fprintf('For E0_{nr} = %.2f V/m, L2 Detuning = %.2f MHz:\n', ...
                    obj.E0_nr, obj.detuning_L2/(2*pi*1e6));
            fprintf('  W/2pi: %.5e Hz\n', best_fit_params(1)/(2*pi));
            fprintf('  a0:    %.5e\n', best_fit_params(2));
            fprintf('  a1:    %.5e 1/Hz\n\n', best_fit_params(3));
            
            % Plot fit if requested
            if plotFlag
                Delta_fit_full = linspace(min(detuning_fit), max(detuning_fit), 1000);
                asym_fit_curve = fitFunc(best_fit_params, Delta_fit_full);
                figure;
                plot(detuning_Hz/1000, obj.asymmetry, 'ko', 'MarkerFaceColor','k', 'DisplayName','All Data');
                hold on;
                plot(detuning_fit/(1000*2*pi), asymmetry_fit, 'bo', 'MarkerFaceColor','b', 'DisplayName','Fit Data');
                plot(Delta_fit_full/(1000*2*pi), asym_fit_curve, 'r-', 'LineWidth',2, 'DisplayName','Fit Curve');
                xlabel('\Delta/2\pi [kHz]');
                ylabel('Asymmetry');
                title('Asymmetry vs. Detuning');
                legend('show');
                grid on;
                set(gca, 'FontName','Times New Roman','FontSize',15);
            end
            
            fitResults = best_fit_params;
        end
        
        %% Plot Population Evolution
        function plotPopulationEvolution(obj, idx, method, includeL2)
            % plotPopulationEvolution Plot the state populations over time.
            %
            %   plotPopulationEvolution(obj, idx, method, includeL2) plots the time
            %   evolution of the populations for states |a>, |b>, and |c> for a given
            %   simulation run.
            %
            %   Inputs:
            %       idx      - Index of the simulation result to plot.
            %       method   - 'OBE' or 'Schrodinger' (default is 'OBE').
            %       includeL2- Flag ('yes' or 'no') indicating whether to include the
            %                  L2 laser depletion effects in the additional TDPT term.
            
            if nargin < 4
                method = 'OBE';
            end
            
            time = obj.tspan;
            timeDiff = (time(end) - time(1)) * 1e6;  % Convert total time span to microseconds
            data = obj.rhoSolutions{idx};
            
            % Extract populations depending on the simulation method
            if strcmpi(method, 'Schrodinger')
                pop_a = abs(data(:,1)).^2;
                pop_b = abs(data(:,2)).^2;
                pop_c = abs(data(:,3)).^2;
            elseif strcmpi(method, 'OBE')
                pop_a = data(:,1);
                pop_b = data(:,5);
                pop_c = data(:,9);
            else
                error('Unknown method. Use "Schrodinger" or "OBE".');
            end
            
            % Plot the populations
            figure;
            plot(time*1e6, pop_a, 'g', 'LineWidth',1.5); hold on;
            plot(time*1e6, pop_b, 'r', 'LineWidth',1.5);
            plot(time*1e6, pop_c, 'b', 'LineWidth',1.5);
            
            % Retrieve additional parameters for TDPT calculation
            Delta_local = obj.detuning_range(idx);
            a = obj.v / obj.sigma_u;
            amp_a_L2 = data(round(length(pop_a)*(54.6+51.1)/timeDiff), 1);
            amp_b_L2 = data(round(length(pop_b)*(54.6+51.1)/timeDiff), 2);
            amp_c_L2 = data(round(length(pop_c)*(54.6+51.1)/timeDiff), 3);
            amp_a_end = data(end, 1);
            pop_a_end = abs(amp_a_end)^2;

            % disp(['amp. |a> at L2: ', num2str(amp_a_L2)]);
            % disp(['amp. |b> at L2: ', num2str(amp_b_L2)]);
            % disp(['amp. |c> at L2: ', num2str(amp_c_L2)]);
            % disp(['actual |a> amp. at end : ', num2str(amp_a_end)]);
            disp(['actual |a> pop. at end : ', num2str(pop_a_end)]);

            T_max = max(time + 51.1e-6);  % Effective upper bound (shifted)
            tauGrid = linspace(0, T_max, obj.timesteps);
            
            % Define weighting and integration terms depending on inclusion of L2
            if strcmpi(includeL2, 'no')
                % Without additional L2 depletion weighting
                WTerm = @(t) (1i * obj.W / Delta_local) .* (exp(-1i * Delta_local * t) - 1);
                integrand_nr = @(tau) exp( -((a*(tau - obj.t0 - 51.1e-6)).^2)/(2*sqrt(2)) - 1i * Delta_local * tau );
            else  % includeL2 == 'yes'
                % Include additional weighting from the L2 depletion pulse
                weight = @(tau) interp1(time + 51.1e-6, data(:,1), tau, 'linear', 'extrap');
                intergrand_W = @(tau) weight(tau) .* exp(-1i * Delta_local * tau);
                I_cum_W = cumtrapz(tauGrid, intergrand_W(tauGrid));
                I_interp_W = @(tt) interp1(tauGrid, I_cum_W, tt, 'linear', 'extrap');
                WTerm = @(t) obj.W * I_interp_W(t);
                
                integrand_nr = @(tau) weight(tau) .* exp( -((a*(tau - obj.t0 - 51.1e-6)).^2)/(2*sqrt(2)) - 1i * Delta_local * tau );
            end
            
            % Compute the non-reversing term (nrTerm) via cumulative integration
            I_cum = cumtrapz(tauGrid, integrand_nr(tauGrid));
            I_interp = @(tt) interp1(tauGrid, I_cum, tt, 'linear', 'extrap');
            nrTerm = @(t) -(1i * obj.d12 * obj.E0_nr) * I_interp(t);
            
            % Stark term (common for both cases)
            starkTerm = @(t) (t >= obj.T_f1 & t <= obj.T_f1 + obj.T_e) .* ( ...
                                obj.d12 * obj.E0_stark .* exp(-1i * Delta_local * (t - obj.T_f1)) ./ (obj.omega_stark^2 - Delta_local^2) .* ...
         ( - Delta_local .* sin(obj.omega_stark * (t - obj.T_f1)) + 1i * obj.omega_stark .* cos(obj.omega_stark * (t - obj.T_f1)) ...
         - 1i * obj.omega_stark .* exp(1i * Delta_local * (t - obj.T_f1)) )).* exp(-1i* Delta_local * (obj.T_f1)) + ... 
                              (t > obj.T_f1 + obj.T_e) .* 2 * obj.omega_stark * obj.d12 * obj.E0_stark .* ...
                                   exp(-1i * Delta_local * obj.T_e / 2) .* ...
                                   sin(Delta_local * obj.T_e / 2) ./ ...
                                   (obj.omega_stark^2 - Delta_local^2).* exp(-1i* Delta_local * (obj.T_f1));
            
            % Combine all contributions to form the total TDPT amplitude function
            totalFunction = @(t) WTerm(t) + starkTerm(t) + nrTerm(t);
            fVals = totalFunction(time + 51.1e-6);
            fValsMagSq = abs(fVals).^2;
            
            
            % Plot the squared magnitude of the total function (interpreted as a population)
            plot(time*1e6, fValsMagSq, 'k-', 'LineWidth', 1);
            xlabel('Time [\mus]');
            ylabel('Population');
            title(['Population Evolution (' method '), no fix, E_{nr}=-3V/m; E_{stark}=10V/m; 10% Laser']);
            legend('|a>', '|b>', '|c>', 'TDPT','Location','Best');
            grid on;
            set(gca, 'FontName','Times New Roman','FontSize',15);
            % ylim([0 0.1]);
        end
        
        %% Calculate analytic asymmetry vs. detuning
        function obj = calculateAnalyticAsymmetry(obj)
            % calculateAnalyticAsymmetry computes the analytic asymmetry for each
            % detuning value stored in obj.detuning_range.
            %
            %   asymmetry_arr = calculateAnalyticAsymmetry(obj) returns an array of
            %   asymmetry values computed using the analytic expression.
            
            a = obj.v^2 / (2 * sqrt(2) * obj.sigma_u^2);
            t_diff = obj.t0 - obj.t_center_L2;
            N_det = length(obj.detuning_range);
            
            % Derived constants
            K1 = (2 * obj.d12 * obj.E0_stark) / obj.omega_stark;
            K3 = (- obj.d13^2 * obj.E0_L2^2 * sqrt(pi) * obj.sigma_L2) / 4;
            K4 = obj.d12 * sqrt(pi / a);
            T  = obj.T_e + obj.T_f2 + t_diff;
            phi = obj.T_e / 2;
            % disp(a);
            
            % Initialize the array for asymmetry values
            asymmetry_arr = zeros(1, N_det);
            
            % Loop over each detuning value to compute the asymmetry
            for j = 1:N_det
                Delta_local = obj.detuning_range(j);
                asymmetry_arr(j) = (-2 * K1 * K4 * obj.E0_nr * sin(Delta_local * phi) * exp( - (Delta_local)^2/(4*a) ) * ...
                    exp( (K3 * obj.gamma_c/2) / ( obj.detuning_L2^2 + (obj.gamma_c/2)^2 ) ) * ...
                    sin( Delta_local * ( T - phi ) - (K3 * obj.detuning_L2) / ( obj.detuning_L2^2 + (obj.gamma_c/2)^2 ))) / ...
                    (K1^2 * ( sin(Delta_local * phi) )^2 + obj.E0_nr^2 * K4^2 * exp( - (Delta_local)^2/(2*a) ) * ...
                    exp( (K3 * obj.gamma_c) / ( obj.detuning_L2^2 + (obj.gamma_c/2)^2 ) ) );
                % approximated form
                % gamma = -2 * K4 * obj.E0_nr * exp( (K3 * obj.gamma_c/2) / ( obj.detuning_L2^2 + (obj.gamma_c/2)^2 ) )/ K1 ;
                % beta = (K3 * obj.detuning_L2) / ( obj.detuning_L2^2 + (obj.gamma_c/2)^2 );
                % asymmetry_arr(j) = (gamma * (1 - (Delta_local)^2/(4*a) ) * ...
                %     sin( Delta_local * ( T - phi ) - beta)) / ...
                %     sin(Delta_local * phi) ;
                % asymmetry_arr(j) = -2 * gamma * sin(beta) / ( Delta_local * obj.T_e) + ...
                %     Delta_local * gamma * ((1/(a * obj.T_e)-obj.T_e/12) + (2*t_diff + obj.T_e + 2 * obj.T_f1)^2 / (4 * obj.T_e)) * sin(beta);
                % asymmetry_arr(j) = (-50*2*pi/Delta_local) .* (obj.omega_stark./(-0.5*obj.d12*obj.E0_stark)) .* ((obj.T_e + obj.T_f2) ./ obj.T_e) + ...
                %     1.1228e-05 * Delta_local;
                % asymmetry_arr(j) = -(50*2*pi/Delta_local) .* ((obj.omega_stark^2-Delta_local^2)./((-obj.d12*obj.E0_stark*obj.omega_stark))) .* ((sin(Delta_local*(obj.T_e/2+obj.T_f2/2)) ./ sin(Delta_local*obj.T_e/2)).* cos(Delta_local*obj.T_f1/2) ) + 1.1228e-05.*Delta_local;
                %   asymmetry_arr(j) = (2*50*2*pi./Delta_local) .* ...
                % ((obj.omega_stark ) / (obj.d12*obj.E0_stark)) .* ...
                % (sin((Delta_local/2)*(obj.T_e + obj.T_f1 + obj.T_f2)) ./ sin((Delta_local/2)*obj.T_e)) .* ...
                % cos((Delta_local/2)*(obj.T_f1 - obj.T_f2)) + 1.1228e-05*Delta_local;
            end
            obj.analytic_asymmetry = asymmetry_arr;
            disp(obj.analytic_asymmetry);
            % disp((K3 * obj.detuning_L2) / ( obj.detuning_L2^2 + (obj.gamma_c/2)^2 ));
            % disp(-2* K4 * exp( (K3 * obj.gamma_c) / ( obj.detuning_L2^2 + (obj.gamma_c/2)^2 )) * obj.E0_nr / K1);
        end

        %% Plot Asymmetry vs. Detuning
        function plotAsymmetry(obj)
            % plotAsymmetry Plot the asymmetry as a function of detuning.
            obj.analytic_asymmetry = obj.calculateAnalyticAsymmetry();
            
            % Convert detuning to kHz
            detuning_kHz = obj.detuning_range / (2*pi*1000);
            
            % Plot the computed asymmetry vs. detuning
            figure;
            hold on;
            plot(detuning_kHz, obj.asymmetry, 'k-', 'LineWidth', 1.5);
            plot(detuning_kHz, obj.analytic_asymmetry, 'r-', 'LineWidth', 1.5);
            xlabel('\Delta/2\pi [kHz]');
            ylabel('Asymmetry');
            legend('OBE','Analytic');
            title_str = sprintf('E_{nr} = %.0f V/m, \\delta_{L2} = %.0f MHz', ...
                    obj.E0_nr, obj.detuning_L2/(2*pi*1e6));
            title(title_str, 'Interpreter', 'tex'); 
            grid on;
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
        end

        %% Plot Electric Field vs Time
        function plotElectricField(obj)
            % plotElectricField Plot the time-dependent electric fields E(t)
            %   This function plots both the forward (E_stark + E_nr) and
            %   reversed (-E_stark + E_nr) electric fields over the simulation time.
            
            % Use the same time span as in the simulation
            time = obj.tspan-50e-6;
            
            % Define the Stark and non-reversing fields (using the same expressions)
            T = 43.7e-6;
            E_stark = (abs(time) <= T) .* (obj.E0_stark * sin(obj.omega_stark * (time + 43.7e-6)));
            % E_nr    = obj.E0_nr * sech(obj.v * (time - obj.t0) / obj.sigma_u);
            Rabi    = 300 * exp(-((time - obj.t_center_L2).^2) / (2*obj.sigma_L2^2)) * obj.E0_L2 / 8.514e2;

            % Compute the original E_nr using the sech function
            E_nr_original = obj.E0_nr * sech(obj.v * (time - obj.t0) / obj.sigma_u);
            % flatVal = obj.E0_nr * sech(obj.v * (obj.t_center_L2 - obj.t0) / obj.sigma_u);
            % delta = 1e-6;  % Adjust delta to make the transition sharper or smoother
            % 
            % % Create smooth transition masks: S1 goes from 0 to 1 around 65 µs, and S2 goes from 1 to 0 around 75 µs.
            % S1 = 0.5 * (1 + tanh((time - obj.t_center_L2 + 5e-6) / delta));
            % S2 = 0.5 * (1 + tanh((obj.t_center_L2 + 5e-6 - time) / delta));
            % mask = S1 .* S2;  % mask is ~1 between 65 and 75 µs, 0 outside, with smooth transitions
            % 
            % % Blend the original E_nr with the flat value using the mask
            % E_nr = mask * flatVal + (1 - mask) .* E_nr_original;
            

            % E_nr_original = obj.E0_nr * sech(obj.v * (time - obj.t0) / obj.sigma_u) .* (0.5*(1 - tanh((-time + obj.t_center_L2)/1e-6)));
            
            % Calculate forward and reversed total fields
            E_forward = E_stark + E_nr_original;
            
            % Plot both fields vs. time (convert time to microseconds for clarity)
            figure('Position', [100, 100, 600, 300]);
            hold on;
            plot(time*1e6, E_stark, 'Color', [0.6 0.6 0.6], 'LineStyle', ':', 'LineWidth', 1.5);
            plot(time*1e6, E_nr_original, 'Color', [0.6 0.6 0.6], 'LineStyle', ':', 'LineWidth', 1.5);
            plot(time*1e6, E_forward, 'k-', 'LineWidth', 1.5);

            plot(time*1e6, Rabi, 'r-', 'LineWidth', 1);
            xlabel('Time (\mu s)');
            ylabel('Electric Field Strength (arb. unit)');
            % set(gca, 'XTick', [], 'YTick', []);
            % title('Electric Field vs Time');
            % legend({'$\mathcal{E}_{nr} + \mathcal{E}_{\mathrm{Stark}}$', ...
            %     '$\mathcal{E}_{L2}$'}, ...
            %     'Interpreter', 'latex', 'FontSize', 16);
            % legend boxoff
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
            grid off;
            box on;
        end

        %% Analytic Prediction
        function analyticSol(obj, order)
            Delta_local = obj.detuning_range(1);
            % First, compute the "even" amplitude from the Stark field:
            stark_term_even = 2 * obj.d12 * obj.omega_stark * obj.E0_stark ...
                * exp(-1i * Delta_local * obj.T_e/2) / (obj.omega_stark^2 - Delta_local^2) ...
                * sin(Delta_local * obj.T_e/2) .* exp(-1i * Delta_local * (obj.T_f1));
            % disp(['analytic |a> at L2: ', num2str(stark_term_even)]);
        
            % Next, the "odd" amplitude factor that includes the second-order Stark correction:
            stark_term_odd = (1 + (obj.d12^2 * obj.E0_stark^2 * ((exp(1i * obj.T_e * Delta_local) - 1) ...
                * obj.omega_stark^3 - 1i*pi*(Delta_local^3 - Delta_local * obj.omega_stark^2))) ...
                / (obj.omega_stark*(Delta_local^2 - obj.omega_stark^2)^2)) ...
                 .* exp(-1i* Delta_local * (obj.T_f1 + obj.T_e + obj.T_f2 + 2e-6)); % 2mu s after L2 center
        
            % Approximate up to the Nth order of the dynamic phase expansion (L2 field)
            N = order;
            z = obj.detuning_L2 + 1i*(obj.gamma_c / 2); %unkown why it's plus instead of minus !!!!!!
            phi_dyn = 0;
            for n = 1:N
                % Generalized binomial coefficient binom(1/2, n) = Gamma(1/2+1)/(Gamma(n+1)*Gamma(1/2-n+1))
                binom = gamma(1/2 + 1) / ( gamma(n + 1) * gamma(1/2 - n + 1) );
                term = binom * ((obj.d13 * obj.E0_L2)^(2*n) * obj.sigma_L2) ...
                    / (2 * z^(2*n - 1)) * sqrt(pi / n);
                phi_dyn = phi_dyn + term;
            end
            Dynamic_phase_after_L2_odd = exp(-1i * phi_dyn); 
            % disp(['analytic |b> at L2: ', num2str(stark_term_odd * Dynamic_phase_after_L2_odd)]);
        
            % NR field contribution: compute Gaussian envelope factor
            a = obj.v^2 / (2 * sqrt(2) * obj.sigma_u^2);
            t_diff = obj.t0 - obj.t_center_L2;
        
            NR_factor = (exp(-1i * Delta_local * t_diff)) * exp(- Delta_local^2 / (4*a)) ...
                * sqrt(pi / a) * obj.d12 * obj.E0_nr;
        
            % Combine everything for c_+^(1,(0+2),N)(t2)
            c_plus = stark_term_even ...
                - 1i * stark_term_odd * Dynamic_phase_after_L2_odd * NR_factor; % * exp( -1i * 2.4);
            prob = abs(c_plus)^2;

            % disp(['analytic |a> at the end: ', num2str(c_plus)]);
            fprintf('analytic prob. (Even) = %.8f\n', prob);        
        end
        
        %% Local function: Erf for complex arguments using series expansion
        function y = erf_complex(obj, z)
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
                % Break if the new term is below tolerance for all elements
                if all(abs(term_new(:)) < tol)
                    break;
                end
                k = k + 1;
            end
        end
    end
end





