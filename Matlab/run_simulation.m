%% run_simulation.m
% -------------------------------------------------------------------------
% This script runs a detuning scan for the quantum simulation and fits the
% asymmetry data to extract the parameters: W, a0, and a1.
%
% The simulation uses the QuantumSimulation class to:
%   1. Set the simulation parameters (detunings, field amplitudes, etc.).
%   2. Run the simulation over a range of parameters.
%   3. Fit the asymmetry data using a theoretical model.
%   4. Store the fitted parameters in tables.
%   5. Plot the evolution of the population and asymmetry curves.
%
% Before running the script, ensure that the QuantumSimulation class is
% in the MATLAB path.
%
% -------------------------------------------------------------------------
% clearvars; close all; clc;  % Uncomment to clear workspace and close figures

%% Define Scan Parameters
% Modify these arrays to change the scan range for the L2 detuning multiplier
% and the non-reversing field amplitude (E0_nr).
detuning_L2_multipliers = 3:1:3;          % L2 detuning multpliers [MHz]
% E0_nr_values = -12:3:12;                  % Non-reversing field amplitudes [V/m]
E0_nr_values = 12;                   % Non-reversing field amplitudes [V/m]

% Determine the number of scan points for each parameter
nDet = numel(detuning_L2_multipliers);
nE0  = numel(E0_nr_values);

% Preallocate matrices for the fitted parameters:
%   - W_matrix  : Extracted W/2pi in Hz
%   - a0_matrix : Fitted parameter a0
%   - a1_matrix : Fitted parameter a1 (in 1/Hz)
W_matrix  = zeros(nDet, nE0);
a0_matrix = zeros(nDet, nE0);
a1_matrix = zeros(nDet, nE0);

%% Loop Over Parameter Values
% For each combination of L2 detuning multiplier and E0_nr, create a new
% simulation instance, set the simulation parameters, run the simulation,
% and fit the asymmetry data.
for i = 1:nDet
    for j = 1:nE0
        % Create a new instance of the QuantumSimulation class
        sim = QuantumSimulation();
        
        % -----------------------------------------------------------------
        % Set simulation parameters:
        %   - Delta       : Energy gap for state |b> [rad/s]
        %   - E0_stark    : Stark field amplitude [V/m]
        %   - E0_nr       : Non-reversing field amplitude [V/m]
        %   - E0_L2       : Depletion laser amplitude [V/m]
        %   - detuning_L2 : L2 laser detuning [rad/s]
        % -----------------------------------------------------------------
        % sim.Delta = 2*pi*1e3 * 2; %2*pi*1e3 * 2
        sim.E0_stark = 40; %40
        sim.E0_nr = E0_nr_values(j);
        sim.E0_L2 = 8.514e2 * 0.11;          % 8.514e2*1 is theoretical amplitude for a 8mW laser
        sim.detuning_L2 = 2*pi*1e6 * detuning_L2_multipliers(i);
        
        % -----------------------------------------------------------------
        % Run the simulation:
        % Choose the evolution method by uncommenting the desired option.
        % -----------------------------------------------------------------
        % sim = sim.runVary_Detuning('OBE');        % Use OBE formalism
        sim = sim.runVary_Detuning('Schrodinger');   % Use Schrödinger evolution
        sim = sim.calculateAnalyticAsymmetry();
        fitResults = sim.fitAsymmetry(true);
        
        % Store the fitted parameters:
        %   - fitResults(1) is W in [rad/s], so divide by 2pi for Hz.
        W_matrix(i, j)  = fitResults(1) / (2*pi);
        a0_matrix(i, j) = fitResults(2);
        a1_matrix(i, j) = fitResults(3);
    end
end

%% Create and Display Tables of Fitted Parameters
% Generate table row names (from L2 detuning multipliers in MHz) and column
% names (from E0_nr values in V/m), then create and display the tables.
rowNames = compose('%d', detuning_L2_multipliers);
varNames = compose('%d', E0_nr_values);

W_table  = array2table(W_matrix, 'RowNames', rowNames, 'VariableNames', varNames);
a0_table = array2table(a0_matrix, 'RowNames', rowNames, 'VariableNames', varNames);
a1_table = array2table(a1_matrix, 'RowNames', rowNames, 'VariableNames', varNames);

disp('W/2pi (Hz)');
disp(W_table);
disp('a0');
disp(a0_table);
disp('a1 (1/Hz)');
disp(a1_table);

%% Plot Simulation Results
% Plot the population evolution using the OBE/Schrödinger formalism:
% sim.plotPopulationEvolution(1, 'OBE', 'yes');
% sim.plotPopulationEvolution(1, 'Schrodinger', 'yes');

%% Plot E field
sim.plotElectricField();

%% Output analytic result
% sim.analyticSol(3);

%% plot the asymmetry vs. detuning curve:
sim.plotAsymmetry();


