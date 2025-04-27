%% run_simulation.m
% Main script
detuning_L2_multipliers = -3:1:3;   % L2 detuning multipliers [MHz]
E0_nr_values           = -12:3:12; % Non-reversing field amplitudes [V/m]
nr_pos = 40:10:150;                
L2_pos = 40:10:150;

cfg = SimConfig.default();
cfg.E0_L2 = 1204.1*0.09;
cfg.E0_nr = 3;
cfg.t_nr = 70.1e-6; % 70.1e-6
cfg.t_L2 = 52.6e-6; % 52.6e-6
cfg.detuning_L2 = 2*pi*1e6 * 3; 
cfg.tspan = linspace(-51.1e-6,200e-6,10000);
cfg.detuning_range = linspace(-4000, 4000, 100) * 2*pi;
cfg.fieldMethod = 'default';  % 'before', 'after','modified','default'

% profile = threeLevel.FieldProfile(cfg);
% solver  = threeLevel.SchrodingerSolver(cfg, profile);
% scan    = threeLevel.ParameterScan(solver, cfg);
% results = scan.run();
% 
% fitter  = threeLevel.AsymmetryFitter(results, cfg);
% params  = fitter.fit();

% L2 detuning and Enr strength scan
[W_mat,a0_mat,a1_mat] = threeLevel.ParameterScan.scan_L2_detune_and_Enr(...
    detuning_L2_multipliers, E0_nr_values, cfg, @threeLevel.SchrodingerSolver);
    
% L2 and Enr position scan
% [W_p, a0_p, a1_p] = threeLevel.ParameterScan.scan_L2_Enr_position(...
%     nr_pos, L2_pos, cfg, @threeLevel.SchrodingerSolver);
%%
% PlotUtils.plotAsymmetry(cfg, results.detuning, results.asymmetry, params, false);
% profile = threeLevel.FieldProfile(cfg);
% PlotUtils.plotField(cfg, profile);

% [W_tab_pos,a0_tab_pos,a1_tab_pos] = PlotUtils.createPositionScanTables(...
%     W_p, a0_p, a1_p, nr_pos, L2_pos);
% 
% PlotUtils.plotPositionMap(cfg, W_p, nr_pos, L2_pos, 'time');

% Plot Simuation result for W, a0, a1
[W_tab,a0_tab,a1_tab] = PlotUtils.create_Detuning_strength_scan_Tables(...
    W_mat, a0_mat, a1_mat, detuning_L2_multipliers, E0_nr_values);
PlotUtils.plotXvsEnr(W_mat, E0_nr_values, detuning_L2_multipliers, 'W/2\pi[Hz]');
PlotUtils.plotXvsEnr(a0_mat, E0_nr_values, detuning_L2_multipliers, 'a_0');
PlotUtils.plotXvsEnr(a1_mat, E0_nr_values, detuning_L2_multipliers, 'a_1[1/Hz]');


% Plot analytic result for W, a0, a1
% W_analytic_mat = PlotUtils.W_analytic(cfg, E0_nr_values, detuning_L2_multipliers);
% PlotUtils.plotXvsEnr(W_analytic_mat, E0_nr_values, detuning_L2_multipliers, 'W/2\pi[Hz]');
% a0_analytic_mat = PlotUtils.a0_analytic(cfg, E0_nr_values, detuning_L2_multipliers);
% PlotUtils.plotXvsEnr(a0_analytic_mat, E0_nr_values, detuning_L2_multipliers, 'a_0');
% a1_analytic_mat = PlotUtils.a1_analytic(cfg, E0_nr_values, detuning_L2_multipliers);
% PlotUtils.plotXvsEnr(a1_analytic_mat, E0_nr_values, detuning_L2_multipliers, 'a_1[1/Hz]');
