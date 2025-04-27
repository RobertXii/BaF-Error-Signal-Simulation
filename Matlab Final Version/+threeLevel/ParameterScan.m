%% +threeLevel/ParameterScan.m
classdef ParameterScan
    % Runs parallel scan over detuning values
    properties
        solver
        cfg SimConfig
    end
    methods
        function obj = ParameterScan(sol,cfg)
            obj.solver = sol;
            obj.cfg = cfg;
        end
        
        function res = run(obj)
            N = numel(obj.cfg.detuning_range);
            asym = zeros(1,N);
            sols = cell(1,N);
            parfor k=1:N
                e_detuning = obj.cfg.detuning_range(k);
                [asym(k), sols{k}] = obj.solver.solve(e_detuning);
            end
            res.detuning = obj.cfg.detuning_range;
            res.asymmetry = asym;
            res.solutions = sols;
        end
    end

    methods(Static)
        function [W_mat, a0_mat, a1_mat] = scan_L2_detune_and_Enr(detuning_multipliers, E0nr_vals, cfg, SolverClass)
            nD = numel(detuning_multipliers);
            nE = numel(E0nr_vals);
            W_mat  = zeros(nD, nE);
            a0_mat = zeros(nD, nE);
            a1_mat = zeros(nD, nE);

            for i = 1:nD
                % set L2 detuning in rad/s
                cfg.detuning_L2 = 2*pi*1e6 * detuning_multipliers(i);
                for j = 1:nE
                    cfg.E0_nr = E0nr_vals(j);

                    profile = threeLevel.FieldProfile(cfg);
                    solver  = SolverClass(cfg, profile);
                    scan    = threeLevel.ParameterScan(solver, cfg);
                    results = scan.run();

                    fitter  = threeLevel.AsymmetryFitter(results, cfg);
                    params  = fitter.fit();

                    W_mat(i, j)  = params(1);
                    a0_mat(i, j) = params(2);
                    a1_mat(i, j) = params(3);

                    disp(params);
                end
            end
        end

        function [W_mat_pos, a0_mat_pos, a1_mat_pos] = scan_L2_Enr_position(nr_pos, L2_pos, cfg, SolverClass)
            nnr = numel(nr_pos);
            nL2 = numel(L2_pos);
            W_mat_pos  = zeros(nnr, nL2);
            a0_mat_pos = zeros(nnr, nL2);
            a1_mat_pos = zeros(nnr, nL2);

            cfg.detuning_L2 = 2*pi*1e6 * -3;
            cfg.E0_nr = 12;

            for i = 1:nnr
                cfg.t_nr = nr_pos(i) * 1e-6;            % NR field center time [s]
                for j = 1:nL2
                    cfg.t_L2 = L2_pos(j) * 1e-6; % L2 pulse center time [s]

                    profile = threeLevel.FieldProfile(cfg);
                    solver  = SolverClass(cfg, profile);
                    scan    = threeLevel.ParameterScan(solver, cfg);
                    results = scan.run();

                    fitter  = threeLevel.AsymmetryFitter(results, cfg);
                    params  = fitter.fit();

                    W_mat_pos(i, j)  = params(1);
                    a0_mat_pos(i, j) = params(2);
                    a1_mat_pos(i, j) = params(3);

                    disp(params);
                end
            end
        end
    end
end