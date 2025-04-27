%% +threeLevel/@SchrodingerSolver/SchrodingerSolver.m
classdef SchrodingerSolver
    % Non-Hermitian Schrödinger evolution + asymmetry
    properties
        cfg SimConfig
        profile
    end

    methods
        function obj = SchrodingerSolver(cfg,profile)
            obj.cfg = cfg;
            obj.profile = profile;
        end

        function [asym, sol] = solve(obj, Delta)
            Hf = @(t) obj.buildH(Delta, obj.profile.E_total(t), obj.profile.Rabi(t));
            Hr = @(t) obj.buildH(Delta, obj.profile.E_total_reversed(t), obj.profile.Rabi(t));
            psi_f = obj.evolve(Hf);
            psi_r = obj.evolve(Hr);
            S_p = abs(psi_f(end,1)).^2;
            S_m = abs(psi_r(end,1)).^2;
            asym = (S_p - S_m)/(S_p + S_m);
            sol = psi_f;
        end

        function psi = evolve(obj,Hfun)
            psi0 = [0;1;0];
            P3 = diag([0,0,1]);
            Heff = @(t) Hfun(t) - 1i*(obj.cfg.gamma_c/2)*P3;
            odeF = @(t,psi) -1i*Heff(t)*psi;
            opts = odeset('RelTol',1e-13,'AbsTol',1e-12);
            [~,psi] = ode45(odeF, obj.cfg.tspan, psi0, opts);
        end

        function H = buildH(obj,Delta,Est_nr,R)
            H = [0,         1i*obj.cfg.W+obj.cfg.d12*Est_nr,        0;
                -1i*obj.cfg.W+obj.cfg.d12*Est_nr,       Delta,      R;
                 0,                 R,           -obj.cfg.detuning_L2];
        end
    end
end
