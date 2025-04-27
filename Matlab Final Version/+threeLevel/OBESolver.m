%% +threeLevel/@OBESolver/OBESolver.m
classdef OBESolver
    % Solves OBE evolution and computes asymmetry
    properties
        cfg SimConfig
        profile FieldProfile
    end
    methods
        function obj = OBESolver(cfg,profile)
            obj.cfg = cfg;
            obj.profile = profile;
        end
        function [asym, sol] = solve(obj, Delta)
            Hf = @(t) obj.buildH(Delta, obj.profile.E_stark(t)+obj.profile.E_nr(t), obj.profile.Rabi(t));
            Hr = @(t) obj.buildH(Delta, -obj.profile.E_stark(t)+obj.profile.E_nr(t), obj.profile.Rabi(t));
            opts = odeset('RelTol',1e-13,'AbsTol',1e-12);
            [~, sol_f] = ode45(@(t,r) obj.masterEq(t,r,Hf), obj.cfg.tspan, obj.initState(), opts);
            [~, sol_r] = ode45(@(t,r) obj.masterEq(t,r,Hr), obj.cfg.tspan, obj.initState(), opts);
            rho_f = reshape(sol_f(end,:),3,3);
            rho_r = reshape(sol_r(end,:),3,3);
            S_p = real(rho_f(1,1));
            S_m = real(rho_r(1,1));
            asym = (S_p - S_m)/(S_p + S_m);
            sol = sol_f;
        end
        function dr = masterEq(obj,t,rho_vec,Hfun)
            rho = reshape(rho_vec,3,3);
            H = Hfun(t);
            comm = -1i*(H*rho - rho*H);
            diss = -obj.cfg.gamma_c/2 * [0,0,rho(1,3);0,0,rho(2,3);rho(3,1),rho(3,2),2*rho(3,3)];
            dr = (comm + diss)(:);
        end
        function rho0 = initState(obj)
            rho0 = zeros(3); rho0(1,1)=1; rho0 = rho0(:);
        end
        function H = buildH(obj,Delta,Est_nr,R)
            H = [0,1i*obj.cfg.W+obj.cfg.d12*Est_nr,0;
                -1i*obj.cfg.W+obj.cfg.d12*Est_nr,Delta,R;
                 0,R,-obj.cfg.detuning_L2];
        end
    end
end
