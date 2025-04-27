%% +threeLevel/FieldProfile.m
classdef FieldProfile
    % Defines time-dependent electric fields
    properties
        cfg SimConfig
    end
    methods
        function obj = FieldProfile(cfg)
            obj.cfg = cfg;
        end
        function E = E_stark(obj,t)
            T = 43.7e-6;
            E = (abs(t)<=T) .* (obj.cfg.E0_stark * sin(obj.cfg.omega_stark * (t+T)));
        end
        function E = E_nr(obj,t)
            E = obj.cfg.E0_nr * sech(obj.cfg.v * (t - obj.cfg.t_nr) / obj.cfg.sigma_u);
        end
        function E = E_nr_afterL2(obj,t)
            E = obj.cfg.E0_nr * sech(obj.cfg.v * (t - obj.cfg.t_nr) / obj.cfg.sigma_u) .* (0.5*(1 + tanh((t - obj.cfg.t_L2)/1e-6)));
        end
        function E = E_nr_beforeL2(obj,t)
            E = obj.cfg.E0_nr * sech(obj.cfg.v * (t - obj.cfg.t_nr) / obj.cfg.sigma_u) .* (0.5*(1 + tanh((obj.cfg.t_L2 - t)/1e-6)));
        end
        function E = E_nr_modified(obj,t)
            flatVal = obj.cfg.E0_nr * sech(obj.cfg.v * (obj.cfg.t_L2 - obj.cfg.t_nr) / obj.cfg.sigma_u);
            delta = 1e-6;
            S1 = 0.5 * (1 + tanh((t - obj.cfg.t_L2 + 5e-6) / delta));
            S2 = 0.5 * (1 + tanh((obj.cfg.t_L2 + 5e-6 - t) / delta));
            mask = S1 .* S2;
            E = mask .* flatVal + (1 - mask) .* obj.E_nr(t);
        end
        function E = E_total(obj, t)
            switch obj.cfg.fieldMethod
                case 'default'
                    E = obj.E_stark(t) + obj.E_nr(t);
                case 'before'
                    E = obj.E_stark(t) + obj.E_nr_beforeL2(t);
                case 'after'
                    E = obj.E_stark(t) + obj.E_nr_afterL2(t);
                case 'modified'
                    E = obj.E_stark(t) + obj.E_nr_modified(t);
                otherwise
                    error('Unknown fieldMethod: %s', obj.cfg.fieldMethod);
            end
        end
        function R = Rabi(obj,t)
            half = obj.cfg.E0_L2 * obj.cfg.d13 / 2;
            R = half * exp(-((t - obj.cfg.t_L2).^2) / (2 * obj.cfg.sigma_L2^2));
        end
    end
end
