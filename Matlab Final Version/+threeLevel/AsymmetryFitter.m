%% +threeLevel/AsymmetryFitter.m
classdef AsymmetryFitter
    % Fits asymmetry data to extract W/2π (Hz), a0, and a1
    properties
        results  
        cfg      
    end
    methods
        function obj = AsymmetryFitter(results, cfg)
            obj.results = results;
            obj.cfg = cfg;
        end

        function params = fit(obj)
            % FIT  Fit the computed asymmetry vs. detuning.
            %   params = fit(obj) returns [W/2π (Hz), a0, a1].

            det_Hz = obj.results.detuning / (2*pi);
            mask   = abs(det_Hz) > 1000;          % |Δ| > 1 kHz
            Delta_fit    = det_Hz(mask) * 2*pi;       % back to rad/s
            Asym_fit    = obj.results.asymmetry(mask);

            fitFunc = @(p,Delta) ...
              (2*p(1)./Delta) .* ...
              ((obj.cfg.omega_stark^2 - Delta.^2) ...
               ./ (obj.cfg.d12*obj.cfg.E0_stark*obj.cfg.omega_stark)) .* ...
              (sin((Delta/2)*(obj.cfg.T_e + obj.cfg.T_f1 + obj.cfg.T_f2)) ...
               ./ sin((Delta/2)*obj.cfg.T_e)) .* ...
              cos((Delta/2)*(obj.cfg.T_f1 - obj.cfg.T_f2)) ...
              + p(2) + p(3).*Delta;

            errFun = @(p) sum((Asym_fit - fitFunc(p, Delta_fit)).^2);
            opts   = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);
            init   = [1, 0.01, 1e-8];              % [W (rad/s), a0, a1]
            bestP  = fminsearch(errFun, init, opts);

            params = [ bestP(1)/(2*pi), bestP(2), bestP(3) ];
        end
    end
end