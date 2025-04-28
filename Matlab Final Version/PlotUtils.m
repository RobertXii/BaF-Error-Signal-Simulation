%% +threeLevel/PlotUtils.m
classdef PlotUtils
    methods(Static)
        function plotAsymmetry(cfg, det, asym, params, plotFlag)
            if nargin<5, plotFlag = false; end

            W_rad = params(1) * 2*pi;  % convert from Hz to rad/s
            a0    = params(2);
            a1    = params(3);

            fitFunc = @(p,Delta) ...
              (2*p(1)./Delta) .* ...
              ((cfg.omega_stark^2 - Delta.^2) ...
               ./ (cfg.d12*cfg.E0_stark*cfg.omega_stark)) .* ...
              (sin((Delta/2)*(cfg.T_e + cfg.T_f1 + cfg.T_f2)) ...
               ./ sin((Delta/2)*cfg.T_e)) .* ...
              cos((Delta/2)*(cfg.T_f1 - cfg.T_f2)) ...
              + p(2) + p(3).*Delta;
            fit_curve = fitFunc([W_rad, a0, a1], det);

            % Mask out the region |Δ| ≤ 1000 Hz to avoid divergence
            det_Hz = det/(2*pi);
            mask = abs(det_Hz) <= 1000;
            fit_curve(mask) = NaN;
            det_kHz = det/(1000*2*pi);

            Analytic_asym = PlotUtils.computeAnalyticAsymmetry(cfg);

            figure; hold on;
            plot(det_kHz, asym, 'k-', 'LineWidth', 1.5, 'MarkerFaceColor','k','DisplayName','Simulated Data');
            plot(det_kHz, Analytic_asym,'r-', 'LineWidth', 1.5, 'DisplayName','Analytic Result');
            if plotFlag
                plot(det_kHz, fit_curve,'r-', 'LineWidth', 2, 'DisplayName','Fit Curve');
            end
            xlabel('Δ/2π [kHz]'); ylabel('Asymmetry');
            legend('show','FontSize', 18);
            set(gca,'FontName','Times New Roman','FontSize',15);
            set(gcf,'Units','centimeters','Position',[10,10,14,11])
            grid off; box on;legend boxoff;
        end

        function A_analytic = computeAnalyticAsymmetry(cfg)
                num = -2 .* cfg.K1 .* cfg.K4 .* cfg.E0_nr .* ...
                      sin(cfg.detuning_range .* cfg.phi) .* ...
                      exp(-cfg.detuning_range.^2 ./ (4 .* cfg.a)) .* ...
                      cfg.alpha .* ...
                      sin( cfg.detuning_range .* (cfg.T_all - cfg.T_f1 - cfg.phi) - cfg.beta );
            
                den =   cfg.K1^2 .* sin(cfg.detuning_range .* cfg.phi).^2  + ...
                      (cfg.E0_nr^2) .* (cfg.K4^2) .* exp(-cfg.detuning_range.^2 ./ (2 .* cfg.a)) .* cfg.alpha.^2;

                A_analytic = num ./ den;
        end

        function W_analytic_mat = W_analytic(cfg, E0_nr_values, detuning_L2_multipliers)
            detuning_L2 = detuning_L2_multipliers * 1e6 * 2*pi;   
        
            denom   = detuning_L2.^2 + (cfg.gamma_c/2)^2;
            expTerm = exp((cfg.K3*(cfg.gamma_c/2)) ./ denom);
            sinTerm = sin((cfg.K3 .* detuning_L2) ./ denom);
        
            prefac_base = (cfg.d12 / cfg.T_before_L2) * sqrt(pi ./ cfg.a);
            prefac_vec  = prefac_base .* E0_nr_values(:);  
        
            W_analytic_mat_transpose = prefac_vec * (expTerm .* sinTerm) ./ (2*pi);
            W_analytic_mat = W_analytic_mat_transpose.';
        end

        function a0_analytic_mat = a0_analytic(cfg, E0_nr_values, detuning_L2_multipliers)
            detuning_L2 = detuning_L2_multipliers * 1e6 * 2*pi;   
        
            denom   = detuning_L2.^2 + (cfg.gamma_c/2)^2;
            expTerm = exp((cfg.K3*(cfg.gamma_c/2)) ./ denom);
        
            prefac_base = (cfg.omega_stark .* (2.*cfg.T_diff + 2.*cfg.T_all - cfg.T_e) / (cfg.E0_stark * cfg.T_e)) * sqrt(pi ./ cfg.a);
            prefac_vec  = prefac_base .* E0_nr_values(:);  
        
            a0_analytic_mat_transpose = prefac_vec * expTerm;
            a0_analytic_mat = a0_analytic_mat_transpose.';
        end

        function a1_analytic_mat = a1_analytic(cfg, E0_nr_values, detuning_L2_multipliers)
            detuning_L2 = detuning_L2_multipliers * 1e6 * 2*pi;   
        
            denom   = detuning_L2.^2 + (cfg.gamma_c/2)^2;
            expTerm = exp((cfg.K3*(cfg.gamma_c/2)) ./ denom);
            sinTerm = sin((cfg.K3 .* detuning_L2) ./ denom);
        
            timeTerm = ( 1./(cfg.a * cfg.T_e) ...
                  - cfg.T_e/12 ) ...
                + ((2*cfg.T_diff + 2*cfg.T_all - cfg.T_e).^2) ...
                  ./ (4*cfg.T_e);

            prefac_base = -(cfg.omega_stark / cfg.E0_stark) ...
                  * sqrt(pi ./ cfg.a) ...
                  * timeTerm;            
            prefac_vec  = prefac_base .* E0_nr_values(:);  
        
            a1_analytic_mat_transpose = prefac_vec * (expTerm .* sinTerm);
            a1_analytic_mat = a1_analytic_mat_transpose.';
        end

        function plotField(cfg,profile)
            t = cfg.tspan;
            E_tot = profile.E_total(t);
            figure('Position', [100, 100, 600, 300]); hold on;
            plot(t*1e6, profile.E_nr(t), 'k:', 'LineWidth',1);
            plot(t*1e6, profile.E_stark(t), 'k:', 'LineWidth',1);
            plot(t*1e6, 30 * profile.Rabi(t)/(cfg.E0_L2 * cfg.d13 / 2), 'r-', 'LineWidth',1.5);
            plot(t*1e6, E_tot, 'k-', 'LineWidth',1.5);
            xlabel('Time [\mus]'); xlim([-60,210]);
            ylabel('E(t) (arb. unit)'); ylim([-50, 50]);
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
            grid off; box on;
        end

        function plotXvsEnr(X, E_nr_vals, detuning_multipliers, yLabel)
            % This plot X vs E_nr as function of L2 detuning, with input X
            % a nL2 by nr matrix
            nDet = numel(detuning_multipliers);
            data_labels = {'$\delta_{L2} = -3$ MHz', '$\delta_{L2} = -2$ MHz', ...
               '$\delta_{L2} = -1$ MHz', '$\delta_{L2} = 0$ MHz', ...
               '$\delta_{L2} = +1$ MHz', '$\delta_{L2} = +2$ MHz', ...
               '$\delta_{L2} = +3$ MHz'};

            E_vals = E_nr_vals * 10;
            cmap = jet(nDet);
            figure; hold on;
            for i = 1:nDet
                plot(E_vals, X(i,:), 'o', 'Color', cmap(i,:), ...
                     'MarkerFaceColor', cmap(i,:), 'LineWidth',1.5);
                p = polyfit(E_vals, X(i,:), 1);
                E_fit = linspace(min(E_vals), max(E_vals), 100);
                W_fit = polyval(p, E_fit);
                plot(E_fit, W_fit, '--', 'Color', cmap(i,:), 'LineWidth',1.5, 'HandleVisibility','off');
            end
            xlabel('E_{nr0} [mV/cm]');
            ylabel(yLabel);
            yMin = 1.2 * min(X(:));
            yMax = 1.2 * max(X(:));
            ylim([yMin, yMax]);
            legend(data_labels, 'Interpreter', 'latex', 'Location', 'eastoutside');
            grid off; box on;
            set(gca,'FontName','Times New Roman','FontSize',15);
            set(gcf,'Units','centimeters','Position',[5,5,25,8]);
        end

        function [W_table, a0_table, a1_table] = create_Detuning_strength_scan_Tables(...
        W_mat, a0_mat, a1_mat, detuning_L2_multipliers, E0_nr_values)
            rowNames = compose('%+dMHz', detuning_L2_multipliers);
            varNames = compose('%+dV/m', E0_nr_values);
        
            W_table  = array2table(W_mat,  'RowNames',rowNames,'VariableNames',varNames);
            a0_table = array2table(a0_mat, 'RowNames',rowNames,'VariableNames',varNames);
            a1_table = array2table(a1_mat, 'RowNames',rowNames,'VariableNames',varNames);
        
            disp('W/2π [Hz]'); disp(W_table);
            disp('a0');        disp(a0_table);
            disp('a1 [1/H   z]');disp(a1_table);
        end 

        function [W_table, a0_table, a1_table] = createPositionScanTables(W_mat, a0_mat, a1_mat, nr_pos, L2_pos)
            rowNames = compose('nr=%dμs', nr_pos);
            varNames = compose('L2=%dμs', L2_pos);
            W_table  = array2table(W_mat,  'RowNames',rowNames,'VariableNames',varNames);
            a0_table = array2table(a0_mat, 'RowNames',rowNames,'VariableNames',varNames);
            a1_table = array2table(a1_mat, 'RowNames',rowNames,'VariableNames',varNames);
            disp('W/2π [Hz]');     disp(W_table);
            disp('a0');            disp(a0_table);
            disp('a1 [1/Hz]');     disp(a1_table);
        end

        function plotPositionMap(cfg, dataMat, nr_pos, L2_pos, unit)
            %   unit    : 'time' (default) or 'distance'
            if nargin<5 || isempty(unit), unit = 'time'; end
        
            switch lower(unit)
              case 'time'
                xVals = nr_pos;       xLabel = 'NR position [\mus]';
                yVals = L2_pos;       yLabel = 'L2 position [\mus]';
                cross_x = 70.1;
                cross_y = 52.6;
                x_line = (cfg.T_e/2)*1e6;
              case 'distance'
                xVals = nr_pos*1e-6*cfg.v*1e3;  xLabel = 'NR position [mm]';
                yVals = L2_pos*1e-6*cfg.v*1e3;  yLabel = 'L2 position [mm]';
                cross_x = 70.1*cfg.v*1e3;
                cross_y = 52.6*cfg.v*1e3;
                x_line = (cfg.T_e/2)*cfg.v*1e3;
              otherwise
                error('Unit must be ''time'' or ''distance''.');
            end
        
            xi = linspace(min(xVals), max(xVals), max(2,numel(xVals))*10);
            yi = linspace(min(yVals), max(yVals), max(2,numel(yVals))*10);
            [Xq,Yq] = meshgrid(xi, yi);
        
            Zq = interp2(xVals, yVals, dataMat.', Xq, Yq, 'cubic');
        
            figure;
            imagesc(xi, yi, Zq);
            set(gca,'YDir','normal');
            hold on;
            % yl = ylim;
            % xl = xlim;
            % plot([x_line x_line], yl, 'w--', 'LineWidth',1);
            % plot(xl, [x_line x_line], 'w--', 'LineWidth',1);
            % plot([-x_line -x_line], yl, 'w--', 'LineWidth',1);
            % plot(xl, [-x_line -x_line], 'w--', 'LineWidth',1);
            plot(cross_x, cross_y, 'wx', 'MarkerSize',10, 'LineWidth', 3);

            width_cm = 12;
            height_cm = width_cm * (numel(yVals)/numel(xVals));
            set(gcf,'Units','centimeters','Position',[10,10,width_cm+2,height_cm]);

            colormap jet; 
            colorbar; clim([0 60]);
            xlabel(xLabel);
            ylabel(yLabel);
            set(gca,'FontName','Times New Roman','FontSize',15, 'TickDir','out');
        end
    end
end


