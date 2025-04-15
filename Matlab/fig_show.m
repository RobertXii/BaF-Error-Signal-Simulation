%% E_nr vs.W, a_0, a_1
E_nr_values_numeric = [-120,	-90,	-60,	-30,	0,	30,	60,	90,	120];
% Data = [
%     -32.693,-26.451,-18.68,-9.6835,-7.08E-14,9.6835,18.68,26.451,32.693;
% 0.29486,0.2305,0.15857,0.08086,-3.99E-17,-0.08086,-0.15857,-0.2305,-0.29486;
% -5.32E-06,-4.47E-06,-3.25E-06,-1.71E-06,-3.32E-20,1.71E-06,3.25E-06,4.47E-06,5.32E-06;
% -20.01,-15.794,-10.935,-5.5964,-7.08E-14,5.5964,10.935,15.794,20.01;
% 0.24047,0.18597,0.12685,0.064336,-3.99E-17,-0.064336,-0.12685,-0.18597,-0.24047;
% -3.29E-06,-2.67E-06,-1.89E-06,-9.78E-07,-3.32E-20,9.78E-07,1.89E-06,2.67E-06,3.29E-06;
% -4.4166,-3.4821,-2.4082,-1.2314,-7.08E-14,1.2314,2.4082,3.4821,4.4166;
% 0.2571,0.19939,0.13633,0.069243,-3.99E-17,-0.069243,-0.13633,-0.19939,-0.2571;
% -6.82E-07,-5.56E-07,-3.93E-07,-2.04E-07,-3.32E-20,2.04E-07,3.93E-07,5.56E-07,6.82E-07;
% -4.45E-09,3.86E-09,-1.87E-09,3.65E-09,-7.08E-14,-1.43E-09,-2.26E-09,5.28E-09,1.35E-08;
% 0.28833,0.22471,0.15424,0.078544,-3.99E-17,-0.078544,-0.15424,-0.22471,-0.28833;
% -2.54E-15,1.56E-15,-3.47E-16,4.86E-16,-3.32E-20,-2.31E-16,-1.22E-15,1.13E-15,5.08E-15;
% 4.4166,3.4821,2.4082,1.2314,-7.08E-14,-1.2314,-2.4082,-3.4821,-4.4166;
% 0.2571,0.19939,0.13633,0.069243,-3.99E-17,-0.069243,-0.13633,-0.19939,-0.2571;
% 6.82E-07,5.56E-07,3.93E-07,2.04E-07,-3.32E-20,-2.04E-07,-3.93E-07,-5.56E-07,-6.82E-07;
% 20.01,15.794,10.935,5.5964,-7.08E-14,-5.5964,-10.935,-15.794,-20.01;
% 0.24047,0.18597,0.12685,0.064336,-3.99E-17,-0.064336,-0.12685,-0.18597,-0.24047;
% 3.29E-06,2.67E-06,1.89E-06,9.78E-07,-3.32E-20,-9.78E-07,-1.89E-06,-2.67E-06,-3.29E-06;
% 32.693,26.451,18.68,9.6835,-7.08E-14,-9.6835,-18.68,-26.451,-32.693;
% 0.29486,0.2305,0.15857,0.08086,-3.99E-17,-0.08086,-0.15857,-0.2305,-0.29486;
% 5.32E-06,4.47E-06,3.25E-06,1.71E-06,-3.32E-20,-1.71E-06,-3.25E-06,-4.47E-06,-5.32E-06;
% ];
% 
% Data labels for the legend
% data_labels = {'$\delta_{L2} = -3$ MHz', '$\delta_{L2} = -2$ MHz', '$\delta_{L2} = -1.5$ MHz',...
%                '$\delta_{L2} = -1$ MHz', '$\delta_{L2} = 0$ MHz', ...
%                '$\delta_{L2} = +1$ MHz', '$\delta_{L2} = +1.5$ MHz','$\delta_{L2} = +2$ MHz', ...
%                '$\delta_{L2} = +3$ MHz'};
data_labels = {'$\delta_{L2} = -3$ MHz', '$\delta_{L2} = -2$ MHz', ...
               '$\delta_{L2} = -1$ MHz', '$\delta_{L2} = 0$ MHz', ...
               '$\delta_{L2} = +1$ MHz', '$\delta_{L2} = +2$ MHz', ...
               '$\delta_{L2} = +3$ MHz'};
len = length(data_labels);
% 
% % Number of groups
% numGroups = size(Data, 1) / 3; % Number of 3-row groups (21 rows total)
% 
% % Preallocate for groups
% W = zeros(numGroups, size(Data, 2));   % W (1+3N rows)
% a_0 = zeros(numGroups, size(Data, 2)); % a_0 (2+3N rows)
% a_1 = zeros(numGroups, size(Data, 2)); % a_1 (3N rows)
% 
% % Extract groups
% for i = 1:numGroups
%     W(i, :) = Data(3*(i-1) + 1, :);   % 1 + 3N row
%     a_0(i, :) = Data(3*(i-1) + 2, :); % 2 + 3N row
%     a_1(i, :) = Data(3*(i-1) + 3, :); % 3N row
% end

W = [51.596	43.023	31.228	16.509	1.72E-14	-16.509	-31.228	-43.023	-51.596
29.809	23.177	15.871	8.0679	1.72E-14	-8.0679	-15.871	-23.177	-29.809
7.3556	5.5286	3.6915	1.8474	1.72E-14	-1.8474	-3.6915	-5.5286	-7.3556
-2.90E-11	-1.18E-09	-5.12E-10	-1.01E-10	1.72E-14	2.24E-10	3.91E-10	-1.33E-10	2.25E-09
-7.3556	-5.5286	-3.6915	-1.8474	1.72E-14	1.8474	3.6915	5.5286	7.3556
-29.809	-23.177	-15.871	-8.0679	1.72E-14	8.0679	15.871	23.177	29.809
-51.596	-43.023	-31.228	-16.509	1.72E-14	16.509	31.228	43.023	51.596];  

a_0 = [-0.025487	-0.020104	-0.013942	-0.0071504	4.92E-16	0.0071504	0.013942	0.020104	0.025487
-0.11329	-0.086381	-0.058293	-0.029365	4.92E-16	0.029365	0.058293	0.086381	0.11329
-0.031269	-0.023475	-0.015661	-0.0078337	4.92E-16	0.0078337	0.015661	0.023475	0.031269
0.013782	0.010337	0.0068916	0.0034459	4.92E-16	-0.0034459	-0.0068916	-0.010337	-0.013782
-0.031269	-0.023475	-0.015661	-0.0078337	4.92E-16	0.0078337	0.015661	0.023475	0.031269
-0.11329	-0.086381	-0.058293	-0.029365	4.92E-16	0.029365	0.058293	0.086381	0.11329
-0.025487	-0.020104	-0.013942	-0.0071504	4.92E-16	0.0071504	0.013942	0.020104	0.025487]; 

a_1 = [1.01E-05	8.78E-06	6.56E-06	3.54E-06	2.24E-20	-3.54E-06	-6.56E-06	-8.78E-06	-1.01E-05
6.24E-06	4.91E-06	3.40E-06	1.74E-06	2.24E-20	-1.74E-06	-3.40E-06	-4.91E-06	-6.24E-06
1.58E-06	1.19E-06	7.95E-07	3.98E-07	2.24E-20	-3.98E-07	-7.95E-07	-1.19E-06	-1.58E-06
5.97E-16	-1.21E-16	1.59E-16	1.96E-17	2.24E-20	-1.29E-16	-1.61E-16	-1.37E-16	1.44E-15
-1.58E-06	-1.19E-06	-7.95E-07	-3.98E-07	2.24E-20	3.98E-07	7.95E-07	1.19E-06	1.58E-06
-6.24E-06	-4.91E-06	-3.40E-06	-1.74E-06	2.24E-20	1.74E-06	3.40E-06	4.91E-06	6.24E-06
-1.01E-05	-8.78E-06	-6.56E-06	-3.54E-06	2.24E-20	3.54E-06	6.56E-06	8.78E-06	1.01E-05];

% Define a colormap with 7 colors (one for each delta value)
colors = jet(len); % From blue (low delta) to red (high delta)

% Plot W vs E_nr_values_numeric
figure;
hold on;
for i = 1:7 % We have 7 delta values
    plot(E_nr_values_numeric, W(i, :), 'o', 'LineWidth', 1.5, ...
         'Color', colors(i, :), 'MarkerFaceColor', colors(i, :));
end
hold off;
xlabel('E_{nr0}[mV/cm]'); 
ylabel('W/2\pi[Hz]');
grid on;
legend(data_labels, 'Interpreter', 'latex', 'Location', 'eastoutside');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
set(gcf, 'Units', 'centimeters', 'Position', [5, 5, 25, 8]); 

figure;
hold on;
fits = cell(1, len); % Store fit objects if needed later
slopes = zeros(1, len); % Store the slopes

for i = 1:len % We have "len" delta values
    % Scatter plot of data points with legend
    plot(E_nr_values_numeric, W(i, :), 'o', 'LineWidth', 1.5, ...
         'Color', colors(i, :), 'MarkerFaceColor', colors(i, :), ...
         'DisplayName', data_labels{i});
    
    % Fit the data with a linear function
    fit_result = polyfit(E_nr_values_numeric, W(i, :), 1);
    fits{i} = fit_result; % Store the fit result if needed
    
    % Extract slope
    slopes(i) = fit_result(1);
    
    % Generate fitted values for plotting
    E_fit = linspace(min(E_nr_values_numeric), max(E_nr_values_numeric), 100);
    W_fit = polyval(fit_result, E_fit);
    
    % Plot the linear fit as a dotted line (without legend)
    plot(E_fit, W_fit, '--', 'LineWidth', 1.5, 'Color', colors(i, :), 'HandleVisibility', 'off');
end

hold off;
xlabel('E_{nr0}[mV/cm]');
ylabel('W/2\pi[Hz]');
grid on;
legend('Interpreter', 'latex', 'Location', 'eastoutside'); % Legend only for data points
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
set(gcf, 'Units', 'centimeters', 'Position', [5, 5, 25, 8]);

% Print out the slopes in order
disp('Slopes of the linear fits:');
fprintf('%0.4g\n', slopes);


% Plot a_0 vs E_nr_values_numeric
figure;
hold on;
for i = 1:len
    plot(E_nr_values_numeric, a_0(i, :), '-o', 'LineWidth', 1.5, ...
         'Color', colors(i, :), 'MarkerFaceColor', colors(i, :));
end
hold off;
xlabel('E_{nr0}[mV/cm]');
ylabel('a_0');
grid on;
legend(data_labels, 'Interpreter', 'latex', 'Location', 'eastoutside');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
set(gcf, 'Units', 'centimeters', 'Position', [5, 5, 25, 8]); 


% Plot a_1 vs E_nr_values_numeric
figure;
hold on;
for i = 1:len
    plot(E_nr_values_numeric, a_1(i, :), '-o', 'LineWidth', 1.5, ...
         'Color', colors(i, :), 'MarkerFaceColor', colors(i, :));
end
hold off;
xlabel('E_{nr0}[mV/cm]');
ylabel('a_1[1/Hz]');
grid on;
legend(data_labels, 'Interpreter', 'latex', 'Location', 'eastoutside');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
set(gcf, 'Units', 'centimeters', 'Position', [5, 5, 25, 8]); 



%%% a_1 vs W
x = W(:);   % Flatten W into a column vector (x-axis)
y = a_1(:); % Flatten a_1 into a column vector (y-axis)
p = polyfit(x, y, 1); % Linear fit (1st order polynomial)

% Generate fit line data
x_fit = linspace(min(x), max(x), 100);
y_fit = polyval(p, x_fit);

figure;
plot(x_fit, y_fit, 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('Fit: y = (%.3fe-7)x', p(1)*10^7)); % Fit line
hold on;
scatter(x, y, 'ko', 'filled', 'HandleVisibility', 'off'); % Scatter plot with blue filled circles

% Label axes
xlabel('W/2\pi[Hz]');
ylabel('a_1[1/Hz]');
grid on;
legend('Location', 'Best');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
hold off;


%%% W vs E_nr*delta_L2
% Define the row and column vectors
row_values = [-120, -90, -60, -30, 0, 30, 60, 90, 120]; % 1×9
col_values = [-3; -2; -1; 0; 1; 2; 3]; % 9×1

product = col_values * row_values; 
product = product(:);
% x = a_1(:);         % Flatten a1 into a column vector (x-axis)
x = W(:);         % Flatten W into a column vector (x-axis)
p = polyfit(product, x, 1); % Linear fit (1st order polynomial)

x_fit = linspace(min(product), max(product), 100);
y_fit = polyval(p, x_fit);

figure;
plot(x_fit, y_fit, 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('Fit: y = (%.3fe-9)x', p(1)*10^9)); % Fit line
hold on;
scatter(product, x, 'ko', 'filled', 'HandleVisibility', 'off'); % Scatter plot with blue filled circles

% Label axes
% ylabel('W/2\pi[Hz]');
ylabel('a_1[1/Hz]');
xlabel('E_{nr0}*\delta_{L2}(mV/cm*MHz)');
grid on;
legend('Location', 'Best');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
hold off;

%% OPR vs. L2 Detuning

% MATLAB Script to Plot OPR vs. L2 Detuning
clc; clear; close all;

% Load the data from Excel
filename = 'Simulation Results.xlsx';
sheetname = 'OPR vs. L2 detuning';

% Read the data
data = readmatrix(filename, 'Sheet', sheetname);

% Extract L2 detuning values (column headers, skipping first column)
detuning = readmatrix(filename, 'Sheet', sheetname, 'Range', 'B1:GS1');

% Extract laser intensity values (first column)
intensities = data(:,1);

% Extract population values (remaining columns)
population = data(:, 2:end);

% Select six intensities to plot (adjust indices as needed)
intensity_indices = [1, 2, 3, 4, 5, 6]; % Modify for different selections
selected_intensities = intensities(1+intensity_indices);
selected_population = population(1+intensity_indices, :);

% Plot
figure;
hold on;
colors = lines(length(intensity_indices)); % Get different colors for each line

for i = 1:length(intensity_indices)
    plot(detuning, selected_population(i, :), 'LineWidth', 2, 'Color', colors(i, :), ...
         'DisplayName', sprintf(' %.1f of E_{L2}', selected_intensities(i)));
end

% Labels and Legend
xlabel('\delta_{L2}/2\pi (MHz)');
ylabel('OPR (%Population)');
title('OPR vs. \delta_{L2} at Different Laser Intensities');
legend('show');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
hold off;

%% Linear fit
% Given data
x = [-3, -2, -1, 0, 1, 2, 3]';
y = [0.2857,
0.1721,
0.03796,
3.9e-11,
-0.03796,
-0.1721,
-0.2857
]';

% Perform linear fitting
p = polyfit(x, y, 1); % Linear fit (1st order polynomial)

% Generate fit line data
x_fit = linspace(min(x), max(x), 100);
y_fit = polyval(p, x_fit);

% Plot the data
figure;
plot(x, y, 'ko', 'MarkerSize', 8, 'MarkerFaceColor','k','HandleVisibility', 'off'); % Data points
hold on;
plot(x_fit, y_fit, 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('Fit: y = %.4fx', p(1))); % Fit line

% Label axes and title
xlabel('\delta_{L2} [MHz]');
ylabel('Slope S = W/E_{nr}');
legend('Location', 'Best');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
grid on;

% Display fit parameters
disp(['Slope (m): ', num2str(p(1))]);
disp(['Intercept (b): ', num2str(p(2))]);


%% Position scan

data = [-28.26001733	-12.70119415	-5.277725032	-2.262556129	-0.991074141	-0.43526115	-0.18948946	-0.081251749	-0.034166979	-0.014032936	-0.005601183
-51.44000495	-44.07678226	-28.181824	-12.85161611	-5.352792105	-2.27128596	-0.976938659	-0.419176169	-0.177610072	-0.073868813	-0.030018535
-52.63617641	-51.62636238	-49.00673719	-41.86940822	-26.9181004	-12.36355984	-5.120582895	-2.135562743	-0.896296705	-0.373348467	-0.152888099
-48.45861377	-48.24930235	-47.77917964	-46.6955369	-44.06657863	-37.43434427	-24.09721149	-11.06988054	-4.524472211	-1.841610259	-0.74879458
-40.69689331	-40.64869242	-40.5426593	-40.31044535	-39.80344623	-38.69118169	-36.19479339	-30.4647862	-19.5570301	-8.919210393	-3.562482414
-29.6207341	-29.60983487	-29.58572866	-29.53283385	-29.41804948	-29.17250756	-28.65619579	-27.58988426	-25.4182511	-21.04829142	-13.40858271
-15.77477198	-15.77240673	-15.7671257	-15.75542399	-15.72977784	-15.67436644	-15.55685795	-15.31371215	-14.82716597	-13.89982281	-12.26421057
-0.153114304	-0.152622118	-0.151512741	-0.149030018	-0.143530101	-0.13150107	-0.105616439	-0.051062806	0.060841616	0.281924139	0.694258274
15.89890096	15.89899943	15.89922399	15.89973065	15.90086362	15.90336654	15.90881727	15.92047798	15.94486994	15.99445395	16.09138834
30.89399546	30.89401498	30.89405928	30.89415927	30.89438513	30.89488792	30.89599132	30.8983779	30.90343519	30.91390029	30.93489634
43.45222472	43.45222845	43.45223688	43.45225622	43.45230011	43.45239852	43.45261567	43.45308753	43.45409655	43.4562056	43.46049742
52.53246942	52.53246985	52.5324716	52.53247567	52.53248357	52.53250264	52.53254405	52.53263538	52.53283108	52.53324196	52.53408472
57.59797874	57.59797856	57.59797894	57.59797996	57.59798125	57.59798518	57.59799282	57.59801012	57.59804755	57.59812589	57.59828606
58.6836341	58.68363425	58.68363384	58.68363475	58.68363489	58.68363545	58.68363681	58.68363941	58.68364706	58.68366137	58.68369173
56.35405454	56.3540545	56.35405474	56.35405453	56.35405522	56.3540549	56.35405507	56.35405585	56.35405709	56.35405987	56.35406498
51.56468606	51.56468612	51.56468612	51.56468613	51.56468586	51.56468599	51.564686	51.56468616	51.56468673	51.56468648	51.56468814
45.45777593	45.45777528	45.45777513	45.45777561	45.45777565	45.45777538	45.45777605	45.45777558	45.45777552	45.45777559	45.45777582
39.13758373	39.13758368	39.13758369	39.13758356	39.13758362	39.13758381	39.13758375	39.13758357	39.13758376	39.13758408	39.1375837
36.18198783	36.18198738	36.18198784	36.18198807	36.18198808	36.18198804	36.18198755	36.18198796	36.18198747	36.18198779	36.18198779
31.05421636	31.05421664	31.05421636	31.0542167	31.05421634	31.0542167	31.0542165	31.05421652	31.05421647	31.05421639	31.05421616
27.16663301	27.16663343	27.16663296	27.16663339	27.16663326	27.16663349	27.16663352	27.16663324	27.16663304	27.16663358	27.16663309
24.40208232	24.40208196	24.40208235	24.40208218	24.40208206	24.40208156	24.40208226	24.40208181	24.40208228	24.40208185	24.40208179
22.32640078	22.3264004	22.3264009	22.32640083	22.32640055	22.32640061	22.32640074	22.32640029	22.32640027	22.32640048	22.32640085
20.33685242	20.33685305	20.33685299	20.33685313	20.33685297	20.33685254	20.33685288	20.33685293	20.33685303	20.33685296	20.33685298
17.83486163	17.8348613	17.83486144	17.83486159	17.83486192	17.83486149	17.83486192	17.83486126	17.83486168	17.83486176	17.83486163
14.38312165	14.38312141	14.38312104	14.38312185	14.38312161	14.38312177	14.38312142	14.38312135	14.38312169	14.38312148	14.38312152
9.813052377	9.813052843	9.813052586	9.813052253	9.813052391	9.813052996	9.813052248	9.813052756	9.813052783	9.813052833	9.813052725
4.261805249	4.261805931	4.26180515	4.261805194	4.261805448	4.261805309	4.261804596	4.261805156	4.261805435	4.261805398	4.261805417
-1.86409109	-1.864090867	-1.864091176	-1.864090699	-1.864091478	-1.864090678	-1.864091364	-1.864090679	-1.864091007	-1.864091023	-1.864091007];


% Define original x and y coordinate values (assuming they range from 1 to 16 for y and 1 to 21 for x)
y = linspace(1,11,11);
x = linspace(1,29,29);

% Create a finer grid for interpolation
[xq, yq] = meshgrid(linspace(min(x), max(x), 600), linspace(min(y), max(y), 100));

% Interpolate the data on the finer grid (note: data' is used to match dimensions)
data_interp = interp2(x, y, data', xq, yq, 'cubic'); % 'cubic' interpolation for smoothness

% Plot the interpolated data
figure;
imagesc(xq(1, :), yq(:, 1), data_interp); % Use interpolated x, y values
% clim([-54 0]);
colorbar;            % Add color scale
colormap jet;        % Choose colormap
% axis equal tight;    % Keep aspect ratio correct

% Label the axes
ylabel('L2 Position (\mus)');
xlabel('non-reversing field Position (\mus)');
title('0 represents the center position of Stark field');

% Set x and y axis ticks to match original coordinate values and convert to physical coordinates:
% x_phys = (x-6)*10, y_phys = (y-1)*10
xticks(x);
yticks(y);
xticklabels(string(30 + x*20));
yticklabels(string((y+4)*10));

% Ensure proper axis direction (by default, imagesc flips it)
set(gca, 'YDir', 'normal');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);

% Convert physical positions to plot coordinate system:
% x_line1 = (-43.7/10) + 6; % corresponds to x_phys = -43.7 µs
% x_line2 = (43.7/10) + 6;  % corresponds to x_phys = 43.7 µs
% y_line = (43.7/10) + 1;   % corresponds to y_phys = 43.7 µs
x_red = (70.1/10) -4;    % x coordinate in plot space
y_red = (52.6/10) -4;    % y coordinate in plot space

hold on;
% Get current axis limits
x_limits = xlim;
y_limits = ylim;

% % Draw vertical white dotted lines at the specified x positions
% line([x_line1, x_line1], y_limits, 'Color', 'w', 'LineStyle', '--', 'LineWidth', 1);
% line([x_line2, x_line2], y_limits, 'Color', 'w', 'LineStyle', '--', 'LineWidth', 1);

% Draw a horizontal white dotted line at the specified y position
% line(x_limits, [y_line, y_line], 'Color', 'w', 'LineStyle', '--', 'LineWidth', 1);

% Add a red dot at the specified location
plot(x_red, y_red, 'wo', 'MarkerSize', 8, 'MarkerFaceColor', 'w');


%% two data scatter plot
% W vs. L2 & NR field location
x = [60,70,80,90,100]';
y1 = [-0.32979, -0.33079,  -0.3295, -0.32516, -0.31704]';
y2 = [-0.32962, -0.33007,-0.32934, -0.32583, -0.3174]';

% Plot the data
figure;
plot(x, y1, 'ko', 'MarkerSize', 8, 'LineWidth', 2); % Circle markers ('o') with thicker edge
hold on;
plot(x, y2, 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor','k'); % Cross markers ('x') with thicker lines

% Label axes and title
xlabel('Position [\mus]');
ylabel('W/2\pi [Hz]');
legend('Original NR field', 'flat-top NR field');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
xlim([55, 105]);    % Adjust x-axis limits as desired
ylim([-0.335, -0.315]);  % Adjust y-axis limits as desired

grid off;

%% analytical W 

%% Constants and Parameters

% Given constants
E0 = 40;                        % Stark Field Strength [V/m]
E_nr0_values = [-12, -9, -6, -3, 0, 3, 6, 9, 12];  % NR Field Strength values [V/m]
E_L2 = 851.4 * 0.11;            % L2 Field Strength [V/m]
d12 = -33.6 * 2*pi;             % Dipole matrix element d_{12} [rad/(s·V/m)]
d13 = -2.15e4 * 2*pi;           % Dipole matrix element d_{13} [rad/(s·V/m)]
Gamma = 2.7e6 * 2*pi;           % Decay rate [rad/s]
omega = 11.44e3 * 2*pi;         % Stark field angular frequency [rad/s]
v = 616;                      % Beam velocity [m/s]
sigma_u = 7.6e-3;             % NR field width [m]
sigma = 0.52e-6;              % L2 width [s]
Te = 87.4e-6;                 % Interaction time [s]
Tf2 = 8.9e-6;                 % Free evolution time [s]
t0 = 70.1e-6 - 52.6e-6;

% Calculate auxiliary constants
K3 = - (d13^2 * E_L2^2 * sqrt(pi) * sigma) / 4;
a  = v^2 / (2 * sqrt(2) * sigma_u^2);

% Create vector of Delta2 values:
% Delta2 is given in MHz in the labels: -3 MHz, -2 MHz, ..., +3 MHz.
% To convert MHz to rad/s: 1 MHz = 1e6 Hz, then multiply by 2*pi.
Delta2_values = (-3:1:3) * 1e6 * 2*pi;  % 7 values [rad/s]

% Preallocate W matrix:
% Rows: length(Delta2_values) (7), Columns: length(E_nr0_values) (9)
W = zeros(length(Delta2_values), length(E_nr0_values));
a0_mat = zeros(length(Delta2_values), length(E_nr0_values));
a1_mat = zeros(length(Delta2_values), length(E_nr0_values));

for i = 1:length(Delta2_values)
    currDelta2 = Delta2_values(i);
    for j = 1:length(E_nr0_values)
        currE_nr0 = E_nr0_values(j);
        % Compute W using the given formula:
        W(i, j) = (1/(2*pi))*(d12 * currE_nr0) / (Te + Tf2) * sqrt(pi/a) * ...
                  exp( (K3 * Gamma/2) / (currDelta2^2 + (Gamma/2)^2) ) * ...
                  sin( (K3 * currDelta2) / (currDelta2^2 + (Gamma/2)^2) );
        a0_mat(i, j) = (-omega * currE_nr0 * (2*t0 + Te + 2*Tf2) / (E0 * Te)) * sqrt(pi/a) * ...
                        exp((K3 * Gamma/2) / (currDelta2^2 + (Gamma/2)^2)) * ...
                        cos((K3 * currDelta2) / (currDelta2^2 + (Gamma/2)^2));
        a1_mat(i, j) = (-omega * currE_nr0 / E0) * sqrt(pi/a) * ...
                        ( (1/(a*Te) - Te/12) + ((2*t0 + Te + 2*Tf2)^2)/(4*Te) ) * ...
                        exp((K3 * Gamma/2) / (currDelta2^2 + (Gamma/2)^2)) * ...
                        sin((K3 * currDelta2) / (currDelta2^2 + (Gamma/2)^2));
    end
end


% Convert E_nr0 values to mV/cm for plotting
% 1 V/m = 0.1 mV/cm.
E_nr0_numeric = E_nr0_values * 10;  % now in mV/cm

% Create data labels for each Delta2 value
data_labels = {'$\delta_{L2} = -3$ MHz', '$\delta_{L2} = -2$ MHz', ...
               '$\delta_{L2} = -1$ MHz', '$\delta_{L2} = 0$ MHz', ...
               '$\delta_{L2} = +1$ MHz', '$\delta_{L2} = +2$ MHz', ...
               '$\delta_{L2} = +3$ MHz'};

% Use distinct colors for each curve
len = length(data_labels);
colors = jet(len); 

figure;
hold on;
for i = 1:length(Delta2_values)
    % Plot the i-th row of W versus E_nr0_numeric using markers 'o'
    plot(E_nr0_numeric, W(i, :), 'o:', 'LineWidth', 1.5, ...
         'Color', colors(i, :), 'MarkerFaceColor', colors(i, :));
end
hold off;
xlabel('E_{nr0} [mV/cm]');
ylabel('W/2\pi [Hz]');
grid on;
legend(data_labels, 'Interpreter', 'latex', 'Location', 'eastoutside');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
set(gcf, 'Units', 'centimeters', 'Position', [5, 5, 25, 8]);

% Plot a0: 
figure;
hold on;
for i = 1:length(Delta2_values)
    % Plot the i-th row of a0_mat versus E_nr0_numeric with markers and dotted lines
    plot(E_nr0_numeric, a0_mat(i, :), 'o:', 'LineWidth', 1.5, ...
         'Color', colors(i, :), 'MarkerFaceColor', colors(i, :));
end
hold off;
xlabel('E_{nr0} [mV/cm]');
ylabel('a_0');
grid on;
legend(data_labels, 'Interpreter', 'latex', 'Location', 'eastoutside');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
set(gcf, 'Units', 'centimeters', 'Position', [5, 5, 25, 8]);

% Plot a1:
figure;
hold on;
for i = 1:length(Delta2_values)
    % Plot the i-th row of a1_mat versus E_nr0_numeric with markers and dotted lines
    plot(E_nr0_numeric, a1_mat(i, :), 'o:', 'LineWidth', 1.5, ...
         'Color', colors(i, :), 'MarkerFaceColor', colors(i, :));
end
hold off;
xlabel('E_{nr0} [mV/cm]');
ylabel('a_1 [1/Hz]');
grid on;
legend(data_labels, 'Interpreter', 'latex', 'Location', 'eastoutside');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
set(gcf, 'Units', 'centimeters', 'Position', [5, 5, 25, 8]);

