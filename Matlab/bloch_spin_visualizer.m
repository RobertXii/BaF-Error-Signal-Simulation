function bloch_spin_visualizer()
% Visualize spin-1/2 dynamics on the Bloch sphere under
% H(t) = [ 0, d E(t) + i W ; d E(t) - i W, Delta ],
% E(t) = E0 * sin(omega * t).
% Also shows a second trajectory for reversed field E_rev(t) = -E0 * sin(omega * t).
%
% Bloch decomposition (ħ=1 here): H = (1/2) Ω·σ + (Δ/2) I
% with Ω = ( 2 d E,  -2 W,  -Δ ).

%% ========== User-adjustable parameters ==========
E0     = 40;
omega  = 2*pi*11.44e3;
d      = 2*pi*-33.6;
W      = 2*pi*5;
Delta  = 2*pi*3.5e3;

t0     = 0;                 % start time [s]
t1     = 87.4e-6;           % end time [s]
Nt     = 1000;              % number of samples

psi0   = [0; 1];            % initial state (spin down here)
normalize_initial = true;

%% ========== Derived functions ==========
% E_of_t      = @(t)  E0 * sin(omega*t)^2;     % original field
% E_of_t_rev  = @(t) -E0 * sin(omega*t)^2;     % reversed field
% 
% E_of_t      = @(t)  E0 ;     % original field
% E_of_t_rev  = @(t) -E0 ;     % reversed field

% T = 43.7e-6; % pulse half-duration
% E_of_t = @(t) ( (t > 0) & (t <= T) )   *  E0  + ...
%               ( (t >= T) & (t <= 2*T) )  * (-E0);
% E_of_t_rev = @(t) ( (t > 0) & (t <= T) )   * (-E0) + ...
%                   ( (t >= T) & (t <= 2*T) )  *  E0;

Ttot  = 87.4e-6;           % total duration
tc    = Ttot/2;            % center of the window
sep   = Ttot/4;            % +/- quarter-window separation (one sine period over Ttot)
sigma = Ttot/8;            % width; increase to make it more sinusoidal, decrease for sharper lobes
% 
% % Optional hard window so the field is zero outside [0, Ttot]
win   = @(t) (t >= 0) & (t <= Ttot);

% Two-sech “sin-like” field: +lobe then −lobe, centered at t = Ttot/2
E_of_t = @(t) win(t) .* ( ...
            E0 * sech( (t - (tc)) / sigma ) ...
          - E0 * sech( (t - (tc + sep)) / sigma )*0 );

% Reversed version (flip the sign of both lobes)
E_of_t_rev = @(t) win(t) .* ( ...
             -E0 * sech( (t - (tc)) / sigma ) ...
             +E0 * sech( (t - (tc + sep)) / sigma )*0 );

% Hamiltonians (Hermitian)
H_of_t = @(t) [ 0,               d*E_of_t(t)     + 1i*W; ...
                d*E_of_t(t)     - 1i*W,                 Delta ];

H_of_t_rev = @(t) [ 0,           d*E_of_t_rev(t) + 1i*W; ...
                    d*E_of_t_rev(t) - 1i*W,             Delta ];

% Ω(t) (ħ = 1 here)
Omega_of_t     = @(t) [ 2*d*E_of_t(t);      -2*W; -Delta ];
Omega_of_t_rev = @(t) [ 2*d*E_of_t_rev(t);  -2*W; -Delta ];

% Schrödinger RHS: dpsi/dt = -i H(t) psi   (ħ = 1)
schrodinger     = @(t, psi) (-1i) * (H_of_t(t)     * psi);
schrodinger_rev = @(t, psi) (-1i) * (H_of_t_rev(t) * psi);

% Bloch vector from state psi = [a;b]
bloch_vec = @(psi) [ 2*real(conj(psi(1))*psi(2)); ...
                     2*imag(conj(psi(1))*psi(2)); ...
                     abs(psi(1))^2 - abs(psi(2))^2 ];

if normalize_initial, psi0 = psi0 / norm(psi0); end

%% ========== Integrate dynamics (both cases) ==========
tgrid = linspace(t0, t1, Nt).';
opts  = odeset('RelTol',1e-8,'AbsTol',1e-10);

[ts, psis1] = ode45(@(t,psi) schrodinger(t,psi),      tgrid, psi0, opts);
[~,  psis2] = ode45(@(t,psi) schrodinger_rev(t,psi),  tgrid, psi0, opts);

psis1 = psis1.';  % 2 x Nt
psis2 = psis2.';

% Bloch trajectories
S1 = zeros(3, Nt);
S2 = zeros(3, Nt);
for k = 1:Nt
    S1(:,k) = bloch_vec(psis1(:,k));
    S2(:,k) = bloch_vec(psis2(:,k));
end

% Ω(t) for both and common scaling
Omega1 = zeros(3, Nt);
Omega2 = zeros(3, Nt);
for k = 1:Nt
    Omega1(:,k) = Omega_of_t(ts(k));
    Omega2(:,k) = Omega_of_t_rev(ts(k));
end
Omega_mag_max = max([vecnorm(Omega1,2,1), vecnorm(Omega2,2,1)]) + eps;
Omega1_scaled = Omega1 ./ Omega_mag_max;
Omega2_scaled = Omega2 ./ Omega_mag_max;

%% ========== Build UI ==========
fig = figure('Name','Bloch Sphere Spin Dynamics (Original vs Reversed Field)',...
             'Color','w','Position',[100 100 1200 560]);

tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');

% ---- Left: Bloch sphere ----
nexttile(1);
ax3d = gca; hold(ax3d,'on'); axis(ax3d,'equal');
xlabel(ax3d,'x'); ylabel(ax3d,'y'); zlabel(ax3d,'z');
title(ax3d,'Bloch Sphere: s(t) & \Omega(t) (both cases, scaled equally)');

% Gray sphere (thicker dotted trajectories)
[NX,NY,NZ] = sphere(60);
surf(ax3d, NX, NY, NZ, 'FaceAlpha',0.1, 'EdgeAlpha',0.08, 'FaceColor',[0.7 0.7 0.7]);
plot3(ax3d, [-1 1],[0 0],[0 0],'k-','LineWidth',0.6);
plot3(ax3d, [0 0],[-1 1],[0 0],'k-','LineWidth',0.6);
plot3(ax3d, [0 0],[0 0],[-1 1],'k-','LineWidth',0.6);
grid(ax3d,'on'); view(ax3d, 35, 20); xlim(ax3d,[-1 1]); ylim(ax3d,[-1 1]); zlim(ax3d,[-1 1]);

% Trajectories (thicker dotted)
traj1 = plot3(ax3d, S1(1,:), S1(2,:), S1(3,:), ':', 'LineWidth', 2.2, 'Color', [0.3 0.3 0.3]);
traj2 = plot3(ax3d, S2(1,:), S2(2,:), S2(3,:), ':', 'LineWidth', 2.2, 'Color', [0.1 0.4 0.9]);

% State markers at current time
state1 = plot3(ax3d, S1(1,1), S1(2,1), S1(3,1), 'o', 'MarkerSize', 8, 'MarkerFaceColor',[0.85 0.2 0.2], 'MarkerEdgeColor','k');
state2 = plot3(ax3d, S2(1,1), S2(2,1), S2(3,1), 'o', 'MarkerSize', 8, 'MarkerFaceColor',[0.2 0.4 0.9], 'MarkerEdgeColor','k');

% Two Omega arrows (scaled by common max)
omg1 = quiver3(ax3d, 0,0,0, Omega1_scaled(1,1), Omega1_scaled(2,1), Omega1_scaled(3,1), ...
               'LineWidth',2.0,'MaxHeadSize',0.55,'Color',[0.85 0.2 0.2]);
omg2 = quiver3(ax3d, 0,0,0, Omega2_scaled(1,1), Omega2_scaled(2,1), Omega2_scaled(3,1), ...
               'LineWidth',2.0,'MaxHeadSize',0.55,'Color',[0.1 0.4 0.9]);

legend(ax3d, [traj1 traj2 state1 state2 omg1 omg2], ...
       {'traj: E(t)=+E_0 sin','traj: E(t)=-E_0 sin','state (+)', 'state (-)', '\Omega (+)', '\Omega (-)'}, ...
       'Location','southoutside','NumColumns',3);

% Readout
omegaText = text(ax3d, 0.02, 0.02, 1.04, '', 'Units','normalized', 'FontSize',10, ...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom', 'Color',[0.1 0.1 0.1]);

% ---- Right: Populations vs time ----
nexttile(2);
axPop = gca; hold(axPop,'on'); box(axPop,'on');
P1_up   = abs(psis1(1,:)).^2;  P1_down = abs(psis1(2,:)).^2;
P2_up   = abs(psis2(1,:)).^2;  P2_down = abs(psis2(2,:)).^2;

p1 = plot(axPop, ts, P1_up,   'LineWidth',1.7, 'DisplayName','P_{\uparrow} (+)');
p2 = plot(axPop, ts, P1_down, 'LineWidth',1.7, 'DisplayName','P_{\downarrow} (+)');
p3 = plot(axPop, ts, P2_up,   '--', 'LineWidth',1.8, 'DisplayName','P_{\uparrow} (-)');
p4 = plot(axPop, ts, P2_down, '--', 'LineWidth',1.8, 'DisplayName','P_{\downarrow} (-)');

tMarker = plot(axPop, ts(1)*[1 1], [0 1], ':', 'Color',[0.2 0.2 0.2], 'DisplayName','t');
legend(axPop,'Location','southoutside','Orientation','horizontal');
xlabel(axPop,'t [s]'); ylabel(axPop,'Population'); ylim(axPop,[-0.02 1.02]);
title(axPop, 'Populations vs Time (solid: +E, dashed: -E)');

% ---- Slider for time ----
uic = uicontrol('Style','slider', 'Units','normalized', ...
    'Position',[0.10 0.02 0.80 0.04], 'Min', ts(1), 'Max', ts(end), 'Value', ts(1));
uic.SliderStep = [1/(Nt-1), 10/(Nt-1)];

% Time label
tLabel = uicontrol('Style','text','Units','normalized','BackgroundColor','w', ...
    'Position',[0.91 0.02 0.08 0.04],'String',sprintf('t = %.3e s', ts(1)));

% Callback and initial update
uic.Callback = @(src,evt) updateFrame(src.Value);
updateFrame(ts(1));

%% ========== Nested updater ==========
    function updateFrame(tvalue)
        [~, idx] = min(abs(ts - tvalue));

        % Update state markers
        set(state1,'XData',S1(1,idx),'YData',S1(2,idx),'ZData',S1(3,idx));
        set(state2,'XData',S2(1,idx),'YData',S2(2,idx),'ZData',S2(3,idx));

        % Update Omega arrows
        set(omg1,'UData',Omega1_scaled(1,idx),'VData',Omega1_scaled(2,idx),'WData',Omega1_scaled(3,idx));
        set(omg2,'UData',Omega2_scaled(1,idx),'VData',Omega2_scaled(2,idx),'WData',Omega2_scaled(3,idx));

        % Update time cursor
        set(tMarker,'XData', ts(idx)*[1 1]);

        % Readout (true |Ω| values shown; scaled on sphere)
        % Om1   = Omega1(:,idx);   Om2   = Omega2(:,idx);
        % Om1Hz = Om1/(2*pi);      Om2Hz = Om2/(2*pi);
        % Om1Mag = norm(Om1);      Om2Mag = norm(Om2);
        % txt = sprintf(['\\Omega_{+}(t)=[%.2e, %.2e, %.2e] rad/s  (|\\Omega_{+}|=%.2e rad/s = %.2e Hz)\n' ...
        %                '\\Omega_{-}(t)=[%.2e, %.2e, %.2e] rad/s  (|\\Omega_{-}|=%.2e rad/s = %.2e Hz)\n' ...
        %                'scale on sphere = 1 / max_t |\\Omega| = 1 / %.2e'], ...
        %                Om1(1), Om1(2), Om1(3), Om1Mag, Om1Mag/(2*pi), ...
        %                Om2(1), Om2(2), Om2(3), Om2Mag, Om2Mag/(2*pi), ...
        %                Omega_mag_max);
        % set(omegaText,'String',txt);

        % Time label
        tLabel.String = sprintf('t = %.3e s', ts(idx));

        drawnow limitrate;
    end
end