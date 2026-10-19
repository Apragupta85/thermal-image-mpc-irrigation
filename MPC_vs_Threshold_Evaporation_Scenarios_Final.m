%% MPC vs fixed-dose threshold: scenarios 2-5 only
% Purpose: test whether forecast-aware MPC adapts irrigation to changing
% evaporation demand while respecting the moisture bounds.
% Simulation: 48 h; irrigation only at 06:00, 07:00, 08:00 and 09:00.
% IMPORTANT: weather is synthetic, not live forecast data. A shared,
% intermittent synthetic rainfall sequence is retained for both controllers.
%
% Requires MATLAB Optimization Toolbox (quadprog).
% Run from a writable folder. Outputs CSV, editable MATLAB .fig files,
% editable vector SVG files, vector PDFs, and 600-dpi PNGs.

clear; clc; close all;
rng(1,'twister');

%% Parameters
B = 0.009;                 % irrigation effectiveness
xref = 0.30;               % target normalised soil moisture
xmin = 0.15;               % minimum moisture bound
xmax = 0.40;               % maximum moisture bound
umax = 10;                 % maximum hourly irrigation input
N = 24;                    % 24-hour prediction horizon
nSteps = 48;               % 48-hour independent simulation
x0 = 0.30;                 % same initial moisture for all cases
threshold = xref;          % fixed threshold trigger
fixedDose = 8;             % fixed dose per trigger, all scenarios
kEScenarios = [0.18 0.15 0.12 0.09];
waterPenalty = 1e-3;        % small penalty on input to discourage excess use

% Synthetic hourly ambient temperature and relative humidity (same in all
% scenarios). Rainfall is intermittent rather than added every hour, so
% moisture rises outside irrigation hours only when a rainfall event occurs.
% This sequence is identical for MPC and threshold in every scenario.
hourOfDay = (0:nSteps-1)';
TA = 20 + 5*sin(2*pi*(hourOfDay-6)/24);       % deg C
RH = 60 - 20*cos(2*pi*(hourOfDay-6)/24);      % percent
rain = zeros(nSteps,1);
rainEventHours = [13 14 37 38]; % zero-based hours; two short rain events
for h = rainEventHours
    rain(h+1) = 0.02 + 0.01*rand;
end
Es = 0.6108 .* exp((17.27 .* TA) ./ (TA + 237.3));

nScen = numel(kEScenarios);
varNames = {'kE','MaxHourlyEvaporation','MPC_TotalInput','Threshold_TotalInput', ...
    'MPC_RMSE','Threshold_RMSE','MPC_MinMoisture','Threshold_MinMoisture', ...
    'MPC_LowerViolations','Threshold_LowerViolations','MPC_UpperViolations', ...
    'Threshold_UpperViolations','MPC_QP_Failures','MPC_FallbackSteps'};
Results = array2table(nan(nScen,numel(varNames)),'VariableNames',varNames);
X_MPC = zeros(nSteps+1,nScen); X_TH = zeros(nSteps+1,nScen);
U_MPC = zeros(nSteps,nScen); U_TH = zeros(nSteps,nScen);
E_all = zeros(nSteps,nScen);

qpOpts = optimoptions('quadprog','Display','off');

for s = 1:nScen
    kE = kEScenarios(s);
    E = kE .* (1 - RH./100) .* Es;
    % Keep the same formula as the manuscript. These values are treated as
    % hourly loss fractions in the state equation; inspect magnitudes before
    % making paper claims because the formula's units/scaling are critical.
    E_all(:,s) = E;

    xM = zeros(nSteps+1,1); xM(1)=x0;
    xT = zeros(nSteps+1,1); xT(1)=x0;
    uM = zeros(nSteps,1); uT = zeros(nSteps,1);
    qpFailures = 0; fallbackSteps = 0;

    for t = 1:nSteps
        % The controller receives the current state. There is no added
        % measurement noise in this controlled comparison.
        xNow = xM(t);
        xTnow = xT(t);

        % Build 24-hour prediction model from the current time. Every future
        % input is a decision variable, but inputs outside 06:00-10:00 are
        % constrained to zero using the absolute clock hour.
        idx = mod((t-1)+(0:N-1),nSteps)+1;
        Ep = E(idx); rp = rain(idx);
        S = zeros(N,1); G = zeros(N,N); D = zeros(N,1);
        for k = 1:N
            S(k) = prod(1-Ep(1:k));
            for j = 1:k
                if j == k, decay = 1; else, decay = prod(1-Ep(j+1:k)); end
                D(k) = D(k) + rp(j)*decay;
                G(k,j) = B*decay;
            end
        end

        % Only future intervals beginning 06:00, 07:00, 08:00 or 09:00
        % can receive irrigation. This enforces the 06:00-10:00 window.
        allowedFuture = false(N,1);
        for j = 1:N
            absoluteHour = mod((t-1)+(j-1),24);
            allowedFuture(j) = absoluteHour >= 6 && absoluteHour < 10;
        end
        lb = zeros(N,1); ub = zeros(N,1); ub(allowedFuture)=umax;

        % Objective: track xref over the full horizon, with a small input
        % penalty. Hard state bounds are enforced in the QP.
        Q = eye(N);
        H = 2*(G'*Q*G + waterPenalty*eye(N));
        f = 2*G'*Q*(S*xNow + D - xref);
        A = [G; -G];
        b = [xmax*ones(N,1) - S*xNow - D; ...
             -xmin*ones(N,1) + S*xNow + D];

        [uPlan,~,exitflag] = quadprog(H,f,A,b,[],[],lb,ub,[],qpOpts);
        hr = mod(t-1,24);
        if exitflag > 0 && ~isempty(uPlan)
            uM(t) = uPlan(1);
        else
            % Do not silently label solver failure as a water-saving decision.
            qpFailures = qpFailures + 1;
            uM(t) = 0;
            fallbackSteps = fallbackSteps + 1;
        end

        % Minimal safety correction: intervene only if zero irrigation would
        % make the next moisture state fall below xmin. Do not force moisture
        % back to xref; that would create artificial over-irrigation.
        noWaterNext = (1-E(t))*xNow + rain(t);
        if hr >= 6 && hr < 10 && uM(t) <= 1e-8 && noWaterNext < xmin
            need = (xmin - noWaterNext)/B;
            uM(t) = max(0,min(umax,need));
        end

        % Fixed-dose threshold controller: no forecast and no dose adjustment.
        hr = mod(t-1,24);
        if hr >= 6 && hr < 10 && xTnow < threshold
            uT(t) = min(fixedDose,umax);
        end

        % Same evaporation and rainfall for both controllers.
        xM(t+1) = (1-E(t))*xM(t) + B*uM(t) + rain(t);
        xT(t+1) = (1-E(t))*xT(t) + B*uT(t) + rain(t);
    end

    X_MPC(:,s)=xM; X_TH(:,s)=xT; U_MPC(:,s)=uM; U_TH(:,s)=uT;
    Results.kE(s)=kE;
    Results.MaxHourlyEvaporation(s)=max(E);
    Results.MPC_TotalInput(s)=sum(uM);
    Results.Threshold_TotalInput(s)=sum(uT);
    Results.MPC_RMSE(s)=sqrt(mean((xM(1:nSteps)-xref).^2));
    Results.Threshold_RMSE(s)=sqrt(mean((xT(1:nSteps)-xref).^2));
    Results.MPC_MinMoisture(s)=min(xM);
    Results.Threshold_MinMoisture(s)=min(xT);
    Results.MPC_LowerViolations(s)=sum(xM < xmin);
    Results.Threshold_LowerViolations(s)=sum(xT < xmin);
    Results.MPC_UpperViolations(s)=sum(xM > xmax);
    Results.Threshold_UpperViolations(s)=sum(xT > xmax);
    Results.MPC_QP_Failures(s)=qpFailures;
    Results.MPC_FallbackSteps(s)=fallbackSteps;
end

% Save numerical outputs and explicit evaporation audit.
writetable(Results,'MPC_vs_Threshold_Evaporation_Results.csv');
EvaporationAudit = table((0:nSteps-1)',TA,RH,Es,rain, ...
    'VariableNames',{'Hour','Temperature_C','RelativeHumidity_pct','SaturationVapourPressure_kPa','Rainfall'});
writetable(EvaporationAudit,'MPC_Evaporation_Weather_Audit.csv');
disp('Scenario comparison results:'); disp(Results);
disp('Shared synthetic rainfall events (zero-based hours):'); disp(rainEventHours);
disp('Evaporation range by scenario (hourly loss fraction):');
disp(table(kEScenarios(:),max(E_all,[],1)',mean(E_all,1)', ...
    'VariableNames',{'kE','MaxHourlyE','MeanHourlyE'}));

%% Publication-quality plotting defaults
set(groot,'defaultFigureColor','w','defaultAxesFontName','Arial', ...
    'defaultAxesFontSize',12,'defaultAxesFontWeight','bold', ...
    'defaultAxesLineWidth',1.4,'defaultLineLineWidth',2.2);

% Figure 1: four vertically stacked moisture panels for original scenarios 2-5.
fig2=figure('Name','Editable - soil moisture trajectories','Color','w', ...
    'Units','inches','Position',[1 1 12.5 9.0]);
tiled=tiledlayout(fig2,6,1,'TileSpacing','compact','Padding','compact');
tt=0:nSteps;
for s=1:nScen
    ax=nexttile(tiled);
    plot(ax,tt,X_MPC(:,s),'-','LineWidth',2.0); hold(ax,'on');
    plot(ax,tt,X_TH(:,s),'--','LineWidth',2.0);
    yline(ax,xmax,'--','Maximum bound','LineWidth',1.2,'FontSize',9,'LabelHorizontalAlignment','left');
    yline(ax,xref,':','Reference','LineWidth',1.2,'FontSize',9,'LabelHorizontalAlignment','left');
    yline(ax,xmin,'--','Minimum bound','LineWidth',1.2,'FontSize',9,'LabelHorizontalAlignment','left');
    xlim(ax,[0 nSteps]); ylim(ax,[0 0.45]);
    set(ax,'FontSize',10,'FontWeight','bold','LineWidth',1.2,'Box','on','TickDir','out');
    title(ax,sprintf('Scenario %d: k_E = %.2f',s+1,kEScenarios(s)), ...
        'FontSize',11,'FontWeight','bold');
    if s==nScen, xlabel(ax,'Time (hours)','FontSize',11,'FontWeight','bold'); end
    ylabel(ax,'Moisture','FontSize',10,'FontWeight','bold');
    grid(ax,'on'); ax.GridAlpha=0.18;
    if s==1
        legend(ax,{'MPC','Fixed-dose threshold'},'Location','southoutside', ...
            'Orientation','horizontal','FontSize',10,'Box','off');
    end
end
title(tiled,'Soil-moisture response under evaporation scenarios 2-5','FontSize',14,'FontWeight','bold');
savefig(fig2,'MPC_vs_Threshold_Evaporation_Moisture.fig');
exportgraphics(fig2,'MPC_vs_Threshold_Evaporation_Moisture.svg','ContentType','vector');
exportgraphics(fig2,'MPC_vs_Threshold_Evaporation_Moisture.pdf','ContentType','vector');
exportgraphics(fig2,'MPC_vs_Threshold_Evaporation_Moisture.png','Resolution',600);

% Figure 2: grouped cumulative-input bars, matching the original figure.
fig1=figure('Name','Editable - irrigation input comparison','Color','w', ...
    'Units','inches','Position',[1 1 12.0 6.0]);
ax1=axes(fig1);
bar(ax1,[Results.MPC_TotalInput Results.Threshold_TotalInput],0.78);
set(ax1,'XTick',1:nScen,'XTickLabel',compose('%.2f',kEScenarios), ...
    'FontSize',12,'FontWeight','bold','LineWidth',1.5,'Box','on','TickDir','out');
grid(ax1,'on'); ax1.GridAlpha=0.18;
xlabel(ax1,'Evaporation coefficient, k_E','FontSize',13,'FontWeight','bold');
ylabel(ax1,'Cumulative irrigation input (model units)','FontSize',13,'FontWeight','bold');
legend(ax1,{'MPC','Fixed-dose threshold'},'Location','northwest','FontSize',11,'Box','on');
title(ax1,sprintf('Water input under evaporation scenarios 2-5 (fixed dose = %g)',fixedDose), ...
    'FontSize',13,'FontWeight','bold');
savefig(fig1,'MPC_vs_Threshold_Evaporation_Input.fig');
exportgraphics(fig1,'MPC_vs_Threshold_Evaporation_Input.svg','ContentType','vector');
exportgraphics(fig1,'MPC_vs_Threshold_Evaporation_Input.pdf','ContentType','vector');
exportgraphics(fig1,'MPC_vs_Threshold_Evaporation_Input.png','Resolution',600);

fprintf('\nFixed threshold dose = %.2f; irrigation window = 06:00-10:00.\n',fixedDose);
fprintf('MPC QP failures and fallback steps are explicitly reported per scenario.\n');
fprintf('Synthetic rainfall events occur at hours: %s (same sequence in every scenario).\n', mat2str(rainEventHours));
fprintf('Before paper use, verify hourly evaporation scaling and check all bound violations.\n');
