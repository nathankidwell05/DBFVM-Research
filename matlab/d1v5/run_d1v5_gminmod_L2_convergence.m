function convergenceTable = run_d1v5_gminmod_L2_convergence
%% L2 Grid-Convergence Test for the KT-D1V5 Generalized-Minmod Solver %%
% This function repeatedly runs the existing generalized-minmod Sod solver.
% Only Nx and therefore dx change between cases.
% The physical model, tau, CFL, tEnd, theta, and initial conditions stay fixed.
%
% Run this file from the MATLAB Command Window with:
%     convergenceTable = run_d1v5_gminmod_L2_convergence;

    clc;                                               % make the convergence output easy to find
    close all;                                         % close old solver figures

    %% Convergence-Test Settings
    NxValues = [6250;12500;25000;50000];              % each grid doubles the previous resolution
    thetaValue = 1.20;                                 % final generalized-minmod setting
    tauValue = 5.0e-5;                                 % same BGK relaxation time in every run
    numberOfRuns = numel(NxValues);                    % number of resolutions in the study

    thisFolder = fileparts(mfilename('fullpath'));     % folder containing this convergence function
    solverFile = fullfile(thisFolder,'d1v5_sodshock_gminmod_rk3.m'); % current generalized-minmod solver

    if ~isfile(solverFile)                             % stop clearly if the working solver was moved
        error('Generalized-minmod solver was not found at: %s',solverFile);
    end

    %% Storage for Errors and Run Information
    dxValues = zeros(numberOfRuns,1);                  % grid spacing for every resolution
    runTime = zeros(numberOfRuns,1);                   % wall-clock time for every resolution

    L2rho = zeros(numberOfRuns,1);                     % global density L2 error
    L2u = zeros(numberOfRuns,1);                       % global velocity L2 error
    L2p = zeros(numberOfRuns,1);                       % global pressure L2 error
    L2T = zeros(numberOfRuns,1);                       % global temperature L2 error

    waveL2rho = zeros(numberOfRuns,1);                 % density L2 error around the waves
    waveL2u = zeros(numberOfRuns,1);                   % velocity L2 error around the waves
    waveL2p = zeros(numberOfRuns,1);                   % pressure L2 error around the waves
    waveL2T = zeros(numberOfRuns,1);                   % temperature L2 error around the waves
    TpeakExcess = zeros(numberOfRuns,1);               % numerical peak T minus exact peak T

    %% Run the Same Solver on Successively Finer Grids
    fprintf('KT-D1V5 generalized-minmod L2 convergence test\n');
    fprintf('theta = %.2f, tau = %.3e\n',thetaValue,tauValue);

    for runNumber = 1:numberOfRuns
        requestedNx = NxValues(runNumber);             % grid resolution for this convergence run
        fprintf('\n============================================================\n');
        fprintf('Resolution %d/%d: Nx = %d\n',runNumber,numberOfRuns,requestedNx);

        metrics = run_one_resolution(requestedNx,thetaValue,tauValue,solverFile); % run working solver once

        dxValues(runNumber) = metrics.dx;              % store actual grid spacing
        runTime(runNumber) = metrics.runTime;          % store solver runtime
        L2rho(runNumber) = metrics.L2rho;              % store global density error
        L2u(runNumber) = metrics.L2u;                  % store global velocity error
        L2p(runNumber) = metrics.L2p;                  % store global pressure error
        L2T(runNumber) = metrics.L2T;                  % store global temperature error
        waveL2rho(runNumber) = metrics.waveL2rho;      % store wave-region density error
        waveL2u(runNumber) = metrics.waveL2u;          % store wave-region velocity error
        waveL2p(runNumber) = metrics.waveL2p;          % store wave-region pressure error
        waveL2T(runNumber) = metrics.waveL2T;          % store wave-region temperature error
        TpeakExcess(runNumber) = metrics.TpeakExcess;  % store temperature peak excess

        fprintf('Global L2: rho=%.6e, u=%.6e, p=%.6e, T=%.6e\n', ...
            metrics.L2rho,metrics.L2u,metrics.L2p,metrics.L2T);
        fprintf('Wave L2:   rho=%.6e, u=%.6e, p=%.6e, T=%.6e\n', ...
            metrics.waveL2rho,metrics.waveL2u,metrics.waveL2p,metrics.waveL2T);
        fprintf('T peak excess = %.6e\n',metrics.TpeakExcess);
    end

    %% Compute Observed Convergence Orders
    orderRho = observed_order(L2rho,dxValues);         % order from successive global density errors
    orderU = observed_order(L2u,dxValues);             % order from successive global velocity errors
    orderP = observed_order(L2p,dxValues);             % order from successive global pressure errors
    orderT = observed_order(L2T,dxValues);             % order from successive global temperature errors

    waveOrderRho = observed_order(waveL2rho,dxValues); % order based only on the wave region
    waveOrderU = observed_order(waveL2u,dxValues);     % wave-region velocity order
    waveOrderP = observed_order(waveL2p,dxValues);     % wave-region pressure order
    waveOrderT = observed_order(waveL2T,dxValues);     % wave-region temperature order

    %% Assemble and Print the Complete Convergence Table
    convergenceTable = table(NxValues,dxValues,runTime, ...
        L2rho,orderRho,L2u,orderU,L2p,orderP,L2T,orderT, ...
        waveL2rho,waveOrderRho,waveL2u,waveOrderU, ...
        waveL2p,waveOrderP,waveL2T,waveOrderT,TpeakExcess);

    fprintf('\n============================================================\n');
    fprintf('FINAL L2 CONVERGENCE TABLE\n');
    disp(convergenceTable);

    %% Large Log-Log L2 Convergence Plot
    convergenceFigure = figure('Name','KT-D1V5 Generalized-Minmod L2 Convergence', ...
        'Color','w','Position',[80 80 1500 900]);       % large readable figure
    layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    loglog(dxValues,L2rho,'o-','LineWidth',2.2,'MarkerSize',8); hold on;
    loglog(dxValues,L2u,'s-','LineWidth',2.2,'MarkerSize',8);
    loglog(dxValues,L2p,'d-','LineWidth',2.2,'MarkerSize',8);
    loglog(dxValues,L2T,'^-','LineWidth',2.2,'MarkerSize',8);
    set(gca,'XDir','reverse','FontSize',14,'LineWidth',1.1);
    xlabel('Grid spacing, \Deltax','FontSize',16,'FontWeight','bold');
    ylabel('Global L2 error','FontSize',16,'FontWeight','bold');
    title('Whole Domain','FontSize',18,'FontWeight','bold');
    legend('\rho','u','p','T','Location','best','FontSize',13);
    grid on; box on;

    nexttile;
    loglog(dxValues,waveL2rho,'o-','LineWidth',2.2,'MarkerSize',8); hold on;
    loglog(dxValues,waveL2u,'s-','LineWidth',2.2,'MarkerSize',8);
    loglog(dxValues,waveL2p,'d-','LineWidth',2.2,'MarkerSize',8);
    loglog(dxValues,waveL2T,'^-','LineWidth',2.2,'MarkerSize',8);
    set(gca,'XDir','reverse','FontSize',14,'LineWidth',1.1);
    xlabel('Grid spacing, \Deltax','FontSize',16,'FontWeight','bold');
    ylabel('Wave-region L2 error','FontSize',16,'FontWeight','bold');
    title('Rarefaction, Contact, and Shock Region','FontSize',18,'FontWeight','bold');
    legend('\rho','u','p','T','Location','best','FontSize',13);
    grid on; box on;

    title(layout,sprintf('KT-D1V5 FVDBM L2 Convergence: Generalized Minmod \\theta=%.2f',thetaValue), ...
        'FontSize',22,'FontWeight','bold');

    %% Save Numerical Results
    projectFolder = fileparts(fileparts(thisFolder));   % project directory above matlab/d1v5
    resultsFolder = fullfile(projectFolder,'results');  % central folder for compact validation outputs
    if ~isfolder(resultsFolder); mkdir(resultsFolder); end % create it if this is a fresh clone
    matFile = fullfile(resultsFolder,'d1v5_gminmod_L2_convergence_results.mat');
    csvFile = fullfile(resultsFolder,'d1v5_gminmod_L2_convergence_results.csv');
    pngFile = fullfile(resultsFolder,'d1v5_gminmod_L2_convergence_plot.png');
    save(matFile,'convergenceTable','NxValues','dxValues','thetaValue','tauValue');
    writetable(convergenceTable,csvFile);
    exportgraphics(convergenceFigure,pngFile,'Resolution',200); % preserve the convergence graph for reports

    fprintf('Saved MAT results to:\n%s\n',matFile);
    fprintf('Saved CSV table to:\n%s\n',csvFile);
    fprintf('Saved convergence plot to:\n%s\n',pngFile);

    setenv('DBM_NX','');                              % remove convergence-test grid override
    setenv('DBM_THETA','');                           % remove convergence-test theta override
    setenv('DBM_TAU','');                             % remove convergence-test relaxation-time override
end

function metrics = run_one_resolution(nxValue,thetaValue,tauValue,solverFile)
%RUN_ONE_RESOLUTION passes one grid to the existing generalized-minmod solver.
% Environment variables survive the solver script's initial clear command.

    setenv('DBM_NX',num2str(nxValue));                 % request this grid resolution
    setenv('DBM_THETA',num2str(thetaValue,'%.15g'));   % hold theta fixed
    setenv('DBM_TAU',num2str(tauValue,'%.15g'));       % hold tau fixed

    run(solverFile);                                   % execute the existing validated FVDBM solver

    waveMask = x >= (x0-0.30) & x <= (x0+0.40);      % same physical wave window at every resolution

    metrics.dx = dx;                                   % actual grid spacing used by solver
    metrics.runTime = elapsedTime;                     % solver runtime
    metrics.L2rho = sqrt(mean((rho-exact.rho).^2));   % global density L2 norm
    metrics.L2u = sqrt(mean((u-exact.u).^2));         % global velocity L2 norm
    metrics.L2p = sqrt(mean((p-exact.p).^2));         % global pressure L2 norm
    metrics.L2T = sqrt(mean((T-exact.T).^2));         % global temperature L2 norm

    metrics.waveL2rho = sqrt(mean((rho(waveMask)-exact.rho(waveMask)).^2)); % wave-region density L2
    metrics.waveL2u = sqrt(mean((u(waveMask)-exact.u(waveMask)).^2));       % wave-region velocity L2
    metrics.waveL2p = sqrt(mean((p(waveMask)-exact.p(waveMask)).^2));       % wave-region pressure L2
    metrics.waveL2T = sqrt(mean((T(waveMask)-exact.T(waveMask)).^2));       % wave-region temperature L2
    metrics.TpeakExcess = max(T(waveMask))-max(exact.T(waveMask));          % positive means overshoot

    close all;                                         % close this resolution's comparison figure
end

function order = observed_order(errorValues,dxValues)
%OBSERVED_ORDER computes p=log(Ecoarse/Efine)/log(dxcoarse/dxfine).

    order = NaN(size(errorValues));                    % first grid has no coarser partner

    for k = 2:length(errorValues)
        order(k) = log(errorValues(k-1)/errorValues(k)) / ...
            log(dxValues(k-1)/dxValues(k));            % measured order between two adjacent grids
    end
end
