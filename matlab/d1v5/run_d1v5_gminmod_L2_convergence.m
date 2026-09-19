function convergenceTable = run_d1v5_gminmod_L2_convergence(onlyNx)
%% L2 Grid-Convergence Test for the KT-D1V5 Generalized-Minmod Solver %%
% This function calls the existing generalized-minmod Sod solver five times.
% Only Nx and dx change. The physical problem, tau, CFL, tEnd, theta, and
% initial conditions remain fixed.
%
% Run from the MATLAB Command Window with:
%     convergenceTable = run_d1v5_gminmod_L2_convergence;
%
% Optional: run (and cache) only one of the grids listed in NxValues, e.g.
%     run_d1v5_gminmod_L2_convergence(6250);
% This lets several MATLAB sessions compute different grids at the same time.
% Calling the function again with no input then assembles the table and plot
% from the cached runs without repeating any solver run.
%
% Both axes of the convergence plot are logarithmic. For
%     L2 error = C*(dx^p),
% the slope of log(L2 error) against log(dx) is the convergence order p.
% The logarithm base does not affect the measured slope.
%
% Revision notes (September 18, 2026):
%   * Grids are 3125, 6250, 12500, 25000, 50000 (the accidental 400000 entry
%     was removed). tau = 5e-6, theta = 1.20, CFL = 0.025.
%   * CFL is now passed to the solver through DBM_CFL so the driver, not the
%     solver default, fixes the value used in the study.
%   * Every finished grid is cached in results/runs/ with its settings in the
%     filename. A cached grid is reused only if Nx, tau, theta, CFL, and tEnd
%     all match.
%   * Output files carry theta, tau, and CFL in their names so one study can no
%     longer overwrite another.
%   * L1 errors, dt, effective CFL, retries, and fallbacks are recorded.
%   * onCleanup clears the DBM_* environment variables even after an error or
%     a Ctrl+C interruption.

    if nargin < 1; onlyNx = []; end                   % default: assemble the complete five-grid study

    clc;                                               % clear old Command Window text
    close all;                                         % close old figures

    %% Convergence-Test Settings
    NxValues = [3125;6250;12500;25000;50000];          % five grids; each refinement doubles Nx
    thetaValue = 1.20;                                 % selected generalized-minmod setting
    tauValue = 5.0e-6;                                 % fixed BGK relaxation time
    cflValue = 0.025;                                  % fixed transport CFL number
    numberOfRuns = numel(NxValues);                    % number of resolutions

    thisFolder = fileparts(mfilename('fullpath'));     % folder containing this function
    solverFile = fullfile(thisFolder,'d1v5_sodshock_muscl_MC_RK3.m'); % live OneDrive solver name
    if ~isfile(solverFile)
        solverFile = fullfile(thisFolder,'d1v5_sodshock_gminmod_rk3.m'); % GitHub repository solver name
    end
    if ~isfile(solverFile)
        error('The generalized-minmod solver was not found in: %s',thisFolder);
    end

    resultsFolder = fullfile(thisFolder,'results');    % keep convergence output inside the research folder
    runFolder = fullfile(resultsFolder,'runs');        % one cached MAT file per finished grid
    if ~isfolder(runFolder); mkdir(runFolder); end

    cleanupEnvironment = onCleanup(@clear_dbm_environment); %#ok<NASGU> % always restore standalone behavior

    %% Optional Single-Grid Mode
    if ~isempty(onlyNx)                                % run only one requested grid and stop
        if ~any(NxValues == onlyNx)
            error('Nx = %d is not one of the convergence grids.',onlyNx);
        end
        get_or_run_resolution(onlyNx,thetaValue,tauValue,cflValue,solverFile,runFolder);
        convergenceTable = table();                    % nothing to tabulate in single-grid mode
        return;
    end

    %% Storage
    dxValues = zeros(numberOfRuns,1);                  % grid spacing in each run
    dtValues = zeros(numberOfRuns,1);                  % base time step actually used
    effectiveCFL = zeros(numberOfRuns,1);              % max|c|*dt/dx actually used
    runTime = zeros(numberOfRuns,1);                   % wall-clock runtime
    retries = zeros(numberOfRuns,1);                   % RK-stage retries (should be zero)
    fallbacks = zeros(numberOfRuns,1);                 % limiter fallbacks (should be zero)
    L1rho = zeros(numberOfRuns,1);                     % global density L1 error
    L1u = zeros(numberOfRuns,1);                       % global velocity L1 error
    L1p = zeros(numberOfRuns,1);                       % global pressure L1 error
    L1T = zeros(numberOfRuns,1);                       % global temperature L1 error
    L2rho = zeros(numberOfRuns,1);                     % global density L2 error
    L2u = zeros(numberOfRuns,1);                       % global velocity L2 error
    L2p = zeros(numberOfRuns,1);                       % global pressure L2 error
    L2T = zeros(numberOfRuns,1);                       % global temperature L2 error
    waveL2rho = zeros(numberOfRuns,1);                 % wave-region density L2 error
    waveL2u = zeros(numberOfRuns,1);                   % wave-region velocity L2 error
    waveL2p = zeros(numberOfRuns,1);                   % wave-region pressure L2 error
    waveL2T = zeros(numberOfRuns,1);                   % wave-region temperature L2 error
    TpeakExcess = zeros(numberOfRuns,1);               % numerical peak T minus exact peak T

    fprintf('KT-D1V5 generalized-minmod L2 convergence test\n');
    fprintf('Five resolutions, theta = %.2f, tau = %.3e, CFL = %.3f\n',thetaValue,tauValue,cflValue);

    %% Run (or Reload) the Existing Solver on Five Grids
    for runNumber = 1:numberOfRuns
        requestedNx = NxValues(runNumber);             % resolution for this run
        fprintf('\n============================================================\n');
        fprintf('Resolution %d/%d: Nx = %d\n',runNumber,numberOfRuns,requestedNx);

        runData = get_or_run_resolution(requestedNx,thetaValue,tauValue,cflValue,solverFile,runFolder);
        metrics = compute_metrics(runData);            % identical error definitions on every grid

        dxValues(runNumber) = runData.dx;              % save actual grid spacing
        dtValues(runNumber) = runData.dtBase;          % save actual base time step
        effectiveCFL(runNumber) = runData.maxSpeed*runData.dtBase/runData.dx; % transport CFL actually used
        runTime(runNumber) = runData.runTime;          % save runtime
        retries(runNumber) = runData.retryCountTotal;  % save retry count
        fallbacks(runNumber) = runData.fallbackCountTotal; % save fallback count
        L1rho(runNumber) = metrics.L1rho;              % save density L1 error
        L1u(runNumber) = metrics.L1u;                  % save velocity L1 error
        L1p(runNumber) = metrics.L1p;                  % save pressure L1 error
        L1T(runNumber) = metrics.L1T;                  % save temperature L1 error
        L2rho(runNumber) = metrics.L2rho;              % save density error
        L2u(runNumber) = metrics.L2u;                  % save velocity error
        L2p(runNumber) = metrics.L2p;                  % save pressure error
        L2T(runNumber) = metrics.L2T;                  % save temperature error
        waveL2rho(runNumber) = metrics.waveL2rho;      % save wave density error
        waveL2u(runNumber) = metrics.waveL2u;          % save wave velocity error
        waveL2p(runNumber) = metrics.waveL2p;          % save wave pressure error
        waveL2T(runNumber) = metrics.waveL2T;          % save wave temperature error
        TpeakExcess(runNumber) = metrics.TpeakExcess;  % save peak-temperature difference

        fprintf('dt = %.4e, effective CFL = %.4e, retries = %d, fallbacks = %d\n', ...
            dtValues(runNumber),effectiveCFL(runNumber),retries(runNumber),fallbacks(runNumber));
        fprintf('Global L2: rho=%.6e, u=%.6e, p=%.6e, T=%.6e\n', ...
            metrics.L2rho,metrics.L2u,metrics.L2p,metrics.L2T);
        fprintf('Wave L2:   rho=%.6e, u=%.6e, p=%.6e, T=%.6e\n', ...
            metrics.waveL2rho,metrics.waveL2u,metrics.waveL2p,metrics.waveL2T);
    end

    %% Pairwise and Overall Convergence Orders
    orderRho = observed_order(L2rho,dxValues);         % density order after each refinement
    orderU = observed_order(L2u,dxValues);             % velocity order after each refinement
    orderP = observed_order(L2p,dxValues);             % pressure order after each refinement
    orderT = observed_order(L2T,dxValues);             % temperature order after each refinement
    waveOrderRho = observed_order(waveL2rho,dxValues); % wave-region density order
    waveOrderU = observed_order(waveL2u,dxValues);     % wave-region velocity order
    waveOrderP = observed_order(waveL2p,dxValues);     % wave-region pressure order
    waveOrderT = observed_order(waveL2T,dxValues);     % wave-region temperature order

    overallOrderRho = loglog_slope(dxValues,L2rho);    % best-fit density order using all grids
    overallOrderU = loglog_slope(dxValues,L2u);        % best-fit velocity order using all grids
    overallOrderP = loglog_slope(dxValues,L2p);        % best-fit pressure order using all grids
    overallOrderT = loglog_slope(dxValues,L2T);        % best-fit temperature order using all grids

    fprintf('\nOverall log-log slopes using all five grids:\n');
    fprintf('rho = %.4f, u = %.4f, p = %.4f, T = %.4f\n', ...
        overallOrderRho,overallOrderU,overallOrderP,overallOrderT);

    %% Results Table
    convergenceTable = table(NxValues,dxValues,dtValues,effectiveCFL,runTime,retries,fallbacks, ...
        L1rho,L1u,L1p,L1T, ...
        L2rho,orderRho,L2u,orderU,L2p,orderP,L2T,orderT, ...
        waveL2rho,waveOrderRho,waveL2u,waveOrderU, ...
        waveL2p,waveOrderP,waveL2T,waveOrderT,TpeakExcess);
    disp(convergenceTable);

    %% Log-Log Convergence Plot
    convergenceFigure = figure('Name','KT-D1V5 Generalized-Minmod L2 Convergence', ...
        'Color','w','Position',[80 80 1500 900]);       % large readable figure
    layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    loglog(dxValues,L2rho,'o-','LineWidth',2.2,'MarkerSize',8); hold on;
    loglog(dxValues,L2u,'s-','LineWidth',2.2,'MarkerSize',8);
    loglog(dxValues,L2p,'d-','LineWidth',2.2,'MarkerSize',8);
    loglog(dxValues,L2T,'^-','LineWidth',2.2,'MarkerSize',8);
    set(gca,'XDir','reverse','FontSize',14,'LineWidth',1.1); % refinement moves left to right
    xlabel('Grid spacing, \Deltax (log scale)','FontSize',16,'FontWeight','bold');
    ylabel('Global L2 error (log scale)','FontSize',16,'FontWeight','bold');
    title('Whole Domain','FontSize',18,'FontWeight','bold');
    legend('\rho','u','p','T','Location','best','FontSize',13);
    grid on; box on;

    nexttile;
    loglog(dxValues,waveL2rho,'o-','LineWidth',2.2,'MarkerSize',8); hold on;
    loglog(dxValues,waveL2u,'s-','LineWidth',2.2,'MarkerSize',8);
    loglog(dxValues,waveL2p,'d-','LineWidth',2.2,'MarkerSize',8);
    loglog(dxValues,waveL2T,'^-','LineWidth',2.2,'MarkerSize',8);
    set(gca,'XDir','reverse','FontSize',14,'LineWidth',1.1);
    xlabel('Grid spacing, \Deltax (log scale)','FontSize',16,'FontWeight','bold');
    ylabel('Wave-region L2 error (log scale)','FontSize',16,'FontWeight','bold');
    title('Rarefaction, Contact, and Shock Region','FontSize',18,'FontWeight','bold');
    legend('\rho','u','p','T','Location','best','FontSize',13);
    grid on; box on;

    title(layout,sprintf('KT-D1V5 FVDBM L2 Convergence: Generalized Minmod, \\theta=%.2f, \\tau=%.1e, CFL=%.3f', ...
        thetaValue,tauValue,cflValue),'FontSize',20,'FontWeight','bold');

    %% Save Table and Plot
    % Settings are part of every filename so studies with different settings
    % cannot overwrite each other.
    tag = settings_tag(thetaValue,tauValue,cflValue);  % e.g. theta1p20_tau5e-06_CFL0p025
    matFile = fullfile(resultsFolder,['d1v5_gminmod_L2_convergence_' tag '_results.mat']);
    csvFile = fullfile(resultsFolder,['d1v5_gminmod_L2_convergence_' tag '_results.csv']);
    pngFile = fullfile(resultsFolder,['d1v5_gminmod_L2_convergence_' tag '_plot.png']);

    save(matFile,'convergenceTable','NxValues','dxValues','dtValues','effectiveCFL', ...
        'thetaValue','tauValue','cflValue', ...
        'overallOrderRho','overallOrderU','overallOrderP','overallOrderT');
    writetable(convergenceTable,csvFile);
    exportgraphics(convergenceFigure,pngFile,'Resolution',200);

    fprintf('Saved MAT results to:\n%s\n',matFile);
    fprintf('Saved CSV results to:\n%s\n',csvFile);
    fprintf('Saved plot to:\n%s\n',pngFile);
end

function runData = get_or_run_resolution(nxValue,thetaValue,tauValue,cflValue,solverFile,runFolder)
%GET_OR_RUN_RESOLUTION reloads a cached grid or runs the solver and caches it.

    cacheFile = fullfile(runFolder,sprintf('d1v5_gminmod_%s_Nx%d.mat', ...
        settings_tag(thetaValue,tauValue,cflValue),nxValue)); % one file per grid and setting

    if isfile(cacheFile)                               % a finished run with these settings may already exist
        cached = load(cacheFile,'runData');
        runData = cached.runData;
        settingsMatch = runData.Nx == nxValue ...
            && abs(runData.tau-tauValue) <= 1e-12*tauValue ...
            && abs(runData.theta-thetaValue) <= 1e-12 ...
            && abs(runData.CFL-cflValue) <= 1e-12 ...
            && abs(runData.tEnd-0.15) <= 1e-12;
        if settingsMatch
            fprintf('Reusing cached run: %s\n',cacheFile);
            return;
        end
        warning('Cached file %s does not match the requested settings; rerunning.',cacheFile);
    end

    runData = run_one_resolution(nxValue,thetaValue,tauValue,cflValue,solverFile);
    save(cacheFile,'runData','-v7.3');                 % -v7.3 handles large grids
    fprintf('Cached run saved to:\n%s\n',cacheFile);
end

function runData = run_one_resolution(nxValue,thetaValue,tauValue,cflValue,solverFile)
%RUN_ONE_RESOLUTION calls the existing solver with one controlled grid.
% Environment variables survive the solver script's initial clear command.

    setenv('DBM_NX',num2str(nxValue));                 % request this number of cells
    setenv('DBM_THETA',num2str(thetaValue,'%.15g'));   % hold theta at 1.20
    setenv('DBM_TAU',num2str(tauValue,'%.15g'));       % hold tau fixed
    setenv('DBM_CFL',num2str(cflValue,'%.15g'));       % hold CFL fixed
    setenv('DBM_DISABLE_LIVE_PLOT','1');               % animation is unnecessary and fragile during automated grid runs

    run(solverFile);                                   % execute the existing generalized-minmod solver

    % Everything below comes from the solver script's own workspace, so it
    % records the settings that were actually used, not the requested ones.
    % eval looks up each name at run time, so a solver variable such as
    % "step" is read as the variable rather than a MATLAB function.
    solverNames = {'Nx','Lx','dx','x0','tEnd','tau','CFL','limiterTheta','mainLimiter', ...
        'dtBase','maxSpeed','step','elapsedTime','retryCountTotal','fallbackCountTotal', ...
        'repairCountTotal','x','rho','u','p','T','exact'};
    runData = struct();                                % collected solver output
    for k = 1:numel(solverNames)
        if exist(solverNames{k},'var')                 % not every solver version defines every counter
            runData.(solverNames{k}) = eval(solverNames{k}); % copy one solver variable
        else
            runData.(solverNames{k}) = NaN;            % record "not reported" instead of failing
        end
    end
    runData.theta = runData.limiterTheta;              % short names used by the cache check
    runData.steps = runData.step;
    runData.runTime = runData.elapsedTime;
    runData.mainLimiter = char(runData.mainLimiter);

    close all;                                         % close this run's solver figure
end

function metrics = compute_metrics(runData)
%COMPUTE_METRICS evaluates the same error norms on every grid.

    x = runData.x; exact = runData.exact;              % unpack for readability
    waveMask = x >= (runData.x0-0.30) & x <= (runData.x0+0.40); % same physical wave window on every grid

    metrics.L1rho = mean(abs(runData.rho-exact.rho));  % global density L1
    metrics.L1u = mean(abs(runData.u-exact.u));        % global velocity L1
    metrics.L1p = mean(abs(runData.p-exact.p));        % global pressure L1
    metrics.L1T = mean(abs(runData.T-exact.T));        % global temperature L1

    % RMS-form discrete L2 errors. Since Lx is fixed, this produces the
    % same convergence order as sqrt(sum(error.^2)*dx).
    metrics.L2rho = sqrt(mean((runData.rho-exact.rho).^2)); % global density L2
    metrics.L2u = sqrt(mean((runData.u-exact.u).^2));  % global velocity L2
    metrics.L2p = sqrt(mean((runData.p-exact.p).^2));  % global pressure L2
    metrics.L2T = sqrt(mean((runData.T-exact.T).^2));  % global temperature L2
    metrics.waveL2rho = sqrt(mean((runData.rho(waveMask)-exact.rho(waveMask)).^2));
    metrics.waveL2u = sqrt(mean((runData.u(waveMask)-exact.u(waveMask)).^2));
    metrics.waveL2p = sqrt(mean((runData.p(waveMask)-exact.p(waveMask)).^2));
    metrics.waveL2T = sqrt(mean((runData.T(waveMask)-exact.T(waveMask)).^2));
    metrics.TpeakExcess = max(runData.T(waveMask))-max(exact.T(waveMask));
end

function tag = settings_tag(thetaValue,tauValue,cflValue)
%SETTINGS_TAG builds a filename-safe label such as theta1p20_tau5e-06_CFL0p025.

    tag = sprintf('theta%s_tau%s_CFL%s', ...
        strrep(sprintf('%.2f',thetaValue),'.','p'), ...
        strrep(sprintf('%.0e',tauValue),'+',''), ...
        strrep(sprintf('%.3f',cflValue),'.','p'));
end

function clear_dbm_environment
%CLEAR_DBM_ENVIRONMENT removes every override so manual runs use file defaults.

    setenv('DBM_NX','');                               % remove resolution override
    setenv('DBM_THETA','');                            % remove theta override
    setenv('DBM_TAU','');                              % remove relaxation-time override
    setenv('DBM_CFL','');                              % remove CFL override
    setenv('DBM_DISABLE_LIVE_PLOT','');                % restore normal standalone live plotting
end

function order = observed_order(errorValues,dxValues)
%OBSERVED_ORDER computes the order between consecutive grids.

    order = NaN(size(errorValues));                    % first grid has no coarser comparison
    for k = 2:length(errorValues)
        order(k) = log(errorValues(k-1)/errorValues(k)) / ...
            log(dxValues(k-1)/dxValues(k));            % pairwise order p
    end
end

function slope = loglog_slope(dxValues,errorValues)
%LOGLOG_SLOPE fits log(error)=log(C)+p*log(dx) across all five grids.

    coefficients = polyfit(log(dxValues),log(errorValues),1); % natural-log fit
    slope = coefficients(1);                          % slope equals overall order p
end
