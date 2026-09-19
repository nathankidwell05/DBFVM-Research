function run_d1v5_limiter_grid(limiterName,nxValue)
%% One DBM Run for the Method Comparison %%
% Runs the live KT-D1V5 solver with a chosen limiter on one grid and saves
% the result in results/method_runs/. Everything else stays at the live
% settings: tau = 5e-6, CFL = 0.025, theta = 1.20, tEnd = 0.15, Lx = 20.
%
% Examples (several MATLAB sessions can run different cases at once):
%     run_d1v5_limiter_grid('firstorder',6250)
%     run_d1v5_limiter_grid('minmod',12500)
%     run_d1v5_limiter_grid('mc',25000)
% The generalized-minmod runs come from the convergence study cache.

    [projectFolder,solverFile] = project_paths();              % solver folder (OneDrive or repository layout)
    outFolder = fullfile(projectFolder,'results','method_runs'); % method-comparison output folder
    if ~isfolder(outFolder); mkdir(outFolder); end
    outFile = fullfile(outFolder,sprintf('d1v5_%s_tau5e-06_CFL0p025_Nx%d.mat',limiterName,nxValue));
    if isfile(outFile)
        fprintf('Already computed: %s\n',outFile);
        return;
    end

    cleanupEnvironment = onCleanup(@clear_dbm_environment); %#ok<NASGU> % restore standalone behavior
    setenv('DBM_NX',num2str(nxValue));                          % requested grid
    setenv('DBM_TAU','5e-06');                                  % live relaxation time
    setenv('DBM_CFL','0.025');                                  % live CFL number
    setenv('DBM_THETA','1.2');                                  % only used by gminmod
    setenv('DBM_LIMITER',limiterName);                          % the one quantity being varied
    setenv('DBM_DISABLE_LIVE_PLOT','1');                        % no animation in batch runs

    runData = run_solver_once(solverFile);                      % solver runs in its own workspace
    save(outFile,'runData','-v7.3');
    fprintf('Saved %s\n',outFile);
    close all;
end

function runData = run_solver_once(solverFile)
%RUN_SOLVER_ONCE runs the solver in a separate workspace so its "clear"
% cannot delete the onCleanup object of the calling function.

    run(solverFile);                                            % execute the unchanged solver
    names = {'Nx','dx','x0','tEnd','tau','CFL','mainLimiter','limiterTheta','dtBase', ...
        'maxSpeed','step','elapsedTime','retryCountTotal','fallbackCountTotal','x','rho','u','p','T','exact'};
    runData = struct();
    for k = 1:numel(names)
        runData.(names{k}) = eval(names{k});                    % read solver variables at run time
    end
    runData.mainLimiter = char(runData.mainLimiter);
end

function clear_dbm_environment
%CLEAR_DBM_ENVIRONMENT removes every override so manual runs use file defaults.
    setenv('DBM_NX',''); setenv('DBM_TAU',''); setenv('DBM_CFL','');
    setenv('DBM_THETA',''); setenv('DBM_LIMITER',''); setenv('DBM_DISABLE_LIVE_PLOT','');
end
