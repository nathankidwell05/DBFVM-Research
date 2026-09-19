function run_d1v5_theta_sweep(thetaValue)
%% Generalized-Minmod Theta Sweep for the KT-D1V5 Solver %%
% Runs the live solver (d1v5_sodshock_muscl_MC_RK3.m) once with a requested
% generalized-minmod theta. Everything else is held at the live settings:
%     Nx = 6250, tau = 5e-6, CFL = 0.025, tEnd = 0.15, Lx = 20.
% theta = 1 reproduces minmod and theta = 2 reproduces MC, so a sweep over
% theta traces the full diffusion-versus-overshoot tradeoff curve.
%
% Several MATLAB sessions can run different theta values at the same time:
%     run_d1v5_theta_sweep(1.4)
% Each result is saved to results/theta_sweep/ with theta in the filename.

    [projectFolder,solverFile] = project_paths();              % solver folder (OneDrive or repository layout)
    outFolder = fullfile(projectFolder,'results','theta_sweep'); % sweep output folder
    if ~isfolder(outFolder); mkdir(outFolder); end

    cleanupEnvironment = onCleanup(@clear_dbm_environment); %#ok<NASGU> % restore standalone behavior

    setenv('DBM_NX','6250');                                    % matched comparison grid
    setenv('DBM_TAU','5e-06');                                  % live relaxation time
    setenv('DBM_CFL','0.025');                                  % live CFL number
    setenv('DBM_THETA',num2str(thetaValue,'%.15g'));            % the one quantity being varied
    setenv('DBM_DISABLE_LIVE_PLOT','1');                        % no animation in batch runs

    sweep = run_solver_once(solverFile);                        % solver runs in its own workspace

    outFile = fullfile(outFolder,sprintf('d1v5_gminmod_Nx6250_tau5e-06_CFL0p025_theta%s.mat', ...
        strrep(sprintf('%.2f',thetaValue),'.','p')));           % e.g. ..._theta1p40.mat
    save(outFile,'sweep');
    fprintf('Theta-sweep result saved to:\n%s\n',outFile);
    close all;
end

function sweep = run_solver_once(solverFile)
%RUN_SOLVER_ONCE runs the solver in a separate workspace.
% The solver starts with "clear", which would also delete the onCleanup
% object if both lived in the same workspace and wipe the overrides early.

    run(solverFile);                                            % execute the unchanged solver

    names = {'Nx','dx','x0','tEnd','tau','CFL','limiterTheta','dtBase','step', ...
        'elapsedTime','retryCountTotal','fallbackCountTotal','x','rho','u','p','T','exact'};
    sweep = struct();                                           % collected solver output
    for k = 1:numel(names)
        sweep.(names{k}) = eval(names{k});                      % read solver variables at run time
    end
end

function clear_dbm_environment
%CLEAR_DBM_ENVIRONMENT removes every override so manual runs use file defaults.
    setenv('DBM_NX',''); setenv('DBM_TAU',''); setenv('DBM_CFL','');
    setenv('DBM_THETA',''); setenv('DBM_DISABLE_LIVE_PLOT','');
end
