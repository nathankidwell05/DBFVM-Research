function run_cost_benchmark
%% Accuracy-Versus-Cost Benchmark %%
% 1. Runs the conventional Euler solver (first order and MUSCL, HLLC flux,
%    CFL 0.5) on the five study grids and stores errors and wall-clock time
%    (best of three repetitions).
% 2. Times the KT-D1V5 solver alone on the two coarsest grids (generalized
%    minmod and first order), so that its cost per cell per time step can be
%    measured without other MATLAB sessions competing for the processor.
% Run this with nothing else running on the computer.
% Results: results/method_runs/euler_*.mat and results/method_runs/timing_benchmark.mat

    [projectFolder,solverFile] = project_paths();
    outFolder = fullfile(projectFolder,'results','method_runs');
    if ~isfolder(outFolder); mkdir(outFolder); end
    grids = [3125 6250 12500 25000 50000];

    %% Euler solver: all grids, both schemes
    schemes = {'firstorder','muscl'};
    for s = 1:numel(schemes)
        for Nx = grids
            best = inf;
            for rep = 1:3                                   % best of three removes start-up noise
                out = euler_sod_solver(Nx,schemes{s},1.2,0.5);
                best = min(best,out.runTime);
            end
            out.runTime = best;
            save(fullfile(outFolder,sprintf('euler_%s_Nx%d.mat',schemes{s},Nx)),'out');
            fprintf('Euler %-10s Nx=%6d  steps=%5d  time=%.3f s\n',schemes{s},Nx,out.steps,best);
        end
    end

    %% DBM solver: clean timing on the two coarsest grids
    cases = {'gminmod',3125; 'gminmod',6250; 'firstorder',3125};
    dbmTiming = struct('limiter',{},'Nx',{},'steps',{},'runTime',{},'costPerCellStep',{});
    for k = 1:size(cases,1)
        cleanupEnvironment = onCleanup(@clear_dbm_environment); %#ok<NASGU>
        setenv('DBM_NX',num2str(cases{k,2})); setenv('DBM_TAU','5e-06'); setenv('DBM_CFL','0.025');
        setenv('DBM_THETA','1.2'); setenv('DBM_LIMITER',cases{k,1}); setenv('DBM_DISABLE_LIVE_PLOT','1');
        r = run_solver_timed(solverFile);
        clear_dbm_environment;
        dbmTiming(k).limiter = cases{k,1};
        dbmTiming(k).Nx = cases{k,2};
        dbmTiming(k).steps = r.step;
        dbmTiming(k).runTime = r.elapsedTime;
        dbmTiming(k).costPerCellStep = r.elapsedTime/(r.step*cases{k,2});
        fprintf('DBM %-10s Nx=%6d  steps=%6d  time=%.1f s  cost per cell-step=%.3e s\n', ...
            cases{k,1},cases{k,2},r.step,r.elapsedTime,dbmTiming(k).costPerCellStep);
    end
    save(fullfile(outFolder,'timing_benchmark.mat'),'dbmTiming');
end

function r = run_solver_timed(solverFile)
    run(solverFile);
    r.step = eval('step'); r.elapsedTime = eval('elapsedTime');
    close all;
end

function clear_dbm_environment
    setenv('DBM_NX',''); setenv('DBM_TAU',''); setenv('DBM_CFL','');
    setenv('DBM_THETA',''); setenv('DBM_LIMITER',''); setenv('DBM_DISABLE_LIVE_PLOT','');
end
