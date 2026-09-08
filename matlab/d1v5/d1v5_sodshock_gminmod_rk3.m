%% One Dimensional Sod Shock Tube Using KT-D1V5, Generalized Minmod, and SSP-RK3 %%
% This program solves a one-dimensional Sod shock tube using a finite-volume DBM.
% The kinetic model is the Kataoka-Tsutahara D1V5 compressible model.
% The transport term is solved using a finite-volume upwind flux.
% The face values are reconstructed with MUSCL.
% The slope limiter uses generalized minmod to balance minmod stability and MC sharpness.
% The time stepping has been upgraded from Forward Euler to SSP-RK3.
% The purpose of these upgrades is to reduce diffusion error and shock/contact smearing.
% This file is intentionally comment-heavy so each command can be followed while learning.

clear;                 % clears variables from the MATLAB workspace
clc;                   % clears the command window
close all;             % closes old plots so this run starts clean
set(0,'DefaultFigureVisible','on');  % show the final comparison and temperature peak-excess marker

%% Physical Parameters
gamma = 1.4;           % ratio of specific heats for an air-like ideal gas
b = 2/(gamma - 1);     % internal degrees of freedom parameter; gamma=1.4 gives b=5
tau = 5.0e-5;          % BGK relaxation time used in the current reference setup
CFL = 0.10;            % CFL number; smaller is safer, larger is faster but can be unstable
tEnd = 0.15;           % final simulation time; keep short enough that waves do not hit boundaries
mainLimiter = "gminmod";  % generalized minmod provides a controlled compromise between minmod and standard MC
limiterTheta = 1.2;       % theta=1 is minmod; theta=2 is MC; 1.2 adds sharpness while staying conservative
thetaText = getenv('DBM_THETA'); % read an optional theta supplied by a controlled study
if ~isempty(thetaText); limiterTheta = str2double(thetaText); end % override the default only when requested
tauText = getenv('DBM_TAU'); % read an optional relaxation time supplied by a controlled study
if ~isempty(tauText); tau = str2double(tauText); end % override the default only when requested

%% Domain
Nx = 50000;               % high-resolution grid used in the current reference setup
nxText = getenv('DBM_NX'); % optionally override the resolution for a controlled study
if ~isempty(nxText); Nx = str2double(nxText); end % apply the requested resolution
Lx = 20.0;              % length of the shock tube
dx = Lx/Nx;            % cell width
x = ((1:Nx)' - 0.5)*dx; % cell-center locations; column vector from left to right
x0 = 0.5*Lx;           % initial diaphragm location

%% Initial Conditions
rhoL = 1.0;            % left density
uL = 0.0;              % left velocity
pL = 1.0;              % left pressure

rhoR = 0.125;          % right density
uR = 0.0;              % right velocity
pR = 0.1;              % right pressure for classic Sod shock tube

rho = rhoR*ones(Nx,1);  % storage for density in each cell; start by filling with right density
u = uR*ones(Nx,1);      % storage for velocity in each cell; start by filling with right velocity
p = pR*ones(Nx,1);      % storage for pressure in each cell; start by filling with right pressure
leftCells = x < x0;                       % logical mask for cells left of the diaphragm
rho(leftCells) = rhoL;                    % assign left density to all left cells at once
u(leftCells) = uL;                        % assign left velocity to all left cells at once
p(leftCells) = pL;                        % assign left pressure to all left cells at once
T = p ./ rho;                             % nondimensional ideal gas relation: p = rho*T

%% KT-D1V5 Velocity Set
% D1 means one physical dimension.
% V5 means five discrete molecular velocities.
% This is not a spatial lattice-streaming method.
% The finite-volume scheme still moves the populations through physical cell faces.

v1 = 1.0;                                 % slow molecular speed used by KT D1V5
v2 = 3.0;                                 % fast molecular speed used by KT D1V5
eta0 = 1.75;                              % internal-energy parameter attached to the rest population
c = [-v2, -v1, 0, v1, v2];                % velocity ordering: fast-left, slow-left, rest, slow-right, fast-right
h2 = [0, 0, eta0^2, 0, 0];               % internal-energy backpack; only the rest population carries eta0^2
Q = length(c);                            % number of discrete velocity populations

%% Build Initial Equilibrium
% The solver evolves f, not rho/u/p directly.
% At the start, we set f equal to feq so the initial state is clean and physical.

f = build_KT_D1V5_equilibrium(rho,u,T,b,v1,v2,eta0,c); % initialize the five populations at equilibrium

%% Time Step
maxSpeed = max(abs(c));                      % largest molecular speed in the D1V5 velocity set
dt_adv = CFL*dx/maxSpeed;                    % advection time-step restriction from finite-volume transport
dt_col = 0.25*tau;                           % explicit collision restriction; keeps BGK relaxation stable
dtBase = min(dt_adv,dt_col);                 % base time step used by the solver

time = 0.0;                                  % current physical time
step = 0;                                    % time-step counter
plotEvery = 50;                              % update live plot every this many steps
printEvery = 500;                            % print progress every this many steps
makeLivePlot = Nx <= 5000;                   % automatically disable live plotting on expensive high-resolution runs
plotSample = unique(round(linspace(1,Nx,min(Nx,2000)))); % sample points used for large-grid plotting
maxRetry = 12;                                % maximum number of times to shrink dt if an RK stage fails
minDt = 1.0e-12;                              % minimum allowed dt before stopping the run
timeTolerance = max(10*eps(max(1,tEnd)),1.0e-8*dtBase); % prevent a roundoff-sized extra final step
retryCountTotal = 0;                          % count how many time-step retries happened
fallbackCountTotal = 0;                       % count how many times the solver had to use a safer limiter

fprintf('Running KT-D1V5 + generalized minmod + SSP-RK3 Sod shock tube...\n'); % print solver description
fprintf('Main limiter: %s\n',char(mainLimiter));                      % print the normal limiter choice
fprintf('Generalized-minmod theta: %.3f\n',limiterTheta);             % print the slope-control value used in this run
fprintf('Positivity repair: disabled\n');                             % reject nonphysical states instead of altering them
fprintf('Grid: Nx=%d, dx=%.3e\n',Nx,dx);                             % print grid information
fprintf('dt_adv=%.3e, dt_col=%.3e, dt=%.3e\n',dt_adv,dt_col,dtBase); % print time-step information
fprintf('Estimated steps: %.0f\n',ceil(tEnd/dtBase));                 % estimate number of time steps before starting

%% Create Live Plot Handles
if makeLivePlot                                  % only create live plots when the grid is not gigantic
    figure('Name','KT-D1V5 Generalized-Minmod SSP-RK3 Sod Shock Tube'); % create live figure window

    subplot(4,1,1);                              % first subplot for density
    hRho = plot(x(plotSample),rho(plotSample),'b-','LineWidth',1.5); % create density plot handle
    ylabel('\rho');                              % label density axis
    grid on;                                     % show grid

    subplot(4,1,2);                              % second subplot for velocity
    hU = plot(x(plotSample),u(plotSample),'b-','LineWidth',1.5); % create velocity plot handle
    ylabel('u');                                 % label velocity axis
    grid on;                                     % show grid

    subplot(4,1,3);                              % third subplot for pressure
    hP = plot(x(plotSample),p(plotSample),'b-','LineWidth',1.5); % create pressure plot handle
    ylabel('p');                                 % label pressure axis
    grid on;                                     % show grid

    subplot(4,1,4);                              % fourth subplot for temperature
    hT = plot(x(plotSample),T(plotSample),'b-','LineWidth',1.5); % create temperature plot handle
    ylabel('T');                                 % label temperature axis
    xlabel('x');                                 % label x axis
    grid on;                                     % show grid
end                                              % end live plot setup

%% Main Time Loop
tic;                                            % start timing the solver

while time < tEnd-timeTolerance                % continue until final time is reached within roundoff

    dt = dtBase;                               % reset dt to the stable base value at the start of each step

    if time + dt > tEnd                        % check whether this step would overshoot final time
        dt = tEnd - time;                      % shrink only the final step so time lands exactly on tEnd
    end                                        % end final-step adjustment

    stepAccepted = false;                      % reset flag saying whether this time step worked
    retryCount = 0;                            % reset retry counter for this time step

    while ~stepAccepted                        % keep trying until the RK3 step is physically acceptable

        if retryCount < 4                                            % first few tries use the sharp limiter
            activeLimiter = mainLimiter;                             % use generalized minmod during normal operation
        elseif retryCount < 8                                        % if generalized minmod keeps failing, soften the reconstruction
            activeLimiter = "minmod";                                % minmod is more diffusive but more stable
            fallbackCountTotal = fallbackCountTotal + 1;             % count that a fallback limiter was used
        else                                                         % if even minmod is unsafe, use the safest transport
            activeLimiter = "firstorder";                            % first-order upwind shuts off MUSCL slopes completely
            fallbackCountTotal = fallbackCountTotal + 1;             % count that a stronger fallback was used
        end                                                          % end limiter fallback choice

        %% SSP-RK3 Time Step Upgrade
        % This section replaces the older Forward Euler update:
        %     f_new = f + dt*rhs
        % Forward Euler only checks the RHS once, then jumps forward.
        % SSP-RK3 checks the RHS three times and blends the results.
        % This is more accurate in time and usually more stable with MUSCL limiters.

        rhs1 = compute_rhs_D1V5(f,b,tau,v1,v2,eta0,c,h2,dx,activeLimiter,limiterTheta); % stage 1 RHS using the current f
        f1 = f + dt*rhs1;                                           % stage 1 trial update
        f1 = apply_transmissive_bc(f1);                             % apply boundary conditions to stage 1 state
        % Do not repair here before the check.
        % If stage 1 creates bad T/rho, we want the retry/fallback logic to see it first.

        if ~is_physical_D1V5(f1,b,c,h2)                             % check whether stage 1 created bad rho or T
            retryCount = retryCount + 1;                            % count this failed try
            retryCountTotal = retryCountTotal + 1;                  % count total failed tries
            dt = 0.5*dt;                                            % cut dt in half and retry
            fprintf('Retry at step %d after stage 1 failed, limiter=%s, new dt = %.3e\n',step+1,char(activeLimiter),dt); % explain retry
            if retryCount > maxRetry || dt < minDt                  % stop if dt became absurdly small
                diagnose_state_D1V5(f1,b,c,h2,'stage 1 failed');    % print detailed physical-state diagnostics
                error('Time step failed after too many retries');    % stop because adaptive retry could not fix it
            end                                                     % end retry safety check
            continue;                                               % restart RK3 attempt with smaller dt
        end                                                         % end stage 1 physical check

        rhs2 = compute_rhs_D1V5(f1,b,tau,v1,v2,eta0,c,h2,dx,activeLimiter,limiterTheta); % stage 2 RHS using the stage 1 state
        f2 = 0.75*f + 0.25*(f1 + dt*rhs2);                          % SSP-RK3 blend for stage 2
        f2 = apply_transmissive_bc(f2);                             % apply boundary conditions to stage 2 state
        % Do not repair here before the check.
        % If stage 2 creates bad T/rho, we want the retry/fallback logic to see it first.

        if ~is_physical_D1V5(f2,b,c,h2)                             % check whether stage 2 created bad rho or T
            retryCount = retryCount + 1;                            % count this failed try
            retryCountTotal = retryCountTotal + 1;                  % count total failed tries
            dt = 0.5*dt;                                            % cut dt in half and retry
            fprintf('Retry at step %d after stage 2 failed, limiter=%s, new dt = %.3e\n',step+1,char(activeLimiter),dt); % explain retry
            if retryCount > maxRetry || dt < minDt                  % stop if dt became absurdly small
                diagnose_state_D1V5(f2,b,c,h2,'stage 2 failed');    % print detailed physical-state diagnostics
                error('Time step failed after too many retries');    % stop because adaptive retry could not fix it
            end                                                     % end retry safety check
            continue;                                               % restart RK3 attempt with smaller dt
        end                                                         % end stage 2 physical check

        rhs3 = compute_rhs_D1V5(f2,b,tau,v1,v2,eta0,c,h2,dx,activeLimiter,limiterTheta); % stage 3 RHS using the stage 2 state
        f_new = (1/3)*f + (2/3)*(f2 + dt*rhs3);                     % final SSP-RK3 update
        f_new = apply_transmissive_bc(f_new);                       % apply boundary conditions to final updated state
        % Do not repair here before the check.
        % If the final RK3 state is bad, the time step should be rejected instead of patched.

        if ~is_physical_D1V5(f_new,b,c,h2)                          % check whether the final updated state is physical
            retryCount = retryCount + 1;                            % count this failed try
            retryCountTotal = retryCountTotal + 1;                  % count total failed tries
            dt = 0.5*dt;                                            % cut dt in half and retry
            fprintf('Retry at step %d after final stage failed, limiter=%s, new dt = %.3e\n',step+1,char(activeLimiter),dt); % explain retry
            if retryCount > maxRetry || dt < minDt                  % stop if dt became absurdly small
                diagnose_state_D1V5(f_new,b,c,h2,'final stage failed'); % print detailed physical-state diagnostics
                error('Time step failed after too many retries');    % stop because adaptive retry could not fix it
            end                                                     % end retry safety check
            continue;                                               % restart RK3 attempt with smaller dt
        end                                                         % end final physical check

        stepAccepted = true;                                        % mark this time step as successful

    end                                                            % end adaptive RK3 retry loop

    step = step + 1;                           % count this accepted time step

    f = f_new;                                                  % accept the SSP-RK3 update
    time = time + dt;                                           % advance physical time

    %% Live Visualization
    if makeLivePlot && (mod(step,plotEvery) == 0 || time >= tEnd) % only update plot occasionally to save speed
        [rhoPlot,uPlot,TPlot,pPlot] = recover_macros_D1V5(f,b,c,h2); % recover macroscopic fields for plotting
        set(hRho,'YData',rhoPlot(plotSample));                  % update density line using sampled points
        set(hU,'YData',uPlot(plotSample));                      % update velocity line using sampled points
        set(hP,'YData',pPlot(plotSample));                      % update pressure line using sampled points
        set(hT,'YData',TPlot(plotSample));                      % update temperature line using sampled points
        sgtitle(sprintf('KT-D1V5 + Generalized Minmod (theta=%.2f) + SSP-RK3: t = %.4f, step = %d',limiterTheta,time,step)); % update title
        drawnow limitrate;                                      % redraw efficiently without slowing MATLAB too much
    end                                                         % end live plot block

    if mod(step,printEvery) == 0 || time >= tEnd                % print progress occasionally
        [rhoCheck,uCheck,TCheck,pCheck] = recover_macros_D1V5(f,b,c,h2); % recover fields so we can monitor stability
        fprintf('step %d, time %.5f / %.5f\n',step,time,tEnd);  % show current step and time
        fprintf('    min(rho)=%.3e, min(T)=%.3e, min(p)=%.3e, max(|u|)=%.3e, retries=%d, fallbacks=%d\n', ...
            min(rhoCheck),min(TCheck),min(pCheck),max(abs(uCheck)),retryCountTotal,fallbackCountTotal); % print physical sanity checks
    end                                                         % end print block

end                                                             % end main time loop

time = tEnd;                                                     % report the requested final time after the roundoff-safe loop

elapsedTime = toc;                                               % stop timing the solver
fprintf('Main solver finished in %.2f seconds.\n',elapsedTime);   % print runtime

%% Final Macroscopic Recovery
[rho,u,T,p] = recover_macros_D1V5(f,b,c,h2);                    % recover final rho, u, T, and p from final f

%% Exact Sod Comparison
exact = sod_exact_solution(x,tEnd,x0,rhoL,uL,pL,rhoR,uR,pR,gamma); % exact Euler Sod solution for comparison

L1_rho = mean(abs(rho - exact.rho));                            % average density error
L1_u = mean(abs(u - exact.u));                                  % average velocity error
L1_p = mean(abs(p - exact.p));                                  % average pressure error
L1_T = mean(abs(T - exact.T));                                  % average temperature error

fprintf('\nValidation Summary\n');                             % print validation header
fprintf('L1(rho) = %.6e\n',L1_rho);                            % print density error
fprintf('L1(u)   = %.6e\n',L1_u);                              % print velocity error
fprintf('L1(p)   = %.6e\n',L1_p);                              % print pressure error
fprintf('L1(T)   = %.6e\n',L1_T);                              % print temperature error

errorRegion = x >= (x0-0.30) & x <= (x0+0.40);               % focus measurement on the moving-wave region
waveL1_rho = mean(abs(rho(errorRegion)-exact.rho(errorRegion))); % wave-region density error
waveL1_u = mean(abs(u(errorRegion)-exact.u(errorRegion)));       % wave-region velocity error
waveL1_p = mean(abs(p(errorRegion)-exact.p(errorRegion)));       % wave-region pressure error
waveL1_T = mean(abs(T(errorRegion)-exact.T(errorRegion)));       % wave-region temperature error
waveLinf_T = max(abs(T(errorRegion)-exact.T(errorRegion)));      % largest local temperature error
[T_peak_DBM,localPeakIndex] = max(T(errorRegion));               % highest numerical temperature inside the wave region
[T_peak_exact,localExactPeakIndex] = max(exact.T(errorRegion));  % highest exact temperature inside the same wave region
waveIndices = find(errorRegion);                                 % convert wave-region positions back to full-grid indices
T_peak_x = x(waveIndices(localPeakIndex));                       % physical location of the numerical temperature peak
T_exact_peak_x = x(waveIndices(localExactPeakIndex));            % first physical location where the exact peak occurs
T_peak_excess = T_peak_DBM-T_peak_exact;                         % signed excess: positive means overshoot; negative means under-prediction
T_peak_excess_positive = max(0,T_peak_excess);                   % overshoot-only value; zero means the numerical peak did not exceed exact
T_peak_excess_percent = 100*T_peak_excess/T_peak_exact;          % signed peak excess as a percentage of the exact peak
fprintf('RUN_RESULT theta=%.4f tau=%.6e Nx=%d rho=%.8e u=%.8e p=%.8e T=%.8e LinfT=%.8e peakT=%.8e\n', ...
    limiterTheta,tau,Nx,waveL1_rho,waveL1_u,waveL1_p,waveL1_T,waveLinf_T,T_peak_excess); % parseable summary

fprintf('\nTemperature Peak Diagnostic (wave region)\n');          % make the peak result easy to find in the Command Window
fprintf('T peak, DBM       = %.8e at x = %.8f\n',T_peak_DBM,T_peak_x); % report numerical peak value and location
fprintf('T peak, exact     = %.8e at x = %.8f\n',T_peak_exact,T_exact_peak_x); % report exact peak value and first peak location
fprintf('T peak excess     = %.8e\n',T_peak_excess);                % positive means DBM overshoot; negative means DBM peak is too low
fprintf('T positive excess = %.8e\n',T_peak_excess_positive);       % overshoot magnitude with negative values clipped to zero
fprintf('T peak excess (%%) = %.6f %%\n',T_peak_excess_percent);     % normalize peak difference by the exact temperature peak

resultsFolder = fullfile(fileparts(mfilename('fullpath')),'results'); % keep generated files beside the solver in one ignored folder
if ~isfolder(resultsFolder); mkdir(resultsFolder); end               % create the output folder on the first run
resultFile = fullfile(resultsFolder,'d1v5_sodshock_gminmod_rk3_latest.mat'); % descriptive output name
save(resultFile, ...
    'x','rho','u','p','T','exact','Nx','dx','dtBase','tEnd','step','tau','CFL','L1_rho','L1_u','L1_p','L1_T', ...
    'waveL1_rho','waveL1_u','waveL1_p','waveL1_T','waveLinf_T','T_peak_DBM','T_peak_exact','T_peak_x','T_exact_peak_x', ...
    'T_peak_excess','T_peak_excess_positive','T_peak_excess_percent','retryCountTotal','fallbackCountTotal', ...
    'mainLimiter','limiterTheta'); % save the error and peak diagnostics for post-processing

%% Final Comparison Plot
figure('Name','KT-D1V5 Generalized-Minmod SSP-RK3 Compared With Exact Sod'); % create final figure

subplot(4,1,1);                                                % density subplot
plot(x(plotSample),rho(plotSample),'g--','LineWidth',1.8);     % generalized-minmod result is dashed green
hold on;                                                       % hold figure for exact solution
plot(x(plotSample),exact.rho(plotSample),'k-','LineWidth',1.8); % exact solution is solid black
ylabel('\rho');                                                % label density axis
legend('DBM','Exact');                                         % label curves
grid on;                                                       % show grid

subplot(4,1,2);                                                % velocity subplot
plot(x(plotSample),u(plotSample),'g--','LineWidth',1.8);       % generalized-minmod result is dashed green
hold on;                                                       % hold figure for exact solution
plot(x(plotSample),exact.u(plotSample),'k-','LineWidth',1.8);  % exact solution is solid black
ylabel('u');                                                   % label velocity axis
legend('DBM','Exact');                                         % label curves
grid on;                                                       % show grid

subplot(4,1,3);                                                % pressure subplot
plot(x(plotSample),p(plotSample),'g--','LineWidth',1.8);       % generalized-minmod result is dashed green
hold on;                                                       % hold figure for exact solution
plot(x(plotSample),exact.p(plotSample),'k-','LineWidth',1.8);  % exact solution is solid black
ylabel('p');                                                   % label pressure axis
legend('DBM','Exact');                                         % label curves
grid on;                                                       % show grid

subplot(4,1,4);                                                % temperature subplot
plot(x(plotSample),T(plotSample),'g--','LineWidth',1.8);       % generalized-minmod result is dashed green
hold on;                                                       % hold figure for exact solution
plot(x(plotSample),exact.T(plotSample),'k-','LineWidth',1.8);  % exact solution is solid black
plot(T_peak_x,T_peak_DBM,'go','MarkerFaceColor','g','MarkerSize',5); % mark the numerical temperature peak being measured
yline(T_peak_exact,'k:','LineWidth',1.0);                       % show the exact peak level so excess is visible vertically
ylabel('T');                                                   % label temperature axis
xlabel('x');                                                   % label x axis
legend('DBM','Exact','DBM peak','Exact peak level','Location','best'); % label curves and peak diagnostic
title(sprintf('Peak excess = %.3e (%.3f%%)',T_peak_excess,T_peak_excess_percent)); % display signed excess above temperature panel
grid on;                                                       % show grid

sgtitle(sprintf('KT-D1V5 FVDBM with Generalized Minmod (theta=%.2f) and SSP-RK3',limiterTheta)); % final plot title

%% Local Functions
% MATLAB requires local functions to stay at the bottom of the script.
% Do not put normal script commands below this point.

function rhs = compute_rhs_D1V5(f,b,tau,v1,v2,eta0,c,h2,dx,limiterType,limiterTheta)
%COMPUTE_RHS_D1V5 builds df/dt for the KT-D1V5 finite-volume DBM equation.
% The RHS has two pieces:
%   1. transport: - d(c_q*f_q)/dx
%   2. collision: (feq_q - f_q)/tau

    [rho,u,T,p] = recover_macros_D1V5(f,b,c,h2);       % recover macroscopic fields from current f

    if any(rho(:) <= 0) || any(T(:) <= 0)              % check for nonphysical density or temperature
        diagnose_state_D1V5(f,b,c,h2,'RHS input failed'); % print detailed diagnostic information
        error('Nonphysical state detected while building RHS'); % stop immediately if the model collapses
    end                                                % end physical-state check

    feq = build_KT_D1V5_equilibrium(rho,u,T,b,v1,v2,eta0,c); % build equilibrium populations from current macros
    collision = (feq - f)/tau;                         % BGK collision term pushing f toward feq
    if limiterType == "hybridmc"                       % only the older binary hybrid method needs a separate troubled-cell mask
        troubledCell = build_macroscopic_troubled_cell_mask(rho,p,T); % locate shocks/contacts for the older hybrid option
    else                                               % generalized minmod controls its slope continuously through theta
        troubledCell = false(size(rho));               % unused placeholder keeps the flux-function inputs consistent
    end                                                % end troubled-cell sensor choice
    fluxDiv = build_MUSCL_flux_divergence(f,c,dx,limiterType,troubledCell,limiterTheta); % finite-volume transport divergence using selected limiter
    rhs = -fluxDiv + collision;                        % combine negative transport divergence with collision source

end                                                    % end RHS function

function fluxDiv = build_MUSCL_flux_divergence(f,c,dx,limiterType,troubledCell,limiterTheta)
%BUILD_MUSCL_FLUX_DIVERGENCE computes d(c_q*f_q)/dx for every population.
% This is the finite-volume transport part of the solver.
% MUSCL reconstructs left/right face values.
% MC is sharp, minmod is safer, and firstorder turns off the MUSCL slope.

    [Nx,Q] = size(f);                                  % get number of cells and number of populations
    fluxDiv = zeros(Nx,Q);                             % storage for flux divergence in every cell and population

    for q = 1:Q                                        % loop through each discrete velocity population

        fq = f(:,q);                                   % pull out one population across all cells
        sigma = zeros(Nx,1);                           % storage for the limited slope in each cell

        %% MUSCL + MC Limiter Upgrade
        % This section replaces the older minmod slope.
        % Minmod chooses the smallest safe slope, which is stable but smeary.
        % MC means monotonized central.
        % MC still blocks dangerous slopes near shocks, but it allows larger slopes in smooth regions.
        % Bigger safe slopes mean sharper shocks and less contact-discontinuity smearing.

        if limiterType == "firstorder"                 % first-order upwind uses no MUSCL slope
            sigma(:) = 0.0;                            % zero slopes make the face value equal to the cell value
        else                                           % otherwise build a limited MUSCL slope
            dL = fq(2:Nx-1) - fq(1:Nx-2);              % left differences for all interior cells at once
            dR = fq(3:Nx) - fq(2:Nx-1);                % right differences for all interior cells at once

            if limiterType == "mc"                     % MC limiter is sharper and less diffusive
                sigma(2:Nx-1) = MC_limiter_vector(dL,dR); % MC-limited slopes for all interior cells at once
            elseif limiterType == "gminmod"            % generalized minmod allows controlled sharpness between minmod and MC
                sigma(2:Nx-1) = generalized_minmod_limiter_vector(dL,dR,limiterTheta); % theta controls the permitted MUSCL slope
            elseif limiterType == "hybridmc"           % use MC in smooth cells and minmod close to macroscopic discontinuities
                slopeMC = MC_limiter_vector(dL,dR);     % calculate sharp MC slopes for this population
                slopeMinmod = minmod_limiter_vector(dL,dR); % calculate safe minmod slopes for this population
                slopeMC(troubledCell(2:Nx-1)) = slopeMinmod(troubledCell(2:Nx-1)); % share the same limiter decision across all five populations
                sigma(2:Nx-1) = slopeMC;                % store the completed hybrid slope for this population
            elseif limiterType == "minmod"             % minmod limiter is safer but more diffusive
                sigma(2:Nx-1) = minmod_limiter_vector(dL,dR); % minmod-limited slopes for all interior cells at once
            else                                       % unknown limiter name means the setup is wrong
                error('Unknown limiter type: %s',char(limiterType)); % stop and report bad limiter name
            end                                        % end limiter type choice
        end                                            % end first-order/MUSCL choice

        Fface = zeros(Nx+1,1);                         % storage for face fluxes; Nx cells have Nx+1 faces

        interiorFaces = 2:Nx;                          % interior faces only; first and last faces are boundaries
        leftCells = 1:Nx-1;                            % cell index on the left side of each interior face
        rightCells = 2:Nx;                             % cell index on the right side of each interior face

        fLeftFace = fq(leftCells) + 0.5*sigma(leftCells); % value at face reconstructed from left cell
        fRightFace = fq(rightCells) - 0.5*sigma(rightCells); % value at face reconstructed from right cell

        if c(q) > 0                                    % positive molecular velocity moves left to right
            Fface(interiorFaces) = c(q)*fLeftFace;     % upwind value comes from left cell
        elseif c(q) < 0                                % negative molecular velocity moves right to left
            Fface(interiorFaces) = c(q)*fRightFace;    % upwind value comes from right cell
        else                                           % zero molecular velocity does not move through faces
            Fface(interiorFaces) = 0.0;                % rest population has no transport flux
        end                                            % end upwind decision

        Fface(1) = Fface(2);                           % left boundary face copies nearest interior face
        Fface(end) = Fface(end-1);                     % right boundary face copies nearest interior face

        fluxDiv(:,q) = (Fface(2:end) - Fface(1:end-1))/dx; % right-face flux minus left-face flux divided by dx

    end                                                % end population loop

end                                                    % end MUSCL flux divergence function

function feq = build_KT_D1V5_equilibrium(rho,u,T,b,v1,v2,eta0,c)
%BUILD_KT_D1V5_EQUILIBRIUM computes the KT-D1V5 equilibrium populations.
% These equations come from the KT D1V5 compressible equilibrium form.
% The coefficients are designed so equilibrium moments recover rho, rho*u, and energy.

    Nx = length(rho);                                  % number of cells
    Q = length(c);                                     % number of populations
    feq = zeros(Nx,Q);                                 % storage for equilibrium populations

    u2 = u.^2;                                         % squared macroscopic velocity

    A_rest = ((b-1)/eta0^2).*T;                        % rest-population coefficient carrying internal energy
    A_slow = (-v2^2 + (((b-1)*v2^2/eta0^2)+1).*T + u2) ./ (2*(v1^2-v2^2)); % coefficient for slow populations
    A_fast = (-v1^2 + (((b-1)*v1^2/eta0^2)+1).*T + u2) ./ (2*(v2^2-v1^2)); % coefficient for fast populations

    B_slow = (-v2^2 + (b+2).*T + u2) ./ (2*v1^2*(v1^2-v2^2)); % directional shift for slow populations
    B_fast = (-v1^2 + (b+2).*T + u2) ./ (2*v2^2*(v2^2-v1^2)); % directional shift for fast populations

    for q = 1:Q                                        % loop through each population
        if q == 3                                      % q=3 is the rest population for ordering [-v2,-v1,0,v1,v2]
            feq(:,q) = rho .* A_rest;                  % rest population has no u*c directional shift
        elseif abs(c(q)) == v1                         % slow populations have speed magnitude v1
            feq(:,q) = rho .* (A_slow + B_slow.*u.*c(q)); % slow KT equilibrium population
        else                                           % remaining nonzero populations are fast populations with speed v2
            feq(:,q) = rho .* (A_fast + B_fast.*u.*c(q)); % fast KT equilibrium population
        end                                            % end population type decision
    end                                                % end equilibrium loop

end                                                    % end KT equilibrium function

function [rho,u,T,p] = recover_macros_D1V5(f,b,c,h2)
%RECOVER_MACROS_D1V5 recovers rho, velocity, temperature, and pressure from f.
% This is moment recovery.
% Density is the sum of populations.
% Momentum is the sum of f*c.
% Temperature is recovered from the energy moment.

    rho = sum(f,2);                                    % density is zeroth moment over q
    momentum = f*c';                                   % momentum is first moment over q
    u = momentum ./ rho;                               % velocity is momentum divided by density
    energy = f*(c.^2 + h2)';                           % energy moment includes molecular speed squared plus h_i^2
    T = (energy ./ rho - u.^2)/b;                      % temperature from energy after subtracting macroscopic kinetic part
    p = rho .* T;                                      % pressure from nondimensional ideal-gas relation

end                                                    % end macro recovery function

function f = apply_transmissive_bc(f)
%APPLY_TRANSMISSIVE_BC applies simple zero-gradient boundary conditions.
% The boundary cell copies the nearest interior cell.
% This is simple and works as long as waves do not strongly hit the boundary before tEnd.

    f(1,:) = f(2,:);                                   % left boundary copies cell 2
    f(end,:) = f(end-1,:);                             % right boundary copies cell Nx-1

end                                                    % end boundary function

function s = MC_limiter_vector(dL,dR)
%MC_LIMITER_VECTOR computes the MC limiter for many cells at once.
% This is the same limiter as MC_limiter, but vectorized for high-resolution runs.
% Vectorized means MATLAB works on whole arrays instead of stepping cell-by-cell.

    a = 2*dL;                                          % first MC candidate: twice the left slope
    b = 0.5*(dL+dR);                                   % second MC candidate: centered slope
    c = 2*dR;                                          % third MC candidate: twice the right slope

    s = zeros(size(dL));                               % default to zero slope wherever the limiter is not comfortable
    positive = (a > 0) & (b > 0) & (c > 0);            % true where all candidate slopes are positive
    negative = (a < 0) & (b < 0) & (c < 0);            % true where all candidate slopes are negative

    s(positive) = min([a(positive),b(positive),c(positive)],[],2); % choose smallest positive safe slope
    s(negative) = max([a(negative),b(negative),c(negative)],[],2); % choose negative slope closest to zero

end                                                    % end vectorized MC limiter function

function s = minmod_limiter_vector(dL,dR)
%MINMOD_LIMITER_VECTOR computes the minmod limiter for many cells at once.
% Minmod is more diffusive than MC.
% We use it as a safety fallback when MC creates negative temperature near a shock/contact.

    s = zeros(size(dL));                               % default to zero slope where signs disagree
    positive = (dL > 0) & (dR > 0);                    % true where both neighboring slopes are positive
    negative = (dL < 0) & (dR < 0);                    % true where both neighboring slopes are negative

    s(positive) = min(dL(positive),dR(positive));      % choose smaller positive slope
    s(negative) = max(dL(negative),dR(negative));      % choose negative slope closest to zero

end                                                    % end vectorized minmod limiter function

function s = generalized_minmod_limiter_vector(dL,dR,theta)
%GENERALIZED_MINMOD_LIMITER_VECTOR controls MUSCL sharpness with theta.
% theta=1 is the conservative end of this limiter family.
% theta=2 gives the standard monotonized-central (MC) limiter.
% Intermediate values reduce minmod diffusion without immediately using full MC.

    if theta < 1.0 || theta > 2.0                      % generalized-minmod TVD parameter should stay between one and two
        error('limiterTheta must remain between 1.0 and 2.0'); % stop rather than silently running an invalid test value
    end                                                % end theta-range check

    leftCandidate = theta*dL;                          % permitted slope based on the change to the left
    centerCandidate = 0.5*(dL+dR);                    % centered slope using information from both sides
    rightCandidate = theta*dR;                        % permitted slope based on the change to the right

    s = zeros(size(dL));                               % default to zero wherever the candidates disagree in sign
    allPositive = leftCandidate > 0 & centerCandidate > 0 & rightCandidate > 0; % smooth positive slope candidates
    allNegative = leftCandidate < 0 & centerCandidate < 0 & rightCandidate < 0; % smooth negative slope candidates

    s(allPositive) = min([leftCandidate(allPositive), ...
        centerCandidate(allPositive),rightCandidate(allPositive)],[],2); % retain the smallest safe positive slope
    s(allNegative) = max([leftCandidate(allNegative), ...
        centerCandidate(allNegative),rightCandidate(allNegative)],[],2); % retain the negative slope closest to zero

end                                                    % end generalized-minmod limiter function

function troubledCell = build_macroscopic_troubled_cell_mask(rho,p,T)
%BUILD_MACROSCOPIC_TROUBLED_CELL_MASK identifies cells close to a shock or contact.
% The old hybrid method examined each population separately.
% That allowed different populations in the same cell to use different limiters,
% which could disturb the combined energy moment and create a temperature spike.
% This function makes one decision from rho, p, and T and shares it across all populations.

    Nx = length(rho);                                  % number of physical cells
    troubledCell = false(Nx,1);                        % begin by assuming every cell is smooth
    fields = {rho,p,T};                                % density sees contacts, pressure sees shocks, and temperature sees thermal jumps
    smallNumber = 1.0e-14;                             % prevents division by zero in constant regions

    for fieldNumber = 1:length(fields)                 % examine each macroscopic field separately
        phi = fields{fieldNumber};                     % select rho, pressure, or temperature
        dL = phi(2:Nx-1) - phi(1:Nx-2);                % change entering each interior cell from its left side
        dR = phi(3:Nx) - phi(2:Nx-1);                  % change leaving each interior cell on its right side

        localScale = abs(phi(2:Nx-1)) + smallNumber;   % local magnitude used to make the jump sensor dimensionless
        relativeJump = max(abs(dL),abs(dR)) ./ localScale; % large value means a strong local jump compared with the field itself
        roughness = abs(dR-dL) ./ ...
            (abs(dR)+abs(dL)+smallNumber);             % near zero for a smooth slope and near one for a sudden change

        fieldTrouble = (roughness > 0.20 & relativeJump > 0.0025) ...
            | relativeJump > 0.08;                     % flag a curved moderate jump or any very strong jump
        troubledCell(2:Nx-1) = troubledCell(2:Nx-1) | fieldTrouble; % combine warnings from rho, pressure, and temperature
    end                                                % end macroscopic-field loop

    originalTrouble = troubledCell;                    % save the initially detected cells before adding a buffer
    troubledCell(2:end) = troubledCell(2:end) | originalTrouble(1:end-1); % also protect one cell to the right
    troubledCell(1:end-1) = troubledCell(1:end-1) | originalTrouble(2:end); % also protect one cell to the left

end                                                    % end shared macroscopic troubled-cell sensor

function ok = is_physical_D1V5(f,b,c,h2)
%IS_PHYSICAL_D1V5 checks whether the recovered macroscopic state is physical.
% The solver can tolerate some negative individual f_i values.
% What cannot go negative are recovered density and recovered temperature.

    [rho,~,T,~] = recover_macros_D1V5(f,b,c,h2);       % recover density and temperature from f
    ok = all(isfinite(rho(:))) && all(isfinite(T(:))) && all(rho(:) > 0) && all(T(:) > 0); % true if state is usable

end                                                    % end physical-state check function

function diagnose_state_D1V5(f,b,c,h2,labelText)
%DIAGNOSE_STATE_D1V5 prints where and how the state became nonphysical.
% This makes MATLAB tell us the actual failure instead of only saying the run crashed.

    [rho,u,T,p] = recover_macros_D1V5(f,b,c,h2);       % recover macroscopic variables from f
    [minRho,idxRho] = min(rho);                        % smallest density and its cell index
    [minT,idxT] = min(T);                              % smallest temperature and its cell index
    [minP,idxP] = min(p);                              % smallest pressure and its cell index

    fprintf('\nPhysical-state diagnostic: %s\n',labelText); % print label for where the failure happened
    fprintf('    min(rho)=%.6e at cell %d\n',minRho,idxRho); % print density diagnostic
    fprintf('    min(T)  =%.6e at cell %d\n',minT,idxT);     % print temperature diagnostic
    fprintf('    min(p)  =%.6e at cell %d\n',minP,idxP);     % print pressure diagnostic
    fprintf('    max(|u|)=%.6e\n',max(abs(u)));              % print largest velocity magnitude

end                                                    % end diagnostic function

function exact = sod_exact_solution(x,t,x0,rhoL,uL,pL,rhoR,uR,pR,gamma)
%SOD_EXACT_SOLUTION computes the exact Euler solution for the Sod shock tube.
% This is only used for plotting and error checks.
% It does not influence the DBM solver.

    if t <= 0                                         % if time is zero
        exact.rho = rhoR*ones(size(x));               % initialize density with right state
        exact.u = uR*ones(size(x));                   % initialize velocity with right state
        exact.p = pR*ones(size(x));                   % initialize pressure with right state
        exact.rho(x < x0) = rhoL;                     % assign left density
        exact.u(x < x0) = uL;                         % assign left velocity
        exact.p(x < x0) = pL;                         % assign left pressure
        exact.T = exact.p ./ exact.rho;               % compute temperature
        return;                                       % leave exact solver
    end                                               % end initial-time case

    aL = sqrt(gamma*pL/rhoL);                         % left sound speed
    aR = sqrt(gamma*pR/rhoR);                         % right sound speed
    pStar = max(0.5*(pL+pR),1e-10);                   % starting guess for star pressure

    for iter = 1:50                                   % Newton iteration for star pressure
        [fL,dfL] = pressure_function(pStar,rhoL,pL,aL,gamma); % left pressure function
        [fR,dfR] = pressure_function(pStar,rhoR,pR,aR,gamma); % right pressure function
        pNew = pStar - (fL + fR + uR - uL)/(dfL + dfR);       % Newton update
        pNew = max(pNew,1e-10);                       % keep pressure positive
        if abs(pNew-pStar)/(pNew+pStar) < 1e-10       % convergence check
            pStar = pNew;                             % accept converged pressure
            break;                                    % leave Newton loop
        end                                           % end convergence check
        pStar = pNew;                                 % update pressure guess
    end                                               % end Newton loop

    [fL,~] = pressure_function(pStar,rhoL,pL,aL,gamma); % left function at star pressure
    [fR,~] = pressure_function(pStar,rhoR,pR,aR,gamma); % right function at star pressure
    uStar = 0.5*(uL + uR + fR - fL);                  % velocity in star region

    xiAll = (x - x0)/t;                               % self-similar coordinate xi=(x-x0)/t
    rhoExact = zeros(size(x));                        % exact density storage
    uExact = zeros(size(x));                          % exact velocity storage
    pExact = zeros(size(x));                          % exact pressure storage

    for k = 1:length(x)                               % loop through x points
        xi = xiAll(k);                                % current self-similar coordinate

        if xi <= uStar                                % point is left of contact discontinuity
            if pStar > pL                             % left wave is a shock
                SL = uL - aL*sqrt((gamma+1)/(2*gamma)*pStar/pL + (gamma-1)/(2*gamma)); % left shock speed
                if xi <= SL                           % undisturbed left state
                    rhoExact(k) = rhoL;               % left density
                    uExact(k) = uL;                   % left velocity
                    pExact(k) = pL;                   % left pressure
                else                                  % left star state
                    rhoExact(k) = rhoL*((pStar/pL + (gamma-1)/(gamma+1))/((gamma-1)/(gamma+1)*pStar/pL + 1)); % left star density
                    uExact(k) = uStar;                % star velocity
                    pExact(k) = pStar;                % star pressure
                end                                   % end left shock decision
            else                                      % left wave is a rarefaction
                aStarL = aL*(pStar/pL)^((gamma-1)/(2*gamma)); % left star sound speed
                SHL = uL - aL;                        % head of left rarefaction
                STL = uStar - aStarL;                 % tail of left rarefaction
                if xi <= SHL                          % undisturbed left state
                    rhoExact(k) = rhoL;               % left density
                    uExact(k) = uL;                   % left velocity
                    pExact(k) = pL;                   % left pressure
                elseif xi > STL                       % left star state
                    rhoExact(k) = rhoL*(pStar/pL)^(1/gamma); % left star density
                    uExact(k) = uStar;                % star velocity
                    pExact(k) = pStar;                % star pressure
                else                                  % inside rarefaction fan
                    uFan = 2/(gamma+1)*(aL + 0.5*(gamma-1)*uL + xi); % fan velocity
                    aFan = 2/(gamma+1)*(aL + 0.5*(gamma-1)*(uL - xi)); % fan sound speed
                    rhoExact(k) = rhoL*(aFan/aL)^(2/(gamma-1)); % fan density
                    uExact(k) = uFan;                 % fan velocity
                    pExact(k) = pL*(aFan/aL)^(2*gamma/(gamma-1)); % fan pressure
                end                                   % end left rarefaction decision
            end                                       % end left wave type decision
        else                                          % point is right of contact discontinuity
            if pStar > pR                             % right wave is a shock
                SR = uR + aR*sqrt((gamma+1)/(2*gamma)*pStar/pR + (gamma-1)/(2*gamma)); % right shock speed
                if xi >= SR                           % undisturbed right state
                    rhoExact(k) = rhoR;               % right density
                    uExact(k) = uR;                   % right velocity
                    pExact(k) = pR;                   % right pressure
                else                                  % right star state
                    rhoExact(k) = rhoR*((pStar/pR + (gamma-1)/(gamma+1))/((gamma-1)/(gamma+1)*pStar/pR + 1)); % right star density
                    uExact(k) = uStar;                % star velocity
                    pExact(k) = pStar;                % star pressure
                end                                   % end right shock decision
            else                                      % right wave is a rarefaction
                aStarR = aR*(pStar/pR)^((gamma-1)/(2*gamma)); % right star sound speed
                SHR = uR + aR;                        % head of right rarefaction
                STR = uStar + aStarR;                 % tail of right rarefaction
                if xi >= SHR                          % undisturbed right state
                    rhoExact(k) = rhoR;               % right density
                    uExact(k) = uR;                   % right velocity
                    pExact(k) = pR;                   % right pressure
                elseif xi <= STR                      % right star state
                    rhoExact(k) = rhoR*(pStar/pR)^(1/gamma); % right star density
                    uExact(k) = uStar;                % star velocity
                    pExact(k) = pStar;                % star pressure
                else                                  % inside rarefaction fan
                    uFan = 2/(gamma+1)*(-aR + 0.5*(gamma-1)*uR + xi); % fan velocity
                    aFan = 2/(gamma+1)*(aR - 0.5*(gamma-1)*(uR - xi)); % fan sound speed
                    rhoExact(k) = rhoR*(aFan/aR)^(2/(gamma-1)); % fan density
                    uExact(k) = uFan;                 % fan velocity
                    pExact(k) = pR*(aFan/aR)^(2*gamma/(gamma-1)); % fan pressure
                end                                   % end right rarefaction decision
            end                                       % end right wave type decision
        end                                           % end contact-side decision
    end                                               % end exact solution loop

    exact.rho = rhoExact;                             % store exact density
    exact.u = uExact;                                 % store exact velocity
    exact.p = pExact;                                 % store exact pressure
    exact.T = pExact ./ rhoExact;                     % store exact temperature

end                                                   % end exact Sod solver

function [f,df] = pressure_function(p,rhoK,pK,aK,gamma)
%PRESSURE_FUNCTION is a helper used by the exact Sod solver.
% It returns the shock or rarefaction pressure function and derivative.

    if p > pK                                         % shock branch
        AK = 2/((gamma+1)*rhoK);                      % shock coefficient A
        BK = (gamma-1)/(gamma+1)*pK;                  % shock coefficient B
        f = (p-pK)*sqrt(AK/(p+BK));                   % shock pressure function
        df = sqrt(AK/(p+BK))*(1 - 0.5*(p-pK)/(p+BK)); % derivative for Newton solve
    else                                              % rarefaction branch
        f = 2*aK/(gamma-1)*((p/pK)^((gamma-1)/(2*gamma)) - 1); % rarefaction pressure function
        df = (1/(rhoK*aK))*(p/pK)^(-(gamma+1)/(2*gamma));      % derivative for Newton solve
    end                                               % end branch decision

end                                                   % end pressure function
