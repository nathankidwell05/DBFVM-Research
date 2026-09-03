%% KT-D1V5 FVDBM Limiter Development Comparison %%
% This cleaner script compares the limiter versions developed during the project.
% Every case uses the same KT-D1V5 equilibrium, finite-volume transport,
% BGK collision, SSP-RK3 time stepping, grid, and Sod initial condition.
% Therefore, the plotted differences come from the spatial limiter.
%
% Limiter history:
%   1. Minmod: stable but more diffusive.
%   2. MC: sharper, but produced a temperature overshoot.
%   3. Population hybrid: each f_i independently chose MC or minmod.
%   4. Macroscopic hybrid: all f_i shared one troubled-cell decision.
%   5. Generalized minmod: one continuous theta value between minmod and MC.

clear;                                              % remove old variables
clc;                                                % clear the Command Window
close all;                                          % close old figures

%% Physical and Numerical Parameters
gamma = 1.4;                                        % specific-heat ratio for an air-like gas
b = 2/(gamma-1);                                    % internal degrees-of-freedom parameter; gamma=1.4 gives b=5
tau = 5.0e-5;                                       % BGK relaxation time
CFL = 0.10;                                         % transport stability factor
tEnd = 0.15;                                        % final comparison time
thetaFinal = 1.20;                                  % final generalized-minmod setting

Nx = 50000;                                         % production resolution; use 10000 for a faster comparison
Lx = 20.0;                                          % long domain keeps waves away from the boundaries
dx = Lx/Nx;                                         % finite-volume cell width
x = ((1:Nx)'-0.5)*dx;                               % cell-center positions
x0 = 0.5*Lx;                                        % initial diaphragm location

%% Sod Initial Condition
rhoL = 1.0;  uL = 0.0;  pL = 1.0;                  % state left of the diaphragm
rhoR = 0.125; uR = 0.0; pR = 0.1;                  % state right of the diaphragm

rho0 = rhoR*ones(Nx,1);                             % begin with the right density everywhere
u0 = uR*ones(Nx,1);                                 % begin with the right velocity everywhere
p0 = pR*ones(Nx,1);                                 % begin with the right pressure everywhere
leftCells = x < x0;                                 % identify cells left of the diaphragm
rho0(leftCells) = rhoL;                             % apply left density
u0(leftCells) = uL;                                 % apply left velocity
p0(leftCells) = pL;                                 % apply left pressure
T0 = p0./rho0;                                      % ideal-gas temperature, p=rho*T

%% KT-D1V5 Discrete-Velocity Model
v1 = 1.0;                                           % slow discrete molecular speed
v2 = 3.0;                                           % fast discrete molecular speed
eta0 = 1.75;                                        % internal-energy parameter on the rest population
c = [-v2,-v1,0,v1,v2];                              % five discrete velocities
h2 = [0,0,eta0^2,0,0];                              % squared internal-energy values
fInitial = build_KT_equilibrium(rho0,u0,T0,b,v1,v2,eta0,c); % start all populations at equilibrium

maxSpeed = max(abs(c));                              % fastest discrete population
dtAdvection = CFL*dx/maxSpeed;                       % finite-volume CFL limit
dtCollision = 0.25*tau;                              % explicit BGK collision limit
dtBase = min(dtAdvection,dtCollision);               % common time step used by every case

%% Limiter Cases in Historical Order
limiterTypes = ["minmod","mc","populationHybrid","macroscopicHybrid","gminmod"];
limiterLabels = {'Minmod','MC','Population hybrid','Macroscopic hybrid', ...
    sprintf('Generalized minmod, \\theta=%.2f',thetaFinal)};
numberOfCases = numel(limiterTypes);                 % total solver versions being compared
results = repmat(struct('rho',[],'u',[],'p',[],'T',[],'steps',0,'runtime',0),numberOfCases,1);

fprintf('KT-D1V5 finite-volume limiter comparison\n');
fprintf('Nx=%d, Lx=%.2f, dx=%.3e, dt=%.3e, tEnd=%.3f\n',Nx,Lx,dx,dtBase,tEnd);

for caseNumber = 1:numberOfCases                     % run every historical limiter with identical settings
    fprintf('\nRunning %d/%d: %s\n',caseNumber,numberOfCases,limiterLabels{caseNumber});
    results(caseNumber) = solve_D1V5_case(fInitial,b,tau,v1,v2,eta0,c,h2,dx, ...
        dtBase,tEnd,limiterTypes(caseNumber),thetaFinal); % solve one complete Sod case
    fprintf('Finished in %.2f s using %d steps.\n', ...
        results(caseNumber).runtime,results(caseNumber).steps); % concise progress report
end

%% Exact Sod Solution and Quantitative Comparison
exact = sod_exact_solution(x,tEnd,x0,rhoL,uL,pL,rhoR,uR,pR,gamma); % exact Euler benchmark
waveMask = x >= (x0-0.30) & x <= (x0+0.40);        % region containing rarefaction, contact, and shock

globalL1rho = zeros(numberOfCases,1);                % whole-domain average density error
globalL1u = zeros(numberOfCases,1);                  % whole-domain average velocity error
globalL1p = zeros(numberOfCases,1);                  % whole-domain average pressure error
globalL1T = zeros(numberOfCases,1);                  % whole-domain average temperature error
waveL1rho = zeros(numberOfCases,1);                  % density error only around the waves
waveL1u = zeros(numberOfCases,1);                    % velocity error only around the waves
waveL1p = zeros(numberOfCases,1);                    % pressure error only around the waves
waveL1T = zeros(numberOfCases,1);                    % temperature error only around the waves
temperaturePeakExcess = zeros(numberOfCases,1);      % numerical maximum T minus exact maximum T

exactTemperaturePeak = max(exact.T(waveMask));       % exact temperature peak in the comparison region

for caseNumber = 1:numberOfCases                     % calculate identical metrics for every limiter
    globalL1rho(caseNumber) = mean(abs(results(caseNumber).rho-exact.rho));
    globalL1u(caseNumber) = mean(abs(results(caseNumber).u-exact.u));
    globalL1p(caseNumber) = mean(abs(results(caseNumber).p-exact.p));
    globalL1T(caseNumber) = mean(abs(results(caseNumber).T-exact.T));

    waveL1rho(caseNumber) = mean(abs(results(caseNumber).rho(waveMask)-exact.rho(waveMask)));
    waveL1u(caseNumber) = mean(abs(results(caseNumber).u(waveMask)-exact.u(waveMask)));
    waveL1p(caseNumber) = mean(abs(results(caseNumber).p(waveMask)-exact.p(waveMask)));
    waveL1T(caseNumber) = mean(abs(results(caseNumber).T(waveMask)-exact.T(waveMask)));
    temperaturePeakExcess(caseNumber) = max(results(caseNumber).T(waveMask))-exactTemperaturePeak;
end

errorTable = table(string(limiterLabels(:)),globalL1rho,globalL1u,globalL1p,globalL1T, ...
    waveL1rho,waveL1u,waveL1p,waveL1T,temperaturePeakExcess, ...
    'VariableNames',{'Limiter','GlobalL1_rho','GlobalL1_u','GlobalL1_p','GlobalL1_T', ...
    'WaveL1_rho','WaveL1_u','WaveL1_p','WaveL1_T','T_peak_excess'}); % collect comparison metrics

fprintf('\nLimiter error comparison\n');
disp(errorTable);                                     % print all errors in one readable table

%% Large Four-Panel Comparison Figure
colors = [0.0000 0.4470 0.7410; ...                 % blue: minmod
          0.8500 0.1500 0.1500; ...                 % red: MC
          0.9290 0.6940 0.1250; ...                 % orange: population hybrid
          0.4940 0.1840 0.5560; ...                 % purple: macroscopic hybrid
          0.0000 0.5500 0.1500];                    % green: generalized minmod
lineStyles = {'-', '-', '-.', ':', '--'};            % generalized minmod is dashed as requested
lineWidths = [1.4,1.4,1.7,2.0,2.6];                 % emphasize later improvements
fieldNames = {'rho','u','p','T'};                    % fields placed in the four panels
yLabels = {'Density, \rho','Velocity, u','Pressure, p','Temperature, T'};

figure('Name','KT-D1V5 Limiter History','Color','w','Position',[60 50 1500 1050]);
layout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact'); % large readable 2-by-2 layout

for panel = 1:4                                      % build density, velocity, pressure, and temperature panels
    ax = nexttile(layout);                            % select the next large panel
    hold(ax,'on');                                    % allow every limiter and exact solution on the same axes

    exactHandle = plot(ax,x,exact.(fieldNames{panel}),'k-','LineWidth',2.8); % exact solution is solid black
    caseHandles = gobjects(numberOfCases,1);          % store handles so legend order is controlled

    for caseNumber = 1:numberOfCases                 % overlay the five limiter histories
        caseHandles(caseNumber) = plot(ax,x,results(caseNumber).(fieldNames{panel}), ...
            'Color',colors(caseNumber,:),'LineStyle',lineStyles{caseNumber}, ...
            'LineWidth',lineWidths(caseNumber));      % generalized minmod appears as dashed green
    end

    xlim(ax,[x0-0.30,x0+0.40]);                      % zoom around the moving waves for easier comparison
    ylabel(ax,yLabels{panel},'FontSize',15,'FontWeight','bold');
    xlabel(ax,'x','FontSize',15,'FontWeight','bold');
    title(ax,yLabels{panel},'FontSize',17,'FontWeight','bold');
    grid(ax,'on');                                    % make wave locations easier to compare
    box(ax,'on');                                     % draw a clear frame around each panel
    ax.FontSize = 14;                                 % enlarge tick labels
    ax.LineWidth = 1.1;                               % make axes easier to see

    if panel == 1                                    % use one shared legend to avoid covering every panel
        legend(ax,[exactHandle;caseHandles],[{'Exact'},limiterLabels], ...
            'Location','best','FontSize',12,'Box','on');
    end
end

title(layout,'KT-D1V5 FVDBM: Limiter Development History','FontSize',22,'FontWeight','bold');

%% Enlarged Temperature-Overshoot Figure
figure('Name','Temperature Overshoot Comparison','Color','w','Position',[90 100 1500 720]);
hold on;
exactTemperatureHandle = plot(x,exact.T,'k-','LineWidth',3.2); % solid black exact temperature
temperatureHandles = gobjects(numberOfCases,1);                % handles for limiter temperature curves

for caseNumber = 1:numberOfCases                               % plot every historical temperature result
    temperatureHandles(caseNumber) = plot(x,results(caseNumber).T, ...
        'Color',colors(caseNumber,:),'LineStyle',lineStyles{caseNumber}, ...
        'LineWidth',lineWidths(caseNumber)+0.4);
end

xlim([x0-0.30,x0+0.40]);                                      % focus on the shock/contact region
yline(exactTemperaturePeak,'k:','Exact peak','LineWidth',1.6,'FontSize',13); % reference level for overshoot
xlabel('x','FontSize',17,'FontWeight','bold');
ylabel('Temperature, T','FontSize',17,'FontWeight','bold');
title(sprintf('Temperature Comparison: Generalized Minmod \\theta = %.2f',thetaFinal), ...
    'FontSize',22,'FontWeight','bold');
legend([exactTemperatureHandle;temperatureHandles],[{'Exact'},limiterLabels], ...
    'Location','best','FontSize',13,'Box','on');
grid on;
box on;
set(gca,'FontSize',15,'LineWidth',1.2);                         % enlarge axes and tick labels

%% Save Comparison Data
outputFile = fullfile(fileparts(mfilename('fullpath')),'d1v5_limiter_history_results.mat');
save(outputFile,'x','x0','results','exact','errorTable','waveMask','Nx','Lx','dx', ...
    'tau','CFL','tEnd','thetaFinal','dtBase');                   % save everything needed for later plots
fprintf('Comparison data saved to:\n%s\n',outputFile);

%% Local Solver Functions
function result = solve_D1V5_case(fInitial,b,tau,v1,v2,eta0,c,h2,dx,dtBase,tEnd,limiterType,theta)
%SOLVE_D1V5_CASE advances one limiter case without retries, repairs, or live plotting.

    f = fInitial;                                       % give this case its own copy of the initial populations
    time = 0.0;                                         % initial physical time
    step = 0;                                           % accepted time-step count
    caseTimer = tic;                                    % measure runtime of this limiter

    while time < tEnd                                   % continue until the requested physical time
        dt = min(dtBase,tEnd-time);                     % use the stable step, shortened only at the end

        rhs1 = compute_rhs(f,b,tau,v1,v2,eta0,c,h2,dx,limiterType,theta); % RK3 stage 1
        f1 = apply_transmissive_bc(f+dt*rhs1);          % first trial state

        rhs2 = compute_rhs(f1,b,tau,v1,v2,eta0,c,h2,dx,limiterType,theta); % RK3 stage 2
        f2 = apply_transmissive_bc(0.75*f+0.25*(f1+dt*rhs2)); % second blended state

        rhs3 = compute_rhs(f2,b,tau,v1,v2,eta0,c,h2,dx,limiterType,theta); % RK3 stage 3
        f = apply_transmissive_bc((1/3)*f+(2/3)*(f2+dt*rhs3)); % final SSP-RK3 state

        time = time+dt;                                 % advance physical time
        step = step+1;                                  % count the completed step
    end

    [rho,u,T,p] = recover_macros(f,b,c,h2);             % recover final macroscopic solution

    if any(~isfinite(rho)) || any(~isfinite(T)) || any(rho <= 0) || any(T <= 0)
        error('%s produced a nonphysical final state.',char(limiterType)); % one final honest validity check
    end

    result.rho = rho;                                   % store density
    result.u = u;                                       % store velocity
    result.p = p;                                       % store pressure
    result.T = T;                                       % store temperature
    result.steps = step;                                % store time-step count
    result.runtime = toc(caseTimer);                    % store wall-clock runtime
end

function rhs = compute_rhs(f,b,tau,v1,v2,eta0,c,h2,dx,limiterType,theta)
%COMPUTE_RHS combines finite-volume transport and BGK collision.

    [rho,u,T,p] = recover_macros(f,b,c,h2);             % recover state needed by equilibrium and macro sensor
    feq = build_KT_equilibrium(rho,u,T,b,v1,v2,eta0,c); % local KT-D1V5 equilibrium
    collision = (feq-f)/tau;                            % BGK relaxation toward equilibrium

    if limiterType == "macroscopicHybrid"              % only this historical case uses the shared sensor
        troubledCells = macroscopic_troubled_cells(rho,p,T); % one mask shared by all populations
    else
        troubledCells = false(size(rho));               % unused placeholder for the other limiters
    end

    fluxDivergence = build_flux_divergence(f,c,dx,limiterType,troubledCells,theta); % FVM transport
    rhs = -fluxDivergence+collision;                    % discrete Boltzmann equation right-hand side
end

function fluxDivergence = build_flux_divergence(f,c,dx,limiterType,troubledCells,theta)
%BUILD_FLUX_DIVERGENCE reconstructs population values and computes face fluxes.

    [Nx,Q] = size(f);                                   % number of physical cells and populations
    fluxDivergence = zeros(Nx,Q);                       % transport contribution for each f_i

    for q = 1:Q                                        % reconstruct and transport one population at a time
        fq = f(:,q);                                    % current population across the domain
        slope = zeros(Nx,1);                            % limited cell slope
        dL = fq(2:Nx-1)-fq(1:Nx-2);                    % change from the left neighbor
        dR = fq(3:Nx)-fq(2:Nx-1);                      % change toward the right neighbor

        if limiterType == "minmod"
            slope(2:Nx-1) = minmod_limiter(dL,dR);     % conservative, diffusive baseline
        elseif limiterType == "mc"
            slope(2:Nx-1) = generalized_minmod(dL,dR,2.0); % full MC limiter
        elseif limiterType == "populationHybrid"
            slope(2:Nx-1) = population_hybrid_limiter(dL,dR); % old independent f_i switching
        elseif limiterType == "macroscopicHybrid"
            limitedSlope = generalized_minmod(dL,dR,2.0); % begin with MC slopes
            safeSlope = minmod_limiter(dL,dR);         % separately build minmod slopes
            mask = troubledCells(2:Nx-1);              % shared physical-cell warning mask
            limitedSlope(mask) = safeSlope(mask);      % all populations switch together in warned cells
            slope(2:Nx-1) = limitedSlope;              % store shared-hybrid result
        elseif limiterType == "gminmod"
            slope(2:Nx-1) = generalized_minmod(dL,dR,theta); % final continuously controlled limiter
        else
            error('Unknown limiter: %s',char(limiterType));
        end

        leftCell = 1:Nx-1;                             % cell immediately left of each interior face
        rightCell = 2:Nx;                              % cell immediately right of each interior face
        faceFromLeft = fq(leftCell)+0.5*slope(leftCell); % MUSCL state reaching face from left
        faceFromRight = fq(rightCell)-0.5*slope(rightCell); % MUSCL state reaching face from right
        faceFlux = zeros(Nx+1,1);                      % Nx cells have Nx+1 faces

        if c(q) > 0                                    % positive population travels left to right
            faceFlux(2:Nx) = c(q)*faceFromLeft;        % upwind state comes from left cell
        elseif c(q) < 0                                % negative population travels right to left
            faceFlux(2:Nx) = c(q)*faceFromRight;       % upwind state comes from right cell
        end

        faceFlux(1) = faceFlux(2);                     % transmissive left flux boundary
        faceFlux(end) = faceFlux(end-1);               % transmissive right flux boundary
        fluxDivergence(:,q) = (faceFlux(2:end)-faceFlux(1:end-1))/dx; % net outgoing face flux
    end
end

function slope = minmod_limiter(dL,dR)
%MINMOD_LIMITER keeps the smallest one-sided slope when both sides agree.

    slope = zeros(size(dL));                           % zero slope where signs disagree
    positive = dL > 0 & dR > 0;                       % monotonically increasing cells
    negative = dL < 0 & dR < 0;                       % monotonically decreasing cells
    slope(positive) = min(dL(positive),dR(positive));  % smallest safe positive slope
    slope(negative) = max(dL(negative),dR(negative));  % negative slope closest to zero
end

function slope = generalized_minmod(dL,dR,theta)
%GENERALIZED_MINMOD moves continuously from minmod (theta=1) to MC (theta=2).

    leftCandidate = theta*dL;                          % slope allowed by left difference
    middleCandidate = 0.5*(dL+dR);                    % centered slope using both neighbors
    rightCandidate = theta*dR;                         % slope allowed by right difference
    slope = minmod_three(leftCandidate,middleCandidate,rightCandidate); % choose safe candidate
end

function slope = minmod_three(a,b,candidateC)
%MINMOD_THREE returns the smallest-magnitude candidate when all signs agree.

    slope = zeros(size(a));                            % zero prevents a new extremum when signs disagree
    positive = a > 0 & b > 0 & candidateC > 0;        % all candidates rise
    negative = a < 0 & b < 0 & candidateC < 0;        % all candidates fall
    slope(positive) = min([a(positive),b(positive),candidateC(positive)],[],2);
    slope(negative) = max([a(negative),b(negative),candidateC(negative)],[],2);
end

function slope = population_hybrid_limiter(dL,dR)
%POPULATION_HYBRID_LIMITER reproduces the first unsuccessful MC/minmod hybrid.

    slopeMC = generalized_minmod(dL,dR,2.0);           % sharp MC choice
    slopeMinmod = minmod_limiter(dL,dR);               % safer minmod choice
    small = 1.0e-14;                                   % avoid division by zero in flat regions
    curvatureRatio = abs(dR-dL)./(abs(dR)+abs(dL)+small); % large when neighboring slopes disagree strongly
    slopeSizeRatio = max(abs(dL),abs(dR))./(min(abs(dL),abs(dR))+small); % large for one-sided jumps
    useMinmod = curvatureRatio > 0.35 | slopeSizeRatio > 3.0; % original per-population warning rule
    slope = slopeMC;                                   % MC in cells judged smooth
    slope(useMinmod) = slopeMinmod(useMinmod);         % minmod in warned cells
end

function troubled = macroscopic_troubled_cells(rho,p,T)
%MACROSCOPIC_TROUBLED_CELLS reproduces the later shared rho/p/T sensor.

    Nx = length(rho);                                  % physical-cell count
    troubled = false(Nx,1);                            % begin with no warned cells
    fields = {rho,p,T};                                % density sees contact; pressure sees shock; T sees thermal jump
    small = 1.0e-14;                                   % avoid division by zero

    for fieldNumber = 1:numel(fields)                  % combine warnings from all three fields
        phi = fields{fieldNumber};                     % current macroscopic field
        dL = phi(2:Nx-1)-phi(1:Nx-2);                  % left change
        dR = phi(3:Nx)-phi(2:Nx-1);                    % right change
        relativeJump = max(abs(dL),abs(dR))./(abs(phi(2:Nx-1))+small); % dimensionless jump size
        roughness = abs(dR-dL)./(abs(dR)+abs(dL)+small); % sudden change in slope
        warning = (roughness > 0.20 & relativeJump > 0.0025) | relativeJump > 0.08;
        troubled(2:Nx-1) = troubled(2:Nx-1) | warning; % share the warning among all populations
    end

    original = troubled;                               % preserve sensor result before adding neighbor protection
    troubled(2:end) = troubled(2:end) | original(1:end-1); % protect one cell to the right
    troubled(1:end-1) = troubled(1:end-1) | original(2:end); % protect one cell to the left
end

function feq = build_KT_equilibrium(rho,u,T,b,v1,v2,eta0,c)
%BUILD_KT_EQUILIBRIUM evaluates the published KT-D1V5 equilibrium.

    feq = zeros(length(rho),length(c));                 % equilibrium population storage
    u2 = u.^2;                                         % squared flow velocity
    Arest = ((b-1)/eta0^2).*T;                         % rest-population coefficient
    Aslow = (-v2^2+(((b-1)*v2^2/eta0^2)+1).*T+u2)./(2*(v1^2-v2^2));
    Afast = (-v1^2+(((b-1)*v1^2/eta0^2)+1).*T+u2)./(2*(v2^2-v1^2));
    Bslow = (-v2^2+(b+2).*T+u2)./(2*v1^2*(v1^2-v2^2));
    Bfast = (-v1^2+(b+2).*T+u2)./(2*v2^2*(v2^2-v1^2));

    for q = 1:length(c)                                % construct the five equilibrium populations
        if c(q) == 0
            feq(:,q) = rho.*Arest;                     % stationary population carries internal energy
        elseif abs(c(q)) == v1
            feq(:,q) = rho.*(Aslow+Bslow.*u.*c(q));    % slow left/right populations
        else
            feq(:,q) = rho.*(Afast+Bfast.*u.*c(q));    % fast left/right populations
        end
    end
end

function [rho,u,T,p] = recover_macros(f,b,c,h2)
%RECOVER_MACROS obtains physical fields from moments of the five populations.

    rho = sum(f,2);                                    % zeroth moment: density
    u = (f*c')./rho;                                   % first moment divided by density: velocity
    energy = f*(c.^2+h2)';                             % kinetic plus internal-energy moment
    T = (energy./rho-u.^2)/b;                          % temperature recovered from total energy
    p = rho.*T;                                        % nondimensional ideal-gas pressure
end

function f = apply_transmissive_bc(f)
%APPLY_TRANSMISSIVE_BC copies the nearest interior populations at each end.

    f(1,:) = f(2,:);                                   % zero-gradient left boundary
    f(end,:) = f(end-1,:);                             % zero-gradient right boundary
end

function exact = sod_exact_solution(x,t,x0,rhoL,uL,pL,rhoR,uR,pR,gamma)
%SOD_EXACT_SOLUTION computes the exact Euler Riemann solution used only for validation.

    if t <= 0
        exact.rho = rhoR*ones(size(x));
        exact.u = uR*ones(size(x));
        exact.p = pR*ones(size(x));
        exact.rho(x < x0) = rhoL;
        exact.u(x < x0) = uL;
        exact.p(x < x0) = pL;
        exact.T = exact.p./exact.rho;
        return;
    end

    aL = sqrt(gamma*pL/rhoL);                          % left sound speed
    aR = sqrt(gamma*pR/rhoR);                          % right sound speed
    pStar = max(0.5*(pL+pR),1.0e-10);                 % initial star-pressure estimate

    for iteration = 1:50                              % Newton iteration for star pressure
        [fL,dfL] = pressure_function(pStar,rhoL,pL,aL,gamma);
        [fR,dfR] = pressure_function(pStar,rhoR,pR,aR,gamma);
        pNew = max(pStar-(fL+fR+uR-uL)/(dfL+dfR),1.0e-10);
        if abs(pNew-pStar)/(pNew+pStar) < 1.0e-10
            pStar = pNew;
            break;
        end
        pStar = pNew;
    end

    [fL,~] = pressure_function(pStar,rhoL,pL,aL,gamma);
    [fR,~] = pressure_function(pStar,rhoR,pR,aR,gamma);
    uStar = 0.5*(uL+uR+fR-fL);                        % star-region velocity
    xiAll = (x-x0)/t;                                  % similarity coordinate
    rhoExact = zeros(size(x));
    uExact = zeros(size(x));
    pExact = zeros(size(x));

    for k = 1:length(x)
        xi = xiAll(k);

        if xi <= uStar
            if pStar > pL                              % left shock
                SL = uL-aL*sqrt((gamma+1)/(2*gamma)*pStar/pL+(gamma-1)/(2*gamma));
                if xi <= SL
                    rhoExact(k)=rhoL; uExact(k)=uL; pExact(k)=pL;
                else
                    rhoExact(k)=rhoL*((pStar/pL+(gamma-1)/(gamma+1))/((gamma-1)/(gamma+1)*pStar/pL+1));
                    uExact(k)=uStar; pExact(k)=pStar;
                end
            else                                       % left rarefaction
                aStarL = aL*(pStar/pL)^((gamma-1)/(2*gamma));
                head = uL-aL;
                tail = uStar-aStarL;
                if xi <= head
                    rhoExact(k)=rhoL; uExact(k)=uL; pExact(k)=pL;
                elseif xi > tail
                    rhoExact(k)=rhoL*(pStar/pL)^(1/gamma); uExact(k)=uStar; pExact(k)=pStar;
                else
                    uFan=2/(gamma+1)*(aL+0.5*(gamma-1)*uL+xi);
                    aFan=2/(gamma+1)*(aL+0.5*(gamma-1)*(uL-xi));
                    rhoExact(k)=rhoL*(aFan/aL)^(2/(gamma-1));
                    uExact(k)=uFan;
                    pExact(k)=pL*(aFan/aL)^(2*gamma/(gamma-1));
                end
            end
        else
            if pStar > pR                              % right shock
                SR = uR+aR*sqrt((gamma+1)/(2*gamma)*pStar/pR+(gamma-1)/(2*gamma));
                if xi >= SR
                    rhoExact(k)=rhoR; uExact(k)=uR; pExact(k)=pR;
                else
                    rhoExact(k)=rhoR*((pStar/pR+(gamma-1)/(gamma+1))/((gamma-1)/(gamma+1)*pStar/pR+1));
                    uExact(k)=uStar; pExact(k)=pStar;
                end
            else                                       % right rarefaction
                aStarR = aR*(pStar/pR)^((gamma-1)/(2*gamma));
                head = uR+aR;
                tail = uStar+aStarR;
                if xi >= head
                    rhoExact(k)=rhoR; uExact(k)=uR; pExact(k)=pR;
                elseif xi <= tail
                    rhoExact(k)=rhoR*(pStar/pR)^(1/gamma); uExact(k)=uStar; pExact(k)=pStar;
                else
                    uFan=2/(gamma+1)*(-aR+0.5*(gamma-1)*uR+xi);
                    aFan=2/(gamma+1)*(aR-0.5*(gamma-1)*(uR-xi));
                    rhoExact(k)=rhoR*(aFan/aR)^(2/(gamma-1));
                    uExact(k)=uFan;
                    pExact(k)=pR*(aFan/aR)^(2*gamma/(gamma-1));
                end
            end
        end
    end

    exact.rho = rhoExact;
    exact.u = uExact;
    exact.p = pExact;
    exact.T = pExact./rhoExact;
end

function [pressureValue,derivative] = pressure_function(p,rhoK,pK,aK,gamma)
%PRESSURE_FUNCTION supplies the shock/rarefaction function for the exact solver.

    if p > pK                                          % shock branch
        A = 2/((gamma+1)*rhoK);
        B = (gamma-1)/(gamma+1)*pK;
        pressureValue = (p-pK)*sqrt(A/(p+B));
        derivative = sqrt(A/(p+B))*(1-0.5*(p-pK)/(p+B));
    else                                               % rarefaction branch
        pressureValue = 2*aK/(gamma-1)*((p/pK)^((gamma-1)/(2*gamma))-1);
        derivative = (1/(rhoK*aK))*(p/pK)^(-(gamma+1)/(2*gamma));
    end
end
