%% Two Dimensional Sod Shock Tube Using KT-D2V9, FVM, MUSCL, and Minmod %%
% This program is a teaching/proof-of-concept solver for compressible flow.
% The physical problem is a Sod shock tube placed in a 2D domain.
% The flow is initialized so it only varies in x, but the solver itself is 2D.
% This means f is stored as f(i,j,q), where i is x, j is y, and q is velocity population.
% The kinetic model used here is the published Kataoka-Tsutahara style D2V9 model.
% The transport/advection is solved with a finite-volume method instead of lattice streaming.
% The face values are reconstructed with MUSCL.
% The slopes used by MUSCL are limited with minmod to reduce nonphysical oscillations near shocks.
% The goal is not maximum speed yet.
% The goal is to make every command readable and connect it back to the finite-volume DBM idea.

clear;                 % clears variables from the workspace so old values do not contaminate the run
clc;                   % clears the command window so the printed output is easier to read
close all;             % closes old figures so only this run's plots show up

%% Physical Parameters
gamma = 1.4;           % ratio of specific heats for air-like ideal gas
b = 2/(gamma-1);       % internal degrees of freedom parameter used by KT; gamma=1.4 gives b=5
tau = 1.0e-4;          % BGK relaxation time; smaller tau means faster relaxation toward equilibrium
CFL = 0.35;            % CFL number controlling the explicit finite-volume time step
tEnd = 0.15;           % final physical time; keep short so waves do not hit boundaries too strongly

%% Domain
Nx = 400;              % number of finite-volume cells in the x direction
Ny = 40;               % number of finite-volume cells in the y direction
Lx = 1.0;              % length of the shock tube in the x direction
Ly = 0.10;             % height of the 2D domain in the y direction
dx = Lx/Nx;            % width of each cell in x
dy = Ly/Ny;            % height of each cell in y
x = ((1:Nx)'-0.5)*dx;  % x-location of each cell center as a column vector
y = ((1:Ny)-0.5)*dy;   % y-location of each cell center as a row vector
x0 = 0.5*Lx;           % diaphragm location where the initial left/right states meet

%% Initial Conditions
rhoL = 1.0;            % left density
uxL = 0.0;             % left x-velocity
uyL = 0.0;             % left y-velocity
pL = 1.0;              % left pressure

rhoR = 0.125;          % right density
uxR = 0.0;             % right x-velocity
uyR = 0.0;             % right y-velocity
pR = 0.1;              % right pressure for the classic Sod shock tube

rho = zeros(Nx,Ny);    % storage for density in every 2D cell
ux = zeros(Nx,Ny);     % storage for x-velocity in every 2D cell
uy = zeros(Nx,Ny);     % storage for y-velocity in every 2D cell
p = zeros(Nx,Ny);      % storage for pressure in every 2D cell
T = zeros(Nx,Ny);      % storage for temperature in every 2D cell

for i = 1:Nx                                      % loop over all x cells
    for j = 1:Ny                                  % loop over all y cells
        if x(i) < x0                              % cells left of diaphragm get the left state
            rho(i,j) = rhoL;                      % assign left density
            ux(i,j) = uxL;                        % assign left x-velocity
            uy(i,j) = uyL;                        % assign left y-velocity
            p(i,j) = pL;                          % assign left pressure
        else                                      % cells right of diaphragm get the right state
            rho(i,j) = rhoR;                      % assign right density
            ux(i,j) = uxR;                        % assign right x-velocity
            uy(i,j) = uyR;                        % assign right y-velocity
            p(i,j) = pR;                          % assign right pressure
        end                                       % end left/right initial condition decision
        T(i,j) = p(i,j)/rho(i,j);                 % ideal-gas temperature using p=rho*T in nondimensional units
    end                                           % end y-cell loop
end                                               % end x-cell loop

%% KT-D2V9 Velocity Set
% D2 means the physical velocity has x and y components.
% V9 means there are nine discrete molecular velocities.
% This is not lattice streaming.
% These velocities are used by the DBM model while finite volume moves f through physical cell faces.

v1 = 1.0;                                         % slow KT speed level
v2 = 3.0;                                         % fast KT speed level
eta0 = 2.0;                                      % internal-energy parameter attached to the rest population

cx = [0, v1, -v1, 0, 0, v2/sqrt(2), -v2/sqrt(2), -v2/sqrt(2), v2/sqrt(2)];  % x velocity for each population
cy = [0, 0, 0, v1, -v1, v2/sqrt(2), v2/sqrt(2), -v2/sqrt(2), -v2/sqrt(2)];  % y velocity for each population
eta2 = [eta0^2, 0, 0, 0, 0, 0, 0, 0, 0];                                  % internal-energy backpack squared
Q = length(cx);                                    % number of discrete velocity populations

%% Build Initial Equilibrium
% The simulation evolves f.
% At t=0, we start f at equilibrium so the initial macroscopic state is exactly represented.

f = build_KT_D2V9_equilibrium(rho,ux,uy,T,b,v1,v2,eta0,cx,cy);     % initialize the nine populations at equilibrium

%% Time Step
maxSpeedX = max(abs(cx));                            % largest discrete speed in the x direction
maxSpeedY = max(abs(cy));                            % largest discrete speed in the y direction
dt_adv_x = CFL*dx/maxSpeedX;                         % x-direction advection time-step restriction
dt_adv_y = CFL*dy/maxSpeedY;                         % y-direction advection time-step restriction
dt_adv = min(dt_adv_x,dt_adv_y);                     % most restrictive advection time step
dt_col = 0.25*tau;                                   % explicit collision time-step restriction
dt = min(dt_adv,dt_col);                             % actual time step must satisfy both transport and collision limits

time = 0.0;                                          % current simulation time
step = 0;                                            % counter for the number of time steps taken
plotEvery = 50;                                      % update the live plot every this many steps

fprintf('Running KT-D2V9 finite-volume DBM Sod shock tube...\n');     % print run description
fprintf('Grid: Nx=%d, Ny=%d, dx=%.3e, dy=%.3e\n',Nx,Ny,dx,dy);        % print grid information
fprintf('dt_adv=%.3e, dt_col=%.3e, dt=%.3e\n',dt_adv,dt_col,dt);      % print time-step information

%% Create Live Plot Handles
figure('Name','KT-D2V9 MUSCL-Minmod Sod Shock Tube');                 % create one figure for live visualization
iyMid = ceil(Ny/2);                                                    % centerline index in y used for 1D plotting

subplot(4,1,1);                                                        % first subplot for density
hRho = plot(x,rho(:,iyMid),'b-','LineWidth',1.5);                      % create density line handle
ylabel('\rho');                                                        % label density axis
grid on;                                                              % turn grid on

subplot(4,1,2);                                                        % second subplot for x velocity
hUx = plot(x,ux(:,iyMid),'b-','LineWidth',1.5);                        % create velocity line handle
ylabel('u_x');                                                         % label velocity axis
grid on;                                                              % turn grid on

subplot(4,1,3);                                                        % third subplot for pressure
hP = plot(x,p(:,iyMid),'b-','LineWidth',1.5);                          % create pressure line handle
ylabel('p');                                                           % label pressure axis
grid on;                                                              % turn grid on

subplot(4,1,4);                                                        % fourth subplot for temperature
hT = plot(x,T(:,iyMid),'b-','LineWidth',1.5);                          % create temperature line handle
ylabel('T');                                                           % label temperature axis
xlabel('x');                                                           % label x axis
grid on;                                                              % turn grid on

%% Main Time Loop
while time < tEnd                                      % continue updating until final time is reached

    if time + dt > tEnd                                % check whether the next step would go past final time
        dt = tEnd - time;                              % shrink final step so simulation stops exactly on tEnd
    end                                                % end final-step adjustment

    step = step + 1;                                   % count the current attempted time step

    %% Recover Macroscopic Variables
    rho = sum(f,3);                                    % density is the zeroth moment: sum over all populations
    mx = zeros(Nx,Ny);                                 % storage for x momentum
    my = zeros(Nx,Ny);                                 % storage for y momentum
    energy = zeros(Nx,Ny);                             % storage for total kinetic/internal energy moment

    for q = 1:Q                                        % loop through each population
        mx = mx + f(:,:,q)*cx(q);                      % add q contribution to x momentum
        my = my + f(:,:,q)*cy(q);                      % add q contribution to y momentum
        energy = energy + f(:,:,q)*(cx(q)^2 + cy(q)^2 + eta2(q));  % add q contribution to energy moment
    end                                                % end moment recovery population loop

    ux = mx ./ rho;                                    % x velocity equals x momentum divided by density
    uy = my ./ rho;                                    % y velocity equals y momentum divided by density
    u2 = ux.^2 + uy.^2;                                % squared macroscopic speed
    T = (energy ./ rho - u2) / b;                      % temperature recovered from energy moment
    if any(rho(:) <= 0) || any(T(:) <= 0)              % check for nonphysical density or temperature
        error('Nonphysical state at step %d, time %.6e',step,time);  % stop if the solver has collapsed
    end                                                % end physical-state check

    %% Build Equilibrium
    feq = build_KT_D2V9_equilibrium(rho,ux,uy,T,b,v1,v2,eta0,cx,cy);  % build local KT equilibrium from current macros

    %% Collision
    collision = (feq - f) / tau;                       % BGK collision relaxes f toward feq over time scale tau

    %% Build FVM Flux With MUSCL-Minmod Reconstruction
    fluxDiv = zeros(Nx,Ny,Q);                          % storage for d(cx*f)/dx + d(cy*f)/dy for every cell and q

    for q = 1:Q                                        % loop over each discrete velocity population

        fq = f(:,:,q);                                 % pull out one population over the whole 2D grid
        slopeX = zeros(Nx,Ny);                         % storage for MUSCL slopes in the x direction
        slopeY = zeros(Nx,Ny);                         % storage for MUSCL slopes in the y direction

        %% Compute X Slopes
        for i = 2:Nx-1                                 % only interior x cells have both left and right neighbors
            for j = 1:Ny                               % loop through every y row
                dL = fq(i,j) - fq(i-1,j);              % left difference in x
                dR = fq(i+1,j) - fq(i,j);              % right difference in x
                slopeX(i,j) = minmod(dL,dR);           % safe x slope chosen by minmod limiter
            end                                        % end y loop for x slopes
        end                                            % end x loop for x slopes

        %% Compute Y Slopes
        for i = 1:Nx                                   % loop through every x column
            for j = 2:Ny-1                             % only interior y cells have both south and north neighbors
                dS = fq(i,j) - fq(i,j-1);              % south difference in y
                dN = fq(i,j+1) - fq(i,j);              % north difference in y
                slopeY(i,j) = minmod(dS,dN);           % safe y slope chosen by minmod limiter
            end                                        % end y loop for y slopes
        end                                            % end x loop for y slopes

        %% East/West Face Fluxes
        FxFace = zeros(Nx+1,Ny);                       % x-face fluxes; Nx cells have Nx+1 vertical faces

        for face = 2:Nx                                % loop over interior x faces
            iL = face - 1;                             % index of cell on the left side of the face
            iR = face;                                 % index of cell on the right side of the face

            for j = 1:Ny                               % compute this x-face flux for every y row
                fLeftFace = fq(iL,j) + 0.5*slopeX(iL,j);   % left cell's reconstructed value at its east face
                fRightFace = fq(iR,j) - 0.5*slopeX(iR,j);  % right cell's reconstructed value at its west face

                if cx(q) > 0                           % positive x velocity means population moves west to east
                    FxFace(face,j) = cx(q)*fLeftFace;  % upwind side is the left cell
                elseif cx(q) < 0                       % negative x velocity means population moves east to west
                    FxFace(face,j) = cx(q)*fRightFace; % upwind side is the right cell
                else                                   % zero x velocity means no x-direction flux
                    FxFace(face,j) = 0.0;              % x flux is zero for this population
                end                                    % end x upwind decision
            end                                        % end y row loop for x faces
        end                                            % end interior x-face loop

        FxFace(1,:) = FxFace(2,:);                     % left boundary x-face copies nearest interior x-face
        FxFace(end,:) = FxFace(end-1,:);               % right boundary x-face copies nearest interior x-face

        %% North/South Face Fluxes
        FyFace = zeros(Nx,Ny+1);                       % y-face fluxes; Ny cells have Ny+1 horizontal faces

        for face = 2:Ny                                % loop over interior y faces
            jS = face - 1;                             % index of cell on the south side of the face
            jN = face;                                 % index of cell on the north side of the face

            for i = 1:Nx                               % compute this y-face flux for every x column
                fSouthFace = fq(i,jS) + 0.5*slopeY(i,jS);  % south cell's reconstructed value at its north face
                fNorthFace = fq(i,jN) - 0.5*slopeY(i,jN);  % north cell's reconstructed value at its south face

                if cy(q) > 0                           % positive y velocity means population moves south to north
                    FyFace(i,face) = cy(q)*fSouthFace; % upwind side is the south cell
                elseif cy(q) < 0                       % negative y velocity means population moves north to south
                    FyFace(i,face) = cy(q)*fNorthFace; % upwind side is the north cell
                else                                   % zero y velocity means no y-direction flux
                    FyFace(i,face) = 0.0;              % y flux is zero for this population
                end                                    % end y upwind decision
            end                                        % end x column loop for y faces
        end                                            % end interior y-face loop

        FyFace(:,1) = FyFace(:,2);                     % south boundary y-face copies nearest interior y-face
        FyFace(:,end) = FyFace(:,end-1);               % north boundary y-face copies nearest interior y-face

        %% Flux Divergence
        for i = 1:Nx                                   % loop over all x cells
            for j = 1:Ny                               % loop over all y cells
                divX = (FxFace(i+1,j) - FxFace(i,j)) / dx;     % east flux minus west flux divided by dx
                divY = (FyFace(i,j+1) - FyFace(i,j)) / dy;     % north flux minus south flux divided by dy
                fluxDiv(i,j,q) = divX + divY;          % total finite-volume transport divergence for this q
            end                                        % end y-cell divergence loop
        end                                            % end x-cell divergence loop

    end                                                % end q population loop

    %% Combine Transport and Collision
    rhs = -fluxDiv + collision;                        % DBM RHS: negative transport divergence plus collision source

    %% Update Distribution Function
    f_new = f + dt*rhs;                                % forward-Euler update of all populations

    %% Transmissive Boundary Conditions on f
    f_new(1,:,:) = f_new(2,:,:);                       % left boundary copies nearest interior cell
    f_new(end,:,:) = f_new(end-1,:,:);                 % right boundary copies nearest interior cell
    f_new(:,1,:) = f_new(:,2,:);                       % south boundary copies nearest interior cell
    f_new(:,end,:) = f_new(:,end-1,:);                 % north boundary copies nearest interior cell

    %% Accept Update
    f = f_new;                                         % overwrite old populations with updated populations
    time = time + dt;                                  % advance physical time

    %% Live Visualization
    if mod(step,plotEvery) == 0 || time >= tEnd        % only update plot occasionally so MATLAB does not slow down badly
        rhoPlot = sum(f,3);                            % recover density for plotting
        mxPlot = zeros(Nx,Ny);                         % storage for plotted x momentum
        myPlot = zeros(Nx,Ny);                         % storage for plotted y momentum
        energyPlot = zeros(Nx,Ny);                     % storage for plotted energy moment

        for qPlot = 1:Q                                % loop over populations for plotting recovery
            mxPlot = mxPlot + f(:,:,qPlot)*cx(qPlot);  % add plotted x momentum contribution
            myPlot = myPlot + f(:,:,qPlot)*cy(qPlot);  % add plotted y momentum contribution
            energyPlot = energyPlot + f(:,:,qPlot)*(cx(qPlot)^2 + cy(qPlot)^2 + eta2(qPlot)); % add plotted energy
        end                                            % end plotting recovery loop

        uxPlot = mxPlot ./ rhoPlot;                    % recover plotted x velocity
        uyPlot = myPlot ./ rhoPlot;                    % recover plotted y velocity
        TPlot = (energyPlot ./ rhoPlot - (uxPlot.^2 + uyPlot.^2)) / b;  % recover plotted temperature
        pPlot = rhoPlot .* TPlot;                      % recover plotted pressure

        set(hRho,'YData',rhoPlot(:,iyMid));            % update density line without rebuilding the plot
        set(hUx,'YData',uxPlot(:,iyMid));              % update x-velocity line without rebuilding the plot
        set(hP,'YData',pPlot(:,iyMid));                % update pressure line without rebuilding the plot
        set(hT,'YData',TPlot(:,iyMid));                % update temperature line without rebuilding the plot
        sgtitle(sprintf('KT-D2V9 MUSCL-Minmod Sod Shock Tube: t = %.4f, step = %d',time,step)); % update figure title
        drawnow limitrate;                             % redraw efficiently without forcing MATLAB to waste time
    end                                                % end live visualization block

    if mod(step,100) == 0 || time >= tEnd              % print progress occasionally
        fprintf('step %d, time %.5f / %.5f\n',step,time,tEnd);  % print current step and time
    end                                                % end progress printing block

end                                                    % end main time loop

%% Final Macroscopic Recovery
rho = sum(f,3);                                        % final density from zeroth moment
mx = zeros(Nx,Ny);                                     % final x momentum storage
my = zeros(Nx,Ny);                                     % final y momentum storage
energy = zeros(Nx,Ny);                                 % final energy moment storage

for q = 1:Q                                            % loop over all populations
    mx = mx + f(:,:,q)*cx(q);                          % final x momentum contribution
    my = my + f(:,:,q)*cy(q);                          % final y momentum contribution
    energy = energy + f(:,:,q)*(cx(q)^2 + cy(q)^2 + eta2(q));  % final energy contribution
end                                                    % end final recovery population loop

ux = mx ./ rho;                                        % final x velocity
uy = my ./ rho;                                        % final y velocity
T = (energy ./ rho - (ux.^2 + uy.^2)) / b;             % final temperature
p = rho .* T;                                          % final pressure

%% Exact Sod Comparison on Centerline
exact = sod_exact_solution(x,tEnd,x0,rhoL,uxL,pL,rhoR,uxR,pR,gamma);  % compute exact 1D Euler Sod solution

L1_rho = mean(abs(rho(:,iyMid) - exact.rho));          % L1 density error on centerline
L1_ux = mean(abs(ux(:,iyMid) - exact.u));              % L1 x-velocity error on centerline
L1_p = mean(abs(p(:,iyMid) - exact.p));                % L1 pressure error on centerline
L1_T = mean(abs(T(:,iyMid) - exact.T));                % L1 temperature error on centerline
yVariation = max(std(rho,0,2));                        % y-variation check; should be near zero for planar Sod

fprintf('\nValidation Summary\n');                     % print validation header
fprintf('L1(rho) = %.6e\n',L1_rho);                    % print density error
fprintf('L1(ux)  = %.6e\n',L1_ux);                     % print velocity error
fprintf('L1(p)   = %.6e\n',L1_p);                      % print pressure error
fprintf('L1(T)   = %.6e\n',L1_T);                      % print temperature error
fprintf('max y-row rho std = %.6e\n',yVariation);      % print planar symmetry/y-uniformity diagnostic

%% Final Comparison Plot
figure('Name','KT-D2V9 DBM Compared with Exact Sod Solution');        % create final comparison figure

subplot(4,1,1);                                                        % density comparison subplot
plot(x,rho(:,iyMid),'b-','LineWidth',1.6);                             % plot DBM density
hold on;                                                              % keep exact curve on same axes
plot(x,exact.rho,'k--','LineWidth',1.4);                               % plot exact density
ylabel('\rho');                                                        % label density axis
legend('DBM','Exact');                                                 % label curves
grid on;                                                              % turn grid on

subplot(4,1,2);                                                        % velocity comparison subplot
plot(x,ux(:,iyMid),'b-','LineWidth',1.6);                              % plot DBM x velocity
hold on;                                                              % keep exact curve on same axes
plot(x,exact.u,'k--','LineWidth',1.4);                                 % plot exact velocity
ylabel('u_x');                                                         % label velocity axis
legend('DBM','Exact');                                                 % label curves
grid on;                                                              % turn grid on

subplot(4,1,3);                                                        % pressure comparison subplot
plot(x,p(:,iyMid),'b-','LineWidth',1.6);                               % plot DBM pressure
hold on;                                                              % keep exact curve on same axes
plot(x,exact.p,'k--','LineWidth',1.4);                                 % plot exact pressure
ylabel('p');                                                           % label pressure axis
legend('DBM','Exact');                                                 % label curves
grid on;                                                              % turn grid on

subplot(4,1,4);                                                        % temperature comparison subplot
plot(x,T(:,iyMid),'b-','LineWidth',1.6);                               % plot DBM temperature
hold on;                                                              % keep exact curve on same axes
plot(x,exact.T,'k--','LineWidth',1.4);                                 % plot exact temperature
ylabel('T');                                                           % label temperature axis
xlabel('x');                                                           % label x axis
legend('DBM','Exact');                                                 % label curves
grid on;                                                              % turn grid on

sgtitle('KT-D2V9 Finite-Volume DBM Sod Shock Tube with MUSCL-Minmod'); % title for final comparison figure

%% Two Dimensional Field Plot
figure('Name','2D Field View');                                        % create 2D field figure

subplot(2,2,1);                                                        % density field subplot
imagesc(x,y,rho');                                                     % plot density field; transpose because imagesc uses rows as y
axis xy;                                                              % make y increase upward
colorbar;                                                             % show color scale
title('\rho');                                                        % title density field
xlabel('x');                                                          % label x axis
ylabel('y');                                                          % label y axis

subplot(2,2,2);                                                        % pressure field subplot
imagesc(x,y,p');                                                       % plot pressure field
axis xy;                                                              % make y increase upward
colorbar;                                                             % show color scale
title('p');                                                           % title pressure field
xlabel('x');                                                          % label x axis
ylabel('y');                                                          % label y axis

subplot(2,2,3);                                                        % x velocity field subplot
imagesc(x,y,ux');                                                      % plot x velocity field
axis xy;                                                              % make y increase upward
colorbar;                                                             % show color scale
title('u_x');                                                         % title x velocity field
xlabel('x');                                                          % label x axis
ylabel('y');                                                          % label y axis

subplot(2,2,4);                                                        % temperature field subplot
imagesc(x,y,T');                                                       % plot temperature field
axis xy;                                                              % make y increase upward
colorbar;                                                             % show color scale
title('T');                                                           % title temperature field
xlabel('x');                                                          % label x axis
ylabel('y');                                                          % label y axis

sgtitle('2D Planar Sod Shock Tube Fields');                            % title for 2D field view

%% Local Functions
% MATLAB local functions must stay at the bottom of the file.
% Do not put normal script commands after this point.

function feq = build_KT_D2V9_equilibrium(rho,ux,uy,T,b,v1,v2,eta0,cx,cy)
%BUILD_KT_D2V9_EQUILIBRIUM computes the KT-D2V9 equilibrium populations.
% The output feq has the same shape as f: Nx by Ny by Q.
% This is the part that encodes the Kataoka-Tsutahara compressible equilibrium.
% The coefficients are chosen so equilibrium moments recover density, momentum, and energy.

    [Nx,Ny] = size(rho);                              % get number of x and y cells from rho
    Q = length(cx);                                   % get number of populations from velocity list
    feq = zeros(Nx,Ny,Q);                             % create equilibrium storage

    u2 = ux.^2 + uy.^2;                               % squared macroscopic velocity magnitude

    A_rest = ((b-2)/eta0^2).*T;                       % rest-population coefficient; carries internal energy through eta0
    A_card = (-v2^2 + (((b-2)*v2^2/eta0^2)+2).*T + (v2^2/v1^2).*u2) ./ (4*(v1^2-v2^2));  % coefficient for axial slow-speed populations
    A_diag = (-v1^2 + (((b-2)*v1^2/eta0^2)+2).*T + (v1^2/v2^2).*u2) ./ (4*(v2^2-v1^2));  % coefficient for diagonal fast-speed populations

    B_card = (-v2^2 + (b+2).*T + u2) ./ (2*v1^2*(v1^2-v2^2));  % momentum-shifting coefficient for axial populations
    B_diag = (-v1^2 + (b+2).*T + u2) ./ (2*v2^2*(v2^2-v1^2));  % momentum-shifting coefficient for diagonal populations

    D_card = 1/(2*v1^4);                              % quadratic velocity coefficient for axial populations
    D_diag = 1/(2*v2^4);                              % quadratic velocity coefficient for diagonal populations

    for q = 1:Q                                       % loop through each population
        uDotC = ux*cx(q) + uy*cy(q);                  % projection of fluid velocity onto molecular velocity q

        if q == 1                                     % q=1 is the rest population
            feq(:,:,q) = rho .* A_rest;               % rest population has no directional B or D shift
        elseif q >= 2 && q <= 5                       % q=2..5 are axial/cardinal populations
            feq(:,:,q) = rho .* (A_card + B_card.*uDotC + D_card.*uDotC.^2);  % axial KT equilibrium
        else                                          % q=6..9 are diagonal populations
            feq(:,:,q) = rho .* (A_diag + B_diag.*uDotC + D_diag.*uDotC.^2);  % diagonal KT equilibrium
        end                                           % end population type decision
    end                                               % end equilibrium population loop

end                                                   % end equilibrium function

function s = minmod(a,b)
%MINMOD is the slope limiter used by MUSCL reconstruction.
% If both candidate slopes have the same sign, keep the smaller magnitude.
% If they have opposite signs, return zero because the cell may be near a shock or extremum.
% This prevents MUSCL from creating new overshoots and undershoots near discontinuities.

    if a*b > 0                                        % same sign means both sides agree on increasing/decreasing behavior
        s = sign(a)*min(abs(a),abs(b));               % keep smaller slope to stay conservative and stable
    else                                              % opposite signs mean suspicious local behavior
        s = 0.0;                                      % shut off slope and fall back to first-order locally
    end                                               % end minmod decision

end                                                   % end minmod function

function exact = sod_exact_solution(x,t,x0,rhoL,uL,pL,rhoR,uR,pR,gamma)
%SOD_EXACT_SOLUTION computes the exact Euler solution for the Sod shock tube.
% This is used only for validation of the centerline.
% It does not affect the DBM solver.

    if t <= 0                                         % if time is zero or negative
        exact.rho = rhoR*ones(size(x));               % initialize exact density with right state
        exact.u = uR*ones(size(x));                   % initialize exact velocity with right state
        exact.p = pR*ones(size(x));                   % initialize exact pressure with right state
        exact.rho(x < x0) = rhoL;                     % assign left density left of diaphragm
        exact.u(x < x0) = uL;                         % assign left velocity left of diaphragm
        exact.p(x < x0) = pL;                         % assign left pressure left of diaphragm
        exact.T = exact.p ./ exact.rho;               % compute exact temperature
        return;                                       % leave function because initial exact solution is finished
    end                                               % end t=0 case

    aL = sqrt(gamma*pL/rhoL);                         % left sound speed
    aR = sqrt(gamma*pR/rhoR);                         % right sound speed

    pGuess = 0.5*(pL+pR);                             % starting guess for star-region pressure
    pStar = max(pGuess,1e-8);                         % keep pressure guess positive

    for iter = 1:50                                   % Newton iteration for star pressure
        [fL,dfL] = pressure_function(pStar,rhoL,pL,aL,gamma);  % left wave pressure function and derivative
        [fR,dfR] = pressure_function(pStar,rhoR,pR,aR,gamma);  % right wave pressure function and derivative
        pNew = pStar - (fL + fR + uR - uL)/(dfL + dfR);         % Newton update for pStar
        pNew = max(pNew,1e-10);                       % prevent negative pressure during iteration
        if abs(pNew-pStar)/(pNew+pStar) < 1e-10       % check convergence of pressure solve
            pStar = pNew;                             % accept converged pressure
            break;                                    % leave Newton loop
        end                                           % end convergence check
        pStar = pNew;                                 % update pressure guess
    end                                               % end Newton loop

    [fL,~] = pressure_function(pStar,rhoL,pL,aL,gamma); % left function value at converged pressure
    [fR,~] = pressure_function(pStar,rhoR,pR,aR,gamma); % right function value at converged pressure
    uStar = 0.5*(uL + uR + fR - fL);                  % velocity in the star region

    xOverT = (x - x0)/t;                              % self-similar coordinate xi=(x-x0)/t
    rhoExact = zeros(size(x));                        % storage for exact density
    uExact = zeros(size(x));                          % storage for exact velocity
    pExact = zeros(size(x));                          % storage for exact pressure

    for k = 1:length(x)                               % loop over x sample points
        xi = xOverT(k);                               % current self-similar position

        if xi <= uStar                                % point is left of contact discontinuity
            if pStar > pL                             % left wave is a shock if star pressure exceeds left pressure
                SL = uL - aL*sqrt((gamma+1)/(2*gamma)*pStar/pL + (gamma-1)/(2*gamma));  % left shock speed
                if xi <= SL                           % point is still in undisturbed left state
                    rhoExact(k) = rhoL;               % exact left density
                    uExact(k) = uL;                   % exact left velocity
                    pExact(k) = pL;                   % exact left pressure
                else                                  % point is behind left shock in left star state
                    rhoExact(k) = rhoL*((pStar/pL + (gamma-1)/(gamma+1))/((gamma-1)/(gamma+1)*pStar/pL + 1)); % left star density behind shock
                    uExact(k) = uStar;                % star velocity
                    pExact(k) = pStar;                % star pressure
                end                                   % end left shock region decision
            else                                      % left wave is a rarefaction if star pressure is below left pressure
                aStarL = aL*(pStar/pL)^((gamma-1)/(2*gamma));  % left star sound speed
                SHL = uL - aL;                        % head speed of left rarefaction
                STL = uStar - aStarL;                 % tail speed of left rarefaction
                if xi <= SHL                          % point is in undisturbed left state
                    rhoExact(k) = rhoL;               % exact left density
                    uExact(k) = uL;                   % exact left velocity
                    pExact(k) = pL;                   % exact left pressure
                elseif xi > STL                       % point is in left star region
                    rhoExact(k) = rhoL*(pStar/pL)^(1/gamma);  % left star density in rarefaction case
                    uExact(k) = uStar;                % star velocity
                    pExact(k) = pStar;                % star pressure
                else                                  % point is inside the rarefaction fan
                    uFan = 2/(gamma+1)*(aL + 0.5*(gamma-1)*uL + xi); % exact velocity inside fan
                    aFan = 2/(gamma+1)*(aL + 0.5*(gamma-1)*(uL - xi)); % exact sound speed inside fan
                    rhoExact(k) = rhoL*(aFan/aL)^(2/(gamma-1));       % exact fan density
                    uExact(k) = uFan;                 % exact fan velocity
                    pExact(k) = pL*(aFan/aL)^(2*gamma/(gamma-1));     % exact fan pressure
                end                                   % end left rarefaction region decision
            end                                       % end left-wave type decision
        else                                          % point is right of contact discontinuity
            if pStar > pR                             % right wave is a shock if star pressure exceeds right pressure
                SR = uR + aR*sqrt((gamma+1)/(2*gamma)*pStar/pR + (gamma-1)/(2*gamma));  % right shock speed
                if xi >= SR                           % point is still in undisturbed right state
                    rhoExact(k) = rhoR;               % exact right density
                    uExact(k) = uR;                   % exact right velocity
                    pExact(k) = pR;                   % exact right pressure
                else                                  % point is behind right shock in right star state
                    rhoExact(k) = rhoR*((pStar/pR + (gamma-1)/(gamma+1))/((gamma-1)/(gamma+1)*pStar/pR + 1)); % right star density behind shock
                    uExact(k) = uStar;                % star velocity
                    pExact(k) = pStar;                % star pressure
                end                                   % end right shock region decision
            else                                      % right wave is a rarefaction if star pressure is below right pressure
                aStarR = aR*(pStar/pR)^((gamma-1)/(2*gamma));  % right star sound speed
                SHR = uR + aR;                        % head speed of right rarefaction
                STR = uStar + aStarR;                 % tail speed of right rarefaction
                if xi >= SHR                          % point is in undisturbed right state
                    rhoExact(k) = rhoR;               % exact right density
                    uExact(k) = uR;                   % exact right velocity
                    pExact(k) = pR;                   % exact right pressure
                elseif xi <= STR                      % point is in right star region
                    rhoExact(k) = rhoR*(pStar/pR)^(1/gamma);  % right star density in rarefaction case
                    uExact(k) = uStar;                % star velocity
                    pExact(k) = pStar;                % star pressure
                else                                  % point is inside right rarefaction fan
                    uFan = 2/(gamma+1)*(-aR + 0.5*(gamma-1)*uR + xi); % exact velocity inside fan
                    aFan = 2/(gamma+1)*(aR - 0.5*(gamma-1)*(uR - xi)); % exact sound speed inside fan
                    rhoExact(k) = rhoR*(aFan/aR)^(2/(gamma-1));       % exact fan density
                    uExact(k) = uFan;                 % exact fan velocity
                    pExact(k) = pR*(aFan/aR)^(2*gamma/(gamma-1));     % exact fan pressure
                end                                   % end right rarefaction region decision
            end                                       % end right-wave type decision
        end                                           % end contact-side decision
    end                                               % end exact-solution sample loop

    exact.rho = rhoExact;                             % store exact density in output structure
    exact.u = uExact;                                 % store exact velocity in output structure
    exact.p = pExact;                                 % store exact pressure in output structure
    exact.T = pExact ./ rhoExact;                     % store exact temperature in output structure

end                                                   % end exact Sod solution function

function [f,df] = pressure_function(p,rhoK,pK,aK,gamma)
%PRESSURE_FUNCTION returns the shock/rarefaction function used by the exact Sod solver.
% This is a helper for the Newton solve for the star-region pressure.

    if p > pK                                         % shock branch
        AK = 2/((gamma+1)*rhoK);                      % shock coefficient A
        BK = (gamma-1)/(gamma+1)*pK;                  % shock coefficient B
        f = (p-pK)*sqrt(AK/(p+BK));                   % shock pressure function
        df = sqrt(AK/(p+BK))*(1 - 0.5*(p-pK)/(p+BK)); % derivative of shock pressure function
    else                                              % rarefaction branch
        f = 2*aK/(gamma-1)*((p/pK)^((gamma-1)/(2*gamma)) - 1);       % rarefaction pressure function
        df = (1/(rhoK*aK))*(p/pK)^(-(gamma+1)/(2*gamma));            % derivative of rarefaction pressure function
    end                                               % end shock/rarefaction branch

end                                                   % end pressure function
