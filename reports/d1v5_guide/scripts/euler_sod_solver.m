function out = euler_sod_solver(Nx,scheme,theta,CFL)
%% Conventional Finite-Volume Euler Solver for the Sod Shock Tube %%
% Reference method for the accuracy-versus-cost comparison with the KT-D1V5
% discrete Boltzmann solver. It solves the Euler equations directly for
% (rho, rho*u, E) with:
%   * the HLLC approximate Riemann solver (Toro, Spruce and Speares 1994),
%   * first-order (scheme = 'firstorder') or MUSCL reconstruction of the
%     primitive variables with the same generalized-minmod limiter as the
%     DBM solver (scheme = 'muscl', slope parameter theta),
%   * SSP-RK3 time integration with dt = CFL*dx/max(|u|+a).
% The domain, initial condition, final time, and error definitions are the
% same as in the DBM study (Lx = 20, x0 = 10, tEnd = 0.15, gamma = 1.4).
%
% Example:
%     out = euler_sod_solver(6250,'muscl',1.2,0.5);

    if nargin < 3; theta = 1.2; end
    if nargin < 4; CFL = 0.5; end

    gamma = 1.4;                                   % ratio of specific heats
    Lx = 20; x0 = 0.5*Lx; tEnd = 0.15;             % same problem as the DBM study
    dx = Lx/Nx;                                    % cell width
    x = ((1:Nx)'-0.5)*dx;                          % cell centres

    rho = 0.125*ones(Nx,1); u = zeros(Nx,1); p = 0.1*ones(Nx,1);
    left = x < x0;                                 % same diaphragm rule as the DBM solver
    rho(left) = 1; p(left) = 1;
    U = prim_to_cons(rho,u,p,gamma);               % conservative variables [rho, rho*u, E]

    time = 0; steps = 0;
    timer = tic;
    while time < tEnd - 1e-14
        [r,v,pp] = cons_to_prim(U,gamma);
        a = sqrt(gamma*pp./r);
        dt = CFL*dx/max(abs(v)+a);                 % transport-limited time step
        dt = min(dt,tEnd-time);                    % land exactly on the final time

        U1 = U + dt*rhs(U,dx,gamma,scheme,theta);                       % SSP-RK3 stage 1
        U2 = 0.75*U + 0.25*(U1 + dt*rhs(U1,dx,gamma,scheme,theta));     % stage 2
        U  = (1/3)*U + (2/3)*(U2 + dt*rhs(U2,dx,gamma,scheme,theta));   % stage 3

        time = time + dt; steps = steps + 1;
    end
    runTime = toc(timer);

    [rho,u,p] = cons_to_prim(U,gamma);
    out = struct('Nx',Nx,'dx',dx,'x',x,'x0',x0,'tEnd',tEnd,'scheme',scheme,'theta',theta, ...
        'CFL',CFL,'steps',steps,'runTime',runTime,'rho',rho,'u',u,'p',p,'T',p./rho);
    out.exact = sod_exact(x,tEnd,x0,gamma);
end

%% ------------------------------------------------------------------------
function L = rhs(U,dx,gamma,scheme,theta)
% finite-volume right-hand side -dF/dx with two transmissive ghost cells per side
    Ug = [U(1,:); U(1,:); U; U(end,:); U(end,:)];
    [r,v,p] = cons_to_prim(Ug,gamma);
    W = [r v p];                                   % primitive variables
    n = size(W,1);
    sigma = zeros(n,3);
    if strcmp(scheme,'muscl')
        dL = W(2:n-1,:) - W(1:n-2,:);
        dR = W(3:n,:) - W(2:n-1,:);
        sigma(2:n-1,:) = gminmod(dL,dR,theta);
    end
    % faces between padded cells k and k+1 for k = 2..n-2 (Nx+1 faces)
    k = (2:n-2)';
    WL = W(k,:) + 0.5*sigma(k,:);
    WR = W(k+1,:) - 0.5*sigma(k+1,:);
    F = hllc(WL,WR,gamma);
    L = -(F(2:end,:) - F(1:end-1,:))/dx;
end

function s = gminmod(dL,dR,theta)
% generalized minmod: minmod(theta*dL, (dL+dR)/2, theta*dR)
    a = theta*dL; b = 0.5*(dL+dR); c = theta*dR;
    s = zeros(size(dL));
    pos = a>0 & b>0 & c>0; neg = a<0 & b<0 & c<0;
    s(pos) = min(min(a(pos),b(pos)),c(pos));
    s(neg) = max(max(a(neg),b(neg)),c(neg));
end

function F = hllc(WL,WR,gamma)
% HLLC flux (Toro 2009, chapter 10) with Davis wave-speed estimates
    rL = WL(:,1); uL = WL(:,2); pL = WL(:,3);
    rR = WR(:,1); uR = WR(:,2); pR = WR(:,3);
    aL = sqrt(gamma*pL./rL); aR = sqrt(gamma*pR./rR);
    EL = pL/(gamma-1) + 0.5*rL.*uL.^2; ER = pR/(gamma-1) + 0.5*rR.*uR.^2;
    SL = min(uL-aL,uR-aR); SR = max(uL+aL,uR+aR);
    Ss = (pR - pL + rL.*uL.*(SL-uL) - rR.*uR.*(SR-uR)) ./ (rL.*(SL-uL) - rR.*(SR-uR));

    UL = [rL, rL.*uL, EL]; UR = [rR, rR.*uR, ER];
    FL = [rL.*uL, rL.*uL.^2+pL, uL.*(EL+pL)];
    FR = [rR.*uR, rR.*uR.^2+pR, uR.*(ER+pR)];
    cL = rL.*(SL-uL)./(SL-Ss); cR = rR.*(SR-uR)./(SR-Ss);
    UsL = [cL, cL.*Ss, cL.*(EL./rL + (Ss-uL).*(Ss + pL./(rL.*(SL-uL))))];
    UsR = [cR, cR.*Ss, cR.*(ER./rR + (Ss-uR).*(Ss + pR./(rR.*(SR-uR))))];

    F = FL;                                        % supersonic to the right (SL >= 0)
    m = SL<0 & Ss>=0;  F(m,:) = FL(m,:) + SL(m).*(UsL(m,:)-UL(m,:));
    m = Ss<0 & SR>0;   F(m,:) = FR(m,:) + SR(m).*(UsR(m,:)-UR(m,:));
    m = SR<=0;         F(m,:) = FR(m,:);
end

function U = prim_to_cons(r,u,p,gamma)
    U = [r, r.*u, p/(gamma-1) + 0.5*r.*u.^2];
end

function [r,u,p] = cons_to_prim(U,gamma)
    r = U(:,1); u = U(:,2)./r; p = (gamma-1)*(U(:,3) - 0.5*r.*u.^2);
end

function ex = sod_exact(x,t,x0,g)
% exact Sod solution (same constants as the DBM study)
    rL=1; pL=1; rR=0.125; pR=0.1; aL=sqrt(g*pL/rL); aR=sqrt(g*pR/rR);
    fK = @(pp,r,pk,a) (pp>pk).*((pp-pk).*sqrt((2/((g+1)*r))./(pp+(g-1)/(g+1)*pk))) + ...
        (pp<=pk).*(2*a/(g-1)*((pp/pk).^((g-1)/(2*g))-1));
    ps = fzero(@(pp) fK(pp,rL,pL,aL)+fK(pp,rR,pR,aR),[0.1 1]);
    us = 0.5*(fK(ps,rR,pR,aR)-fK(ps,rL,pL,aL));
    rsL = rL*(ps/pL)^(1/g);
    rsR = rR*((ps/pR+(g-1)/(g+1))/((g-1)/(g+1)*ps/pR+1));
    S = aR*sqrt((g+1)/(2*g)*ps/pR+(g-1)/(2*g));
    head = -aL; tail = us-aL*(ps/pL)^((g-1)/(2*g));
    xi = (x-x0)/t;
    ex.rho = rR*ones(size(x)); ex.u = zeros(size(x)); ex.p = pR*ones(size(x));
    m = xi<=head; ex.rho(m)=rL; ex.p(m)=pL;
    m = xi>head & xi<=tail;
    uF = 2/(g+1)*(aL+xi(m)); aF = 2/(g+1)*(aL-0.5*(g-1)*xi(m));
    ex.rho(m) = rL*(aF/aL).^(2/(g-1)); ex.u(m) = uF; ex.p(m) = pL*(aF/aL).^(2*g/(g-1));
    m = xi>tail & xi<=us; ex.rho(m)=rsL; ex.u(m)=us; ex.p(m)=ps;
    m = xi>us & xi<=S;    ex.rho(m)=rsR; ex.u(m)=us; ex.p(m)=ps;
    ex.T = ex.p./ex.rho;
end
