function verify_kt_d1v5_theory
%% Numerical checks of the analytical results used in the report %%
% 1. KT-D1V5 equilibrium moments (constrained and unconstrained).
% 2. Chapman-Enskog O(tau) momentum and energy fluxes, checked against the
%    closed-form expressions derived in the report by differentiating the
%    actual equilibrium numerically.
% 3. Eigenvalues of the collision Jacobian d(feq)/df.
% 4. Exact Sod constants (star state, wave speeds, positions at t=0.15).
% 5. Equilibrium populations and temperature sensitivities in the Sod states.
% Output is written to ../data/theory_checks.txt for the report.

    here = fileparts(mfilename('fullpath'));
    dataFolder = fullfile(fileparts(here),'data');
    if ~isfolder(dataFolder); mkdir(dataFolder); end
    fid = fopen(fullfile(dataFolder,'theory_checks.txt'),'w');
    out = @(varargin) fprintf_both(fid,varargin{:});

    gamma = 1.4; b = 2/(gamma-1); v1 = 1; v2 = 3; eta0 = 1.75;
    c = [-v2,-v1,0,v1,v2]; h2 = [0,0,eta0^2,0,0];
    a0 = -v1^2*v2^2; a1 = v1^2+v2^2+(b-1)*v1^2*v2^2/eta0^2; a2 = v1^2+v2^2;

    %% 1. Moments
    rng(1);
    worst = zeros(1,5);
    for trial = 1:2000
        rho = 0.1+2*rand; u = -0.95+1.9*rand; T = 0.3+1.5*rand;
        f = feq_kt(rho,u,T,b,v1,v2,eta0,c);
        m = [sum(f)-rho, f*c'-rho*u, f*(c.^2+h2)'-rho*(b*T+u^2), ...
             f*(c.^2)'-rho*(T+u^2), f*(c.^3)'-rho*u*((b+2)*T+u^2)];
        worst = max(worst,abs(m));
    end
    out('Max moment residuals over 2000 random states [M0 M1 E M2 M3]: %s\n',mat2str(worst,3));
    rho = 0.7; u = 0.4; T = 1.1; f = feq_kt(rho,u,T,b,v1,v2,eta0,c);
    out('M4 check: sum(f c^4) = %.12f, closed form = %.12f\n',f*(c.^4)', ...
        rho*(a0+a1*T+a2*u^2));
    out('alpha0 = %.6f, alpha1 = %.6f, alpha2 = %.6f\n',a0,a1,a2);

    %% 2. Chapman-Enskog fluxes
    worstPi = 0; worstQ = 0; worstMass = 0;
    for trial = 1:500
        st = [0.2+rand, -0.9+1.8*rand, 0.1+rand];            % rho, u, p
        g = randn(1,3);                                        % rho_x, u_x, p_x
        [r,uu,pp] = deal(st(1),st(2),st(3)); TT = pp/r;
        dt0 = [-(uu*g(1)+r*g(2)), -(uu*g(2)+g(3)/r), -(uu*g(3)+gamma*pp*g(2))];
        J = zeros(5,3); h = 1e-6;                              % Jacobian of M0..M4 wrt (rho,u,p)
        for k = 1:3
            sp = st; sm = st; sp(k) = sp(k)+h; sm(k) = sm(k)-h;
            J(:,k) = (moments(sp)-moments(sm))/(2*h);
        end
        Mt = J*dt0'; Mx = J*g';
        mass1 = -(Mt(2)+Mx(3));                                 % sum c f1 / tau
        Pi1 = -(Mt(3)+Mx(4));                                   % sum c^2 f1 / tau
        q1 = -0.5*(Mt(4)+Mx(5));                                % (1/2) sum c(c^2+eta^2) f1 / tau
        PiForm = -(b-1)*(uu*g(3)+gamma*pp*g(2));
        P = (v1^2-uu^2)*(v2^2-uu^2);
        qForm = 0.5*(P*g(1) - (a1-(b+2)*TT-(b+5)*uu^2)*g(3) ...
            + uu*(4*r*uu^2+(b+2)*(gamma+1)*pp-2*a2*r)*g(2));
        worstMass = max(worstMass,abs(mass1));
        worstPi = max(worstPi,abs(Pi1-PiForm)/max(1,abs(PiForm)));
        worstQ = max(worstQ,abs(q1-qForm)/max(1,abs(qForm)));
    end
    out('CE check (500 random states): max|mass flux| = %.2e, rel err Pi = %.2e, rel err q = %.2e\n', ...
        worstMass,worstPi,worstQ);

    %% 3. Collision Jacobian
    st = [0.4,0.9,0.3]; f0 = feq_kt(st(1),st(2),st(3)/st(1),b,v1,v2,eta0,c);
    Jc = zeros(5); h = 1e-7;
    for k = 1:5
        fp = f0; fm = f0; fp(k) = fp(k)+h; fm(k) = fm(k)-h;
        Jc(:,k) = (feq_of_f(fp)-feq_of_f(fm))'/(2*h);
    end
    out('Eigenvalues of d(feq)/df at (rho,u,p)=(0.4,0.9,0.3): %s\n',mat2str(sort(real(eig(Jc)))',4));

    %% 4. Exact Sod constants
    [ps,us,rsL,rsR,S,head,tail] = sod_constants(gamma);
    t = 0.15;
    out('Sod: p* = %.6f, u* = %.6f, rho*L = %.6f, rho*R = %.6f\n',ps,us,rsL,rsR);
    out('Sod: T*L = %.6f, T*R = %.6f, TL = 1, TR = 0.8\n',ps/rsL,ps/rsR);
    out('Sod speeds: head = %.6f, tail = %.6f, contact = %.6f, shock = %.6f\n',head,tail,us,S);
    out('Sod positions at t=0.15 relative to x0: %.6f, %.6f, %.6f, %.6f\n',head*t,tail*t,us*t,S*t);
    aR = sqrt(gamma*0.1/0.125); out('Shock Mach number Ms = %.6f\n',S/aR);
    out('Contact P(u*) = (v1^2-u*^2)(v2^2-u*^2) = %.6f\n',(v1^2-us^2)*(v2^2-us^2));

    %% 5. Populations and temperature sensitivity in the Sod states
    states = {'Left (1,0,1)',1,0,1; 'Right (0.125,0,0.1)',0.125,0,0.1; ...
        'Star-left',rsL,us,ps; 'Star-right',rsR,us,ps};
    for k = 1:size(states,1)
        r = states{k,2}; uu = states{k,3}; TT = states{k,4}/r;
        fk = feq_kt(r,uu,TT,b,v1,v2,eta0,c);
        dTdf = (c.^2+h2-2*uu*c+uu^2-b*TT)/(b*r);
        out('%-22s feq = %s ; dT/df_i = %s\n',states{k,1},mat2str(fk,4),mat2str(dTdf,4));
    end
    out('Slow populations negative when T > (v2^2-u^2)/((b-1)v2^2/eta0^2+1); at u=0: T > %.4f\n', ...
        v2^2/((b-1)*v2^2/eta0^2+1));
    out('Fast populations negative when T < (v1^2-u^2)/((b-1)v1^2/eta0^2+1); at u=0: T < %.4f\n', ...
        v1^2/((b-1)*v1^2/eta0^2+1));

    fclose(fid);

    function M = moments(s)
        ff = feq_kt(s(1),s(2),s(3)/s(1),b,v1,v2,eta0,c);
        M = [sum(ff); ff*c'; ff*(c.^2)'; ff*(c.^3)'; ff*(c.^4)'];
    end
    function fe = feq_of_f(ff)
        r = sum(ff); uu = (ff*c')/r; TT = ((ff*(c.^2+h2)')/r-uu^2)/b;
        fe = feq_kt(r,uu,TT,b,v1,v2,eta0,c);
    end
end

function f = feq_kt(rho,u,T,b,v1,v2,eta0,c)
    Ar = (b-1)/eta0^2*T;
    As = (-v2^2+((b-1)*v2^2/eta0^2+1)*T+u^2)/(2*(v1^2-v2^2));
    Af = (-v1^2+((b-1)*v1^2/eta0^2+1)*T+u^2)/(2*(v2^2-v1^2));
    Bs = (-v2^2+(b+2)*T+u^2)/(2*v1^2*(v1^2-v2^2));
    Bf = (-v1^2+(b+2)*T+u^2)/(2*v2^2*(v2^2-v1^2));
    f = rho*[Af+Bf*u*c(1), As+Bs*u*c(2), Ar, As+Bs*u*c(4), Af+Bf*u*c(5)];
end

function [ps,us,rsL,rsR,S,head,tail] = sod_constants(g)
    rL=1; pL=1; rR=0.125; pR=0.1; aL=sqrt(g*pL/rL); aR=sqrt(g*pR/rR);
    fK = @(p,r,pk,a) (p>pk).*((p-pk).*sqrt((2/((g+1)*r))./(p+(g-1)/(g+1)*pk))) + ...
        (p<=pk).*(2*a/(g-1)*((p/pk).^((g-1)/(2*g))-1));
    ps = fzero(@(p) fK(p,rL,pL,aL)+fK(p,rR,pR,aR),[0.1 1]);
    us = 0.5*(fK(ps,rR,pR,aR)-fK(ps,rL,pL,aL));
    rsL = rL*(ps/pL)^(1/g);
    rsR = rR*((ps/pR+(g-1)/(g+1))/((g-1)/(g+1)*ps/pR+1));
    S = aR*sqrt((g+1)/(2*g)*ps/pR+(g-1)/(2*g));
    head = -aL; tail = us-aL*(ps/pL)^((g-1)/(2*g));
end

function fprintf_both(fid,varargin)
    fprintf(varargin{:}); fprintf(fid,varargin{:});
end
