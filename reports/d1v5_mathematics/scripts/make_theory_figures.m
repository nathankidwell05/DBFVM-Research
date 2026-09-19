function make_theory_figures
%% Figures that depend only on the mathematics (no solver output) %%
%   fig_equilibrium_populations.pdf  KT-D1V5 equilibrium versus T and u
%   fig_sweby_diagram.pdf            limiter functions in the Sweby TVD region
%   fig_rk3_stability.pdf            SSP-RK3 stability region and the solver's eigenvalues
%   fig_exact_sod.pdf                exact Sod profiles and the x-t wave diagram
%   fig_contact_diffusivity.pdf      model contact diffusivity D(U) derived in the report

    here = fileparts(mfilename('fullpath'));
    figFolder = fullfile(fileparts(here),'figures');
    if ~isfolder(figFolder); mkdir(figFolder); end
    S = report_style();

    gamma = 1.4; b = 2/(gamma-1); v1 = 1; v2 = 3; eta0 = 1.75;
    c = [-v2,-v1,0,v1,v2];

    %% Equilibrium populations
    fig = new_figure(S,6.5,2.7);
    t = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    labels = {'$f_1^{eq}$ ($c=-3$)','$f_2^{eq}$ ($c=-1$)','$f_3^{eq}$ ($c=0$)', ...
              '$f_4^{eq}$ ($c=+1$)','$f_5^{eq}$ ($c=+3$)'};
    styles = {'-','--','-','--','-'};
    order = [1 2 3 4 5];
    ax = nexttile(t); hold(ax,'on');
    T = linspace(0.3,1.5,300)';
    F = cell2mat(arrayfun(@(TT) feq_kt(1,0,TT,b,v1,v2,eta0,c),T,'UniformOutput',false));
    patch(ax,[0.711 1.141 1.141 0.711],[-0.6 -0.6 2 2],S.band,'EdgeColor','none','HandleVisibility','off');
    for q = order
        plot(ax,T,F(:,q),styles{q},'Color',S.cat(q,:),'LineWidth',S.lw);
    end
    yline(ax,0,'-','Color',S.axis,'HandleVisibility','off');
    xlabel(ax,'$T$  (with $u=0$)'); ylabel(ax,'$f_i^{eq}/\rho$');
    text(ax,0.93,1.75,'Sod $T$ range','HorizontalAlignment','center','FontSize',S.fs-1,'Color',S.ink2);
    ylim(ax,[-0.6 2]); xlim(ax,[0.3 1.5]); style_axes(ax,S);
    title(ax,'(a) Rest frame, varying temperature','FontWeight','normal');

    ax = nexttile(t); hold(ax,'on');
    U = linspace(-1.2,1.2,300)';
    F = cell2mat(arrayfun(@(uu) feq_kt(1,uu,1,b,v1,v2,eta0,c),U,'UniformOutput',false));
    for q = order
        plot(ax,U,F(:,q),styles{q},'Color',S.cat(q,:),'LineWidth',S.lw);
    end
    yline(ax,0,'-','Color',S.axis,'HandleVisibility','off');
    xline(ax,0.9275,':','Color',S.ink2,'LineWidth',1,'HandleVisibility','off');
    text(ax,0.9,1.75,'$u^*$','HorizontalAlignment','right','FontSize',S.fs-1,'Color',S.ink2);
    xlabel(ax,'$u$  (with $T=1$)'); ylabel(ax,'$f_i^{eq}/\rho$');
    ylim(ax,[-0.6 2]); style_axes(ax,S);
    title(ax,'(b) $T=1$, varying flow velocity','FontWeight','normal');
    lg = legend(ax,labels,'Location','northoutside','NumColumns',5,'Orientation','horizontal');
    lg.Layout.Tile = 'north'; lg.Box = 'off';
    save_figure(fig,fullfile(figFolder,'fig_equilibrium_populations.pdf'));

    %% Sweby diagram
    fig = new_figure(S,6.5,2.9);
    t = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    ax = nexttile(t); hold(ax,'on');
    r = linspace(0,3.5,700);
    tvdUpper = min(2*r,2);
    patch(ax,[r fliplr(r)],[zeros(size(r)) fliplr(tvdUpper)],S.band,'EdgeColor','none');
    phiGM = @(r,th) max(0,min(min(th*r,(1+r)/2),th));
    thetas = [1 1.2 1.5 2];
    names = {'minmod ($\theta=1$)','$\theta=1.2$ (solver)','$\theta=1.5$','MC ($\theta=2$)'};
    cols = S.seq([1 3 4 5],:);
    lws = [S.lw S.lw+0.9 S.lw S.lw];
    for k = 1:numel(thetas)
        plot(ax,r,phiGM(r,thetas(k)),'Color',cols(k,:),'LineWidth',lws(k));
    end
    plot(ax,1,1,'o','MarkerSize',6,'MarkerFaceColor',S.ink,'MarkerEdgeColor','w','HandleVisibility','off');
    text(ax,1.08,0.88,'$\phi(1)=1$','FontSize',S.fs-1,'Color',S.ink);
    xlabel(ax,'$r = \Delta_{j-1/2}/\Delta_{j+1/2}$'); ylabel(ax,'$\phi(r)$');
    xlim(ax,[0 3.5]); ylim(ax,[0 2.2]); style_axes(ax,S);
    legend(ax,[{'Sweby TVD region'},names],'Location','southeast','Box','off','FontSize',S.fs-1);
    title(ax,'(a) Generalized-minmod family','FontWeight','normal');

    ax = nexttile(t); hold(ax,'on');
    % Which candidate is active for theta = 1.2: left (theta r), centre ((1+r)/2), right (theta)
    th = 1.2; rr = linspace(0.001,3.5,2000);
    cand = [th*rr; (1+rr)/2; th*ones(size(rr))];
    [~,active] = min(cand,[],1);
    names2 = {'$\theta\,\Delta_{j-1/2}$ (left) active','$\frac12(\Delta_{j-1/2}+\Delta_{j+1/2})$ (centred) active','$\theta\,\Delta_{j+1/2}$ (right) active'};
    for k = 1:3
        m = active==k;
        plot(ax,rr(m),phiGM(rr(m),th),'.','Color',S.cat(k,:),'MarkerSize',7);
    end
    rA = 1/(2*th-1); rB = 2*th-1;
    xline(ax,rA,':','Color',S.ink2,'HandleVisibility','off'); xline(ax,rB,':','Color',S.ink2,'HandleVisibility','off');
    text(ax,rA,2.05,sprintf('$r=1/(2\\theta-1)=%.3f$',rA),'HorizontalAlignment','center','FontSize',S.fs-2,'Color',S.ink2);
    text(ax,rB,1.85,sprintf('$r=2\\theta-1=%.1f$',rB),'HorizontalAlignment','center','FontSize',S.fs-2,'Color',S.ink2);
    xlabel(ax,'$r$'); ylabel(ax,'$\phi(r)$ for $\theta=1.2$');
    xlim(ax,[0 3.5]); ylim(ax,[0 2.2]); style_axes(ax,S);
    legend(ax,names2,'Location','southeast','Box','off','FontSize',S.fs-2);
    title(ax,'(b) Which slope candidate is selected','FontWeight','normal');
    save_figure(fig,fullfile(figFolder,'fig_sweby_diagram.pdf'));

    %% SSP-RK3 stability region
    fig = new_figure(S,6.5,3.0);
    t = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    ax = nexttile(t); hold(ax,'on');
    [X,Y] = meshgrid(linspace(-3.2,1,600),linspace(-3,3,600));
    Z = X+1i*Y; R = abs(1+Z+Z.^2/2+Z.^3/6);
    contourf(ax,X,Y,double(R<=1),[0.5 0.5],'FaceColor',S.band,'EdgeColor',S.seq(4,:),'LineWidth',S.lw);
    Ref = abs(1+Z); contour(ax,X,Y,Ref,[1 1],'--','Color',S.ink2,'LineWidth',1);
    plot(ax,-0.25,0,'o','MarkerSize',7,'MarkerFaceColor',S.cat(2,:),'MarkerEdgeColor','w');
    plot(ax,-2.5127,0,'s','MarkerSize',6,'MarkerFaceColor',S.ink,'MarkerEdgeColor','w');
    text(ax,0.12,0.45,'$-\Delta t/\tau=-0.25$','HorizontalAlignment','left','FontSize',S.fs-2,'Color',S.ink);
    text(ax,-2.5127,-0.4,'$-2.513$','HorizontalAlignment','center','FontSize',S.fs-2,'Color',S.ink);
    xline(ax,0,'-','Color',S.axis); yline(ax,0,'-','Color',S.axis);
    axis(ax,'equal'); xlim(ax,[-3.2 1.6]); ylim(ax,[-2.6 2.6]); style_axes(ax,S);
    xlabel(ax,'$\mathrm{Re}(z)$'); ylabel(ax,'$\mathrm{Im}(z)$');
    legend(ax,{'SSP-RK3 $|R(z)|\le 1$','Forward Euler','BGK mode (solver)','Real-axis limit'}, ...
        'Location','southoutside','NumColumns',2,'Box','off','FontSize',S.fs-2);
    title(ax,'(a) Stability region, $z=\lambda\Delta t$','FontWeight','normal');

    ax = nexttile(t); hold(ax,'on');
    % Fourier symbol of first-order upwind (circle) and MUSCL (unlimited Fromm-type)
    k = linspace(0,2*pi,400);
    nu = 3*1.25e-6/(20/50000);                     % effective CFL on the finest grid
    zUp = -nu*(1-exp(-1i*k));
    zFromm = -nu*(1-exp(-1i*k)).*(1+0.25*(exp(1i*k)-exp(-1i*k)));
    plot(ax,real(zUp),imag(zUp),'-','Color',S.seq(3,:),'LineWidth',S.lw);
    plot(ax,real(zFromm),imag(zFromm),'-','Color',S.cat(2,:),'LineWidth',S.lw);
    [X,Y] = meshgrid(linspace(-0.05,0.01,400),linspace(-0.03,0.03,400));
    Z = X+1i*Y; R = abs(1+Z+Z.^2/2+Z.^3/6);
    contour(ax,X,Y,R,[1 1],'-','Color',S.ink2,'LineWidth',1);
    xline(ax,0,'-','Color',S.axis); yline(ax,0,'-','Color',S.axis);
    axis(ax,'equal'); style_axes(ax,S);
    xlabel(ax,'$\mathrm{Re}(z)$'); ylabel(ax,'$\mathrm{Im}(z)$');
    legend(ax,{'Upwind, $\nu=0.0094$','MUSCL (unlimited), $\nu=0.0094$','SSP-RK3 boundary'}, ...
        'Location','southoutside','NumColumns',2,'Box','off','FontSize',S.fs-2);
    title(ax,'(b) Transport spectrum, $N_x=50{,}000$','FontWeight','normal');
    save_figure(fig,fullfile(figFolder,'fig_rk3_stability.pdf'));

    %% Exact Sod solution and wave diagram
    [ps,us,rsL,rsR,Ssh,head,tail] = sod_constants(gamma);
    x = linspace(9.6,10.5,4000)'; x0 = 10; tEnd = 0.15;
    ex = sod_exact(x,tEnd,x0,gamma,ps,us,rsL,rsR,Ssh,head,tail);
    fig = new_figure(S,6.5,4.6);
    t = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
    fields = {'rho','u','p','T'}; ylab = {'$\rho$','$u$','$p$','$T$'};
    tiles = [1 2 4 5];
    for kk = 1:4
        ax = nexttile(t,tiles(kk)); hold(ax,'on');
        xl = x0+[head tail us Ssh]*tEnd;
        for m = 1:4; xline(ax,xl(m),':','Color',S.axis); end
        plot(ax,x,ex.(fields{kk}),'-','Color',S.ink,'LineWidth',S.lw);
        xlim(ax,[9.6 10.5]); style_axes(ax,S); ylabel(ax,ylab{kk});
        if kk>2; xlabel(ax,'$x$'); end
    end
    ax = nexttile(t,3,[2 1]); hold(ax,'on');
    tt = linspace(0,tEnd,2);
    for s = linspace(head,tail,7)
        plot(ax,x0+s*tt,tt,'-','Color',S.seq(2,:),'LineWidth',0.9);
    end
    plot(ax,x0+us*tt,tt,'--','Color',S.cat(2,:),'LineWidth',S.lw);
    plot(ax,x0+Ssh*tt,tt,'-','Color',S.cat(8,:),'LineWidth',S.lw+0.4);
    patch(ax,[9.7 10.4 10.4 9.7],[tEnd tEnd tEnd*1.04 tEnd*1.04],S.band,'EdgeColor','none');
    text(ax,10.05,tEnd*1.02,'error window','HorizontalAlignment','center','FontSize',S.fs-2,'Color',S.ink2);
    text(ax,9.62,0.125,'rarefaction','HorizontalAlignment','left','FontSize',S.fs-2,'Color',S.seq(4,:));
    text(ax,x0+us*0.13-0.012,0.13,'contact','HorizontalAlignment','right','FontSize',S.fs-2,'Color',S.cat(2,:));
    text(ax,x0+Ssh*0.06+0.03,0.06,'shock','HorizontalAlignment','left','FontSize',S.fs-2,'Color',S.cat(8,:));
    xlim(ax,[9.6 10.5]); ylim(ax,[0 tEnd*1.04]); style_axes(ax,S);
    xlabel(ax,'$x$'); ylabel(ax,'$t$');
    title(ax,'Wave diagram','FontWeight','normal');
    save_figure(fig,fullfile(figFolder,'fig_exact_sod.pdf'));

    %% Contact diffusivity D(U)
    fig = new_figure(S,6.5,2.5);
    t = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    ax = nexttile(t); hold(ax,'on');
    U = linspace(0,1.3,400);
    P = (v1^2-U.^2).*(v2^2-U.^2);
    plot(ax,U,P/(v1^2*v2^2),'-','Color',S.seq(4,:),'LineWidth',S.lw);
    yline(ax,0,'-','Color',S.axis); xline(ax,us,':','Color',S.ink2);
    plot(ax,us,(v1^2-us^2)*(v2^2-us^2)/9,'o','MarkerFaceColor',S.cat(2,:),'MarkerEdgeColor','w','MarkerSize',7);
    text(ax,us+0.03,0.3,sprintf('$u^*$: %.3f',(v1^2-us^2)*(v2^2-us^2)/9),'FontSize',S.fs-1,'Color',S.ink);
    xlabel(ax,'contact velocity $U$'); ylabel(ax,'$P(U)/P(0)$');
    xlim(ax,[0 1.3]); ylim(ax,[-0.4 1.05]); style_axes(ax,S);
    title(ax,'(a) Galilean dependence of contact diffusion','FontWeight','normal');
    ax = nexttile(t); hold(ax,'on');
    tauList = [5e-6 5e-5]; dxList = 20./[3125 6250 12500 25000 50000];
    Pst = (v1^2-us^2)*(v2^2-us^2);
    for k = 1:2
        D = tauList(k)*Pst./((b+2)*[ps/rsR ps/rsL]);
        w = sqrt(4*pi*D*tEnd);
        patch(ax,[1e-4 1e-2 1e-2 1e-4],[w(1) w(1) w(2) w(2)],S.cat(k,:),'FaceAlpha',0.25,'EdgeColor','none');
    end
    loglog(ax,dxList,dxList,'--','Color',S.ink2);
    set(ax,'XScale','log','YScale','log','XDir','reverse');
    text(ax,4e-3,1.6e-3,'$\tau=5\times10^{-6}$','FontSize',S.fs-1,'Color',S.ink);
    text(ax,4e-3,5.4e-3,'$\tau=5\times10^{-5}$','FontSize',S.fs-1,'Color',S.ink);
    text(ax,1.1e-3,0.8e-3,'$\delta=\Delta x$','FontSize',S.fs-1,'Color',S.ink2);
    xlabel(ax,'$\Delta x$'); ylabel(ax,'predicted contact thickness $\delta_c$');
    xlim(ax,[3e-4 8e-3]); ylim(ax,[3e-4 1e-2]); style_axes(ax,S);
    set(ax,'XTick',fliplr(dxList),'XTickLabel',{'$4\times10^{-4}$','$8\times10^{-4}$','$1.6\times10^{-3}$','$3.2\times10^{-3}$','$6.4\times10^{-3}$'},'XMinorGrid','off');
    title(ax,'(b) Kinetic contact thickness at $t=0.15$','FontWeight','normal');
    save_figure(fig,fullfile(figFolder,'fig_contact_diffusivity.pdf'));
end

%% ---------------------------------------------------------------------
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

function ex = sod_exact(x,t,x0,g,ps,us,rsL,rsR,S,head,tail)
    xi = (x-x0)/t; aL = sqrt(g);
    ex.rho = 0.125*ones(size(x)); ex.u = zeros(size(x)); ex.p = 0.1*ones(size(x));
    m = xi<=head; ex.rho(m)=1; ex.p(m)=1;
    m = xi>head & xi<=tail;
    uF = 2/(g+1)*(aL+xi(m)); aF = 2/(g+1)*(aL-0.5*(g-1)*xi(m));
    ex.rho(m) = (aF/aL).^(2/(g-1)); ex.u(m) = uF; ex.p(m) = (aF/aL).^(2*g/(g-1));
    m = xi>tail & xi<=us; ex.rho(m)=rsL; ex.u(m)=us; ex.p(m)=ps;
    m = xi>us & xi<=S;    ex.rho(m)=rsR; ex.u(m)=us; ex.p(m)=ps;
    ex.T = ex.p./ex.rho;
end
