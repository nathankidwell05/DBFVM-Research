function make_results_figures
%% Figures and LaTeX tables built from solver output %%
% Reads
%   results/runs/d1v5_gminmod_theta1p20_tau5e-06_CFL0p025_Nx*.mat   (convergence grids)
%   results/d1v5_limiter_history_Nx6250_tau5e-06_CFL0p025_results.mat
%   results/theta_sweep/*.mat
% Writes vector PDFs to ../figures and LaTeX table bodies to ../data.

    here = fileparts(mfilename('fullpath'));
    reportFolder = fileparts(here);
    projectFolder = project_paths();                    % solver folder; results live in projectFolder/results
    figFolder = fullfile(reportFolder,'figures');
    dataFolder = fullfile(reportFolder,'data');
    if ~isfolder(dataFolder); mkdir(dataFolder); end
    S = report_style();
    gamma = 1.4; b = 5; v1 = 1; v2 = 3;

    %% Load the convergence grids
    NxList = [3125 6250 12500 25000 50000];
    runs = {};
    for Nx = NxList
        fn = fullfile(projectFolder,'results','runs', ...
            sprintf('d1v5_gminmod_theta1p20_tau5e-06_CFL0p025_Nx%d.mat',Nx));
        if isfile(fn)
            tmp = load(fn,'runData'); runs{end+1} = tmp.runData; %#ok<AGROW>
        end
    end
    nG = numel(runs);
    fprintf('Loaded %d convergence grids\n',nG);
    Nx = cellfun(@(r) r.Nx,runs); dx = cellfun(@(r) r.dx,runs);
    x0 = runs{1}.x0; tEnd = runs{1}.tEnd; tau = runs{1}.tau;
    gridNames = arrayfun(@(n) sprintf('$N_x=%s$',commas(n)),Nx,'UniformOutput',false);
    gridCols = S.seq(end-nG+1:end,:);

    [ps,us,rsL,rsR,Ssh,head,tail] = sod_constants(gamma);

    %% Metrics
    vars = {'rho','u','p','T'};
    waveOf = @(r) r.x >= x0-0.30 & r.x <= x0+0.40;
    M = struct();
    for k = 1:nG
        r = runs{k}; w = waveOf(r);
        for v = 1:4
            e = r.(vars{v}) - r.exact.(vars{v});
            M.L1(k,v) = mean(abs(e));
            M.L2(k,v) = sqrt(mean(e.^2));
            M.wL1(k,v) = mean(abs(e(w)));
            M.wL2(k,v) = sqrt(mean(e(w).^2));
            M.Linf(k,v) = max(abs(e(w)));
        end
        M.Tpeak(k) = max(r.T(w)) - max(r.exact.T(w));
        M.dt(k) = r.dtBase; M.cfl(k) = r.maxSpeed*r.dtBase/r.dx;
        M.retry(k) = r.retryCountTotal; M.fall(k) = r.fallbackCountTotal;
        M.steps(k) = r.steps; M.time(k) = r.runTime;
        % conservation: mass, momentum, energy in the code normalisation
        E = 0.5*r.rho.*(b*r.T + r.u.^2);
        M.mass(k) = sum(r.rho)*r.dx; M.mom(k) = sum(r.rho.*r.u)*r.dx; M.energy(k) = sum(E)*r.dx;
        left0 = r.x < x0;                                  % discrete initial condition (as in the solver)
        M.mass0(k) = (sum(left0)*1 + sum(~left0)*0.125)*r.dx;
        M.energy0(k) = (sum(left0)*2.5 + sum(~left0)*0.25)*r.dx;
        M.leftCells(k) = sum(left0);
        % feature-zone contributions to the global L2 error (add in quadrature)
        zones = {[x0-0.25, x0+0.07],[x0+0.07, x0+0.20],[x0+0.20, x0+0.40]};
        for v = 1:4
            e = r.(vars{v}) - r.exact.(vars{v});
            for z = 1:3
                m = r.x >= zones{z}(1) & r.x < zones{z}(2);
                M.zone(k,v,z) = sqrt(sum(e(m).^2)/numel(r.x));
            end
        end
        % exact solution averaged over each cell (the consistent finite-volume reference)
        exAvg = cell_averaged_exact(r.x,r.dx,x0,tEnd,gamma,ps,us,rsL,rsR,Ssh,head,tail);
        for v = 1:4
            eA = r.(vars{v}) - exAvg.(vars{v});
            M.L2avg(k,v) = sqrt(mean(eA.^2));
            M.L1avg(k,v) = mean(abs(eA));
        end
        xsExact = x0+Ssh*tEnd;
        M.phase(k) = mod(xsExact/r.dx,1);                  % where the exact shock sits inside its cell
        mShock = abs(r.x-xsExact) < 4*r.dx;
        eU = r.u - r.exact.u; eUa = r.u - exAvg.u;
        M.shockU(k) = sqrt(sum(eU(mShock).^2)/numel(r.x));
        M.shockUavg(k) = sqrt(sum(eUa(mShock).^2)/numel(r.x));
        % contact and shock thickness and location from the density profile
        [M.dc(k),M.xc(k)] = thickness(r.x,r.rho,[x0+0.07 x0+0.20],rsL,rsR);
        [M.ds(k),M.xs(k)] = thickness(r.x,r.rho,[x0+0.20 x0+0.32],rsR,0.125);
        [M.dcT(k),~] = thickness(r.x,r.T,[x0+0.07 x0+0.20],ps/rsL,ps/rsR);
    end
    massExact = 10*1 + 10*0.125; momExact = (1-0.1)*tEnd; energyExact = 10*2.5 + 10*0.25;

    orders = @(E) [NaN, log(E(1:end-1)./E(2:end))./log(dx(1:end-1)./dx(2:end))];
    fitSlope = @(E) polyfit(log(dx),log(E),1);

    %% Tables: convergence
    fid = fopen(fullfile(dataFolder,'tab_conv_setup.tex'),'w');
    for k = 1:nG
        fprintf(fid,'%s & %s & %s & %s & %s & %d & %d \\\\\n',commas(Nx(k)),sci(dx(k),2),sci(M.dt(k),3), ...
            sci(M.cfl(k),3),commas(M.steps(k)),M.retry(k),M.fall(k));
    end
    fclose(fid);
    write_err_table(fullfile(dataFolder,'tab_conv_L1.tex'),Nx,M.L1,[]);
    write_err_table(fullfile(dataFolder,'tab_conv_L2.tex'),Nx,M.L2,M.Tpeak);
    write_err_table(fullfile(dataFolder,'tab_conv_waveL2.tex'),Nx,M.wL2,[]);
    fid = fopen(fullfile(dataFolder,'tab_conv_orders.tex'),'w');
    names = {'$\rho$','$u$','$p$','$T$'};
    for v = 1:4
        oL2 = orders(M.L2(:,v)'); oL1 = orders(M.L1(:,v)'); ow = orders(M.wL2(:,v)');
        pL2 = fitSlope(M.L2(:,v)'); pL1 = fitSlope(M.L1(:,v)'); pw = fitSlope(M.wL2(:,v)');
        fprintf(fid,'%s',names{v});
        for k = 2:nG; fprintf(fid,' & %.3f',oL2(k)); end
        fprintf(fid,' & \\textbf{%.3f} & %.3f & %.3f \\\\\n',pL2(1),pL1(1),pw(1));
    end
    fclose(fid);
    write_err_table(fullfile(dataFolder,'tab_conv_L2avg.tex'),Nx,M.L2avg,[]);
    fid = fopen(fullfile(dataFolder,'tab_shock_phase.tex'),'w');
    for k = 1:nG
        fprintf(fid,'%s & %.2f & %s & %.4f & %s & %.4f \\\\\n',commas(Nx(k)),M.phase(k), ...
            sci(M.shockU(k),3),M.shockU(k)/sqrt(dx(k)),sci(M.shockUavg(k),3),M.shockUavg(k)/sqrt(dx(k)));
    end
    fclose(fid);
    fid = fopen(fullfile(dataFolder,'tab_conv_orders_avg.tex'),'w');
    for v = 1:4
        oA = orders(M.L2avg(:,v)'); pA = fitSlope(M.L2avg(:,v)'); pA1 = fitSlope(M.L1avg(:,v)');
        fprintf(fid,'%s',names{v});
        for k = 2:nG; fprintf(fid,' & %.3f',oA(k)); end
        fprintf(fid,' & \\textbf{%.3f} & %.3f \\\\\n',pA(1),pA1(1));
    end
    fclose(fid);
    fid = fopen(fullfile(dataFolder,'tab_conservation.tex'),'w');
    for k = 1:nG
        fprintf(fid,'%s & %s & %s & %s \\\\\n',commas(Nx(k)), sci(M.mass(k)-M.mass0(k),2), ...
            sci(M.mom(k)-momExact,2), sci(M.energy(k)-M.energy0(k),2));
    end
    fclose(fid);
    fid = fopen(fullfile(dataFolder,'tab_features.tex'),'w');
    for k = 1:nG
        fprintf(fid,'%s & %s & %.2f & %s & %.2f & %s & %s \\\\\n',commas(Nx(k)), ...
            sci(M.dc(k),2),M.dc(k)/dx(k),sci(M.ds(k),2),M.ds(k)/dx(k), ...
            sci(M.xc(k)-(x0+us*tEnd),2),sci(M.xs(k)-(x0+Ssh*tEnd),2));
    end
    fclose(fid);
    pc = polyfit(log(dx),log(M.dc),1); pshk = polyfit(log(dx),log(M.ds),1);
    Pst = (v1^2-us^2)*(v2^2-us^2);
    Dc = tau*Pst./((b+2)*[ps/rsR ps/rsL]);
    dKin = sqrt(4*pi*Dc*tEnd);
    % contact-only L2 floor from an erf profile of the kinetic width (density)
    wErf = sqrt(4*Dc*tEnd); Jrho = rsL-rsR;
    floorRho = sqrt(Jrho^2*wErf*(2-sqrt(2))/(2*sqrt(pi))/20);
    fid = fopen(fullfile(dataFolder,'numbers.tex'),'w');
    fprintf(fid,'%% auto-generated by make_results_figures.m\n');
    fprintf(fid,'\\newcommand{\\contactWidthSlope}{%.2f}\n',pc(1));
    fprintf(fid,'\\newcommand{\\shockWidthSlope}{%.2f}\n',pshk(1));
    fprintf(fid,'\\newcommand{\\kinContactMin}{%s}\n',strrep(sci(min(dKin),2),'$',''));
    fprintf(fid,'\\newcommand{\\kinContactMax}{%s}\n',strrep(sci(max(dKin),2),'$',''));
    fprintf(fid,'\\newcommand{\\rhoFloorMin}{%s}\n',strrep(sci(min(floorRho),2),'$',''));
    fprintf(fid,'\\newcommand{\\rhoFloorMax}{%s}\n',strrep(sci(max(floorRho),2),'$',''));
    for v = 1:4
        p2 = fitSlope(M.L2(:,v)');
        fprintf(fid,'\\newcommand{\\fitLtwo%s}{%.3f}\n',char('A'+v-1),p2(1));
    end
    for v = 1:4
        pA = fitSlope(M.L2avg(:,v)');
        fprintf(fid,'\\newcommand{\\fitLtwoAvg%s}{%.3f}\n',char('A'+v-1),pA(1));
    end
    pMatchA = log(M.shockU(1)/M.shockU(4))/log(dx(1)/dx(4));
    pMatchB = log(M.shockU(2)/M.shockU(5))/log(dx(2)/dx(5));
    fprintf(fid,'\\newcommand{\\shockPhaseSlopeA}{%.3f}\n',pMatchA);
    fprintf(fid,'\\newcommand{\\shockPhaseSlopeB}{%.3f}\n',pMatchB);
    for v = [1 4]
        for z = 1:3
            pz = fitSlope(squeeze(M.zone(:,v,z))');
            fprintf(fid,'\\newcommand{\\zoneSlope%s%s}{%.2f}\n',char('A'+v-1),char('A'+z-1),pz(1));
        end
    end
    fclose(fid);
    save(fullfile(dataFolder,'convergence_metrics.mat'),'M','Nx','dx');
    disp(array2table([Nx' dx' M.L2 M.Tpeak'],'VariableNames',{'Nx','dx','rho','u','p','T','Tpk'}));

    %% Figure: finest-grid solution and pointwise error
    rf = runs{end}; w = rf.x >= x0-0.30 & rf.x <= x0+0.40;
    fig = new_figure(S,6.5,5.2);
    t = tiledlayout(fig,4,2,'TileSpacing','tight','Padding','compact');
    ylab = {'$\rho$','$u$','$p$','$T$'};
    for v = 1:4
        ax = nexttile(t,2*v-1); hold(ax,'on');
        plot(ax,rf.x(w),rf.exact.(vars{v})(w),'-','Color',S.ink,'LineWidth',1.0);
        plot(ax,rf.x(w),rf.(vars{v})(w),'--','Color',S.cat(1,:),'LineWidth',S.lw);
        xlim(ax,[x0-0.3 x0+0.4]); ylabel(ax,ylab{v}); style_axes(ax,S);
        if v<4; ax.XTickLabel = []; else; xlabel(ax,'$x$'); end
        if v==1; legend(ax,{'Exact Euler',sprintf('FVDBM, $N_x=%s$',commas(rf.Nx))},'Location','northeast','Box','off'); title(ax,'Solution','FontWeight','normal'); end
        ax = nexttile(t,2*v); hold(ax,'on');
        plot(ax,rf.x(w),rf.(vars{v})(w)-rf.exact.(vars{v})(w),'-','Color',S.cat(2,:),'LineWidth',1.0);
        yline(ax,0,'-','Color',S.axis);
        xlim(ax,[x0-0.3 x0+0.4]); ylabel(ax,['error in ' ylab{v}]); style_axes(ax,S);
        if v<4; ax.XTickLabel = []; else; xlabel(ax,'$x$'); end
        if v==1; title(ax,'Pointwise error (numerical $-$ exact)','FontWeight','normal'); end
    end
    save_figure(fig,fullfile(figFolder,'fig_solution_finest.pdf'));

    %% Figure: grid-refinement zooms
    fig = new_figure(S,6.5,4.6);
    t = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
    zooms = {'rho',[x0-0.02 x0+0.03],'(a) $\rho$ at rarefaction tail'; ...
             'rho',[x0+0.12 x0+0.16],'(b) $\rho$ at contact'; ...
             'rho',[x0+0.245 x0+0.28],'(c) $\rho$ at shock'; ...
             'T',[x0+0.12 x0+0.28],'(d) $T$ across contact and shock'};
    for zI = 1:4
        ax = nexttile(t); hold(ax,'on');
        fld = zooms{zI,1}; xr = zooms{zI,2};
        m = rf.x>=xr(1) & rf.x<=xr(2);
        plot(ax,rf.x(m),rf.exact.(fld)(m),'-','Color',S.ink,'LineWidth',1.0);
        for k = 1:nG
            r = runs{k}; m = r.x>=xr(1)-0.01 & r.x<=xr(2)+0.01;
            plot(ax,r.x(m),r.(fld)(m),'-','Color',gridCols(k,:),'LineWidth',1.2);
        end
        xlim(ax,xr); style_axes(ax,S); title(ax,zooms{zI,3},'FontWeight','normal'); xlabel(ax,'$x$');
        if zI==1; legend(ax,[{'Exact'},gridNames],'Location','northeast','Box','off','FontSize',S.fs-2); end
    end
    save_figure(fig,fullfile(figFolder,'fig_grid_zoom.pdf'));

    %% Figure: convergence
    fig = new_figure(S,6.5,3.0);
    t = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
    markers = {'o','s','d','^'};
    sets = {M.L1,'global $L_1$'; M.L2,'global $L_2$'; M.L2avg,'global $L_2$, cell-averaged exact'};
    for sI = 1:3
        ax = nexttile(t); hold(ax,'on');
        E = sets{sI,1};
        for v = 1:4
            loglog(ax,dx,E(:,v),['-' markers{v}],'Color',S.cat(v,:),'MarkerFaceColor',S.cat(v,:), ...
                'MarkerEdgeColor','w','MarkerSize',5,'LineWidth',1.2);
        end
        ref = min(E(end,:))*0.55;
        for p = [1/3 1/2 1]
            loglog(ax,dx,ref*(dx/dx(end)).^p,':','Color',S.ink2,'LineWidth',0.9,'HandleVisibility','off');
            text(ax,dx(1)*1.08,ref*(dx(1)/dx(end))^p,sprintf('$%s$',slopeText(p)),'FontSize',S.fs-2, ...
                'Color',S.ink2,'HorizontalAlignment','right');
        end
        set(ax,'XScale','log','YScale','log','XDir','reverse'); style_axes(ax,S);
        xlim(ax,[dx(end)*0.8 dx(1)*2.6]); set_dx_ticks(ax,dx);
        xlabel(ax,'$\Delta x\ [\times10^{-3}]$'); title(ax,sets{sI,2},'FontWeight','normal');
        if sI==1; lgC = legend(ax,{'$\rho$','$u$','$p$','$T$'},'Orientation','horizontal','Box','off'); lgC.Layout.Tile = 'south'; end
    end
    save_figure(fig,fullfile(figFolder,'fig_convergence.pdf'));

    %% Figure: error by wave feature
    fig = new_figure(S,6.5,2.6);
    t = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    zoneNames = {'rarefaction zone','contact zone','shock zone'};
    for vv = [1 4]
        ax = nexttile(t); hold(ax,'on');
        for z = 1:3
            loglog(ax,dx,squeeze(M.zone(:,vv,z)),['-' markers{z}],'Color',S.cat(z,:),'MarkerFaceColor',S.cat(z,:), ...
                'MarkerEdgeColor','w','MarkerSize',5,'LineWidth',1.2);
        end
        loglog(ax,dx,M.L2(:,vv),'-','Color',S.ink,'LineWidth',1.6);
        if vv==1
            yline(ax,floorRho(1),'--','Color',S.cat(2,:),'LineWidth',0.9);
            yline(ax,floorRho(2),'--','Color',S.cat(2,:),'LineWidth',0.9,'HandleVisibility','off');
        end
        set(ax,'XScale','log','YScale','log','XDir','reverse'); style_axes(ax,S);
        xlabel(ax,'$\Delta x\ [\times10^{-3}]$'); ylabel(ax,'contribution to global $L_2$');
        title(ax,sprintf('(%s) %s',char('a'+(vv>1)),ylab{vv}),'FontWeight','normal');
        if vv==1
            lgZ = legend(ax,[zoneNames,{'total','finite-$\tau$ contact floor (model, $\rho$ only)'}],'NumColumns',3,'Box','off','FontSize',S.fs-2);
            lgZ.Layout.Tile = 'south';
        end
        set_dx_ticks(ax,dx);
    end
    save_figure(fig,fullfile(figFolder,'fig_feature_errors.pdf'));

    %% Figure: contact and shock thickness
    fig = new_figure(S,4.2,3.0);
    ax = axes(fig); hold(ax,'on');
    loglog(ax,dx,M.dc,'-o','Color',S.cat(2,:),'MarkerFaceColor',S.cat(2,:),'MarkerEdgeColor','w','LineWidth',1.3);
    loglog(ax,dx,M.ds,'-s','Color',S.cat(8,:),'MarkerFaceColor',S.cat(8,:),'MarkerEdgeColor','w','LineWidth',1.3);
    patch(ax,[dx(end)*0.7 dx(1)*1.4 dx(1)*1.4 dx(end)*0.7],[dKin(1) dKin(1) dKin(2) dKin(2)],S.cat(2,:), ...
        'FaceAlpha',0.18,'EdgeColor','none');
    loglog(ax,dx,dx,':','Color',S.ink2);
    set(ax,'XScale','log','YScale','log','XDir','reverse'); style_axes(ax,S);
    xlim(ax,[dx(end)*0.7 dx(1)*1.4]);
    xlabel(ax,'$\Delta x\ [\times10^{-3}]$'); ylabel(ax,'thickness $\delta$');
    legend(ax,{sprintf('contact, slope %.2f',pc(1)),sprintf('shock, slope %.2f',pshk(1)), ...
        'kinetic contact thickness (model)','$\delta=\Delta x$'},'Location','northeast','Box','off','FontSize',S.fs-1);
    set_dx_ticks(ax,dx);
    save_figure(fig,fullfile(figFolder,'fig_feature_widths.pdf'));

    %% Figure: sub-cell shock phase and the velocity error
    fig = new_figure(S,6.5,2.7);
    t = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    ax = nexttile(t); hold(ax,'on');
    sU = 1e3*M.shockU;                                  % plotted in units of 1e-3
    loglog(ax,dx,sU,'-','Color',S.axis,'LineWidth',1.0,'HandleVisibility','off');
    grpA = [1 4]; grpB = [2 5];
    loglog(ax,dx(grpA),sU(grpA),'--o','Color',S.cat(1,:),'MarkerFaceColor',S.cat(1,:),'MarkerEdgeColor','w','LineWidth',1.2);
    loglog(ax,dx(grpB),sU(grpB),'--s','Color',S.cat(2,:),'MarkerFaceColor',S.cat(2,:),'MarkerEdgeColor','w','LineWidth',1.2);
    loglog(ax,dx(3),sU(3),'d','Color',S.cat(3,:),'MarkerFaceColor',S.cat(3,:),'MarkerEdgeColor','w');
    for k = 1:nG
        text(ax,dx(k)*0.93,sU(k),sprintf('%.2f',M.phase(k)),'FontSize',S.fs-2,'Color',S.ink2,'HorizontalAlignment','left');
    end
    ylim(ax,[2 10]);
    set(ax,'XScale','log','YScale','log','XDir','reverse'); style_axes(ax,S);
    xlim(ax,[dx(end)*0.8 dx(1)*1.25]);
    xlabel(ax,'$\Delta x\ [\times10^{-3}]$'); ylabel(ax,'shock-zone $L_2(u)$ $[\times10^{-3}]$'); set_dx_ticks(ax,dx);
    legend(ax,{sprintf('phase $\\approx0.55$: slope %.2f',log(M.shockU(1)/M.shockU(4))/log(8)), ...
        sprintf('phase $\\approx0.1$: slope %.2f',log(M.shockU(2)/M.shockU(5))/log(8)),'phase 0.26'}, ...
        'Location','southwest','Box','off','FontSize',S.fs-2);
    title(ax,'(a) Labels: sub-cell position of exact shock','FontWeight','normal');
    ax = nexttile(t); hold(ax,'on');
    loglog(ax,dx,1e3*M.L2(:,2),'-o','Color',S.cat(2,:),'MarkerFaceColor',S.cat(2,:),'MarkerEdgeColor','w','LineWidth',1.2);
    loglog(ax,dx,1e3*M.L2avg(:,2),'-s','Color',S.cat(1,:),'MarkerFaceColor',S.cat(1,:),'MarkerEdgeColor','w','LineWidth',1.2);
    loglog(ax,dx,1e3*M.L2avg(end,2)*0.8*(dx/dx(end)).^0.5,':','Color',S.ink2);
    set(ax,'XScale','log','YScale','log','XDir','reverse'); style_axes(ax,S);
    xlim(ax,[dx(end)*0.8 dx(1)*1.25]);
    xlabel(ax,'$\Delta x\ [\times10^{-3}]$'); ylabel(ax,'global $L_2(u)$ $[\times10^{-3}]$'); set_dx_ticks(ax,dx);
    legend(ax,{'vs.\ point-sampled exact','vs.\ cell-averaged exact','$\propto\Delta x^{1/2}$'},'Location','northeast','Box','off','FontSize',S.fs-2);
    title(ax,'(b) Choice of exact reference','FontWeight','normal');
    save_figure(fig,fullfile(figFolder,'fig_shock_phase.pdf'));

    %% Limiter comparison
    limFile = fullfile(projectFolder,'results','d1v5_limiter_history_Nx6250_tau5e-06_CFL0p025_results.mat');
    if isfile(limFile)
        L = load(limFile);
        labs = {'Minmod','MC','Population hybrid','Macroscopic hybrid','Gen.\ minmod $\theta=1.2$'};
        fid = fopen(fullfile(dataFolder,'tab_limiters.tex'),'w');
        et = L.errorTable;
        for k = 1:height(et)
            fprintf(fid,'%s & %s & %s & %s & %s & %s & %s \\\\\n',labs{k},sci(et.GlobalL1_rho(k),3), ...
                sci(et.GlobalL1_u(k),3),sci(et.GlobalL1_p(k),3),sci(et.GlobalL1_T(k),3), ...
                sci(et.WaveL1_T(k),3),sci(et.T_peak_excess(k),3));
        end
        fclose(fid);
        fig = new_figure(S,6.5,4.6);
        t = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
        styles = {'-','-','-.',':','--'}; lw = [1.2 1.2 1.3 1.6 1.8];
        zoomsL = {'T',[x0-0.3 x0+0.4],'(a) $T$, wave window'; 'T',[x0+0.125 x0+0.155],'(b) $T$ at contact'; ...
                  'T',[x0+0.245 x0+0.275],'(c) $T$ at shock'; 'rho',[x0+0.125 x0+0.155],'(d) $\rho$ at contact'};
        for zI = 1:4
            ax = nexttile(t); hold(ax,'on');
            fld = zoomsL{zI,1}; xr = zoomsL{zI,2}; m = L.x>=xr(1)-0.005 & L.x<=xr(2)+0.005;
            plot(ax,L.x(m),L.exact.(fld)(m),'-','Color',S.ink,'LineWidth',1.0);
            for k = 1:5
                plot(ax,L.x(m),L.results(k).(fld)(m),styles{k},'Color',S.cat(k,:),'LineWidth',lw(k));
            end
            xlim(ax,xr); style_axes(ax,S); xlabel(ax,'$x$'); title(ax,zoomsL{zI,3},'FontWeight','normal');
            if zI==1; legend(ax,[{'Exact'},labs],'Location','northwest','Box','off','FontSize',S.fs-2); end
        end
        save_figure(fig,fullfile(figFolder,'fig_limiters.pdf'));
    else
        fprintf('Limiter comparison not found yet.\n');
    end

    %% Theta sweep
    thetaFiles = dir(fullfile(projectFolder,'results','theta_sweep','*.mat'));
    if ~isempty(thetaFiles)
        th = []; E1 = []; wT = []; pk = []; E2T = [];
        allRuns = {};
        for k = 1:numel(thetaFiles)
            tmp = load(fullfile(thetaFiles(k).folder,thetaFiles(k).name),'sweep');
            allRuns{end+1} = tmp.sweep; %#ok<AGROW>
        end
        base = runs{Nx==6250};                          % theta = 1.2 from the convergence study
        allRuns{end+1} = struct('limiterTheta',1.2,'x',base.x,'rho',base.rho,'u',base.u,'p',base.p, ...
            'T',base.T,'exact',base.exact,'retryCountTotal',base.retryCountTotal,'fallbackCountTotal',base.fallbackCountTotal);
        for k = 1:numel(allRuns)
            r = allRuns{k}; w = r.x >= x0-0.30 & r.x <= x0+0.40;
            th(k) = r.limiterTheta; %#ok<AGROW>
            for v = 1:4; E1(k,v) = mean(abs(r.(vars{v})-r.exact.(vars{v}))); end %#ok<AGROW>
            wT(k) = mean(abs(r.T(w)-r.exact.T(w))); %#ok<AGROW>
            E2T(k) = sqrt(mean((r.T-r.exact.T).^2)); %#ok<AGROW>
            pk(k) = max(r.T(w))-max(r.exact.T(w)); %#ok<AGROW>
            fb(k) = r.fallbackCountTotal; %#ok<AGROW>
        end
        [th,o] = sort(th); E1 = E1(o,:); wT = wT(o); pk = pk(o); E2T = E2T(o); fb = fb(o);
        fid = fopen(fullfile(dataFolder,'tab_theta.tex'),'w');
        for k = 1:numel(th)
            if abs(th(k)-1.2)<1e-9; pre = '\textbf{'; post = '}'; else; pre=''; post=''; end
            fprintf(fid,'%s%.2f%s & %s & %s & %s & %s & %s & %s & %d \\\\\n',pre,th(k),post,sci(E1(k,1),3),sci(E1(k,2),3), ...
                sci(E1(k,3),3),sci(E1(k,4),3),sci(wT(k),3),sci(pk(k),3),fb(k));
        end
        fclose(fid);
        save(fullfile(dataFolder,'theta_sweep_metrics.mat'),'th','E1','wT','pk','E2T');
        fig = new_figure(S,6.5,2.8);
        t = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
        ax = nexttile(t); hold(ax,'on');
        i1 = find(abs(th-1)<1e-9,1);
        for v = 1:4
            plot(ax,th,E1(:,v)/E1(i1,v),['-' markers{v}],'Color',S.cat(v,:),'MarkerFaceColor',S.cat(v,:), ...
                'MarkerEdgeColor','w','MarkerSize',5,'LineWidth',1.2);
        end
        xline(ax,1.2,':','Color',S.ink2,'HandleVisibility','off');
        style_axes(ax,S); xlabel(ax,'$\theta$'); ylabel(ax,'global $L_1$ / minmod value');
        legend(ax,{'$\rho$','$u$','$p$','$T$'},'Location','northeast','Box','off');
        title(ax,'(a) Average error versus $\theta$','FontWeight','normal');
        ax = nexttile(t); hold(ax,'on');
        px = 1e3*pk; py = 1e4*E1(:,4);
        plot(ax,px,py,'-','Color',S.seq(3,:),'LineWidth',1.2);
        scatter(ax,px,py,28,S.seq(4,:),'filled','MarkerEdgeColor','w');
        for k = 1:numel(th)
            text(ax,px(k),py(k),sprintf('  %.1f',th(k)),'FontSize',S.fs-2,'Color',S.ink2, ...
                'VerticalAlignment','bottom');
        end
        k12 = find(abs(th-1.2)<1e-9,1);
        plot(ax,px(k12),py(k12),'o','MarkerSize',9,'MarkerEdgeColor',S.cat(2,:),'LineWidth',1.4);
        xline(ax,0,'-','Color',S.axis);
        style_axes(ax,S); xlabel(ax,'temperature peak excess $[\times10^{-3}]$'); ylabel(ax,'global $L_1(T)$ $[\times10^{-4}]$');
        title(ax,'(b) Diffusion--overshoot tradeoff (labels: $\theta$)','FontWeight','normal');
        save_figure(fig,fullfile(figFolder,'fig_theta_sweep.pdf'));
    end
end

%% ---------------------------------------------------------------------
function [delta,xMid] = thickness(x,q,xr,qLeft,qRight)
% maximum-gradient thickness |jump|/max|dq/dx| and mid-level crossing location
    m = find(x>=xr(1) & x<=xr(2));
    xs = x(m); qs = q(m);
    g = abs(gradient(qs,xs));
    delta = abs(qLeft-qRight)/max(g);
    level = 0.5*(qLeft+qRight);
    s = sign(qs-level); idx = find(s(1:end-1).*s(2:end) <= 0,1,'last');
    xMid = xs(idx) + (level-qs(idx))*(xs(idx+1)-xs(idx))/(qs(idx+1)-qs(idx));
end

function ex = cell_averaged_exact(x,dx,x0,t,g,ps,us,rsL,rsR,S,head,tail)
% average the exact solution over each cell with 64 sub-samples near the waves
    ex = sod_exact(x,t,x0,g,ps,us,rsL,rsR,S,head,tail);
    m = find(x >= x0+head*t-2*dx & x <= x0+S*t+2*dx);
    nSub = 64; offs = ((1:nSub)-0.5)/nSub - 0.5;
    xs = x(m) + dx*offs;                                % cells-by-subsamples
    sub = sod_exact(xs(:),t,x0,g,ps,us,rsL,rsR,S,head,tail);
    E = 0.5*sub.rho.*(5*sub.T) + 0.5*sub.rho.*sub.u.^2; % average conserved variables
    rhoA = mean(reshape(sub.rho,size(xs)),2);
    momA = mean(reshape(sub.rho.*sub.u,size(xs)),2);
    EA = mean(reshape(E,size(xs)),2);
    ex.rho(m) = rhoA; ex.u(m) = momA./rhoA;
    ex.T(m) = (2*EA./rhoA - ex.u(m).^2)/5; ex.p(m) = ex.rho(m).*ex.T(m);
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

function write_err_table(fn,Nx,E,extra)
    fid = fopen(fn,'w');
    for k = 1:numel(Nx)
        fprintf(fid,'%s',commas(Nx(k)));
        for v = 1:size(E,2); fprintf(fid,' & %s',sci(E(k,v),3)); end
        if ~isempty(extra); fprintf(fid,' & %s',sci(extra(k),3)); end
        fprintf(fid,' \\\\\n');
    end
    fclose(fid);
end

function s = sci(v,d)
% LaTeX scientific notation, e.g. 1.23\times10^{-3}
    if v == 0; s = '$0$'; return; end
    e = floor(log10(abs(v))); m = v/10^e;
    if abs(round(m,d-1)) >= 10; m = m/10; e = e+1; end
    s = sprintf(['$%.' num2str(d-1) 'f\\times10^{%d}$'],m,e);
end

function s = commas(n)
    s = regexprep(sprintf('%d',round(n)),'(\d)(?=(\d{3})+$)','$1,');
end

function set_dx_ticks(ax,dx)
% label the x axis at the actual grid spacings
    labels = arrayfun(@(d) sprintf('%.1f',1e3*d),dx,'UniformOutput',false);   % in units of 1e-3
    [dsorted,o] = sort(dx);
    set(ax,'XTick',dsorted,'XTickLabel',labels(o),'XMinorTick','off','XMinorGrid','off','XTickLabelRotation',0);
    ax.XLabel.String = '$\Delta x\ [\times10^{-3}]$';
end

function s = slopeText(p)
    if abs(p-1/3)<1e-9; s = '\propto\Delta x^{1/3}';
    elseif abs(p-1/2)<1e-9; s = '\propto\Delta x^{1/2}';
    else; s = '\propto\Delta x'; end
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
