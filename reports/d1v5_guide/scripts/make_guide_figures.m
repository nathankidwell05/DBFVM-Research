function make_guide_figures(onlyThese)
%% Figures and Tables for the KT-D1V5 Solver Guide %%
% Every figure is created in MATLAB and saved three ways in the "figures"
% directory next to the report folder (the research folder on OneDrive):
%     <name>.fig  (open and edit in MATLAB)
%     <name>.png  (quick viewing)
%     <name>.pdf  (vector version used by the LaTeX document)
% Table rows for the document are written to ../data/.
%
% Optional input: a cell array of figure names to rebuild only those, e.g.
%     make_guide_figures({'G10_accuracy_vs_cost'})

    if nargin < 1; onlyThese = {}; end
    here = fileparts(mfilename('fullpath'));
    projectFolder = project_paths();
    figFolder = fullfile(fileparts(fileparts(here)),'figures');  % research folder (OneDrive) or reports/ (GitHub)
    dataFolder = fullfile(fileparts(here),'data');
    if ~isfolder(figFolder); mkdir(figFolder); end
    if ~isfolder(dataFolder); mkdir(dataFolder); end
    S = report_style();
    want = @(name) isempty(onlyThese) || any(strcmp(onlyThese,name));

    gamma = 1.4; b = 5; v1 = 1; v2 = 3; eta0 = 1.75; c = [-v2,-v1,0,v1,v2];
    x0 = 10; tEnd = 0.15;
    [ps,us,rsL,rsR,Ssh,head,tail] = sod_constants(gamma);
    grids = [3125 6250 12500 25000 50000];
    vars = {'rho','u','p','T'}; ylab = {'$\rho$','$u$','$p$','$T$'};
    markers = {'o','s','d','^'};

    %% Load data --------------------------------------------------------------
    runs = cell(1,numel(grids));
    for k = 1:numel(grids)
        tmp = load(fullfile(projectFolder,'results','runs', ...
            sprintf('d1v5_gminmod_theta1p20_tau5e-06_CFL0p025_Nx%d.mat',grids(k))),'runData');
        runs{k} = tmp.runData;
    end
    dx = cellfun(@(r) r.dx,runs);
    M = metrics_for(runs,x0);

    %% G01 velocity set ------------------------------------------------------
    if want('G01_velocity_set')
        fig = new_figure(S,6.0,2.1); ax = axes(fig,'Position',[0.02 0.05 0.96 0.9]); hold(ax,'on');
        plot(ax,[-4.2 4.2],[0 0],'-','Color',S.ink2,'LineWidth',1.0);
        for k = 1:5
            plot(ax,c(k),0,'o','MarkerSize',8,'MarkerFaceColor',S.cat(1,:),'MarkerEdgeColor','w');
            if c(k)==0; lbl = sprintf('$c_%d = 0$',k); else; lbl = sprintf('$c_%d = %+d$',k,c(k)); end
            text(ax,c(k),-0.32,lbl,'HorizontalAlignment','center','FontSize',S.fs);
            text(ax,c(k),0.62,sprintf('$f_%d$',k),'HorizontalAlignment','center','FontSize',S.fs+1);
            if c(k) ~= 0
                quiver(ax,c(k),0.25,0.3*c(k),0,0,'Color',S.cat(1,:),'LineWidth',1.6,'MaxHeadSize',0.9);
            end
        end
        plot(ax,0,0,'o','MarkerSize',16,'MarkerEdgeColor',S.cat(2,:),'LineWidth',1.3);
        text(ax,0,-0.72,'rest population: carries internal energy $\eta_0^2$','HorizontalAlignment','center', ...
            'FontSize',S.fs-1,'Color',S.cat(2,:));
        text(ax,-2,1.0,'left-moving','HorizontalAlignment','center','FontSize',S.fs-1,'Color',S.ink2);
        text(ax,2,1.0,'right-moving','HorizontalAlignment','center','FontSize',S.fs-1,'Color',S.ink2);
        axis(ax,'off'); xlim(ax,[-4.3 4.3]); ylim(ax,[-0.95 1.15]);
        save_all(fig,figFolder,'G01_velocity_set');
    end

    %% G02 exact Sod solution --------------------------------------------------
    if want('G02_exact_sod')
        xx = linspace(9.6,10.5,4000)';
        ex = sod_exact(xx,tEnd,x0,gamma);
        fig = new_figure(S,6.5,4.4);
        t = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
        tiles = [1 2 4 5]; xl = x0+[head tail us Ssh]*tEnd;
        for kk = 1:4
            ax = nexttile(t,tiles(kk)); hold(ax,'on');
            for m = 1:4; xline(ax,xl(m),':','Color',S.axis); end
            plot(ax,xx,ex.(vars{kk}),'-','Color',S.ink,'LineWidth',S.lw);
            xlim(ax,[9.6 10.5]); style_axes(ax,S); ylabel(ax,ylab{kk});
            if kk>2; xlabel(ax,'$x$'); end
        end
        ax = nexttile(t,3,[2 1]); hold(ax,'on');
        tt = [0 tEnd];
        for s = linspace(head,tail,7); plot(ax,x0+s*tt,tt,'-','Color',S.seq(2,:),'LineWidth',0.9); end
        plot(ax,x0+us*tt,tt,'--','Color',S.cat(2,:),'LineWidth',S.lw);
        plot(ax,x0+Ssh*tt,tt,'-','Color',S.cat(8,:),'LineWidth',S.lw+0.4);
        patch(ax,[9.7 10.4 10.4 9.7],[tEnd tEnd tEnd*1.04 tEnd*1.04],S.band,'EdgeColor','none');
        text(ax,10.05,tEnd*1.02,'error window','HorizontalAlignment','center','FontSize',S.fs-2,'Color',S.ink2);
        text(ax,9.62,0.125,'rarefaction','FontSize',S.fs-2,'Color',S.seq(4,:));
        text(ax,x0+us*0.13-0.012,0.13,'contact','HorizontalAlignment','right','FontSize',S.fs-2,'Color',S.cat(2,:));
        text(ax,x0+Ssh*0.06+0.03,0.06,'shock','FontSize',S.fs-2,'Color',S.cat(8,:));
        xlim(ax,[9.6 10.5]); ylim(ax,[0 tEnd*1.04]); style_axes(ax,S);
        xlabel(ax,'$x$'); ylabel(ax,'$t$'); title(ax,'Wave paths','FontWeight','normal');
        save_all(fig,figFolder,'G02_exact_sod');
    end

    %% G03 equilibrium populations -------------------------------------------
    if want('G03_equilibrium')
        fig = new_figure(S,6.5,2.7);
        t = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
        labels = {'$f_1^{eq}$ ($c=-3$)','$f_2^{eq}$ ($c=-1$)','$f_3^{eq}$ ($c=0$)','$f_4^{eq}$ ($c=+1$)','$f_5^{eq}$ ($c=+3$)'};
        styles = {'-','--','-','--','-'};
        ax = nexttile(t); hold(ax,'on');
        T = linspace(0.3,1.5,300)';
        F = cell2mat(arrayfun(@(TT) feq_kt(1,0,TT,b,v1,v2,eta0,c),T,'UniformOutput',false));
        patch(ax,[0.711 1.141 1.141 0.711],[-0.6 -0.6 2 2],S.band,'EdgeColor','none','HandleVisibility','off');
        for q = 1:5; plot(ax,T,F(:,q),styles{q},'Color',S.cat(q,:),'LineWidth',S.lw); end
        yline(ax,0,'-','Color',S.axis,'HandleVisibility','off');
        text(ax,0.93,1.75,'Sod range of $T$','HorizontalAlignment','center','FontSize',S.fs-1,'Color',S.ink2);
        xlabel(ax,'$T$ (with $u=0$)'); ylabel(ax,'$f_i^{eq}/\rho$'); ylim(ax,[-0.6 2]); xlim(ax,[0.3 1.5]); style_axes(ax,S);
        title(ax,'(a) Gas at rest, varying temperature','FontWeight','normal');
        ax = nexttile(t); hold(ax,'on');
        U = linspace(-1.2,1.2,300)';
        F = cell2mat(arrayfun(@(uu) feq_kt(1,uu,1,b,v1,v2,eta0,c),U,'UniformOutput',false));
        for q = 1:5; plot(ax,U,F(:,q),styles{q},'Color',S.cat(q,:),'LineWidth',S.lw); end
        yline(ax,0,'-','Color',S.axis,'HandleVisibility','off');
        xline(ax,us,':','Color',S.ink2,'HandleVisibility','off');
        text(ax,us-0.03,1.75,'$u^*$','HorizontalAlignment','right','FontSize',S.fs-1,'Color',S.ink2);
        xlabel(ax,'$u$ (with $T=1$)'); ylabel(ax,'$f_i^{eq}/\rho$'); ylim(ax,[-0.6 2]); style_axes(ax,S);
        title(ax,'(b) $T=1$, varying flow velocity','FontWeight','normal');
        lg = legend(ax,labels,'NumColumns',5,'Box','off'); lg.Layout.Tile = 'north';
        save_all(fig,figFolder,'G03_equilibrium');
    end

    %% G04 MUSCL schematic -----------------------------------------------------
    if want('G04_muscl_schematic')
        fig = new_figure(S,6.0,2.9); ax = axes(fig,'Position',[0.04 0.14 0.92 0.8]); hold(ax,'on');
        avg = [0.6 0.9 1.9 2.2];
        slope = [0.2 0.3 0.36 0.1];
        for j = 1:4
            plot(ax,[j-1 j-1],[0 2.7],'-','Color',S.grid,'LineWidth',1);
            plot(ax,[j-1 j],[avg(j) avg(j)],'-','Color',S.ink,'LineWidth',1.3);
            plot(ax,[j-1 j],avg(j)+slope(j)*[-0.5 0.5],'-','Color',S.cat(1,:),'LineWidth',2.2);
            plot(ax,j-0.5,avg(j),'o','MarkerSize',5,'MarkerFaceColor',S.ink,'MarkerEdgeColor',S.ink);
        end
        plot(ax,[4 4],[0 2.7],'-','Color',S.grid,'LineWidth',1);
        fl = avg(2)+0.5*slope(2); fr = avg(3)-0.5*slope(3);
        plot(ax,2,fl,'o','MarkerSize',8,'MarkerFaceColor',S.cat(2,:),'MarkerEdgeColor','w');
        plot(ax,2,fr,'o','MarkerSize',8,'MarkerFaceColor',S.cat(2,:),'MarkerEdgeColor','w');
        text(ax,2.07,fl-0.12,'$f^{-}_{j+1/2}$ (from cell $j$)','FontSize',S.fs,'Color',S.cat(2,:));
        text(ax,1.93,fr+0.14,'$f^{+}_{j+1/2}$ (from cell $j{+}1$)','FontSize',S.fs,'Color',S.cat(2,:),'HorizontalAlignment','right');
        names = {'$j-1$','$j$','$j+1$','$j+2$'};
        for j = 1:4; text(ax,j-0.5,-0.18,names{j},'HorizontalAlignment','center','FontSize',S.fs); end
        text(ax,2,2.85,'face $x_{j+1/2}$','HorizontalAlignment','center','FontSize',S.fs-1,'Color',S.ink2);
        axis(ax,'off'); xlim(ax,[-0.1 4.1]); ylim(ax,[-0.3 3.0]);
        save_all(fig,figFolder,'G04_muscl_schematic');
    end

    %% G05 limiter functions ---------------------------------------------------
    if want('G05_limiters')
        fig = new_figure(S,6.5,2.9);
        t = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
        ax = nexttile(t); hold(ax,'on');
        r = linspace(0,3.5,700);
        patch(ax,[r fliplr(r)],[zeros(size(r)) fliplr(min(2*r,2))],S.band,'EdgeColor','none');
        phiGM = @(r,th) max(0,min(min(th*r,(1+r)/2),th));
        th = [1 1.2 1.5 2]; cols = S.seq([1 3 4 5],:); lws = [S.lw S.lw+0.9 S.lw S.lw];
        for k = 1:4; plot(ax,r,phiGM(r,th(k)),'Color',cols(k,:),'LineWidth',lws(k)); end
        plot(ax,1,1,'o','MarkerSize',6,'MarkerFaceColor',S.ink,'MarkerEdgeColor','w','HandleVisibility','off');
        xlabel(ax,'smoothness ratio $r$'); ylabel(ax,'slope factor $\phi(r)$');
        xlim(ax,[0 3.5]); ylim(ax,[0 2.2]); style_axes(ax,S);
        legend(ax,{'no-new-oscillation region','minmod ($\theta=1$)','$\theta=1.2$ (this solver)','$\theta=1.5$','MC ($\theta=2$)'}, ...
            'Location','southeast','Box','off','FontSize',S.fs-1);
        title(ax,'(a) Limiter family','FontWeight','normal');
        ax = nexttile(t); hold(ax,'on');
        thv = 1.2; rr = linspace(0.001,3.5,2000);
        [~,active] = min([thv*rr; (1+rr)/2; thv*ones(size(rr))],[],1);
        nm = {'left-sided slope $\theta\,d_L$','centred slope $(d_L+d_R)/2$','right-sided slope $\theta\,d_R$'};
        for k = 1:3; m = active==k; plot(ax,rr(m),phiGM(rr(m),thv),'.','Color',S.cat(k,:),'MarkerSize',7); end
        xline(ax,1/(2*thv-1),':','Color',S.ink2,'HandleVisibility','off'); xline(ax,2*thv-1,':','Color',S.ink2,'HandleVisibility','off');
        xlabel(ax,'smoothness ratio $r$'); ylabel(ax,'$\phi(r)$ for $\theta=1.2$');
        xlim(ax,[0 3.5]); ylim(ax,[0 2.2]); style_axes(ax,S);
        legend(ax,nm,'Location','southeast','Box','off','FontSize',S.fs-2);
        title(ax,'(b) Candidate selected at $\theta=1.2$','FontWeight','normal');
        save_all(fig,figFolder,'G05_limiters');
    end

    %% G06 finest-grid solution ------------------------------------------------
    if want('G06_finest_solution')
        rf = runs{end}; w = rf.x >= x0-0.30 & rf.x <= x0+0.40;
        fig = new_figure(S,6.5,5.2);
        t = tiledlayout(fig,4,2,'TileSpacing','tight','Padding','compact');
        for v = 1:4
            ax = nexttile(t,2*v-1); hold(ax,'on');
            plot(ax,rf.x(w),rf.exact.(vars{v})(w),'-','Color',S.ink,'LineWidth',1.0);
            plot(ax,rf.x(w),rf.(vars{v})(w),'--','Color',S.cat(1,:),'LineWidth',S.lw);
            xlim(ax,[x0-0.3 x0+0.4]); ylabel(ax,ylab{v}); style_axes(ax,S);
            if v<4; ax.XTickLabel = []; else; xlabel(ax,'$x$'); end
            if v==1; legend(ax,{'exact','DBM, $N_x=50{,}000$'},'Location','northeast','Box','off'); title(ax,'Solution','FontWeight','normal'); end
            ax = nexttile(t,2*v); hold(ax,'on');
            plot(ax,rf.x(w),rf.(vars{v})(w)-rf.exact.(vars{v})(w),'-','Color',S.cat(2,:),'LineWidth',1.0);
            yline(ax,0,'-','Color',S.axis);
            xlim(ax,[x0-0.3 x0+0.4]); ylabel(ax,['error in ' ylab{v}]); style_axes(ax,S);
            if v<4; ax.XTickLabel = []; else; xlabel(ax,'$x$'); end
            if v==1; title(ax,'Pointwise error (computed $-$ exact)','FontWeight','normal'); end
        end
        save_all(fig,figFolder,'G06_finest_solution');
    end

    %% G07 grid refinement close-ups ------------------------------------------
    if want('G07_grid_zoom')
        rf = runs{end};
        fig = new_figure(S,6.5,4.6);
        t = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
        zooms = {'rho',[x0-0.02 x0+0.03],'(a) Density at the rarefaction tail'; ...
                 'rho',[x0+0.12 x0+0.16],'(b) Density at the contact'; ...
                 'rho',[x0+0.245 x0+0.28],'(c) Density at the shock'; ...
                 'T',[x0+0.12 x0+0.28],'(d) Temperature, contact to shock'};
        gridNames = arrayfun(@(n) sprintf('$N_x=%s$',commas(n)),grids,'UniformOutput',false);
        for zI = 1:4
            ax = nexttile(t); hold(ax,'on');
            fld = zooms{zI,1}; xr = zooms{zI,2};
            m = rf.x>=xr(1) & rf.x<=xr(2);
            plot(ax,rf.x(m),rf.exact.(fld)(m),'-','Color',S.ink,'LineWidth',1.0);
            for k = 1:numel(runs)
                r = runs{k}; m = r.x>=xr(1)-0.01 & r.x<=xr(2)+0.01;
                plot(ax,r.x(m),r.(fld)(m),'-','Color',S.seq(k,:),'LineWidth',1.2);
            end
            xlim(ax,xr); style_axes(ax,S); title(ax,zooms{zI,3},'FontWeight','normal'); xlabel(ax,'$x$');
            if zI==1; legend(ax,[{'exact'},gridNames],'Location','northeast','Box','off','FontSize',S.fs-2); end
        end
        save_all(fig,figFolder,'G07_grid_zoom');
    end

    %% G08 convergence of the DBM solver ----------------------------------------
    if want('G08_convergence')
        fig = new_figure(S,6.5,3.0);
        t = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
        sets = {M.L1,'(a) Mean absolute error, $L_1$'; M.L2,'(b) Root-mean-square error, $L_2$'};
        for sI = 1:2
            ax = nexttile(t); hold(ax,'on');
            E = sets{sI,1};
            for v = 1:4
                plot(ax,dx,E(:,v),['-' markers{v}],'Color',S.cat(v,:),'MarkerFaceColor',S.cat(v,:), ...
                    'MarkerEdgeColor','w','MarkerSize',5,'LineWidth',1.2);
            end
            set(ax,'XScale','log','YScale','log','XDir','reverse'); style_axes(ax,S);
            xlim(ax,[dx(end)*0.38 dx(1)*1.25]); set_dx_ticks(ax,dx);
            yl = [min(E(:))*0.2 max(E(:))*1.8]; ylim(ax,yl);
            draw_reference_slopes(ax,dx,min(E(1,:))*0.55,[2 1 0.5],S,yl);
            title(ax,sets{sI,2},'FontWeight','normal');
            if sI==1; lg = legend(ax,{'$\rho$','$u$','$p$','$T$'},'Orientation','horizontal','Box','off'); lg.Layout.Tile = 'south'; end
        end
        save_all(fig,figFolder,'G08_convergence');
    end

    %% Method comparison data ------------------------------------------------
    methodsReady = isfile(fullfile(projectFolder,'results','method_runs','timing_benchmark.mat'));
    if methodsReady
        C = method_comparison(projectFolder,runs,grids,x0);
        write_method_tables(C,dataFolder);
    end

    %% G09 order comparison ---------------------------------------------------
    if want('G09_order_comparison') && methodsReady
        fig = new_figure(S,6.5,3.1);
        t = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
        for sI = 1:2
            ax = nexttile(t); hold(ax,'on');
            for m = 1:numel(C)
                E = C(m).L1rho; if sI==2; E = C(m).L2rho; end
                plot(ax,C(m).dx,E,C(m).style,'Color',C(m).color,'MarkerFaceColor',C(m).fill,'MarkerEdgeColor',C(m).color, ...
                    'MarkerSize',5,'LineWidth',C(m).lw);
            end
            allE = cell2mat(arrayfun(@(m) m.L1rho(:)',C,'UniformOutput',false)); if sI==2; allE = cell2mat(arrayfun(@(m) m.L2rho(:)',C,'UniformOutput',false)); end
            set(ax,'XScale','log','YScale','log','XDir','reverse'); style_axes(ax,S);
            xlim(ax,[dx(end)*0.38 dx(1)*1.25]); set_dx_ticks(ax,dx);
            yl = [min(allE)*0.2 max(allE)*1.8]; ylim(ax,yl);
            coarse = arrayfun(@(m) C(m).L1rho(1),1:numel(C)); if sI==2; coarse = arrayfun(@(m) C(m).L2rho(1),1:numel(C)); end
            draw_reference_slopes(ax,dx,min(coarse)*0.55,[1 0.5],S,yl);
            if sI==1; title(ax,'(a) Density, $L_1$','FontWeight','normal'); else; title(ax,'(b) Density, $L_2$','FontWeight','normal'); end
            if sI==1; lg = legend(ax,{C.label},'NumColumns',2,'Box','off','FontSize',S.fs-1); lg.Layout.Tile = 'south'; end
        end
        save_all(fig,figFolder,'G09_order_comparison');
    end

    %% G10 accuracy versus cost ------------------------------------------------
    if want('G10_accuracy_vs_cost') && methodsReady
        fig = new_figure(S,6.5,3.4);
        t = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
        ax = nexttile(t); hold(ax,'on');
        for m = 1:numel(C)
            plot(ax,C(m).time,C(m).L1rho,C(m).style,'Color',C(m).color,'MarkerFaceColor',C(m).fill,'MarkerEdgeColor',C(m).color, ...
                'MarkerSize',5,'LineWidth',C(m).lw);
        end
        set(ax,'XScale','log','YScale','log'); style_axes(ax,S);
        xlabel(ax,'computing time [s]'); ylabel(ax,'density error $L_1(\rho)$');
        title(ax,'(a) Accuracy versus computing time','FontWeight','normal');
        lg = legend(ax,{C.label},'NumColumns',2,'Box','off','FontSize',S.fs-1); lg.Layout.Tile = 'south';
        ax = nexttile(t); hold(ax,'on');
        iD = find(strcmp({C.name},'gminmod') & strcmp({C.kind},'dbm'));
        iE = find(strcmp({C.name},'muscl') & strcmp({C.kind},'euler'));
        plot(ax,C(iD).Nx,C(iD).steps,C(iD).style,'Color',C(iD).color,'MarkerFaceColor',C(iD).fill,'LineWidth',1.6,'MarkerSize',5,'HandleVisibility','off');
        plot(ax,C(iE).Nx,C(iE).steps,C(iE).style,'Color',C(iE).color,'MarkerFaceColor',C(iE).fill,'LineWidth',1.4,'MarkerSize',5,'HandleVisibility','off');
        text(ax,C(iD).Nx(2),C(iD).steps(2)*1.6,'DBM: $\Delta t=\tau/4$ (collision limit)','FontSize',S.fs-1,'Color',C(iD).color);
        text(ax,C(iE).Nx(2),C(iE).steps(2)*0.45,'Euler: $\Delta t$ from CFL $=0.5$','FontSize',S.fs-1,'Color',C(iE).color);
        set(ax,'XScale','log','YScale','log'); style_axes(ax,S);
        xlim(ax,[2500 62000]); ylim(ax,[50 5e5]);
        set(ax,'XTick',[3125 6250 12500 25000 50000],'XTickLabel',{'3,125','6,250','12,500','25,000','50,000'},'XMinorTick','off');
        xlabel(ax,'number of cells $N_x$'); ylabel(ax,'number of time steps');
        title(ax,'(b) Time steps required','FontWeight','normal');
        save_all(fig,figFolder,'G10_accuracy_vs_cost');
    end

    %% G11 shock position within a cell ---------------------------------------
    if want('G11_shock_phase')
        fig = new_figure(S,6.5,2.8);
        t = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
        ax = nexttile(t); hold(ax,'on');
        sU = 1e3*M.shockU;
        plot(ax,dx,sU,'-','Color',S.axis,'LineWidth',1.0,'HandleVisibility','off');
        plot(ax,dx([1 4]),sU([1 4]),'--o','Color',S.cat(1,:),'MarkerFaceColor',S.cat(1,:),'MarkerEdgeColor','w','LineWidth',1.2);
        plot(ax,dx([2 5]),sU([2 5]),'--s','Color',S.cat(2,:),'MarkerFaceColor',S.cat(2,:),'MarkerEdgeColor','w','LineWidth',1.2);
        plot(ax,dx(3),sU(3),'d','Color',S.cat(3,:),'MarkerFaceColor',S.cat(3,:),'MarkerEdgeColor','w');
        for k = 1:numel(dx)
            text(ax,dx(k)*0.93,sU(k),sprintf('%.2f',M.phase(k)),'FontSize',S.fs-2,'Color',S.ink2);
        end
        set(ax,'XScale','log','YScale','log','XDir','reverse'); style_axes(ax,S);
        xlim(ax,[dx(end)*0.8 dx(1)*1.25]); ylim(ax,[2 10]); set_dx_ticks(ax,dx);
        ylabel(ax,'shock-region $L_2(u)$ $[\times10^{-3}]$');
        legend(ax,{sprintf('shock near cell centre: order %.2f',log(M.shockU(1)/M.shockU(4))/log(8)), ...
            sprintf('shock near cell face: order %.2f',log(M.shockU(2)/M.shockU(5))/log(8)),'intermediate position'}, ...
            'Location','southwest','Box','off','FontSize',S.fs-2);
        title(ax,'(a) Labels: shock position within its cell','FontWeight','normal');
        ax = nexttile(t); hold(ax,'on');
        plot(ax,dx,1e3*M.L2(:,2),'-o','Color',S.cat(2,:),'MarkerFaceColor',S.cat(2,:),'MarkerEdgeColor','w','LineWidth',1.2);
        plot(ax,dx,1e4*M.L1(:,2),'-s','Color',S.cat(1,:),'MarkerFaceColor',S.cat(1,:),'MarkerEdgeColor','w','LineWidth',1.2);
        set(ax,'XScale','log','YScale','log','XDir','reverse'); style_axes(ax,S);
        xlim(ax,[dx(end)*0.8 dx(1)*1.25]); set_dx_ticks(ax,dx);
        ylabel(ax,'velocity error');
        legend(ax,{'$L_2(u)$ $[\times10^{-3}]$','$L_1(u)$ $[\times10^{-4}]$'},'Location','southwest','Box','off','FontSize',S.fs-1);
        title(ax,'(b) $L_2$ and $L_1$ velocity errors','FontWeight','normal');
        save_all(fig,figFolder,'G11_shock_phase');
    end

    %% G12 limiter comparison profiles ----------------------------------------
    limFile = fullfile(projectFolder,'results','d1v5_limiter_history_Nx6250_tau5e-06_CFL0p025_results.mat');
    if want('G12_limiter_profiles') && isfile(limFile)
        L = load(limFile);
        labs = {'minmod','MC','population hybrid','macroscopic hybrid','generalized minmod $\theta=1.2$'};
        fig = new_figure(S,6.5,4.6);
        t = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
        styles = {'-','-','-.',':','--'}; lw = [1.2 1.2 1.3 1.6 1.8];
        zoomsL = {'T',[x0-0.3 x0+0.4],'(a) Temperature, full wave region'; 'T',[x0+0.125 x0+0.155],'(b) Temperature at the contact'; ...
                  'T',[x0+0.245 x0+0.275],'(c) Temperature at the shock'; 'rho',[x0+0.125 x0+0.155],'(d) Density at the contact'};
        for zI = 1:4
            ax = nexttile(t); hold(ax,'on');
            fld = zoomsL{zI,1}; xr = zoomsL{zI,2}; m = L.x>=xr(1)-0.005 & L.x<=xr(2)+0.005;
            plot(ax,L.x(m),L.exact.(fld)(m),'-','Color',S.ink,'LineWidth',1.0);
            for k = 1:5; plot(ax,L.x(m),L.results(k).(fld)(m),styles{k},'Color',S.cat(k,:),'LineWidth',lw(k)); end
            xlim(ax,xr); style_axes(ax,S); xlabel(ax,'$x$'); title(ax,zoomsL{zI,3},'FontWeight','normal');
            if zI==1; legend(ax,[{'exact'},labs],'Location','northwest','Box','off','FontSize',S.fs-2); end
        end
        save_all(fig,figFolder,'G12_limiter_profiles');
    end

    %% G13 theta sweep ----------------------------------------------------------
    thetaFiles = dir(fullfile(projectFolder,'results','theta_sweep','*.mat'));
    if want('G13_theta_tradeoff') && ~isempty(thetaFiles)
        th = []; E1 = []; pk = [];
        for k = 1:numel(thetaFiles)
            tmp = load(fullfile(thetaFiles(k).folder,thetaFiles(k).name),'sweep'); r = tmp.sweep;
            w = r.x >= x0-0.30 & r.x <= x0+0.40;
            th(end+1) = r.limiterTheta; %#ok<AGROW>
            E1(end+1,:) = cellfun(@(q) mean(abs(r.(q)-r.exact.(q))),vars); %#ok<AGROW>
            pk(end+1) = max(r.T(w))-max(r.exact.T(w)); %#ok<AGROW>
        end
        r = runs{2}; w = r.x >= x0-0.30 & r.x <= x0+0.40;       % theta = 1.2 on Nx = 6250
        th(end+1) = 1.2; E1(end+1,:) = cellfun(@(q) mean(abs(r.(q)-r.exact.(q))),vars); pk(end+1) = max(r.T(w))-max(r.exact.T(w));
        [th,o] = sort(th); E1 = E1(o,:); pk = pk(o);
        fig = new_figure(S,6.5,2.8);
        t = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
        ax = nexttile(t); hold(ax,'on');
        for v = 1:4
            plot(ax,th,E1(:,v)/E1(1,v),['-' markers{v}],'Color',S.cat(v,:),'MarkerFaceColor',S.cat(v,:),'MarkerEdgeColor','w','MarkerSize',5,'LineWidth',1.2);
        end
        xline(ax,1.2,':','Color',S.ink2,'HandleVisibility','off');
        style_axes(ax,S); xlabel(ax,'$\theta$'); ylabel(ax,'$L_1$ error / minmod value');
        legend(ax,{'$\rho$','$u$','$p$','$T$'},'Location','northeast','Box','off');
        title(ax,'(a) Average error versus $\theta$','FontWeight','normal');
        ax = nexttile(t); hold(ax,'on');
        px = 1e3*pk; py = 1e4*E1(:,4);
        plot(ax,px,py,'-','Color',S.seq(3,:),'LineWidth',1.2);
        scatter(ax,px,py,28,S.seq(4,:),'filled','MarkerEdgeColor','w');
        for k = 1:numel(th); text(ax,px(k),py(k),sprintf('  %.1f',th(k)),'FontSize',S.fs-2,'Color',S.ink2,'VerticalAlignment','bottom'); end
        k12 = find(abs(th-1.2)<1e-9,1);
        plot(ax,px(k12),py(k12),'o','MarkerSize',9,'MarkerEdgeColor',S.cat(2,:),'LineWidth',1.4);
        style_axes(ax,S); xlabel(ax,'temperature overshoot $[\times10^{-3}]$'); ylabel(ax,'$L_1(T)$ $[\times10^{-4}]$');
        title(ax,'(b) Accuracy--overshoot tradeoff (labels: $\theta$)','FontWeight','normal');
        save_all(fig,figFolder,'G13_theta_tradeoff');
    end

    %% G14 contact and shock thickness ----------------------------------------
    if want('G14_feature_thickness')
        Pst = (v1^2-us^2)*(v2^2-us^2);
        Dc = runs{1}.tau*Pst./((b+2)*[ps/rsR ps/rsL]);
        dKin = sqrt(4*pi*Dc*tEnd);
        pc = polyfit(log(dx),log(M.dc),1); pshk = polyfit(log(dx),log(M.ds),1);
        fig = new_figure(S,4.6,3.1); ax = axes(fig); hold(ax,'on');
        plot(ax,dx,M.dc,'-o','Color',S.cat(2,:),'MarkerFaceColor',S.cat(2,:),'MarkerEdgeColor','w','LineWidth',1.3);
        plot(ax,dx,M.ds,'-s','Color',S.cat(8,:),'MarkerFaceColor',S.cat(8,:),'MarkerEdgeColor','w','LineWidth',1.3);
        patch(ax,[dx(end)*0.7 dx(1)*1.4 dx(1)*1.4 dx(end)*0.7],[dKin(1) dKin(1) dKin(2) dKin(2)],S.cat(2,:),'FaceAlpha',0.18,'EdgeColor','none');
        plot(ax,dx,dx,':','Color',S.ink2);
        set(ax,'XScale','log','YScale','log','XDir','reverse'); style_axes(ax,S);
        xlim(ax,[dx(end)*0.7 dx(1)*1.4]); set_dx_ticks(ax,dx);
        ylabel(ax,'thickness');
        legend(ax,{sprintf('contact (slope %.2f)',pc(1)),sprintf('shock (slope %.2f)',pshk(1)), ...
            'contact width set by $\tau$','one cell'},'Location','northeast','Box','off','FontSize',S.fs-1);
        save_all(fig,figFolder,'G14_feature_thickness');
    end
end

%% =========================================================================
function M = metrics_for(runs,x0)
    vars = {'rho','u','p','T'};
    gamma = 1.4; [~,us,rsL,rsR,Ssh] = sod_constants(gamma); tEnd = 0.15;
    for k = 1:numel(runs)
        r = runs{k}; w = r.x >= x0-0.30 & r.x <= x0+0.40;
        for v = 1:4
            e = r.(vars{v}) - r.exact.(vars{v});
            M.L1(k,v) = mean(abs(e)); M.L2(k,v) = sqrt(mean(e.^2));
        end
        M.Tpeak(k) = max(r.T(w)) - max(r.exact.T(w));
        xs = x0+Ssh*tEnd; M.phase(k) = mod(xs/r.dx,1);
        m = abs(r.x-xs) < 4*r.dx; eU = r.u - r.exact.u;
        M.shockU(k) = sqrt(sum(eU(m).^2)/numel(r.x));
        M.dc(k) = thickness(r.x,r.rho,[x0+0.07 x0+0.20],rsL,rsR);
        M.ds(k) = thickness(r.x,r.rho,[x0+0.20 x0+0.32],rsR,0.125);
    end
    M.us = us;
end

function C = method_comparison(projectFolder,gminRuns,grids,x0)
% collect errors, steps and computing time for every method
    S = report_style();
    folder = fullfile(projectFolder,'results','method_runs');
    tb = load(fullfile(folder,'timing_benchmark.mat'));
    cG = mean([tb.dbmTiming(strcmp({tb.dbmTiming.limiter},'gminmod')).costPerCellStep]);
    cF = tb.dbmTiming(strcmp({tb.dbmTiming.limiter},'firstorder')).costPerCellStep;
    defs = {'dbm','firstorder','DBM, first-order upwind',S.seq(2,:),':o',cF;
            'dbm','minmod','DBM, MUSCL minmod',S.seq(3,:),'-.s',cG;
            'dbm','gminmod','DBM, MUSCL $\theta=1.2$ (this solver)',S.seq(5,:),'-d',cG;
            'dbm','mc','DBM, MUSCL MC',S.seq(4,:),'--^',cG;
            'euler','firstorder','Euler, first-order HLLC',S.cat(2,:),':v',NaN;
            'euler','muscl','Euler, MUSCL-HLLC $\theta=1.2$',S.cat(8,:),'-p',NaN};
    C = struct('label',{},'color',{},'fill',{},'style',{},'lw',{},'Nx',{},'dx',{},'L1rho',{},'L2rho',{}, ...
        'L1T',{},'L2T',{},'steps',{},'time',{},'kind',{},'name',{});
    for m = 1:size(defs,1)
        kind = defs{m,1}; name = defs{m,2};
        E = struct('Nx',[],'dx',[],'L1rho',[],'L2rho',[],'L1T',[],'L2T',[],'steps',[],'time',[]);
        for Nx = grids
            if strcmp(kind,'dbm')
                if strcmp(name,'gminmod')
                    r = gminRuns{grids==Nx};
                else
                    f = fullfile(folder,sprintf('d1v5_%s_tau5e-06_CFL0p025_Nx%d.mat',name,Nx));
                    if ~isfile(f); continue; end
                    tmp = load(f,'runData'); r = tmp.runData;
                end
                steps = r.step; time = defs{m,6}*steps*Nx;   % clean cost per cell-step times work
            else
                f = fullfile(folder,sprintf('euler_%s_Nx%d.mat',name,Nx));
                if ~isfile(f); continue; end
                tmp = load(f,'out'); r = tmp.out; steps = r.steps; time = r.runTime;
            end
            E.Nx(end+1) = Nx; E.dx(end+1) = r.dx;
            E.L1rho(end+1) = mean(abs(r.rho-r.exact.rho)); E.L2rho(end+1) = sqrt(mean((r.rho-r.exact.rho).^2));
            E.L1T(end+1) = mean(abs(r.T-r.exact.T)); E.L2T(end+1) = sqrt(mean((r.T-r.exact.T).^2));
            E.steps(end+1) = steps; E.time(end+1) = time;
        end
        sty = defs{m,5}; filled = ~strcmp(name,'firstorder');
        C(m).label = defs{m,3}; C(m).color = defs{m,4}; C(m).style = sty;
        if filled; C(m).fill = defs{m,4}; else; C(m).fill = 'w'; end
        C(m).lw = 1.2 + 0.6*strcmp(name,'gminmod');
        C(m).kind = kind; C(m).name = name;
        fn = fieldnames(E); for q = 1:numel(fn); C(m).(fn{q}) = E.(fn{q}); end
    end
end

function write_method_tables(C,dataFolder)
    fid = fopen(fullfile(dataFolder,'tab_methods.tex'),'w');
    for m = 1:numel(C)
        g = C(m).Nx <= 25000;                                % same four grids for every method
        p1 = polyfit(log(C(m).dx(g)),log(C(m).L1rho(g)),1); p2 = polyfit(log(C(m).dx(g)),log(C(m).L2rho(g)),1);
        k = find(C(m).Nx==25000,1);
        fprintf(fid,'%s & %.2f & %.2f & %s & %s & %s & %s \\\\\n',strrep(C(m).label,' (this solver)',''), ...
            p1(1),p2(1),sci(C(m).L1rho(k),2),sci(C(m).L2rho(k),2),commas(C(m).steps(k)),fmt_time(C(m).time(k)));
    end
    fclose(fid);
    save(fullfile(dataFolder,'method_comparison.mat'),'C');
end

function s = fmt_time(t)
    if t >= 100; s = sprintf('%s',commas(round(t))); elseif t >= 1; s = sprintf('%.1f',t); else; s = sprintf('%.2f',t); end
end

function draw_reference_slopes(ax,dx,ref,orders,S,yl)
% reference lines anchored at the coarsest grid, fanning out toward finer grids
    for p = orders
        y = ref*(dx/dx(1)).^p;
        loglog(ax,dx,y,':','Color',S.ink2,'LineWidth',1.0,'HandleVisibility','off');
        if p == 2; lbl = 'slope 2 (2nd order)'; elseif p == 1; lbl = 'slope 1 (1st order)'; else; lbl = 'slope 1/2'; end
        k = find(y >= yl(1)*1.6,1,'last');                  % last visible point of the line
        if p == 1; va = 'top'; elseif p == 2; va = 'middle'; else; va = 'bottom'; end
        text(ax,dx(k)*0.94,y(k),lbl,'FontSize',S.fs-2,'Color',S.ink2,'HorizontalAlignment','left','VerticalAlignment',va);
    end
end

function set_dx_ticks(ax,dx)
    labels = arrayfun(@(d) sprintf('%.1f',1e3*d),dx,'UniformOutput',false);
    [dsorted,o] = sort(dx);
    set(ax,'XTick',dsorted,'XTickLabel',labels(o),'XMinorTick','off','XMinorGrid','off','XTickLabelRotation',0);
    xlabel(ax,'$\Delta x\ [\times10^{-3}]$');
end

function save_all(fig,folder,name)
% save the figure as an editable MATLAB .fig, a PNG, and a vector PDF
    exportgraphics(fig,fullfile(folder,[name '.pdf']),'ContentType','vector','BackgroundColor','white');
    exportgraphics(fig,fullfile(folder,[name '.png']),'Resolution',220,'BackgroundColor','white');
    set(fig,'Visible','on');                                   % so the .fig opens visibly in MATLAB
    savefig(fig,fullfile(folder,[name '.fig']));
    close(fig);
    fprintf('saved %s (.fig/.png/.pdf)\n',name);
end

function delta = thickness(x,q,xr,qLeft,qRight)
    m = x>=xr(1) & x<=xr(2);
    g = abs(gradient(q(m),x(m)));
    delta = abs(qLeft-qRight)/max(g);
end

function s = sci(v,d)
    e = floor(log10(abs(v))); mm = v/10^e;
    if abs(round(mm,d-1)) >= 10; mm = mm/10; e = e+1; end
    s = sprintf(['$%.' num2str(d-1) 'f\\times10^{%d}$'],mm,e);
end

function s = commas(n)
    s = regexprep(sprintf('%d',round(n)),'(\d)(?=(\d{3})+$)','$1,');
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

function ex = sod_exact(x,t,x0,g)
    [ps,us,rsL,rsR,S,head,tail] = sod_constants(g);
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
