function out = d1v5_order_estimate(varargin)
%D1V5_ORDER_ESTIMATE  Observed order of accuracy from the L2 convergence study.
%
%   d1v5_order_estimate
%       Uses the cached generalized-minmod runs (theta = 1.20, tau = 5e-6,
%       CFL = 0.025) in results/runs, prints the L1 and L2 errors, the
%       order between each pair of grids, the least-squares order over all
%       grids, and a short verdict on what order of accuracy the solver is
%       showing on this problem.
%
%   out = d1v5_order_estimate(...)  also returns the numbers in a struct.
%
%   Options (name-value pairs):
%       'Limiter'  'gminmod' (default), 'minmod', 'mc' or 'firstorder'.
%                  Anything other than gminmod is read from results/method_runs.
%       'Grids'    grids to use (default: every cached grid found).
%       'Variable' 'rho' (default), 'u', 'p', 'T', or 'all' for a summary of each.
%       'Plot'     true to add a log-log plot with the fitted slope (default false).
%
%   Examples:
%       d1v5_order_estimate
%       d1v5_order_estimate('Variable','all')
%       d1v5_order_estimate('Limiter','firstorder')
%       out = d1v5_order_estimate('Grids',[6250 12500 25000 50000],'Plot',true);
%
%   Note on interpretation: the Sod solution contains a shock and a contact
%   discontinuity.  No shock-capturing scheme, whatever its formal order,
%   converges at its formal rate on such a solution.  The useful question is
%   therefore not "is the number 2?" but "is it in the range a second-order
%   scheme produces on this problem?".  The verdict printed below compares the
%   measured orders against reference bands and against the first-order
%   version of this same solver.

    p = inputParser;
    p.addParameter('Limiter','gminmod',@(s) ischar(s) || isstring(s));
    p.addParameter('Grids',[],@isnumeric);
    p.addParameter('Variable','rho',@(s) ischar(s) || isstring(s));
    p.addParameter('Plot',false,@(x) islogical(x) || isnumeric(x));
    p.parse(varargin{:});
    limiter = char(p.Results.Limiter);
    variable = char(p.Results.Variable);

    %% Load the cached runs -------------------------------------------------
    projectFolder = project_paths();
    [dx,Nx,err] = load_errors(projectFolder,limiter,p.Results.Grids);
    if numel(Nx) < 2
        error('Need at least two cached grids for "%s". Run the convergence study first.',limiter);
    end
    vars = {'rho','u','p','T'};
    labels = {'density','velocity','pressure','temperature'};

    %% Report ----------------------------------------------------------------
    fprintf('\n====================================================================\n');
    fprintf(' Observed order of accuracy  -  KT-D1V5 solver, limiter: %s\n',limiter);
    fprintf(' Grids: %s   (Sod shock tube, t = 0.15)\n',strjoin(string(Nx),', '));
    fprintf('====================================================================\n');

    if strcmpi(variable,'all'); list = 1:4; else; list = find(strcmp(vars,lower(variable))); end
    if isempty(list); error('Variable must be rho, u, p, T or all.'); end

    out = struct('Nx',Nx,'dx',dx,'limiter',limiter);
    for v = list
        L1 = err.L1(:,v)'; L2 = err.L2(:,v)';
        o1 = pairwise(L1,dx); o2 = pairwise(L2,dx);
        [f1,r1] = lsq_order(dx,L1); [f2,r2] = lsq_order(dx,L2);
        fprintf('\n--- %s ---\n',upper(labels{v}));
        fprintf('  %8s  %11s %8s   %11s %8s\n','Nx','L1 error','order','L2 error','order');
        for k = 1:numel(Nx)
            s1 = '     -  '; s2 = '     -  ';
            if k > 1; s1 = sprintf('%8.2f',o1(k)); s2 = sprintf('%8.2f',o2(k)); end
            fprintf('  %8d  %11.4e %s   %11.4e %s\n',Nx(k),L1(k),s1,L2(k),s2);
        end
        fprintf('  least-squares fit over all grids:  L1 order = %.2f (R^2 = %.3f),  L2 order = %.2f (R^2 = %.3f)\n', ...
            f1,r1,f2,r2);
        if numel(Nx) >= 4
            fine = numel(Nx)-2:numel(Nx);
            fprintf('  finest three grids only:           L1 order = %.2f,             L2 order = %.2f\n', ...
                lsq_order(dx(fine),L1(fine)),lsq_order(dx(fine),L2(fine)));
        end
        out.(vars{v}) = struct('L1',L1,'L2',L2,'orderL1',o1,'orderL2',o2,'fitL1',f1,'fitL2',f2,'R2L1',r1,'R2L2',r2);
    end

    %% Verdict ---------------------------------------------------------------
    fitL1 = lsq_order(dx,err.L1(:,1)');                 % density drives the verdict
    fitL2 = lsq_order(dx,err.L2(:,1)');
    fitL1T = lsq_order(dx,err.L1(:,4)');
    print_verdict(fitL1,fitL2,fitL1T,limiter);
    out.verdictOrderL1 = fitL1; out.verdictOrderL2 = fitL2;

    %% Optional plot ----------------------------------------------------------
    if p.Results.Plot
        v = list(1);
        figure('Color','w','Name','Observed order of accuracy');
        loglog(dx,err.L1(:,v),'o-','LineWidth',1.4,'MarkerFaceColor','w'); hold on;
        loglog(dx,err.L2(:,v),'s-','LineWidth',1.4,'MarkerFaceColor','w');
        ref = min(err.L2(:,v))*0.4;
        loglog(dx,ref*(dx/dx(1)).^1,'k:','LineWidth',1);
        loglog(dx,ref*(dx/dx(1)).^2,'k--','LineWidth',1);
        set(gca,'XDir','reverse'); grid on; box on;
        xlabel('\Deltax'); ylabel(sprintf('error in %s',labels{v}));
        legend({sprintf('L1 (order %.2f)',lsq_order(dx,err.L1(:,v)')), ...
                sprintf('L2 (order %.2f)',lsq_order(dx,err.L2(:,v)')), ...
                'slope 1 (first order)','slope 2 (second order)'},'Location','southwest');
        title(sprintf('KT-D1V5 observed order, limiter: %s',limiter));
    end

    if nargout == 0; clear out; end                    % print only, unless a result is requested
end

%% ------------------------------------------------------------------------
function [dx,Nx,err] = load_errors(projectFolder,limiter,requestedGrids)
% collect cached runs and evaluate the error norms
    if strcmp(limiter,'gminmod')
        folder = fullfile(projectFolder,'results','runs');
        pattern = 'd1v5_gminmod_theta1p20_tau5e-06_CFL0p025_Nx*.mat';
        field = 'runData';
    else
        folder = fullfile(projectFolder,'results','method_runs');
        pattern = sprintf('d1v5_%s_tau5e-06_CFL0p025_Nx*.mat',limiter);
        field = 'runData';
    end
    files = dir(fullfile(folder,pattern));
    if isempty(files)
        error('No cached runs found in %s matching %s',folder,pattern);
    end
    Nx = []; dx = []; L1 = []; L2 = [];
    vars = {'rho','u','p','T'};
    for k = 1:numel(files)
        S = load(fullfile(files(k).folder,files(k).name),field);
        r = S.(field);
        if ~isempty(requestedGrids) && ~any(requestedGrids == r.Nx); continue; end
        Nx(end+1) = r.Nx; dx(end+1) = r.dx; %#ok<AGROW>
        for v = 1:4
            e = r.(vars{v}) - r.exact.(vars{v});
            L1(numel(Nx),v) = mean(abs(e)); %#ok<AGROW>
            L2(numel(Nx),v) = sqrt(mean(e.^2)); %#ok<AGROW>
        end
    end
    [Nx,order] = sort(Nx); dx = dx(order);
    err.L1 = L1(order,:); err.L2 = L2(order,:);
end

function o = pairwise(E,dx)
% order between successive grids
    o = NaN(size(E));
    for k = 2:numel(E)
        o(k) = log(E(k-1)/E(k))/log(dx(k-1)/dx(k));
    end
end

function [slope,R2] = lsq_order(dx,E)
% least-squares slope of log(error) against log(dx), with R^2
    x = log(dx(:)); y = log(E(:));
    c = polyfit(x,y,1); slope = c(1);
    yhat = polyval(c,x);
    R2 = 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
end

function print_verdict(fitL1,fitL2,fitL1T,limiter)
% plain-language interpretation of the measured orders
    fprintf('\n--- VERDICT -------------------------------------------------------\n');
    fprintf('  Measured on this problem:  L1 order = %.2f,  L2 order = %.2f  (density)\n',fitL1,fitL2);
    fprintf('\n  Reference values for the Sod problem (shock + contact present):\n');
    fprintf('      first-order scheme :  L1 about 0.5-0.7,  L2 about 0.25-0.45\n');
    fprintf('      second-order scheme:  L1 about 0.8-1.0,  L2 about 0.33-0.55\n');
    fprintf('      (this solver with first-order upwind measures L1 = 0.64, L2 = 0.43)\n');
    fprintf('      A formal order of 1 or 2 is only attainable where the solution is smooth.\n');

    if fitL1 >= 0.78
        grade = 'SECOND-ORDER behaviour';
        note  = 'matches a conventional second-order (MUSCL) scheme on this problem';
    elseif fitL1 >= 0.72
        grade = 'between first and second order, closer to second';
        note  = 'the reconstruction is helping, but less than a clean MUSCL scheme would';
    elseif fitL1 >= 0.60
        grade = 'FIRST-ORDER behaviour';
        note  = 'close to the 0.64 measured for first-order upwind on this problem';
    else
        grade = 'below first-order behaviour';
        note  = 'unexpected: check the run for fallbacks, a too-large time step or a bad grid';
    end
    fprintf('\n  => %s (%s).\n',grade,note);
    fprintf('     Temperature gives L1 order = %.2f for comparison.\n',fitL1T);
    if ~strcmp(limiter,'firstorder') && fitL1 < 0.72                 % a limited scheme should do better
        fprintf('     NOTE: this is lower than expected for limiter "%s".\n',limiter);
    end
    fprintf('\n  Caution: velocity orders are noisy on this problem because almost all of\n');
    fprintf('  the velocity error sits in the few cells at the shock, and its size depends\n');
    fprintf('  on where the shock happens to fall inside a cell.  Judge the order from the\n');
    fprintf('  density and temperature, or from L1 rather than L2.\n');
    fprintf('====================================================================\n\n');
end
