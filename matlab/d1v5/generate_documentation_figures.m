%% Generate Extra Documentation Figures From Completed D1V5 Runs %%
% This script does not rerun or modify the solver.
% It reads the MAT files produced by the reference and limiter-history scripts.
% It then makes figures that explain where the error occurs and why the final
% generalized-minmod limiter was selected as a compromise.

clear;                                              % remove old workspace variables
clc;                                                % clear the Command Window
close all;                                          % close figures from an earlier plotting run

%% Locate Saved MATLAB Results
thisFolder = fileparts(mfilename('fullpath'));       % matlab/d1v5 folder containing this script
projectFolder = fileparts(fileparts(thisFolder));   % top-level repository folder
resultsFolder = fullfile(projectFolder,'results');  % versioned figure output folder

referenceFile = fullfile(thisFolder,'results', ...
    'd1v5_sodshock_gminmod_rk3_latest.mat');         % output from the main reference solver
limiterFile = fullfile(thisFolder, ...
    'd1v5_limiter_history_results.mat');             % output from the limiter-history comparison

if ~isfile(referenceFile)                            % stop with a useful instruction if data are absent
    error(['Run d1v5_sodshock_gminmod_rk3 first. Missing file: ' referenceFile]);
end

if ~isfile(limiterFile)                              % stop with a useful instruction if comparison is absent
    error(['Run d1v5_sodshock_limiter_history_comparison first. Missing file: ' limiterFile]);
end

if ~isfolder(resultsFolder)                          % create the central output folder in a fresh clone
    mkdir(resultsFolder);
end

reference = load(referenceFile);                     % load main-solver fields into one named structure
comparison = load(limiterFile);                      % load limiter results into a second named structure

Nx = reference.Nx;                                  % grid resolution used by the saved reference run
x0 = comparison.x0;                                 % initial diaphragm location
waveMask = reference.x >= x0-0.30 & reference.x <= x0+0.40; % common wave comparison window

%% Figure 1: Spatial Absolute Error
errorFigure = figure('Name','D1V5 Spatial Absolute Error','Color','w', ...
    'Position',[70 50 1500 980]);                    % large four-panel output
errorLayout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

fieldNames = {'rho','u','p','T'};                    % fields to compare with the exact solution
fieldLabels = {'Density error, |rho-rho_{exact}|', ...
    'Velocity error, |u-u_{exact}|', ...
    'Pressure error, |p-p_{exact}|', ...
    'Temperature error, |T-T_{exact}|'};             % readable titles and y-axis labels

for panel = 1:4                                      % create one spatial error panel for each field
    ax = nexttile(errorLayout);                      % select the next panel
    numericalValue = reference.(fieldNames{panel});  % final DBM field
    exactValue = reference.exact.(fieldNames{panel});% matching exact Euler field
    absoluteError = abs(numericalValue-exactValue);  % pointwise magnitude of the difference

    plot(ax,reference.x(waveMask),absoluteError(waveMask), ...
        'Color',[0.00 0.55 0.15],'LineWidth',2.1);   % show only the wave-region error
    xlabel(ax,'x','FontSize',14,'FontWeight','bold');% label physical position
    ylabel(ax,'Absolute error','FontSize',14,'FontWeight','bold'); % label error magnitude
    title(ax,fieldLabels{panel},'FontSize',16,'FontWeight','bold'); % describe the field
    grid(ax,'on');                                    % help locate error peaks
    box(ax,'on');                                     % frame the panel
    ax.FontSize = 13;                                 % make tick labels readable
end

title(errorLayout,sprintf('Where the Generalized-Minmod Error Occurs: Nx=%d',Nx), ...
    'FontSize',21,'FontWeight','bold');               % shared explanation above the figure
exportgraphics(errorFigure,fullfile(resultsFolder, ...
    sprintf('d1v5_gminmod_absolute_error_Nx%d.png',Nx)),'Resolution',180); % save figure

%% Figure 2: Global and Wave-Region L1 Bars
errorTable = comparison.errorTable;                  % table produced by matched limiter comparison
limiterLabels = cellstr(errorTable.Limiter);         % labels for x-axis categories
limiterLabels = strrep(limiterLabels,'\theta','theta'); % avoid TeX commands in tick labels

globalErrors = [errorTable.GlobalL1_rho,errorTable.GlobalL1_u, ...
    errorTable.GlobalL1_p,errorTable.GlobalL1_T];     % whole-domain L1 errors by field
waveErrors = [errorTable.WaveL1_rho,errorTable.WaveL1_u, ...
    errorTable.WaveL1_p,errorTable.WaveL1_T];         % wave-only L1 errors by field

barFigure = figure('Name','D1V5 Limiter L1 Error Bars','Color','w', ...
    'Position',[60 40 1600 1050]);                   % large comparison output
barLayout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
shortFieldLabels = {'Density','Velocity','Pressure','Temperature'}; % panel titles

for panel = 1:4                                      % build one bar panel per macroscopic field
    ax = nexttile(barLayout);                        % select the next panel
    bar(ax,[globalErrors(:,panel),waveErrors(:,panel)],'grouped'); % compare two error regions
    ax.XTick = 1:numel(limiterLabels);                % one x position per limiter
    ax.XTickLabel = limiterLabels;                    % identify the limiter under every group
    ax.XTickLabelRotation = 18;                       % prevent long labels from overlapping
    ylabel(ax,'L1 error','FontSize',14,'FontWeight','bold'); % label error metric
    title(ax,shortFieldLabels{panel},'FontSize',17,'FontWeight','bold'); % identify field
    grid(ax,'on');                                    % make heights easier to compare
    box(ax,'on');                                     % frame the panel
    ax.FontSize = 11;                                 % readable category labels
    if panel == 1                                     % use one legend instead of repeating it four times
        legend(ax,{'Whole domain','Wave region'},'Location','best','FontSize',12);
    end
end

title(barLayout,sprintf('Matched Limiter Errors: Nx=%d',Nx), ...
    'FontSize',22,'FontWeight','bold');               % shared figure title
exportgraphics(barFigure,fullfile(resultsFolder, ...
    sprintf('d1v5_limiter_error_bars_Nx%d.png',Nx)),'Resolution',180); % save bars

%% Figure 3: Accuracy Versus Temperature Overshoot
tradeoffFigure = figure('Name','D1V5 Limiter Tradeoff','Color','w', ...
    'Position',[80 80 1500 700]);                    % two-panel decision figure
tradeoffLayout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tradeoffLayout);                      % left panel: direct temperature tradeoff
peakExcess = max(errorTable.T_peak_excess,0);        % plot positive overshoot magnitude only
scatter(ax1,errorTable.WaveL1_T,peakExcess,115,'filled'); % every point is one limiter
hold(ax1,'on');                                      % allow text labels on the same axes
for caseNumber = 1:height(errorTable)                % label each point directly
    labelX = errorTable.WaveL1_T(caseNumber)*1.006;   % normally place the label to the point's right
    labelAlignment = 'left';                         % normal horizontal text alignment
    if errorTable.WaveL1_T(caseNumber) > 0.95*max(errorTable.WaveL1_T)
        labelX = errorTable.WaveL1_T(caseNumber)*0.994; % move right-edge labels to the point's left
        labelAlignment = 'right';                    % keep the full label inside the axes
    end
    text(ax1,labelX,peakExcess(caseNumber)*1.03,limiterLabels{caseNumber}, ...
        'FontSize',11,'HorizontalAlignment',labelAlignment); % identify the limiter
end
xlabel(ax1,'Wave-region L1 temperature error','FontSize',14,'FontWeight','bold');
ylabel(ax1,'Temperature peak excess','FontSize',14,'FontWeight','bold');
title(ax1,'Smaller Is Better in Both Directions','FontSize',17,'FontWeight','bold');
set(ax1,'YScale','log','FontSize',12);                % log scale separates the overshoot magnitudes
grid(ax1,'on');                                      % aid quantitative comparison
box(ax1,'on');                                       % frame the plot

ax2 = nexttile(tradeoffLayout);                      % right panel: normalized all-field errors
normalizedGlobal = 100*globalErrors./globalErrors(1,:); % minmod equals 100 percent in every field
bar(ax2,normalizedGlobal,'grouped');                 % compare relative global error by limiter
ax2.XTick = 1:numel(limiterLabels);                  % one group per limiter
ax2.XTickLabel = limiterLabels;                      % identify each limiter
ax2.XTickLabelRotation = 18;                         % fit the long generalized-minmod label
yline(ax2,100,'k--','Minmod baseline','LineWidth',1.4); % values below line improve on minmod
ylabel(ax2,'Global L1 error (% of minmod)','FontSize',14,'FontWeight','bold');
title(ax2,'Average Error Relative to Minmod','FontSize',17,'FontWeight','bold');
legend(ax2,{'Density','Velocity','Pressure','Temperature'}, ...
    'Location','best','FontSize',11);                % identify grouped field bars
grid(ax2,'on');                                      % help compare percentages
box(ax2,'on');                                       % frame the plot
ax2.FontSize = 11;                                   % readable axes

title(tradeoffLayout,sprintf('Limiter Accuracy and Overshoot Tradeoff: Nx=%d',Nx), ...
    'FontSize',21,'FontWeight','bold');               % shared figure title
exportgraphics(tradeoffFigure,fullfile(resultsFolder, ...
    sprintf('d1v5_limiter_tradeoff_Nx%d.png',Nx)),'Resolution',180); % save tradeoff

fprintf('Documentation figures saved in:\n%s\n',resultsFolder); % confirm output location
