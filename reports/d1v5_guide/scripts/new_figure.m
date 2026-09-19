function fig = new_figure(S,widthIn,heightIn)
%NEW_FIGURE creates an invisible, print-sized white figure with LaTeX text.
    fig = figure('Visible','off','Color','w','Units','inches','Position',[1 1 widthIn heightIn]);
    set(fig,'DefaultTextInterpreter','latex','DefaultLegendInterpreter','latex', ...
        'DefaultAxesTickLabelInterpreter','latex','DefaultAxesFontSize',S.fs, ...
        'DefaultTextFontSize',S.fs,'DefaultLegendFontSize',S.fs-1);
end
