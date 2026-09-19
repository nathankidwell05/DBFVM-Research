function style_axes(ax,S)
%STYLE_AXES recessive grid and axes, primary-ink labels.
    set(ax,'Box','off','TickDir','out','LineWidth',0.6,'XColor',S.ink2,'YColor',S.ink2, ...
        'GridColor',S.grid,'GridAlpha',1,'MinorGridColor',S.grid,'MinorGridAlpha',0.6, ...
        'FontSize',S.fs,'TickLabelInterpreter','latex','Layer','top');
    grid(ax,'on');
    ax.XLabel.Color = S.ink; ax.YLabel.Color = S.ink; ax.Title.Color = S.ink;
    ax.Title.FontSize = S.fs;
end
