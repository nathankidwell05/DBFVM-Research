function save_figure(fig,fileName)
%SAVE_FIGURE writes a vector PDF for LaTeX and closes the figure.
    exportgraphics(fig,fileName,'ContentType','vector','BackgroundColor','white');
    fprintf('saved %s\n',fileName);
    close(fig);
end
