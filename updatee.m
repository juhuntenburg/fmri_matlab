function [points,h_indexedit]=updatee(points,axes_height,pos,to_correct)
    [xi ~]=ginput(5);
    for i=1:5
        xi(i)=round(xi(i));
    end
    points=[points xi'];
    for i=1:length(to_correct)
        line((to_correct(i)*ones(1,2)), axes_height,'Color', 'r');
    end
    for i=1:length(points)
        line((points(i)*ones(1,2)), axes_height,'Color', 'r');
    end
    h_indexedit = uicontrol(gcf, 'Units', 'characters', 'Position', pos,...
        'String', int2str(unique([to_correct points])), 'Style', 'edit', ...
        'HorizontalAlignment', 'left', ...
        'Tag', 'indexedit');
end