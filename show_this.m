function show_this(pack_colors, marker ,alpha, N, data, file_name, xtext, ytext)
    rows = size(data,1);
    for ii=1:rows
        h = plot(N,data(ii,:,1),strcat(pack_colors(ii),marker));
        set(h, 'MarkerFaceColor', get(h, 'Color'));
        hold on
        Legenda{ii} = strcat('\alpha = ',sprintf('%.2f',alpha(ii)));
    end
    legend(Legenda,'Location','best')
    xlabel(xtext)
    ylabel(ytext)
    print(gcf,strcat(file_name,'.png'),'-dpng')
end