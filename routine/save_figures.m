function [ret] = save_figures(fig, fname_cell, type)
% Save figure
% INPUT
%     fig              Figure handler
%     fname_cell       File name
%     type             File type
% 
% Jun 2019 Shoichi Koyama, Gilles Chardon, and Laurent Daudet

num = length(fig);

if strcmp(type,'fig')==1
    for ii=1:num
        saveas(fig(ii),fname_cell{ii});
    end
else
    for ii=1:num
        set(fig(ii),'Units','Inches');
        pos = get(fig(ii),'Position');
        set(fig(ii),'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
        saveas(fig(ii),fname_cell{ii},type);
    end
end

ret = 0;

end