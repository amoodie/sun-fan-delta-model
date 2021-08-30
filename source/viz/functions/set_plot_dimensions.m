function set_plot_dimensions(width,height)
% plotting parameters
set(gcf,'color','w')
set(gca,'outerposition',[0 0 1 1],'units','normalized')
set(gcf,'Units','centimeters','position',[0 0 width height])
%set(gcf,'Units','centimeters','position',[0 0 desired_width pos(4)*desired_width/pos(3)]) % this way scales for proper height, but width is set
set(gcf,'PaperPositionMode','auto')
set(gcf, 'PaperUnits', 'centimeters');
siz=get(gcf,'position');
siz=siz(3:4); % just get width and height
set(gcf, 'PaperSize', siz)
end