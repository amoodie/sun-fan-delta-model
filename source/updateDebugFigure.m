function [fig] = updateDebugFigure(fig, grid)

    % grab the subplots
    sp1 = subplot(2, 3, 1);
    sp2 = subplot(2, 3, 2);
    sp3 = subplot(2, 3, 3);
%     sp4 = subplot(2, 3, 4);
%     sp5 = subplot(2, 3, 5);
%     sp6 = subplot(2, 3, 6);
    sp4 = subplot(2, 3, 4:6);

    subplot(sp1); hold on; cla;
        imagesc(grid.z);
        contourlevel = grid.oceanLevel;
        if isfinite(contourlevel)
            contour3(grid.z,[contourlevel,contourlevel],'k-')
        end
        colormap(sp1, 'summer');
        title('bed elev');
        gridToChannelArrows(grid)
        sp1.YDir = 'reverse';
        xlim([0.5, grid.size(2)+0.5])
        ylim([0.5, grid.size(1)+0.5])
        axis equal

    subplot(sp2); cla;
        imagesc(grid.deltaz);
        colormap(sp2, 'turbo');
        caxis([-1e-4, 1e-4]);
        title('bed change');
        sp1.YDir = 'reverse';
        xlim([0.5, grid.size(2)+0.5])
        ylim([0.5, grid.size(1)+0.5])
        axis equal

    subplot(sp3); cla;
        stage = grid.H;
        upperlim = prctile(stage(:), 90);
        stage(~isfinite(stage)) = NaN; grid.z((~isfinite(stage)));
        imagesc(stage);
%         caxis([min(min(stage)), upperlim])
        colormap(sp3, 'winter');
        title('stage');
        sp1.YDir = 'reverse';
        xlim([0.5, grid.size(2)+0.5])
        ylim([0.5, grid.size(1)+0.5])
        axis equal

%     subplot(sp4); cla;
%         imagesc(grid.S.alongFlow);
%         colormap(sp4);
%         title('slope');
% 
%     subplot(sp5); cla;
%         imagesc(grid.Qw);
%         colormap(sp5, 'cool');
%         title('Qw');
% 
%     subplot(sp6); cla;
%         imagesc(grid.channelFlag);
%         colormap(sp6, 'gray');
%         title('channelFlag');
    
    subplot(sp4); cla; hold on;
        plot(grid.z0(:, 100), 'k-')
        plot(grid.z(:, 100), 'Color', [0.6 0.6 0.6])
        plot(grid.z(:, 102), 'Color', [0.6 0.6 0.6])
        plot(grid.z(:, 98), 'Color', [0.6 0.6 0.6])

    drawnow


end