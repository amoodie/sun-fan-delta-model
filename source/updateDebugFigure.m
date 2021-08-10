function [fig] = updateDebugFigure(fig, grid)

    sp1 = subplot(2, 3, 1);
    sp2 = subplot(2, 3, 2);
    sp3 = subplot(2, 3, 3);
    sp4 = subplot(2, 3, 4);
    sp5 = subplot(2, 3, 5);
    sp6 = subplot(2, 3, 6);

    subplot(sp1)
        imagesc(grid.z);
        colormap(sp1, 'summer');
        xlim([40,60]);
        ylim([0.5,20]);
        title('bed elev');
    subplot(sp2)
        imagesc(grid.deltaz);
        colormap(sp2, 'turbo');
        caxis([-1e-4, 1e-4]);
        xlim([40,60]);
        ylim([0.5,20]);
        title('bed change');
    subplot(sp3)
        imagesc(grid.z + grid.H);
        colormap(sp3, 'winter');
        xlim([40,60]);
        ylim([0.5,20]);
        title('stage');
    subplot(sp4)
        imagesc(grid.S.alongFlow);
        colormap(sp4);
        xlim([40,60]);
        ylim([0.5,20]);
        title('slope');
    subplot(sp5)
        imagesc(grid.Qw);
        colormap(sp5, 'cool');
        xlim([40,60]);
        ylim([0.5,20]);
        title('Qw');
    subplot(sp6)
        imagesc(grid.channelFlag);
        colormap(sp6, 'gray');
        xlim([40,60]);
        ylim([0.5,20]);
        title('channelFlag');
    drawnow


end