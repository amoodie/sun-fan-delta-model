1. Run model using sunFanModel_v4.m
2. Examples to visualize output:

discharge: imagesc(grid.Qw)
elevation: imagesc(grid.z)
channel cells: imshow(grid.channelFlag)
ocean cells: imshow(grid.oceanFlag)
