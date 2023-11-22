clearvars
clc

load('../processed/20231122_active.mat');

%Get average position of background beads

avgBGpos = zeros(trackerBG.NumTracks, 2);
for iBead = 1:trackerBG.NumTracks

    ct = getTrack(trackerBG, iBead);

    avgBGpos(iBead, :) = mean(ct.Centroid, 1);

end

pos = cat(1, particlePos.Centroid);

plot(pos(:, 1), pos(:, 2), '-x')

%% Try clustering by density
idx = dbscan(pos, 4, 2);

gscatter(pos(:, 1), pos(:, 2), idx);
legend('off')
hold on
plot(avgBGpos(:, 1), avgBGpos(:, 2), 'o')
hold off
axis image