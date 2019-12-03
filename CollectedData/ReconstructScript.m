close all
clear
clc

% Parameters.
[tagPosition, rxPosition, freq] = tag_antenna_positions3D_func();
TagNum = size(tagPosition, 1);
RecvNum = size(rxPosition, 1);
FreqNum = length(freq);
c = physconst('LightSpeed');

%% Data Collection.

% Threshold for standard deviation of phase.
phStdTh = 15; % In degrees. A very high value means no rejection.
% fprintf('Using phase standard deviation threshold %d\n', phStdTh);
rssiDiffTh = .5;
RJT = 1;

% Specify data directory and name.
dirName = 'D:\commercial_RFID_imaging_apps\Vicon room results\20190524 (12by12 room, tag regular, ant ceiling, 4_80_6, systematic tests)\';
fileName = 'Pragya_standing_x2y10';

% Load the raw data.
load([dirName, 'data_wo_',fileName,'_v1_to_v2_rm0.mat']); % Looking at 2 minute cases
% data_wo = [chindexlist, tagindexlist, antennalist, rssiimpinjlist, rssiimpinjlist_d, phasedeglist];

% Generate the data matrices.
PhaseDataCollect1 = cell(TagNum, RecvNum, FreqNum);
PhaseCollectStd1 = zeros(TagNum, RecvNum, FreqNum);
PhaseCollectMean1 = zeros(TagNum, RecvNum, FreqNum);
RSSIDataCollect1 = cell(TagNum, RecvNum, FreqNum);
RSSICollectStd1 = zeros(TagNum, RecvNum, FreqNum);
RSSICollectMean1 = zeros(TagNum, RecvNum, FreqNum);
CompNumDataCollect1 = cell(TagNum, RecvNum, FreqNum);
CompNumCollectMean1 = zeros(TagNum, RecvNum, FreqNum);
for j = 1:RecvNum
    phase_t = phasedeglist(antennalist == j);
    rssi_t = rssiimpinjlist_d(antennalist == j);
    cha_t = chindexlist(antennalist == j);
    tag_t = tagindexlist(antennalist == j);
    for k = 1:FreqNum
        phase_tt = phase_t(cha_t == k);
        rssi_tt = rssi_t(cha_t == k);
        tag_tt = tag_t(cha_t == k);
        for i = 1:TagNum
            phase_ttt = phase_tt(tag_tt == i);
            phase_ttt = rad2deg(unwrap(deg2rad(phase_ttt)));
            rssi_ttt = rssi_tt(tag_tt == i);
            if ~isempty(phase_ttt)
                PhaseDataCollect1{i, j, k} = phase_ttt;
                PhaseCollectStd1(i, j, k) = std(phase_ttt);
                RSSIDataCollect1{i, j, k} = rssi_ttt;
                RSSICollectStd1(i, j, k) = std(rssi_ttt);
                CompNumDataCollect1{i, j, k} = rssi_ttt.*cos(deg2rad(phase_ttt)) + sqrt(-1).*rssi_ttt.*sin(deg2rad(phase_ttt));
                if PhaseCollectStd1(i, j, k) > phStdTh
                    % First completely reject the channel using phase
                    % standard deviation threshold.
                    PhaseCollectMean1(i, j, k) = NaN;
                    RSSICollectMean1(i, j, k) = NaN;
                    CompNumCollectMean1(i, j, k) = NaN;
                else
                    PhaseCollectMean1(i, j, k) = mean(phase_ttt);
                    RSSICollectMean1(i, j, k) = mean(rssi_ttt);
                    CompNumCollectMean1(i, j, k) = mean(CompNumDataCollect1{i, j, k});
                end
            else
                % Simply there's no data received for this channel.
                PhaseDataCollect1{i, j, k} = NaN;
                PhaseCollectStd1(i, j, k) = NaN;
                PhaseCollectMean1(i, j, k) = NaN;
                RSSIDataCollect1{i, j, k} = NaN;
                RSSICollectStd1(i, j, k) = NaN;
                RSSICollectMean1(i, j, k) = NaN;
                CompNumDataCollect1{i, j, k} = NaN;
                CompNumCollectMean1(i, j, k) = NaN;
            end                
        end
    end
end
clear chindexlist tagindexlist antennalist rssiimpinjlist rssiimpinjlist_d phasedeglist

% Load the raw data.
load([dirName, '\data_w_',fileName,'_v1-2_rm0.mat']); 
% data_w = [chindexlist, tagindexlist, antennalist, rssiimpinjlist, rssiimpinjlist_d, phasedeglist];

% Generate the data matrices.
PhaseDataCollect2 = cell(TagNum, RecvNum, FreqNum);
PhaseCollectStd2 = zeros(TagNum, RecvNum, FreqNum);
PhaseCollectMean2 = zeros(TagNum, RecvNum, FreqNum);
RSSIDataCollect2 = cell(TagNum, RecvNum, FreqNum);
RSSICollectStd2 = zeros(TagNum, RecvNum, FreqNum);
RSSICollectMean2 = zeros(TagNum, RecvNum, FreqNum);
CompNumDataCollect2 = cell(TagNum, RecvNum, FreqNum);
CompNumCollectMean2 = zeros(TagNum, RecvNum, FreqNum);
for j = 1:RecvNum
    phase_t = phasedeglist(antennalist == j);
    rssi_t = rssiimpinjlist_d(antennalist == j);
    cha_t = chindexlist(antennalist == j);
    tag_t = tagindexlist(antennalist == j);
    for k = 1:FreqNum
        phase_tt = phase_t(cha_t == k);
        rssi_tt = rssi_t(cha_t == k);
        tag_tt = tag_t(cha_t == k);
        for i = 1:TagNum
            phase_ttt = phase_tt(tag_tt == i);
            phase_ttt = rad2deg(unwrap(deg2rad(phase_ttt)));
            rssi_ttt = rssi_tt(tag_tt == i);
            if ~isempty(phase_ttt)
                PhaseDataCollect2{i, j, k} = phase_ttt;
                PhaseCollectStd2(i, j, k) = std(phase_ttt);
                RSSIDataCollect2{i, j, k} = rssi_ttt;
                RSSICollectStd2(i, j, k) = std(rssi_ttt);
                CompNumDataCollect2{i, j, k} = rssi_ttt.*cos(deg2rad(phase_ttt)) + sqrt(-1).*rssi_ttt.*sin(deg2rad(phase_ttt));
                if PhaseCollectStd2(i, j, k) > phStdTh
                    % First completely reject the channel using phase
                    % standard deviation threshold.
                    PhaseCollectMean2(i, j, k) = NaN;
                    RSSICollectMean2(i, j, k) = NaN;
                    CompNumCollectMean2(i, j, k) = NaN;
                else
                    PhaseCollectMean2(i, j, k) = mean(phase_ttt);
                    RSSICollectMean2(i, j, k) = mean(rssi_ttt);
                    CompNumCollectMean2(i, j, k) = mean(CompNumDataCollect2{i, j, k});
                end
            else
                % Simply there's no data received for this channel.
                PhaseDataCollect2{i, j, k} = NaN;
                PhaseCollectStd2(i, j, k) = NaN;
                PhaseCollectMean2(i, j, k) = NaN;
                RSSIDataCollect2{i, j, k} = NaN;
                RSSICollectStd2(i, j, k) = NaN;
                RSSICollectMean2(i, j, k) = NaN;
                CompNumDataCollect2{i, j, k} = NaN;
                CompNumCollectMean2(i, j, k) = NaN;
            end                
        end
    end
end
clear chindexlist tagindexlist antennalist rssiimpinjlist rssiimpinjlist_d phasedeglist

% Reject the channels with blocked RSSI's.
if RJT == 1
    for n = 1:TagNum
        for k = 1:RecvNum
            for m = 1:FreqNum
                if RSSICollectMean2(n, k, m)/RSSICollectMean1(n, k, m) < rssiDiffTh
                    RSSICollectMean1(n, k, m) = NaN;
                    PhaseCollectMean1(n, k, m) = NaN;
                    CompNumCollectMean1(n, k, m) = NaN;
                    RSSICollectMean2(n, k, m) = NaN;
                    PhaseCollectMean2(n, k, m) = NaN;
                    CompNumCollectMean2(n, k, m) = NaN;
                end
            end
        end
    end
end


G_wo_s = CompNumCollectMean1;
G_w_s = CompNumCollectMean2;

% Find the zero values in the calibration matrix, and replace them with the
% maximum value in the matrix. (We can replace those zeros with any number
% we like, actually.)
idx_wo = find(RSSICollectMean1 ~= RSSICollectMean1);
idx_w = find(RSSICollectMean2 ~= RSSICollectMean2);
idx = [idx_wo; idx_w];
idx = unique(idx);
% ConstParam = max(max(max(abs(G_wo_s))));
ConstParam = 0;
% G_wo_s(idx_wo) = ConstParam;
% G_w_s(idx_wo) = ConstParam;
% G_w_s(idx_w) = G_wo_s(idx_w);
G_wo_s(idx) = ConstParam;
G_w_s(idx) = ConstParam;

% How many channels are lost.
NumLost = length(idx);
PercLost = NumLost/length(G_wo_s(:));
fprintf('Percentage Tx-Rx-Freq combination lost: %3.2f%% \n', 100*PercLost);


%% Calibration.
% Calculate the distance from tags to receivers.
r_x = bsxfun(@minus, tagPosition(:, 1), rxPosition(:, 1)');
r_y = bsxfun(@minus, tagPosition(:, 2), rxPosition(:, 2)');
r_z = bsxfun(@minus, tagPosition(:, 3), rxPosition(:, 3)');
r = sqrt(r_x.^2 + r_y.^2 + r_z.^2);

% Perform the calibration.
G_calib = zeros(TagNum, RecvNum, FreqNum);
% Differential receiving: First calculate the ratio of the data
% received by one receiver antenna, to the data received by all the
% other receiver antennas. Perform this for the data from both the
% case when no object is present, and the case when object is
% present. Then calculate the ratio of the above-calculated ratios
% of the case when object is present to the case when no object is
% present. The results are stored in the matrix named G_calib.
for m = 1:FreqNum
    for n = 1:TagNum
        for k = 1:RecvNum
            for k_pair = 1:RecvNum
                if k_pair ~= k
                    if G_w_s(n, k_pair, m) ~= 0 && G_wo_s(n, k_pair, m) ~= 0 ...
                            && G_w_s(n, k, m) ~= 0 && G_wo_s(n, k, m) ~= 0
                        % Both denominators and both numerators
                        % should not be zero.
                        G_wo_s_ratio = G_wo_s(n, k, m)/G_wo_s(n, k_pair, m);
                        G_w_s_ratio = G_w_s(n, k, m)/G_w_s(n, k_pair, m);
                        G_calib(n, k, m) = G_calib(n, k, m) + (G_w_s_ratio/G_wo_s_ratio - 1) ...
                            *exp(-1j*(2*pi*freq(m)/c)*r(n, k));
                    else
                        G_calib(n, k, m) = G_calib(n, k, m);
                    end
                end
            end
        end
    end
end
b = G_calib(:);

%% Reconstruction.
% The voxel coordinates.
x_v=0:0.12:3.6;
y_v=0:0.12:3.6;
z_v=0:0.3:2.4;

NxVoxel = length(x_v);
NyVoxel = length(y_v);
NzVoxel = length(z_v);

% Calculate all the distances needed.
VoxelCoord = combvec(x_v, y_v, z_v);
p_xVoxel = VoxelCoord(1, :);
p_yVoxel = VoxelCoord(2, :);
p_zVoxel = VoxelCoord(3, :);

p_tagx = zeros(TagNum, NxVoxel*NyVoxel*NzVoxel);
p_tagy = zeros(TagNum, NxVoxel*NyVoxel*NzVoxel);
p_tagz = zeros(TagNum, NxVoxel*NyVoxel*NzVoxel);
p_tagx(:, 1) = tagPosition(:, 1);
p_tagy(:, 1) = tagPosition(:, 2);
p_tagz(:, 1) = tagPosition(:, 3);
p_tagx = repmat(p_tagx(:, 1), [1, NxVoxel*NyVoxel*NzVoxel]);
p_tagy = repmat(p_tagy(:, 1), [1, NxVoxel*NyVoxel*NzVoxel]);
p_tagz = repmat(p_tagz(:, 1), [1, NxVoxel*NyVoxel*NzVoxel]);

p_recvx = zeros(RecvNum, NxVoxel*NyVoxel*NzVoxel);
p_recvy = zeros(RecvNum, NxVoxel*NyVoxel*NzVoxel);
p_recvz = zeros(RecvNum, NxVoxel*NyVoxel*NzVoxel);
p_recvx(:, 1) = rxPosition(:, 1);
p_recvy(:, 1) = rxPosition(:, 2);
p_recvz(:, 1) = rxPosition(:, 3);
p_recvx = repmat(p_recvx(:, 1), [1, NxVoxel*NyVoxel*NzVoxel]);
p_recvy = repmat(p_recvy(:, 1), [1, NxVoxel*NyVoxel*NzVoxel]);
p_recvz = repmat(p_recvz(:, 1), [1, NxVoxel*NyVoxel*NzVoxel]);

DistTagVoxel = sqrt((p_tagx - repmat(p_xVoxel, [TagNum, 1])).^2 ...
    +(p_tagy - repmat(p_yVoxel, [TagNum, 1])).^2+(p_tagz - repmat(p_zVoxel, [TagNum, 1])).^2);

DistRecvVoxel = sqrt((p_recvx - repmat(p_xVoxel, [RecvNum, 1])).^2 ...
    +(p_recvy - repmat(p_yVoxel, [RecvNum, 1])).^2+(p_recvz - repmat(p_zVoxel, [RecvNum, 1])).^2);

DistTemp = repmat(DistTagVoxel, RecvNum, 1) + repelem(DistRecvVoxel, TagNum, 1);
DistUplink = repmat(DistTemp, length(freq), 1);

% Calculate the frequency constant.
FreqConst = (repelem(-1j*2*pi*freq/c, TagNum*RecvNum)).';
Freq_ratio = (repelem(freq/min(freq), TagNum*RecvNum).^2).';
% A = Freq_ratio.*exp(FreqConst.*DistUplink);
A = exp(FreqConst.*DistUplink);

% Reconstruct the object reflectivity for each voxel.
tic;
imgComplex = A'*b;
tComp = toc;

%% Post-Processing.
% Apply threshold for reconstructed reflectivity intensity.
imgComplexAbs = abs(imgComplex).^2;
% imgComplexAbs = abs(imgComplex);
ReconsDistMax = max(imgComplexAbs);
ReconsDistMin = min(imgComplexAbs);
ReconsDistNorm = (imgComplexAbs - ReconsDistMin)/(ReconsDistMax - ReconsDistMin);

Threshold = ReconsDistMax*0.9;
% image_3D = imgComplexAbs;
% image_3D(imgComplexAbs < Threshold)=nan;
% image_3D = reshape(image_3D, [length(x_v), length(y_v), length(z_v)]);
imgComplexAbs = reshape(imgComplexAbs, [length(x_v), length(y_v), length(z_v)]);
imgComplexAbs(imgComplexAbs < Threshold) = 0;
imgComplexAbs(imgComplexAbs >= Threshold) = 1;

% Clustering.
roomSize = [x_v(1), x_v(length(x_v)); y_v(1), y_v(length(y_v)); z_v(1), z_v(length(z_v))];
voxelSize = [x_v(2)-x_v(1); y_v(2)-y_v(1); z_v(2)-z_v(1)];
[cDistribution, clusters,centroidDist] = i4block_components(imgComplexAbs, roomSize, voxelSize);

fprintf('Initial cluster number = %d\n',length(clusters.centroid));
for i = 1:size(clusters.centroid,1)
    fprintf('Clusters centroid: [%3.2f, %3.2f, %3.2f] with element number %d\n',clusters.centroid(i,:),clusters.elemNum(i));
end
fprintf('\n');
opts.distTh = 0.2; % distance threshold, clusters with centers closer than this will be combined
opts.XYdistTh = 0.2;
opts.elemNumTh = 0.61; % clusters with element number less than 60% of the maximum will be rejected
opts.minHeightRatio = 0.6; % Minimum height ratio compared to largest object, exact ht depends on voxel size etc.
clusterOut = clusterProcess(clusters,opts);

centroid = clusterOut.centroid;
elemNum = clusterOut.elemNum;
for i = 1:size(centroid,1)
    fprintf('Clusters centroid: [%3.2f, %3.2f, %3.2f] with element number %d\n',centroid(i,:),elemNum(i));
end


%% Plotting.
% Convert the relative/normalized brightness of the reconstructed
% reflectivity vector into a 3D matrix.
ReconsDistNorm = reshape(ReconsDistNorm, [length(x_v), length(y_v), length(z_v)]);

% Plot the relative/normalized brightness of the reconstructed
% reflectivity, in a 3D grid.
[X_V, Y_V, Z_V] = meshgrid(x_v, y_v, z_v);
ReconsDistNorm = permute(ReconsDistNorm, [2 1 3]);
h = slice(X_V, Y_V, Z_V, ReconsDistNorm,x_v,y_v,z_v);
% imgComplexAbs(imgComplexAbs == 0) = NaN;
% h = slice(X_V, Y_V, Z_V, imgComplexAbs,x_v,y_v,z_v);
xlabel('x / m','FontSize',14);
ylabel('y / m','FontSize',14);
zlabel('z / m','FontSize',14);
xlim([x_v(1), x_v(length(x_v))]);
ylim([y_v(1), y_v(length(y_v))]);
zlim([z_v(1), z_v(length(z_v))]);

set(h, 'EdgeColor','none',...
    'FaceColor','interp',...
    'FaceAlpha','interp');
alpha('color');
a = alphamap('rampup',256);
imgThresh = 200;
a(1:imgThresh) = 0;
alphamap(a);
set(gca, 'fontweight', 'bold');

hold on
scatter3(2.4, 1.8, 1.2, 5, 'bo', 'LineWidth', 5);

