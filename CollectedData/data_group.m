function [dr, dp, dr_avg, dp_avg, dc_avg] = data_group(data, channel, tag, ant)

channelind = data(:, 1);
tagind = data(:, 2);
antennaind = data(:, 3);
rssi = data(:, 5);
phase = data(:, 6);

j = 0;
dr_temp = 0;
dp_temp = 0;


for i = 1:length(antennaind)
    if(channelind(i) == channel && tagind(i) == tag && antennaind(i) == ant && ~isnan(tagind(i)))
        j = j + 1;
        if (j == 1) 
            dr_temp = rssi(i);
            dp_temp = phase(i);
        else
            dr_temp = cat(1, dr_temp, rssi(i));
            dp_temp = cat(1, dp_temp, phase(i));
        end
    end
end
dr_sort = sort(dr_temp);
dp_sort = sort(dp_temp);

% threshold for deciding whether a value is too large or too small
dev = 0;    % the percentage of the dataset to be truncated
dr_min = dr_sort(round(length(dr_sort)*dev) + 1);
dr_max = dr_sort(length(dr_sort)-round(length(dr_sort)*dev));
dp_min = dp_sort(round(length(dp_sort)*dev) + 1);
dp_max = dp_sort(length(dp_sort)-round(length(dp_sort)*dev));
% matrices used for recording the values which are too large or too small
rmat = 0;
pmat = 0;
for i = 1:length(dr_temp)
    if (dr_temp(i) < dr_min || dr_temp(i) > dr_max)
        rmat = cat(1, rmat, dr_temp(i));
    end
end
for i = 1:length(dp_temp)
    if (dp_temp(i) < dp_min || dp_temp(i) > dp_max)
        pmat = cat(1, pmat, dp_temp(i));
    end
end
rmat(1) = [];
pmat(1) = [];

% truncate the data too large or too small off the dataset
for i = 1:length(rmat)
    dp_temp(dr_temp == rmat(i)) = [];
    dr_temp(dr_temp == rmat(i)) = [];
end
for i = 1:length(pmat)
    dr_temp(dp_temp == pmat(i)) = [];
    dp_temp(dp_temp == pmat(i)) = [];
end

dr = dr_temp;
dp = dp_temp;
dr_avg = mean(dr);
dp_avg = mean(dp);

dc = dr.*cos(deg2rad(dp)) + sqrt(-1).*dr.*sin(deg2rad(dp));
dc_avg = mean(dc);

disp('Completed for:');
disp(channel);
disp(tag);
disp(ant);

end