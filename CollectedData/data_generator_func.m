function [C1, C2, C3, C4, C5, C6] = data_generator_func(data)

% System parameters

ChaNum = 6;

TagNum = 80;

AntNum = 4;

%% Extract the data.
% The six columns of the data are (from left to right): channel index, tag
% index, antenna ID, tag RSSI/dBm, tag RSSI/nW, and tag phase.

% [~, text, raw] = xlsread(filename, 1, 'A2:F4000');
% data = str2double(raw);


%% Data extraction and grouping for each channel
C = zeros(ChaNum, TagNum, AntNum);

for i = 1:ChaNum
    [~, ~, C(i, :, :)] = data_mat(data, i, TagNum, AntNum);
end

C1 = C(1, :, :);
C2 = C(2, :, :);
C3 = C(3, :, :);
C4 = C(4, :, :);
C5 = C(5, :, :);
C6 = C(6, :, :);

end
