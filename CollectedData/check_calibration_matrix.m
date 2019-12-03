close all
clear
clc

% System parameters

ChaNum = 6;

TagNum = 80;

AntNum = 4;

%% Extract the data.
% [~, text, raw] = xlsread(filename, 1, 'A2:F4000');
% The six columns of the data are (from left to right): channel index, tag
% index, antenna ID, tag RSSI/dBm, tag RSSI/nW, and tag phase.

% data = str2double(raw);

% Load the data.
load('data_wo_me_standing_x9y3_v1_to_v2.mat');
data = [chindexlist, tagindexlist, antennalist, rssiimpinjlist, rssiimpinjlist_d, phasedeglist];

%% Data extraction and grouping for each channel

clist = data(:, 1);
tlist = data(:, 2);
alist = data(:, 3);
rlist = data(:, 4);
rdlist = data(:, 5);
plist = data(:, 6);
R = zeros(ChaNum, TagNum, AntNum);
P = zeros(ChaNum, TagNum, AntNum);
C = zeros(ChaNum, TagNum, AntNum);
for i = 1:AntNum
    % Extract the raw data only for the specified antenna.
    clist_t = clist(alist == i);
    tlist_t = tlist(alist == i);
    alist_t = alist(alist == i);
    rlist_t = rlist(alist == i);
    rdlist_t = rdlist(alist == i);
    plist_t = plist(alist == i);
    data_t = [clist_t, tlist_t, alist_t, rlist_t, rdlist_t, plist_t];
    
    % group the data and get the 3-D data matrix
    for cha = 1:ChaNum
        for j = 1:TagNum
            [~, ~, R(cha, j, i), P(cha, j, i), C(cha, j, i)] = data_group(data_t, cha, j, i);
        end
    end
end

% Check whether the calibration matrix has zeros.
data_wo_containszeros = 0;
empt = find(C == 0);
empt_index = cell(length(empt), 1);

% Store the number of lost data for each antenna.
emptant1 = 0;
emptant2 = 0;
emptant3 = 0;
emptant4 = 0;
if(~isempty(empt))
    disp("Warning: The calibration matrix contains zero(s)! Please calibrate again.");
    data_wo_containszeros = 1;
    % Store the 3-D index of the zero-valued elements of the matrix.
    for i = 1:length(empt)
        dim3 = floor(empt(i)/(ChaNum*TagNum));
        if(empt(i) - dim3*ChaNum*TagNum == 0)
            dim2 = TagNum;
            dim1 = ChaNum;
            empt_index{i, 1} = [dim1, dim2, dim3];
        else
            dim2 = floor((empt(i) - dim3*ChaNum*TagNum)/ChaNum);
            if((empt(i) - dim3*ChaNum*TagNum) - dim2*ChaNum == 0)
                dim1 = ChaNum;
                empt_index{i, 1} = [dim1, dim2, dim3+1];
            else
                dim1 = empt(i) - dim3*ChaNum*TagNum - dim2*ChaNum;
                empt_index{i, 1} = [dim1, dim2+1, dim3+1];
            end
        end
        disp(empt_index{i, 1});
        disp(" ");
        
        if(empt_index{i, 1}(3) == 1)
            emptant1 = emptant1 + 1;
        elseif(empt_index{i, 1}(3) == 2)
            emptant2 = emptant2 + 1;
        elseif(empt_index{i, 1}(3) == 3)
            emptant3 = emptant3 + 1;
        elseif(empt_index{i, 1}(3) == 4)
            emptant4 = emptant4 + 1;
        end
        
%         % Assign non-zero values for zero elements in the calibration
%         % matrix.
%         if(empt(i) - 50 > 0)
%             C(empt(i)) = C(empt(i) - ChaNum*TagNum);
%         else
%             C(empt(i)) = C(empt(i) + ChaNum*TagNum);
%         end
        
    end
    disp(" ");
else
    disp("This calibration matrix is perfect! It doesn't contain zeros.")
    disp(" ");
end

% Display the lost data number for each antenna.
msg1 = strcat('Ant 1 loses: ', num2str(emptant1), ' elements.');
msg2 = strcat('Ant 2 loses: ', num2str(emptant2), ' elements.');
msg3 = strcat('Ant 3 loses: ', num2str(emptant3), ' elements.');
msg4 = strcat('Ant 4 loses: ', num2str(emptant4), ' elements.');
disp(msg1);
disp(msg2);
disp(msg3);
disp(msg4);
