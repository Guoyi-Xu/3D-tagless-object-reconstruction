function [coordinates, location, image_3D, max_intense, x_v, y_v, z_v] = imaging_3D_func(filename_wo, filename_w)
%% Data collection block.
% movefile data_wo.xlsx test3_wo.xlsx
% movefile data_w.xlsx test3_w.xlsx

% Load the data.
load(filename_wo);
data_wo = [chindexlist, tagindexlist, antennalist, rssiimpinjlist, rssiimpinjlist_d, phasedeglist];
load(filename_w);
data_w = [chindexlist, tagindexlist, antennalist, rssiimpinjlist, rssiimpinjlist_d, phasedeglist];
clear chindexlist tagindexlist antennalist rssiimpinjlist rssiimpinjlist_d phasedeglist


% Generate the complex matrix.
[G_wo_s(:, :, 1), G_wo_s(:, :, 2), G_wo_s(:, :, 3), G_wo_s(:, :, 4), G_wo_s(:, :, 5), G_wo_s(:, :, 6)] = data_generator_func(data_wo);
[G_w_s(:, :, 1), G_w_s(:, :, 2), G_w_s(:, :, 3), G_w_s(:, :, 4), G_w_s(:, :, 5), G_w_s(:, :, 6)] = data_generator_func(data_w);

% Find the zero values in the calibration matrix, and replace them with the
% maximum value in the matrix. (We can replace those zeros with any number
% we like, actually.)
empt = find(G_wo_s == 0);
G_wo_s_abs = abs(G_wo_s);
max_mag = max(max(max(G_wo_s_abs)));
index = find(G_wo_s_abs == max_mag);
G_wo_s_new = G_wo_s;
G_wo_s_new(empt) = G_wo_s(index);

% Replace the values of G_w_s at locations where the G_wo_s is equal to
% zero, with the replaced values (the maximum value of G_wo_s_abs) of
% G_wo_s.
G_w_s_new = G_w_s;
G_w_s_new(empt) = G_wo_s_new(empt);
loc_zero = find(G_w_s_new == 0);
G_w_s_new(loc_zero) = G_wo_s_new(loc_zero);

% Create the coefficient matrix and obtain other variables.
[Coeff, Freq, k1, k2, lamda0, r, x_v, y_v, z_v] = coefmat_generator_func();

% load('coef_mat.mat');
%% Differential receiving block.
% Calculate the input data results after applying differential receiving.
for n=1:k2
    for n_pair=1:k2
        G_wo_s_ratio(:,n,:,n_pair) = G_wo_s_new(:,n,:)./G_wo_s_new(:,n_pair,:);
        G_w_s_ratio(:,n,:,n_pair) = G_w_s_new(:,n,:)./G_w_s_new(:,n_pair,:);
    end    
end

G_calib_ratio = (G_w_s_ratio./G_wo_s_ratio-1).*exp(-1i*2*pi./lamda0.*r);
G_calib = sum(G_calib_ratio,4);

% Extend the differential receiving matrix to 6-D for afterwards
% element-wise matrix multiplications.
G_calib_XYZ = zeros(length(x_v),length(y_v),length(z_v),k1,k2,length(Freq));
G_calib_XYZ(1,1,1,:,:,:) = G_calib;
G_calib_XYZ = repmat(G_calib_XYZ(1,1,1,:,:,:),[length(x_v),length(y_v),length(z_v),1,1,1]);

%% Image reconstruction block.
% Set the imaging threshold.
% Threshold = 80e4;

% Calculate the 3-D image. (as well as its 2-D projection)
image_temp = Coeff.* G_calib_XYZ; 
image_tt = sum(image_temp, 6);
image_temp = image_tt;
image_tt = sum(image_temp, 5);
image_temp = image_tt;
image_output = sum(image_temp, 4);
imge2D = sum(image_output, 3);
image_output_intensity = abs(image_output).^2;
image = abs(imge2D);

image_3D = image_output_intensity;
max_intense = max(max(max(image_3D)));
Threshold = max_intense*0.35;
image_3D(image_3D<Threshold)=nan;

location = find(image_3D == max_intense);
coordinates = get_coordinates(location, x_v, y_v, z_v);
disp(coordinates);

% save('image_3D.mat', 'location', 'coordinates', 'image_3D', 'max_intense', 'x_v', 'y_v', 'z_v');


% % Plot the image.
% figure(1);
% x0=600;
% y0=350;
% width=800;
% height=600;
% set(gcf,'position',[x0,y0,width,height]);
% 
% [X_V, Y_V, Z_V] = meshgrid(x_v, y_v, z_v);
% image_3D = permute(image_3D, [2 1 3]);
% h = slice(X_V, Y_V, Z_V, image_3D,[],[],z_v);
% set(h, 'EdgeColor','none', 'FaceColor','interp')
% set(gca, 'FontSize', 30);
% xlabel('x label/m');
% ylabel('y label/m');
% zlabel('z label/m');
% % title('Object location: (1.4 m, 1.64 m)', 'FontSize', 30);
% 
% savefig(gcf, 'figure1.fig');


end
