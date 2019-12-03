close all
clear
clc

filename_wo = 'D:\reader_apps\project29\data_wo_me_standing_x6y6_v1_to_v2.mat';


Num = 2;
file_w = cell(Num, 1);
coordinate_3D = zeros(Num, 3);
location_tot = zeros(Num, 1);
max_mag = zeros(Num, 1);
for i = 1:2
    
    filename = strcat('D:\reader_apps\project29\data_w_me_standing_x6y6_v', num2str(i), '.mat');
    file_w{i} = filename;
    [coordinates, location, image_3D, max_intense, x_v, y_v, z_v] = imaging_3D_func(filename_wo, file_w{i});
    coordinate_3D(i, :) = coordinates;
    location_tot(i) = location;
    max_mag(i) = max_intense;
    save_name = ['image_3D_me_standing_x6y6_v', num2str(i), '.mat'];
    save(save_name, 'location', 'coordinates', 'image_3D', 'max_intense', 'x_v', 'y_v', 'z_v');
    
    % Plot the image.
    figure(1);
    x0=600;
    y0=350;
    width=800;
    height=600;
    set(gcf,'position',[x0,y0,width,height]);

    [X_V, Y_V, Z_V] = meshgrid(x_v, y_v, z_v);
    image_3D = permute(image_3D, [2 1 3]);
    h = slice(X_V, Y_V, Z_V, image_3D,[],[],z_v);
    
    grid on
    set(gca,'xtick',0:0.3:y_v(length(y_v)));
    set(gca,'ytick',0:0.3:x_v(length(y_v)));

    set(h, 'EdgeColor','none', 'FaceColor','interp')
    set(gca, 'FontSize', 18);
    xlabel('x label/m');
    ylabel('y label/m');
    zlabel('z label/m');
    title(['Object location: (', num2str(coordinates(1)), ' m, ', num2str(coordinates(2)), ' m)'], 'FontSize', 22);

    savefig(gcf, ['me_standing_x6y6_v', num2str(i), '.fig']);

    disp(i);
    
end

save('coordinate_3D.mat', 'coordinate_3D');
