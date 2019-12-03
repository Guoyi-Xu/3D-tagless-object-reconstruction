close all
clear
clc

load('image_3D_Pragya_standing_x6y6_rm0.mat');
Threshold = max_intense*0.5;
image_3D(image_3D<Threshold)=nan;

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
set(gca,'xtick',0:1.2:y_v(length(y_v)));
set(gca,'ytick',0:1.2:x_v(length(y_v)));

set(h, 'EdgeColor','none', 'FaceColor','interp')
set(gca, 'FontSize', 18);
xlabel('x / m');
ylabel('y / m');
zlabel('z / m');

hold on
scatter3(1.8, 1.8, 1.2, 200, 'ro', 'LineWidth', 2);
% title(['Object location: (', num2str(coordinates(1)), ' m, ', num2str(coordinates(2)), ' m)'], 'FontSize', 22);

% savefig(gcf, 'me_standing_x6y6_v1-2.fig');



% figure(2);
% imagesc(image)
% xlabel('x label');
% ylabel('y label');

% figure(3);
% imagesc(image_3D(:,:,32))




