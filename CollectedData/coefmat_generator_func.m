function [Coeff, Freq, k1, k2, lamda0, r, x_v, y_v, z_v] = coefmat_generator_func()

% The speed of light.
c = 3e8;

% The voxel coordinates.
x_v=0:0.05:3.6;
y_v=0:0.05:3.6;
z_v=0:0.1:2.8;

% Frequencies used.
Freq = [905.25, 909.25, 913.25, 917.25, 921.25, 925.25]* 1e6;

% Tags positions. (64 tags)
Pos_tag=[121.7-37.2, 0, 122.8; 121.7-28.4, 0, 122.6; 121.7-18.7, 0, 123.1; 121.7-8.4, 0, 123.1;
    121.7+7.3, 0, 123.2; 121.7+18.3, 0, 123.3; 121.7+28.3, 0, 123.4; 121.7+38.6, 0, 123.4;
    365.6, 121.6-37.7, 123.3; 365.6, 121.6-28.6, 123.5; 365.6, 121.6-18.2, 123.3; 365.6, 121.6-8, 123.5;
    365.6, 121.6+8.9, 123.3; 365.6, 121.6+18.5, 123.4; 365.6, 121.6+27.7, 123.6; 365.6, 121.6+36.9, 123.4;
    
    243.7+38.3, 365.4, 120.7; 243.7+29.1, 365.4, 120.7; 243.7+18.2, 365.4, 120.4; 243.7+8.7, 365.4, 120.9;
    243.7-9.1, 365.4, 120.9; 243.7-18.2, 365.4, 120.9; 243.7-27.5, 365.4, 120.7; 243.7-37.8, 365.4, 121.2;
    0, 243.8+37.9, 123.8; 0, 243.8+28.6, 123.6; 0, 243.8+18.7, 123.8; 0, 243.8+8.8, 123.8;
    0, 243.8-8.7, 123.7; 0, 243.8-17.9, 124; 0, 243.8-29, 124.3; 0, 243.8-38.7, 124.3;
    
    243.7-37.2, 0, 123.5; 243.7-27.4, 0, 123.6; 243.7-17.6, 0, 123.6; 243.7-8.4, 0, 123.6;
    243.7+8.4, 0, 123.3; 243.7+18.1, 0, 123.3; 243.7+29.1, 0, 123.5; 243.7+38.6, 0, 123.7;
    365.6, 243.8-37, 123.3; 365.6, 243.8-26.8, 123.4; 365.6, 243.8-17.2, 123.5; 365.6, 243.8-7.9, 123.4;
    365.6, 243.8+9.4, 123; 365.6, 243.8+19.8, 123; 365.6, 243.8+29.5, 123; 365.6, 243.8+38.5, 123.2;
    
    121.7+35.8, 365.4, 124.4; 121.7+27.1, 365.4, 124.2; 121.7+17, 365.4, 124.3; 121.7+7.1, 365.4, 124.3;
    121.7-7, 365.4, 123.4; 121.7-16.1, 365.4, 123.2; 121.7-27.1, 365.4, 123.3; 121.7-38, 365.4, 123.5;
    0, 121.6+36.9, 123.9; 0, 121.6+25.8, 124.2; 0, 121.6+16.9, 123.8; 0, 121.6+7.3, 123.9;
    0, 121.6-7.8, 124; 0, 121.6-17.5, 123.6; 0, 121.6-27.6, 123.6; 0, 121.6-39.1, 123.4;
    
    17.3+36.3*sqrt(2), 15.2-36.3*sqrt(2), 99.9; 17.3+14.9*sqrt(2), 15.2-14.9*sqrt(2), 100.2; 17.3-13.4*sqrt(2), 15.2+13.4*sqrt(2), 102.3; 17.3-35.1*sqrt(2), 15.2+35.1*sqrt(2), 102.1;
    351+35.1*sqrt(2), 14.2+35.1*sqrt(2), 98.1; 351+14*sqrt(2), 14.2+14*sqrt(2), 98; 351-11.8*sqrt(2), 14.2-11.8*sqrt(2), 97.8; 351-35.6*sqrt(2), 14.2-35.6*sqrt(2), 98.1;
    351.6-36*sqrt(2), 349.8+36*sqrt(2), 92.7; 351.6-12.1*sqrt(2), 349.8+12.1*sqrt(2), 93.4; 351.6+9.6*sqrt(2), 349.8-9.6*sqrt(2), 93.7; 351.6+35.2*sqrt(2), 349.8-35.2*sqrt(2), 94.6;
    15.6-36.7*sqrt(2), 351.6-36.7*sqrt(2), 94; 15.6-12.5*sqrt(2), 351.6-12.5*sqrt(2), 94.5; 15.6+11.7*sqrt(2), 351.6+11.7*sqrt(2), 94.8; 15.6+36.1*sqrt(2), 351.6+36.1*sqrt(2), 94.9]*0.01;

% Receiver antennas locations.
Pos_recv=[228, 143.8, 256.6;
    162.4, 145.9, 255.3;
    224.6, 204.2, 256.4;
    163.3, 203.4, 255.9]*0.01;

% Pos_recv=[0, 0, 180.7;
%     273, 0, 176.9;
%     273.5, 272.9, 181.2;
%     0, 276.4, 179]*0.01;

% Number of tags and receiver antennas.
[k1, l1]=size(Pos_tag);
[k2, l2]=size(Pos_recv);

% Variables to be used in imaging scripts.
freq = zeros(k1,k2,length(Freq),k2);
freq(1,1,:,1) = Freq;
freq = repmat(freq(1,1,:,1) ,[k1,k2,1,k2]);
lamda0 = c./freq;
r_x=bsxfun(@minus,Pos_tag(:,1),Pos_recv(:,1)');
r_y=bsxfun(@minus,Pos_tag(:,2),Pos_recv(:,2)');
r_z=bsxfun(@minus,Pos_tag(:,3),Pos_recv(:,3)');
r = sqrt(r_x.^2+r_y.^2+r_z.^2);
r = repmat(r, [1,1,1,k2]);

% The 6th dimension contains the frequency channels.
lamda = zeros(length(x_v),length(y_v),length(z_v),k1,k2,length(Freq));
lamda(1,1,1,1,1,:) = c./Freq;
lamda = repmat(lamda(1,1,1,1,1,:), [length(x_v),length(y_v),length(z_v),k1,k2,1]);
ratio_freq = zeros(length(x_v),length(y_v),length(z_v),k1,k2,length(Freq));
ratio_freq(1,1,1,1,1,:) = Freq/min(Freq);
ratio_freq = repmat(ratio_freq(1,1,1,1,1,:) ,[length(x_v),length(y_v),length(z_v),k1,k2,1]);
ratio_phi = 1; 

% The first three dimensions contain voxel locations.
px = zeros(length(x_v),length(y_v),length(z_v),k1,k2,length(Freq));
py = zeros(length(x_v),length(y_v),length(z_v),k1,k2,length(Freq));
pz = zeros(length(x_v),length(y_v),length(z_v),k1,k2,length(Freq));
px(:,1,1,1,1,1) = x_v;
px = repmat(px(:,1,1,1,1,1), [1,length(y_v),length(z_v),k1,k2,length(Freq)]);
py(1,:,1,1,1,1) = y_v;
py = repmat(py(1,:,1,1,1,1), [length(x_v),1,length(z_v),k1,k2,length(Freq)]);
pz(1,1,:,1,1,1) = z_v;
pz = repmat(pz(1,1,:,1,1,1), [length(x_v),length(y_v),1,k1,k2,length(Freq)]);

% The 4th dimension contains the tag locations.
p_tagx = zeros(length(x_v),length(y_v),length(z_v),k1,k2,length(Freq));
p_tagy = zeros(length(x_v),length(y_v),length(z_v),k1,k2,length(Freq));
p_tagz = zeros(length(x_v),length(y_v),length(z_v),k1,k2,length(Freq));
p_tagx(1,1,1,:,1,1) = Pos_tag(:,1);
p_tagy(1,1,1,:,1,1) = Pos_tag(:,2);
p_tagz(1,1,1,:,1,1) = Pos_tag(:,3);
p_tagx = repmat(p_tagx(1,1,1,:,1,1), [length(x_v),length(y_v),length(z_v),1,k2,length(Freq)]);
p_tagy = repmat(p_tagy(1,1,1,:,1,1), [length(x_v),length(y_v),length(z_v),1,k2,length(Freq)]);
p_tagz = repmat(p_tagz(1,1,1,:,1,1), [length(x_v),length(y_v),length(z_v),1,k2,length(Freq)]);

% The 5th dimension contains the receiver antenna locations.
p_recvx = zeros(length(x_v),length(y_v),length(z_v),k1,k2,length(Freq));
p_recvy = zeros(length(x_v),length(y_v),length(z_v),k1,k2,length(Freq));
p_recvz = zeros(length(x_v),length(y_v),length(z_v),k1,k2,length(Freq));
p_recvx(1,1,1,1,:,1) = Pos_recv(:,1);
p_recvy(1,1,1,1,:,1) = Pos_recv(:,2);
p_recvz(1,1,1,1,:,1) = Pos_recv(:,3);
p_recvx = repmat(p_recvx(1,1,1,1,:,1), [length(x_v),length(y_v),length(z_v),k1,1,length(Freq)]);
p_recvy = repmat(p_recvy(1,1,1,1,:,1), [length(x_v),length(y_v),length(z_v),k1,1,length(Freq)]);
p_recvz = repmat(p_recvz(1,1,1,1,:,1), [length(x_v),length(y_v),length(z_v),k1,1,length(Freq)]);

% Combine 1st with 4th dimension to get tag-to-voxel distance, and combine
% 1st with 5th dimension to get receiver-to-voxel distance.
R1 = sqrt((px-p_tagx).^2+(py-p_tagy).^2+(pz-p_tagz).^2);
% clearvars -except Freq k1 k2 lamda0 r x_v y_v z_v ratio_freq ratio_freq lamda R1 ratio_phi px p_recvx py p_recvy pz p_recvz;
R2 = sqrt((px-p_recvx).^2+(py-p_recvy).^2+(pz-p_recvz).^2);
% clearvars -except Freq k1 k2 lamda0 r x_v y_v z_v ratio_freq ratio_freq lamda R1 R2  ratio_phi;

Coeff = ratio_freq.^2 .* ratio_phi .* exp(1j.*2.*pi./lamda.*(R1+R2));

% clearvars -except Coeff Freq k1 k2 lamda0 r x_v y_v z_v;

end
