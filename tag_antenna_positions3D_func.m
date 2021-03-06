function [tagPosition, rxPosition, freq] = tag_antenna_positions3D_func()

% Tags positions.
tagPosition = [121.7-37.2, 0, 122.8; 121.7-28.4, 0, 122.6; 121.7-18.7, 0, 123.1; 121.7-8.4, 0, 123.1;
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
rxPosition = [228, 143.8, 256.6;
                    162.4, 145.9, 255.3;
                    224.6, 204.2, 256.4;
                    163.3, 203.4, 255.9]*0.01;

% Frequency channels.
freq = [905.25, 909.25, 913.25, 917.25, 921.25, 925.25].* 1e6;

end
