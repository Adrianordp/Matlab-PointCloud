% clear ;close all;clc
% rosshutdown
% rosinit
clear;close all;clc
load('Scans')
pc0 = [scan0.Ranges.*cos(scan0.Angles) scan0.Ranges.*sin(scan0.Angles) ];
pc1 = [scan1.Ranges.*cos(scan1.Angles) scan1.Ranges.*sin(scan1.Angles) ];

x = my_scanmatch(pc0,pc1,[0 0 0]')
transform_scan0 = p_transform(pc0,x)  ;
scan0_transformed = lidarScan(transform_scan0');
%% Plots
close all
plot(scan0)
hold on
plot(scan1)
plot(scan0_transformed)
legend('Posição Inicial','Posição Final','Transformação')





