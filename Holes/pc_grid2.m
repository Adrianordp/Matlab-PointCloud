clear; close all;clc
% xyz = load('pilha_densa.txt');
xyz = load('ensaio_bom_ss.txt');
xyz = xyz(:,1:3);
xy = xyz(:,1:2);
% pc = pointCloud(xyz)
% pcshow(pc)
%% Faz Regressão pra alinhar em X

% Y' = aX
% parametro = a

Y = xyz(:,2);
X = xyz(:,1);

a = inv((X'*X))*Y'*X;


angle = -0.188;% 0.188
R = [cos(angle) -sin(angle);sin(angle) cos(angle)]; %CW
x_c = (max(xyz(:,1)) + min(xyz(:,1)) )/2
y_c = (max(xyz(:,2)) + min(xyz(:,2)) )/2
xy_rotated = (R*(xy' - [x_c;y_c]) + [x_c;y_c])'; %em torno da origem

% xy_rotated = xyz_rotated(:,1:2); %pega XY

% plot(xy_rotated(:,1),xy_rotated(:,2),'r.')
% hold on
% plot(xy(:,1),xy(:,2),'b.')

%% xy to grid
% Translada grid para posição conversível para matrix
x_min = min(xy_rotated(:,1));
y_min = min(xy_rotated(:,2));
map.tfx = fix(x_min);
map.tfy = fix(y_min);

xy_ = xy_rotated - [map.tfx-2 map.tfy-2]; %-2 ?
% plot(xy_(:,1),xy_(:,2),'.')


% Tem que ver isso aqui
precision = 1; %m
res = 1/precision;
xy_ = xy_ * res;

xy_int = fix(xy_);
x_max_i = max(xy_int(:,1));
y_max_i = max(xy_int(:,2));

M = zeros(x_max_i,y_max_i);

n = length(xy_int);
for i=1:n
   M(xy_int(i,1),xy_int(i,2)) = 1 ;    
end

plotmatrix(M);
[row,col] = size(M);

total_cells = row*col;
total_filled = sum(M(:));
pct = total_filled/total_cells





% x-> column
% y -> row




