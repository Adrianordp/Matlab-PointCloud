clear; close all;clc
% xyz = load('pilha_densa.txt');
filename = 'pontos_1574704717_seg.asc'; %AQUI VC SELECIONA A NUVEM JÁ SEGMENTADA
xyz = load(filename); 
xyz = xyz(:,1:3);
xy = xyz(:,1:2);
% pc = pointCloud(xyz)
% pcshow(pc)
%% Faz Regressão pra alinhar em X

% Y' = aX
% parametro = a

% Y = xyz(:,2);
% X = xyz(:,1);

x_c = (max(xyz(:,1)) + min(xyz(:,1)) )/2;
y_c = (max(xyz(:,2)) + min(xyz(:,2)) )/2;

% a = inv((X'*X))*Y'*X;
% TOTCELLS = [];
% ANGLES = -0.180:-0.001:-0.3;
% for angle=ANGLES

%     angle = ANGLES(i);
angle = -0.188;% 0.188 ALINHA
R = [cos(angle) -sin(angle);sin(angle) cos(angle)]; %CW

xy_rotated = (R*(xy' - [x_c;y_c]) + [x_c;y_c])'; %em torno do centro da pilha
% save('xy_rotated.asc','xy_rotated','-ASCII')
% return
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

xy_real = xy_rotated - [map.tfx-2 map.tfy-2]; %-2 ?
% plot(xy_(:,1),xy_(:,2),'.')


% Tem que ver isso aqui
RES = 0.1:0.1:100
total_res = length(RES);
for j=1:total_res
res = RES(j) %m
adjust = 1/res - 1;
xy_ = xy_real / res - [adjust adjust];

xy_int = round(xy_); % 'fix' ou 'round'
x_max_i = max(xy_int(:,1));
y_max_i = max(xy_int(:,2));

% TOTCELLS = [TOTCEL LS;x_max_i*y_max_i];

% end
% plot(TOTCELLS)
% return

M = zeros(x_max_i,y_max_i);

n = length(xy_int);
for i=1:n
   M(xy_int(i,1),xy_int(i,2)) = 1 ;    
end

% plotgrid2(M);
% plotmatrix(M);

[row,col] = size(M);

total_cells = row*col;
total_filled = sum(M(:));
pct(j) = total_filled/total_cells

end

plot(RES,1-pct);
grid on
xlabel('Res [m]');
ylabel('Holes %');
legend(filename)
fig_name = strcat(filename,'.fig');
savefig(fig_name)


% x-> column
% y -> row




