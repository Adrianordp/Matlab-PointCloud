clear; close all;clc
% xyz = load('pilha_densa.txt');
xyz = load('ensaio_bad_ss.txt');
xyz = xyz(:,1:3);
% pc = pointCloud(xyz)
% pcshow(pc)
%% Faz Regressão pra alinhar em X

% Y' = aX
% parametro = a

Y = xyz(:,2);
X = xyz(:,1);

a = inv((X'*X))*Y'*X;

% x0 = -114;
% y0 = a*x0;
% 
% x1 = -500;
% y1 = -500*a;
% angle = atan(a);
ANGLES = 0:0.001:2*pi;
iterations = length(ANGLES);
AREA = zeros(iterations,1);
GDIV = zeros(iterations,1);

parfor i=1:iterations

angle = ANGLES(i);
R = [cos(angle) sin(angle) 0;-sin(angle) cos(angle) 0;0 0 1];
xyz_rotated = (R*xyz')'; %pequena melhoria
xy_rotated = xyz_rotated(:,1:2); %pega XY

% plots para ver.. comentar depois
% plot(xy_rotated(:,1),xy_rotated(:,2),'blue.');
% hold on;
% plot(xyz(:,1),xyz(:,2),'red.');
% return
%% Faz Grid
X = xyz_rotated(:,1);
Y = xyz_rotated(:,2);

% inflation 
inflation_factor = 1;
% max_x = max(X)*inflation_factor;
% min_x = min(X)*inflation_factor;
% max_y = max(Y)*inflation_factor;
% min_y = min(Y)*inflation_factor;

max_x = ceil(max(X)) + inflation_factor;
min_x = floor(min(X)) - inflation_factor;
max_y = ceil(max(Y)) + inflation_factor;
min_y = floor(min(Y)) - inflation_factor;

%largura e altura do grid
dx = max_x - min_x;
dy = max_y - min_y;
%% Problema = encontrar menor quadrilátero com a menor/melhor resolução!
% centro do quadrilátero ?
x_sq = (max_x + min_x) / 2;
y_sq = (max_y + min_y) / 2;

m_factor = 1;
mult_factor = m_factor*dy;
while (mult_factor < dx)
   m_factor = m_factor+1;
   mult_factor = m_factor*dy;
end
err = mult_factor - dx;
new_dx = mult_factor;
new_dy = dy;

%precisa ser td inteiro
min_x_sq = x_sq - new_dx/2;
max_x_sq = x_sq + new_dx/2;
min_y_sq = y_sq - new_dy/2;
max_y_sq = y_sq + new_dy/2;
% return
area = (max_x_sq - min_x_sq)*(max_y_sq - min_y_sq) % <-- Quero minimizar
AREA(i) = area;
%%
g_div = gcd(new_dx,new_dy) % <--- Quero minimizar !
GDIV(i) = g_div;
end
plot(AREA)
figure
plot(GDIV)

% division_factor = 2^(-6); % 0 < div_factor <= 1 %Ajsute da precisão
% cell_length = g_div*division_factor;
% fprintf('Resolução Cell = %d m\n',cell_length);
% 
% % Mesh do 
% xcells = min_x_sq:cell_length:max_x_sq;
% ycells = min_y_sq:cell_length:max_y_sq;
% 
% %Nossa grid foi definida agora.
% nx = length(xcells);
% ny = length(ycells);
% 
% ncells = (nx-1) * (ny-1); 
% % p = [1000 1000]
% x_all = xyz_rotated(:,1);
% y_all = xyz_rotated(:,2);
% plotgrid(xcells,ycells);
% plot(x_all,y_all,'.');hold on;
% ncells
return
