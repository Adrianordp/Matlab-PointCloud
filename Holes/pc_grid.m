clear; close all;clc
% xyz = load('pilha_menor.txt');
xyz = load('pilha_menor.txt');
xyz = xyz(:,1:3);
pc = pointCloud(xyz)
% pcshow(pc)
%% Faz Regressão com os pontos

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
angle = atan(a);
angle = 0;
R = [cos(angle) sin(angle) 0;-sin(angle) cos(angle) 0;0 0 1];
xyz_rotated = (R*xyz')'; %pequena melhorada
%% Faz Grid
X = xyz_rotated(:,1);
Y = xyz_rotated(:,2);

res = 25; %resolução em metros
dx = max(X) - min(X);
dy = max(Y) - min(Y);

% pegamos o maior dentre (dx,dy), e fazemos grid quadrático
if (dx > dy)
    dcell = dx;
else
    dcell = dy;
end

width = ceil(dcell/res);
heigth = ceil(dcell/res);


% Mesh
xcells = linspace(min(X),max(X),width);
ycells = linspace(min(Y),max(Y),heigth);
%Nossa grid foi definida agora.
n_cells = (length(xcells)-1) * (length(ycells)-1);

M = fill_grid([X Y],xcells,ycells);
M_ = flipud(M)

%% plota
close all
plotgrid(xcells,ycells)
plot(X,Y,'.')

total_ums = sum(M(:));
preenchimento = total_ums / (n_cells)

%% KD Tree test































     

