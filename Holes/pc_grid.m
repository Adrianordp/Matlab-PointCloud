clear; close all;clc
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
% angle = 0;
R = [cos(angle) sin(angle) 0;-sin(angle) cos(angle) 0;0 0 1];

%%  Rotate
xyz_rotated = (R*xyz')'; %pequena melhorada

% KD = KDTreeSearcher(xyz_rotated,'distance','euclidean');
KD = KDTreeSearcher(xyz_rotated,'distance','minkowski','P',8);

Y = xyz_rotated(1,:); %pnto de prova

% knnsearch(KD,Y) %pega o mais próximo
% knnsearch(KD,Y,'K',20) %pega K mais próximos

r = 2;
IDX = rangesearch(KD,Y,r) %pega os mais próximos de raio r


return


% save('pilha_rotated.txt','xyz_rotated','-ASCII')

%% Faz Grid
X = xyz_rotated(:,1);
Y = xyz_rotated(:,2);

res = 10; %resolução em metros
dx = max(X) - min(X);
dy = max(Y) - min(Y);
width = ceil(dx/res);
heigth = ceil(dy/res);


% Mesh
xcells = linspace(min(X),max(X),width);
ycells = linspace(min(Y),max(Y),heigth);
%Nossa grid foi definida agora.
n_cells = (length(xcells)-1) * (length(ycells)-1);

%% plota
close all
plotgrid(xcells,ycells)
plot(X,Y,'.')

M = fill_grid([X Y],xcells,ycells);
M_ = flipud(M)


total_ums = sum(M(:));
preenchimento = total_ums / (n_cells);

%% KD Tree test































     

