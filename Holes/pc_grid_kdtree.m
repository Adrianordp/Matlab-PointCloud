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
xyz_rotated = (R*xyz')'; %pequena melhoria
xy_rotated = xyz_rotated(:,1:2); %pega XY

%% Faz Grid
X = xyz_rotated(:,1);
Y = xyz_rotated(:,2);

res = 50; %resolução em metros
dx = max(X) - min(X);
dy = max(Y) - min(Y);

% pegamos o maior dentre (dx,dy), e fazemos grid quadrático
if (dx > dy)
    dcell = dx;
else
    dcell = dy;
end

cell_length = ceil(dcell/res);


% Mesh
xcells = min(X):cell_length:max(X)
ycells = min(Y):cell_length:max(Y)

%Nossa grid foi definida agora.
nx = length(xcells);
ny = length(ycells);

ncells = (nx-1) * (ny-1);
plotgrid(xcells,ycells)

%% separa pontos de query da nuvem
% serão n_cells procuras
% pega pontos medios
xmed = zeros(1,nx-1);
ymed = xmed;
for i=1:nx-1
   xmed(i) = (xcells(i)+xcells(i+1))/2;    
end

for i=1:ny-1
   ymed(i) = (ycells(i)+ycells(i+1))/2;    
end



qpoints = [];
for i=1:nx-1
    for k=1:ny-1
        qpoints = [qpoints;
                  xmed(i) ymed(k)];
    end
end

plot(qpoints(:,1),qpoints(:,2),'.');


%% KD Tree
% KD = KDTreeSearcher(xy_rotated,'BucketSize',10,'distance','minkowski','P',1);
KD = KDTreeSearcher(xy_rotated,'BucketSize',10,'distance','chebychev');


p=qpoints(3,:)

r = (cell_length/2)*0.95; %raio ~ quadrado de lado 20!
IDX = rangesearch(KD,p,r); %pega os mais próximos de raio r
idxs = IDX{1};

%% Plots
x_range = xyz_rotated([idxs],1);
y_range = xyz_rotated([idxs],2);

x_all = xyz_rotated(:,1);
y_all = xyz_rotated(:,2);

plot(x_all,y_all,'.');hold on;
plot(x_range,y_range,'g.')
return


total_ums = sum(M(:));
preenchimento = total_ums / (n_cells)

%% KD Tree test































     

