clear; close all;clc
% xyz = load('pilha_densa.txt');
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

% na verdade, nao precisamos fatiar exatamente as dimensões da pilha. pode
% ser algo mais 'comportado', ou 'quadrado'.

% bolar um algoritmo para que gcd(dx,dy) seja maximizado. o GCD é o tamanho máximo
% da célula. quato maior, melhor

% encontarando máximo multiplo da menor dimensão mais próximo da maior
% dimensão
if (dx > dy)
    minor_dim = dy;
    minor_vec = Y;
    major_dim = dx;
    major_vec = X;
else
    minor_dim = dx; 
    minor_vec = X;
    major_dim = dy;
    major_vec = Y;
end

minor_multiple = minor_dim;
factor = 1;
while (minor_multiple < major_dim)
    factor = factor+1;
    minor_multiple = minor_multiple*factor;
end
err = minor_multiple - major_dim;
%recalcula dx
% testar 
min_x_corrected = min_x - err/2;
max_x_corrected = max_x + err/2;
new_dx = max_x_corrected - min_x_corrected

g_div = gcd(new_dx,dy)
% return
division_factor = 0.1; % 0 < div_factor <= 1
cell_length = g_div*division_factor;
fprintf('Res = %d m\n',cell_length);

% Mesh
xcells = min_x_corrected:cell_length:max_x_corrected;
ycells = min_y:cell_length:max_y;

%Nossa grid foi definida agora.
nx = length(xcells);
ny = length(ycells);

ncells = (nx-1) * (ny-1);
% p = [1000 1000]
x_all = xyz_rotated(:,1);
y_all = xyz_rotated(:,2);
plotgrid(xcells,ycells);
plot(x_all,y_all,'.');hold on;
return
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

% query points
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
KD = KDTreeSearcher(xy_rotated,'BucketSize',50,'distance','chebychev');

r = (cell_length/2); %raio ~ quadrado de lado 2r!

interations = (nx-1)*(ny-1);
M = zeros(1,interations);
for i=1:interations
   i
p = qpoints(i,:);
IDX_ = rangesearch(KD,p,r); 
IDX = IDX_{1};

if(isempty(IDX))
   M(i) = 0; 
else
   M(i) = 1; 
end

%% Plots
x_range = xyz_rotated([IDX],1);
y_range = xyz_rotated([IDX],2);

plot(x_range,y_range,'g.')
drawnow
% pause(.1)
end
total_ums = sum(M(:))
preenchimento = total_ums / (ncells)

%% KD Tree test































     

