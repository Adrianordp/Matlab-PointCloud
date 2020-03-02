clear;close all;clc
%% Load Points
p0 = pcread('p0.pcd');
p1 = pcread('p1.pcd');
% pcshow(p0);
% hold on
% pcshow(p1)
X0 = p0.Location(:,1);
Y0 = p0.Location(:,2);

X1 = p1.Location(:,1);
Y1 = p1.Location(:,2);

P0 = [nonzeros(X0) nonzeros(Y0)];
P1 = [nonzeros(X1) nonzeros(Y1)];
%% Filtros
n = length(P0);
n1 = length(P1);

%% a Few plots
%% Plots
plot(P0(:,1),P0(:,2),'b.') % Original
hold on
plot(P1(:,1),P1(:,2),'r.') % Final
drawnow
%% Point Correspondence
match_iterations = 10;
% Closest point (ñ mt eficiente?)
P0_transformed = P0;
for k=1:match_iterations
    for i=1:n
        p0 = P0_transformed(i,:);
        mindist = 100;
        min_j = 1;
        % Min Dist correspondence
        for j=1:n1
            p1 = P1(j,:);
           dist = sqrt((p0-p1)*(p0-p1)');
        if(dist < mindist)
               min_j = j;
            mindist = dist;           
        end
        end
    P1_cor(i,:) = P1(min_j,:); %correspondences
    end
    % Match / Registration
    T = my_scanmatch(P0_transformed,P1_cor,[0 0 0]');
    P0_transformed = pc_transform(P0_transformed,T);
    plot(P0_transformed(:,1),P0_transformed(:,2),'.')
    drawnow
end









