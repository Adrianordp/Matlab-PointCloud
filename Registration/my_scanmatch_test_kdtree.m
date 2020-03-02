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
%% myICP
match_iterations = 10;
% Closest point (ñ mt eficiente?)

%% KD TREE for closest
P0_transformed = P0;
KD = KDTreeSearcher(P1,'BucketSize',10);
r = 1; %distance
for k=1:match_iterations
    j = 1;
    for i=1:n %em P0
        p0 = P0_transformed(i,:);     
        IDX_ = rangesearch(KD,p0,r); 
        IDX = IDX_{1};
        if(isempty(IDX))
            
        else
        P0_corr(j,:) = p0;
        P1_corr(j,:) = P1(IDX(1),:);
        j = j+1;
        end
    end
    plot(P0_corr(:,1),P0_corr(:,2),'o')
    plot(P1_corr(:,1),P1_corr(:,2),'o')
    % Match / Registration
    T = my_scanmatch(P0_corr,P1_corr,[0 0 0]',5);
    P0_transformed = pc_transform(P0_transformed,T);
    plot(P0_transformed(:,1),P0_transformed(:,2),'.')
    drawnow
end

%% Assumir data assossiation fixa
P1 = pc_transform(P0,[5 5 0]) %test
%Plots
figure
plot(P0(:,1),P0(:,2),'b.') % Original
hold on
plot(P1(:,1),P1(:,2),'r.') % Final
drawnow
P0_transformed = P0;
T = my_scanmatch(P0_transformed,P1,[0 0 0]',10);    
P0_transformed = pc_transform(P0_transformed,T);
   plot(P0_transformed(:,1),P0_transformed(:,2),'.')
   drawnow











