clear;close all;clc
ID = [84.36 87.13 83.31 90.13 89.48 90.22 85.97 88.28 87.99 90.15 88.89]';
Durs = [197 215 187 215 330 335 215 215 330 333 254]';
Vels = [10 8 10 8 4 4 8 8 4 4 6]';
H = [30 30 25 30 25 30 25 30 25 25 27.5]';

n = length(ID);
Y = ID;
% X = Durs';
X = [ones(n,1) Vels];

Beta = inv((X'*X))*X'*Y


disp('Erro médio quadrático')
Err = Y - X*Beta
Errsqr = Err'*Err

%% Cost
MaxID = 100;
beta0 = Beta(1);
beta = Beta(2);
Qd = 1*eye(n); %minimiza ID
Qv = 1*eye(n); %minimiza vel

beta0_ = beta0*ones(n,1);


ID_ = ID;
v_opt = -inv(beta^2*Qd+Qv)*(beta0_'*Qd*beta)'
ID_opt = beta0 + beta*v_opt
% JJ = [];
% min_j = inf;
% for v = -1000:1:1000
%     v_opt = v*ones(n,1);
% J = (beta0_ + beta*v_opt)'*Qd*(beta0_ + beta*v_opt) + v_opt'*Qv*v_opt;
% if(J < min_j)
%     v_min = v_opt;
%     min_j = J;
% end
% JJ = [JJ;J];
% end
% plot(JJ)
