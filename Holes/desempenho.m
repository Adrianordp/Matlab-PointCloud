clear;close all;clc
ID = [84.36 87.13 83.31 90.13 89.48 90.22 85.97 88.28 87.99 90.15 88.89]';
Durs = [197 215 187 215 330 335 215 215 330 333 254]';
Vels = [10 8 10 8 4 4 8 8 4 4 6]';
H = [30 30 25 30 25 30 25 30 25 25 27.5]';

n = length(ID);
Y = ID;
% X = Durs';
X = [ones(n,1) Vels Vels.*Vels H H.*H];

Beta = inv((X'*X))*X'*Y


disp('Erro médio quadrático')
Err = Y - X*Beta
Errsqr = Err'*Err