clear;close all;clc

x = 1:0.1:10;
res = 0.1;
adjust = 1/res - 1;
x_ = x/res - adjust;

x_fix = round(x_);

x_max = max(x_fix);
V = zeros(x_max,1);

for i=1:length(x_fix)
    V(x_fix(i)) = x_fix(i);
end

V