function R = quat2mat(q_)
q = reshape(q_, [4 1]);
s = 1/(q'*q)^2;
R = [1 - 2*s*(q(3)*q(3) + q(4)*q(4)) 2*s*(q(2)*q(3)-q(4)*q(1)) 2*s*(q(2)*q(4)+q(3)*q(1));
     2*s*(q(2)*q(3)+q(4)*q(1)) 1-2*s*(q(2)*q(2)+q(4)*q(4)) 2*s*(q(3)*q(4)-q(2)*q(1));
     2*s*(q(2)*q(4)-q(3)*q(1)) 2*s*(q(3)*q(4)+q(2)*q(1)) 1-2*s*(q(2)*q(2)+q(3)*q(3))];
end