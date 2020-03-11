function tf_endpoint = transform_endpoints(s,pose)
x = pose(1);
y = pose(2);
theta = pose(3);
R = [cos(theta) -sin(theta);sin(theta) cos(theta)];

t = R*s' + [x;y];
tf_endpoint.x = t(1);
tf_endpoint.y = t(2);
end