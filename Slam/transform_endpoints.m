function tf_endpoint = transform_endpoints(s,pose)
R = [cos(pose.theta) -sin(pose.theta);sin(pose.theta) cos(pose.theta)];

t = R*s' + [pose.x;pose.y];
tf_endpoint.x = t(1);
tf_endpoint.y = t(2);
end