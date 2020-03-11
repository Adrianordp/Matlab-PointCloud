clear;close all;clc
% rosinit

%% ROS STUFF

tftree = rostf
transfm = rosmessage('geometry_msgs/TransformStamped');
transfm.ChildFrameId = 'pose';
transfm.Header.FrameId = 'map';
transfm.Transform.Translation.X = 0;
transfm.Transform.Rotation.W = 1;
transfm.Transform.Rotation.X = 0;
transfm.Transform.Rotation.Y = 0;
transfm.Transform.Rotation.Z = 0;



%% Map Specs
map.resolution = 0.05;
map.size = 200; %cells
disp('Map side length') ;
disp(map.resolution*map.size);
map.startx = 0.5; % 0 < x < 1
map.starty = 0.5; % 0 < y < 1

%translada cloud para origem da grid para acessar em forma de matrix
map.tfx = map.startx*map.resolution*map.size;
map.tfy = map.starty*map.resolution*map.size;
map.grid = zeros(map.size);



%% Loop
lidar_read = rossubscriber('/cloud');
first_cloud = false;

pose = [0 0 0]';
update_pose = pose;

% endpois

while(1)
  cloud = receive(lidar_read);
  cloud_time = cloud.Header.Stamp;
  cloud_xyz = cloud.readXYZ;
  cloud_xy = cloud_xyz(:,1:2);
  n = length(cloud_xy);
  
  if(~first_cloud)
      first_cloud = true;
      map = registerCloud(map,cloud_xy);

      plot(cloud_xy(:,1),cloud_xy(:,2),'.');
      plotmatrix(map.grid)
      
%       mapaccess(map,cloud_xy(1,1),cloud_xy(1,2))
      
  else 
      de = 0;
      H = zeros(3);
      dtr = zeros(3,1);
      for i=1:n
      endpoint = cloud_xy(i,:); %sensor reading
      endpoint_tf = transform_endpoints(endpoint,pose);
      %as duas funções abaixo podem ser uma só
      funval = 1 - mapaccess(map,endpoint_tf);
      dm = mapgradient(map,endpoint_tf);
      jac = model_deriv(endpoint,pose);
      
      dtr = dtr + (dm*jac)'*funval; 
      H = H + (dm*jac)'*(dm*jac);
      end
      
      if(H(1,1) ~= 0 && H(2,2) ~= 0)
         searchdir = inv(H)*dtr
         pose = pose + searchdir
         q = eul2quat([pose(3) 0 0]);
         transfm.Transform.Translation.X = pose(1);
         transfm.Transform.Translation.Y = pose(2);
         transfm.Transform.Rotation.W = q(1);
         transfm.Transform.Rotation.X = q(2);
         transfm.Transform.Rotation.Y = q(3);
         transfm.Transform.Rotation.Z = q(4);
         transfm.Header.Stamp = cloud_time;
         tftree.sendTransform(transfm);
         
         updist = sqrt ((update_pose(1:2) - pose(1:2))'*(update_pose(1:2) - pose(1:2)))
         angdist = abs(update_pose(3) - pose(3))
         if(updist > 0.3 || angdist > 0.06)
             cloud_t = transform_cloud(cloud_xy,pose);
             figure(1)
             hold on
             plotcloud(cloud_t);
             map = registerCloud(map,cloud_t);
             update_pose = pose;
         end
         
      end
      
      
%       plotmatrix(map.grid)
%       drawnow
      
      
  end
  
  
  
  
  
   
  
    
end

function jac = model_deriv(endpoint,pose)
theta = pose(3);
jac = [1 0 -sin(theta)*endpoint(1)-cos(theta)*endpoint(2);
       0 1 cos(theta)*endpoint(1)-sin(theta)*endpoint(2)];

end





% function tf = pose2tf(pose)
% tf.R = [cos(pose.phi) -sin(pose.phi)
%     sin(pose.phi) cos(pose.phi)];
% 
% tf.T = [pose.x;pose.y];
% 
% end

%aproximate


%   for(unsigned int i=0;i<size;++i){
% 	if(map_g->data[i] > LOW_LIMIT){ //verifica se ponto eh valido
% 	points++;
% 	column = i%height;
% 	row = i/width;
% //	ROS_INFO("i = %d, row = %d, column = %d",i,row,column);
% 	y = row*resolution;
% 	x = column*resolution;
% 	x = x + p0.position.x;
% 	y = y + p0.position.y;
% //	ROS_INFO("x = %f, y = %f",x,y);
% 	cloud->push_back(pcl::PointXYZ (x,y,0));
% 	  }
%   }




