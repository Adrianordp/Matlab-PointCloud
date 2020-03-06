clear;close all;clc
% rosinit


%% Map Specs
map.resolution = 0.1;
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

pose.x = 0;
pose.y = 0;
pose.theta = 0;

% endpois

while(1)
  cloud = receive(lidar_read);
  cloud_xyz = cloud.readXYZ;
  cloud_xy = cloud_xyz(:,1:2);
  n = length(cloud_xy);
  
  if(~first_cloud)
      first_cloud = true;
      map = registerCloud(map,cloud_xy);

      plot(cloud_xy(:,1),cloud_xy(:,2),'.');
      plotmatrix(map.grid)
      
      mapaccess(map,cloud_xy(1,1),cloud_xy(1,2))
      
  else 
      de = 0;
      for i=1:n
      endpoint = cloud_xy(i,:); %sensor reading
      tf_endp = transform_endpoints(endpoint,pose);
      MG = mapgradient(map,tf_endp.x,tf_endp.y);
      JAC = model_deriv(endpoint,pose);
      H = (MG*JAC)'*(MG*JAC);
      if(det(H) == 0)
          
      else
          MG
          JAC
          (1-mapaccess(map,tf_endp.x,tf_endp.y))
          de = de + inv(H)*(MG*JAC)'*(1-mapaccess(map,tf_endp.x,tf_endp.y))
          
      end
      
      end
      
      
%       plotmatrix(map.grid)
%       drawnow
      
      
  end
  
  
  
  
  
   
  
    
end

function jac = model_deriv(endpoint,pose)

jac = [1 0 -sin(pose.theta)*endpoint(1)-cos(pose.theta)*endpoint(2);
       0 1 cos(pose.theta)*endpoint(1)-sin(pose.theta)*endpoint(2)];

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




