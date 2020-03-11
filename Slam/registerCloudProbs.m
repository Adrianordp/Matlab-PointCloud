function updated_map = registerCloudProbs(map,cloud_xy,pose)
updated_map = map;
n = length(cloud_xy);
% cloud_tfed = (tf.R*cloud_xy' + tf.T)'; %transforma nuvem para coordenada grid
% cloud_xy(:,2) = -cloud_xy(:,2)
cloud_tfed = cloud_xy + [map.tfx map.tfy]; %transforma nuvem para coordenada grid
ix0 = (pose(1) + map.tfx)/ map.resolution ;
iy0 = (pose(2) + map.tfy) /map.resolution ;
for i=1:n
    
    % iy,ix -> ix (mapa ) ~ y (continuo) / iy (mapa) ~ x (continuo). PODE
    % MELHORAR ?
    ix = fix( cloud_tfed(i,1)/map.resolution );
%     iy = fix( cloud_tfed(i,2)/map.resolution );
    iy = fix( cloud_tfed(i,2)/map.resolution );
    if(ix > map.size-1 || iy > map.size-1  || ix < 0 || iy < 0)
    else 
        
        %usar modelo de probabilidade

    updated_map.grid = breseham(updated_map.grid,ix0,iy0,ix,iy);  %free 
    updated_map.grid(ix,iy) = 100; %chamar occ model

    end
end

end


function m = breseham(matrix, x0,y0, x1,y1)
    if abs(y1 - y0) < abs(x1 - x0)
        if x0 > x1
           m  = plotLineLow(matrix,x1, y1, x0, y0);
        else
            m = plotLineLow(matrix,x0, y0, x1, y1);
        end
    else
        if y0 > y1
           m =  plotLineHigh(matrix,x1, y1, x0, y0);
        else
           m =  plotLineHigh(matrix,x0, y0, x1, y1);
        end
    end
end


function m = plotLineLow(matrix,x0,y0, x1,y1)
m = matrix;
    dx = x1 - x0;
    dy = y1 - y0;
    yi = 1;
    if dy < 0
        yi = -1;
        dy = -dy;
    end
    D = 2*dy - dx;
    y = y0;

    for x=x0:1:x1
%         plot(x, y)
        m(x,y) = 0; %chamar free model
        if D > 0
               y = y + yi;
               D = D - 2*dx;
        end
        D = D + 2*dy;
    end
end


function m = plotLineHigh(matrix,x0,y0, x1,y1)
m = matrix;
    dx = x1 - x0;
    dy = y1 - y0;
    xi = 1;
    if dx < 0
        xi = -1;
        dx = -dx;
    end
    D = 2*dx - dy;
    x = x0;

    for y=y0:1:y1
%         plot(x, y,'.')
    m(x,y) = 0; %chamar free model
        if D > 0
               x = x + xi;
               D = D - 2*dy;
        end
        D = D + 2*dx;
    end
end
    
        

% function update = inv_sensor_model()
% 
% 
% end