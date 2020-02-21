% FORÇA BRUTISSIMA O(N^2)
% Alternativas -> 
% * kDSearch inteligente 
% * Varrer NP, e não GRID O(N)
function M = fill_grid(pointcloud,xgrid,ygrid)
grid_width = length(xgrid);
grid_height = length(ygrid);
n_pts = length(pointcloud)
M = zeros(grid_height-1,grid_width-1);

%iguais ?
%rect_area = (xgrid(2)-xgrid(1)) * (ygrid(2)-ygrid(1));


% força menos bruta .. 
found = 0
for k=1:n_pts
    k
for i=1:grid_height-1
    
    for j=1:grid_width-1
        
        if(found)
            found = 0;
           break
        end
        
        rect = get_rect(xgrid,j,ygrid,i);
        rect_area = get_rect_area(rect);
        
        if isInside(pointcloud(k,:),rect,rect_area)
                M(i,j) = 1;
%                 fprintf('Achado em M(%d,%d) \n', i,j);
                found = 1;
        end  
    end
end
    
    
end



end

    function rect = get_rect(xgrid,j,ygrid,i)
    rect = zeros(4,2);
     
     rect(1,:) = [xgrid(j) ygrid(i)];
     rect(2,:) = [xgrid(j+1) ygrid(i)];
     rect(3,:) = [xgrid(j+1) ygrid(i+1)];
     rect(4,:) = [xgrid(j) ygrid(i+1)];
     
    %retorna matrix 4x2 c pontos do retângulo
    end


%P -> ABCD ?
% rect = 4x2
function res = isInside(pt,rect,rect_area) 



APB = [rect(1,:);pt;rect(2,:)];
BPC = [rect(2,:);pt;rect(3,:)];
CPD = [rect(3,:);pt;rect(4,:)];
DPA = [rect(4,:);pt;rect(1,:)];

areas =  [tri_area(APB);
        tri_area(BPC);
        tri_area(CPD);
        tri_area(DPA);];
    
total_area = sum(areas);
if(total_area > rect_area) %ta fora 
    res = 0;
else
    res = 1;
end


    
end

%tri ABC -> 3x2 
function area = tri_area(tri)
AC = [tri(2,1)-tri(1,1) tri(2,2) - tri(1,2) 0];
AB = [tri(3,1)-tri(1,1) tri(3,2) - tri(1,2) 0];
area = norm(cross(AC,AB))/2;

end

%rect 4x2
function area = get_rect_area(rect)
area = (rect(2,1)-rect(1,1))*(rect(3,2)-rect(2,2));
end

