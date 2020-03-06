function dm = mapgradient(map,x,y)
    x = x + map.tfx;
    y = y + map.tfx;

% iy,ix -> ix (mapa ) ~ y (continuo) / iy (mapa) ~ x (continuo)
    ix = fix( x/map.resolution );
    iy = fix( y/map.resolution );
     if(ix > map.size || iy > map.size  || ix < 0 || iy < 0)
        dm = [0 0];
        return;
     end
    
    p00 = map.grid(ix,iy);
    p01 = map.grid(ix+1,iy);
    p10 = map.grid(ix,iy+1);
    p11 = map.grid(ix+1,iy+1);
    
    
    
    dx_x0 = x - map.resolution*ix;
    dy_y0 = y - map.resolution*iy;
    
    dx1_x =map.resolution*(ix+1) - x;
    dy1_y = map.resolution*(iy+1) - y;
    
    dmx = (dy_y0/map.resolution)*(p01 - p00) + (dy1_y/map.resolution)*(p10-p11);
    dmy = (dx_x0/map.resolution)*(p01 - p00) + (dx1_x/map.resolution)*(p10-p11);
    
    dm = [dmx dmy];
    
    


end