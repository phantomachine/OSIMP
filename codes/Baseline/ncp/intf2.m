function z = intf2(xgrid,ygrid,zmat,x,y)
% 2-dim linear interpolation
nx = size(xgrid,1);
ny = size(ygrid,1);

% for x
ix = 1;
if (x<=xgrid(1))
    ix = 1;
elseif (x>=xgrid(nx))
    ix = nx-1;
else
    jx = 2;
    while (jx<=nx)

        if (x<xgrid(jx))
            ix = jx-1;
            jx=nx;
        end

    jx=jx+1;
    end
end

ix = min(max(1,ix),nx-1);

% for y
iy = 1;
if (y<=ygrid(1))
    iy = 1;
elseif (y>=ygrid(ny))
    iy = ny-1;
else
    jy = 2;
    while (jy<=ny)

        if (y<ygrid(jy))
            iy=jy-1;
            jy=ny;
        end

    jy=jy+1;
    end
end

iy = min(max(1,iy),ny-1);

etax = (x-xgrid(ix))/(xgrid(ix+1)-xgrid(ix));
etay = (y-ygrid(iy))/(ygrid(iy+1)-ygrid(iy));
z1 = (1-etay)*zmat(ix,iy)   + etay*zmat(ix,iy+1);
z2 = (1-etay)*zmat(ix+1,iy) + etay*zmat(ix+1,iy+1);
z  = (1-etax)*z1 + etax*z2;
