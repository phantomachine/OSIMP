function w = intf3(xgrid,ygrid,zgrid,wmat,x,y,z)
% 3-dim linear interpolation
nx = size(xgrid,1);
ny = size(ygrid,1);
nz = size(zgrid,1);

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

% for z
iz = 1;
if (z<=zgrid(1))
    iz = 1;
elseif (z>=zgrid(nz))
    iz = nz-1;
else
    jz = 2;
    while (jz<=nz)

        if (z<zgrid(jz))
            iz=jz-1;
            jz=nz;
        end

    jz=jz+1;
    end
end

iz = min(max(1,iz),nz-1);

etax = (x-xgrid(ix))/(xgrid(ix+1)-xgrid(ix));
etay = (y-ygrid(iy))/(ygrid(iy+1)-ygrid(iy));
etaz = (z-zgrid(iz))/(zgrid(iz+1)-zgrid(iz));
w11  = (1-etay)*wmat(ix,iy,iz)   + etay*wmat(ix,iy+1,iz);
w21  = (1-etay)*wmat(ix+1,iy,iz) + etay*wmat(ix+1,iy+1,iz);
w12  = (1-etay)*wmat(ix,iy,iz+1)   + etay*wmat(ix,iy+1,iz+1);
w22  = (1-etay)*wmat(ix+1,iy,iz+1) + etay*wmat(ix+1,iy+1,iz+1);
w1   = (1-etax)*w11 + etax*w21;
w2   = (1-etax)*w12 + etax*w22;
w    = (1-etaz)*w1 + etaz*w2;