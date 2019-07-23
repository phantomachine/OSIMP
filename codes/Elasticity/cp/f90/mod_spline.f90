module mod_spline
      
implicit none

contains

    
function ridiculous(x,y,w) result(z)

real(8), intent(in) :: x, y, w
real(8) z

z  = (w**2)*x + cos(y)

end function ridiculous


function spbas(r,knots) result(T)
! evaluate basis function of cubic spline on r+2 knots
! Feb 4 2010 Takeki Sunakawa

    integer, intent(in) :: r
    real(8), intent(in) :: knots(r+2)
    real(8) T(r,r)
    real(8) dt0, dt1
    integer i, j, k
    
    
    do i = 1,r
        do j = 1,r
            T(i,j) = 0.0d0
        end do
    end do

    do i = 1,r
    
        k = i+1 ! index on knots
        
        dt0 = knots(k)-knots(k-1)
        dt1 = knots(k+1)-knots(k)
        
        if (i==1) then
            
            T(i,i) = 2.0d0*(dt0+dt1)-(dt1-dt0**2/dt1)
            T(i,i+1) = dt0+dt0**2/dt1
            ! BUG: 190605 fixed
            ! T(i,i+1) = dt1+dt0**2/dt1 ! found an error on 11/21/14    

        elseif (i==r) then
            
            T(i,i-1) = dt1+dt1**2/dt0
            T(i,i) = 2.0d0*(dt0+dt1)-(dt0-dt1**2/dt0)
            
        else
            
            ! found an error on 11/21/14 - need to check
            T(i,i-1) = dt0
            T(i,i) = 2.0d0*(dt0+dt1)
            T(i,i+1) = dt1 
        
        end if
        
    end do

end function spbas


function spfit(invT,fval,r,knots) result(c)

    integer, intent(in) :: r
    real(8), intent(in) :: invT(r,r), fval(r+2), knots(r+2)
    real(8) c(4,r+1) 
    real(8) f0, f1, f2, dt0, dt1, df0, df1, a0, a1, a2, b0, b1, b2
    real(8) F(r), s(r+2), SVEC(2)
    integer i, j, k    

    do i = 1,r
        F(i) = 0.0d0
    end do

    do i = 1,r

        k = i+1 ! index on knots

        f0 = fval(k-1)
        f1 = fval(k)
        f2 = fval(k+1)
        dt0 = knots(k)-knots(k-1)
        dt1 = knots(k+1)-knots(k)
        df0 = (f1-f0)/dt0
        df1 = (f2-f1)/dt1

        if (i==1) then

            F(i) = 3.0d0*(dt1*df0+dt0*df1) - 2.0d0*(dt1*df0-dt0**2/dt1*df1)

            a1 = 1.0d0-(dt0/dt1)**2
            a2 = (dt0/dt1)**2
            a0 = 2.0d0*(df0-(dt0/dt1)**2*df1)  

        elseif (i==r) then

            F(i) = 3.0d0*(dt1*df0+dt0*df1) - 2.0d0*(dt0*df1-dt1**2/dt0*df0)

            b1 = 1.0d0-(dt1/dt0)**2
            b2 = (dt1/dt0)**2
            b0 = 2.0d0*(df1-(dt1/dt0)**2*df0)

        else

            F(i) = 3.0d0*(dt1*df0+dt0*df1)

        end if

    end do    

    s(2:r+1) = matmul(invT,F)

    ! not-a-knot
    s(1) = a0 - a1*s(2) + a2*s(3)
    s(r+2) = b0 - b1*s(r+1) + b2*s(r)

    do i = 1,r+1

        k = i+1 ! index on knots

        f0 = fval(k-1)
        f1 = fval(k)
        dt0 = knots(k)-knots(k-1)
        df0 = (f1-f0)/dt0

        c(1,i) = f0
        c(2,i) = s(i)
        c(3,i) = 3.0d0*df0/dt0 - 2.0d0*s(i)/dt0 - s(i+1)/dt0
        c(4,i) = -2.0d0*df0/(dt0**2) + s(i)/(dt0**2) + s(i+1)/(dt0**2)

    end do

end function spfit


subroutine speva(c,x,r,knots,f,df,d2f)

    integer, intent(in) :: r
    real(8), intent(in) :: c(4,r+1), x, knots(r+2)
    real(8), intent(out) :: f, df, d2f
    real(8) t !, x1, x2, f1, f2, s
    integer i, j, k, n

    k = 0
    do j = 1,r+1
        if (knots(j)>x) exit
        k = k+1
    end do

    k = min(k,r+1)
    k = max(k,1)

    t   = knots(k) ! leftknot
    f   = c(1,k) + c(2,k)*(x-t) + c(3,k)*(x-t)**2 + c(4,k)*(x-t)**3
    df  = c(2,k) + 2.0d0*c(3,k)*(x-t) + 3.0d0*c(4,k)*(x-t)**2
    d2f = 2.0d0*c(3,k) + 6.0d0*c(4,k)*(x-t)        

end subroutine speva


function spfit2(invTx,invTy,fmat,rx,ry,xknots,yknots) result(cmat2)

    integer, intent(in) :: rx, ry
    real(8), intent(in) :: invTx(rx,rx), invTy(ry,ry), fmat(rx+2,ry+2), xknots(rx+2), yknots(ry+2)
    real(8) cmat(4,(ry+1),4*(rx+1)), cmat2(16,rx+1,ry+1), dmat(4*(rx+1),ry+2)
    real(8) c(4,ry+1), d(4,rx+1), fvalx(ry+2), fvaly(rx+2)
    integer ix, iy


    do iy = 1,ry+2

        fvaly = fmat(:,iy)
        d = spfit(invTx,fvaly,rx,xknots)
        dmat(:,iy) = reshape(d, (/4*(rx+1)/))

    end do


    do ix = 1,4*(rx+1)

        fvalx = reshape(dmat(ix,:), (/ry+2/))
        c = spfit(invTy,fvalx,ry,yknots) ! 4x(ry+1)
        cmat(:,:,ix) = reshape(c, (/4,ry+1/))

    end do


    do iy = 1,ry+1

        do ix = 1,rx+1

            cmat2(:,ix,iy) = reshape(cmat(:,iy,4*(ix-1)+1:4*ix), (/16/))

        end do

    end do

end function spfit2


subroutine speva2(cmat,x,y,rx,ry,xknots,yknots,f,dfx,dfy)

    integer, intent(in) :: rx, ry
    real(8), intent(in) :: cmat(16,rx+1,ry+1), x, y, xknots(rx+2), yknots(ry+2)
    real(8), intent(out) :: f, dfx, dfy
    real(8) tx, ty, c(16)
    integer jx, jy, kx, ky


    kx = 0

    do jx = 1,rx+1
        if (xknots(jx)>x) exit
        kx = kx+1
    end do

    kx = min(kx,rx+1)
    kx = max(kx,1)        

    tx = xknots(kx) ! leftknot

    ky = 0

    do jy = 1,ry+1
        if (yknots(jy)>y) exit
        ky = ky+1
    end do

    ky = min(ky,ry+1)
    ky = max(ky,1)

    ty = yknots(ky) ! leftknot

    c = reshape(cmat(:,kx,ky), (/16/))

    f = c(1) + c(2)*(y-ty) + c(3)*(y-ty)**2 + c(4)*(y-ty)**3 &
        + (x-tx)*(c(5) + c(6)*(y-ty) + c(7)*(y-ty)**2 + c(8)*(y-ty)**3) &
        + (x-tx)**2*(c(9) + c(10)*(y-ty) + c(11)*(y-ty)**2 + c(12)*(y-ty)**3) &
        + (x-tx)**3*(c(13) + c(14)*(y-ty) + c(15)*(y-ty)**2 + c(16)*(y-ty)**3)
    dfx = c(5) + c(6)*(y-ty) + c(7)*(y-ty)**2 + c(8)*(y-ty)**3 &
        + 2.0d0*(x-tx)*(c(9) + c(10)*(y-ty) + c(11)*(y-ty)**2 + c(12)*(y-ty)**3) &
        + 3.0d0*(x-tx)**2*(c(13) + c(14)*(y-ty) + c(15)*(y-ty)**2 + c(16)*(y-ty)**3)
    dfy = c(2) + c(6)*(x-tx) + c(10)*(x-tx)**2 + c(14)*(x-tx)**3 &
        + 2.0d0*(y-ty)*(c(3) + c(7)*(x-tx) + c(11)*(x-tx)**2 + c(15)*(x-tx)**3) &
        + 3.0d0*(y-ty)**2*(c(4) + c(8)*(x-tx) + c(12)*(x-tx)**2 + c(16)*(x-tx)**3)
!    d2fx = 2.0d0*(c(9) + c(10)*(y-ty) + c(11)*(y-ty)**2 + c(12)*(y-ty)**3) &
!        + 6.0d0*(x-tx)*(c(13) + c(14)*(y-ty) + c(15)*(y-ty)**2 + c(16)*(y-ty)**3)
!    d2fy = 2.0d0*(c(3) + c(7)*(x-tx) + c(11)*(x-tx)**2 + c(15)*(x-tx)**3) &
!        + 6.0d0*(y-ty)*(c(4) + c(8)*(x-tx) + c(12)*(x-tx)**2 + c(16)*(x-tx)**3)
!    dfxy = c(6) + 2.0d0*c(10)*(x-tx) + 3.0d0*c(14)*(x-tx)**2 &
!        + 2.0d0*(y-ty)*(c(7) + 2.0d0*c(11)*(x-tx) + 3.0d0*c(15)*(x-tx)**2) &
!        + 3.0d0*(y-ty)**2*(c(8) + 2.0d0*c(12)*(x-tx) + 3.0d0*c(16)*(x-tx)**2)

end subroutine speva2


function spfit3(invTx,invTy,invTz,fmat,rx,ry,rz,xknots,yknots,zknots) result(cmat3)

integer, intent(in) :: rx, ry, rz
real(8), intent(in) :: invTx(rx,rx), invTy(ry,ry), invTz(rz,rz), fmat(rx+2,ry+2,rz+2), xknots(rx+2), yknots(ry+2), zknots(rz+2)
real(8) cmat(4,(ry+1),4*(rx+1)), cmat2(4,ry+1,4*(rx+1),4*(rz+1)), cmat3(64,rx+1,ry+1,rz+1)
real(8) dmat(4*(rx+1),ry+2), emat(4*(rz+1),rx+2,ry+2), emat2(rx+2,ry+2)
real(8) c(4,ry+1), d(4,rx+1), e(4,(rz+1)), fvalx(rx+2), fvaly(rx+2), fvalz(rz+2)
integer iy, ix, iz, index


! fit (rx+2)x(ry+2) z-splines on f(x,y,:)
do ix = 1,rx+2
    
    do iy = 1,ry+2
        
        fvalz = reshape(fmat(ix,iy,:),(/rz+2/))
        e = spfit(invTz,fvalz,rz,zknots) ! 4x(rz+1)
        emat(:,ix,iy) = reshape(e,(/4*(rz+1)/))
        
    end do
    
end do


do iz = 1,4*(rz+1)
    
    emat2 = reshape(emat(iz,:,:),(/rx+2,ry+2/))

    do iy = 1,ry+2

        fvalx = emat2(:,iy)
        ! fit univariate x-spline
        d = spfit(invTx,fvalx,rx,xknots)
        dmat(:,iy) = reshape(d,(/4*(rx+1)/));

    end do

    do ix = 1,4*(rx+1)

        fvaly = reshape(dmat(ix,:), (/ry+2/))
        ! fit univariate y-spline,
        c = spfit(invTy,fvaly,ry,yknots) ! 4x(ry+1)
        cmat(:,:,ix) = c !reshape(c,(/4*(rx+1)/)) 4x(ry+1)x(rx+1)

    end do
    
    cmat2(:,:,:,iz) = cmat
    
end do


! reshaping cmat2, 64x(rx+1)x(ry+1)x(rz+1)
do iz = 1,rz+1

    do iy = 1,ry+1

        do ix = 1,rx+1

            cmat3(:,ix,iy,iz) = reshape(cmat2(:,iy,4*(ix-1)+1:4*ix,4*(iz-1)+1:4*iz),(/64/))

        end do

    end do
    
end do


end function spfit3


subroutine speva3(cmat,x,y,z,rx,ry,rz,xknots,yknots,zknots,f,dfx,dfy,dfz)

integer, intent(in) :: rx, ry, rz
real(8), intent(in) :: cmat(64,rx+1,ry+1,rz+1), x, y, z, xknots(rx+2), yknots(ry+2), zknots(rz+2)
real(8), intent(out) :: f, dfx, dfy, dfz
real(8) tx, ty, tz, c(64)
integer jx, jy, jz, kx, ky, kz


kx = 0

do jx = 1,rx+1
    if (xknots(jx)>x) exit
    kx = kx+1
end do

kx = min(kx,rx+1)
kx = max(kx,1)        

ky = 0

do jy = 1,ry+1
    if (yknots(jy)>y) exit
    ky = ky+1
end do

ky = min(ky,ry+1)
ky = max(ky,1)

kz = 0

do jz = 1,rz+1
    if (zknots(jz)>z) exit
    kz = kz+1
end do

kz = min(kz,rz+1)
kz = max(kz,1)

c = reshape(cmat(:,kx,ky,kz), (/64/))
tx = xknots(kx) ! leftknot
ty = yknots(ky) ! leftknot
tz = zknots(kz) ! leftknot

f =               c(1) +  c(2)*(y-ty)  + c(3)*(y-ty)**2  + c(4)*(y-ty)**3 &
    + (x-tx)*(     c(5) +  c(6)*(y-ty)  + c(7)*(y-ty)**2  + c(8)*(y-ty)**3) &
    + (x-tx)**2*(   c(9) +  c(10)*(y-ty) + c(11)*(y-ty)**2 + c(12)*(y-ty)**3) &
    + (x-tx)**3*(   c(13) + c(14)*(y-ty) + c(15)*(y-ty)**2 + c(16)*(y-ty)**3) &
    + (z-tz)*(     c(17) + c(18)*(y-ty) + c(19)*(y-ty)**2 + c(20)*(y-ty)**3 &
    + (x-tx)*(     c(21) + c(22)*(y-ty) + c(23)*(y-ty)**2 + c(24)*(y-ty)**3) &
    + (x-tx)**2*(   c(25) + c(26)*(y-ty) + c(27)*(y-ty)**2 + c(28)*(y-ty)**3) &
    + (x-tx)**3*(   c(29) + c(30)*(y-ty) + c(31)*(y-ty)**2 + c(32)*(y-ty)**3)) &
    + (z-tz)**2*(   c(33) + c(34)*(y-ty) + c(35)*(y-ty)**2 + c(36)*(y-ty)**3 &
    + (x-tx)*(     c(37) + c(38)*(y-ty) + c(39)*(y-ty)**2 + c(40)*(y-ty)**3) &
    + (x-tx)**2*(   c(41) + c(42)*(y-ty) + c(43)*(y-ty)**2 + c(44)*(y-ty)**3) &
    + (x-tx)**3*(   c(45) + c(46)*(y-ty) + c(47)*(y-ty)**2 + c(48)*(y-ty)**3)) &
    + (z-tz)**3*(   c(49) + c(50)*(y-ty) + c(51)*(y-ty)**2 + c(52)*(y-ty)**3 &
    + (x-tx)*(     c(53) + c(54)*(y-ty) + c(55)*(y-ty)**2 + c(56)*(y-ty)**3) &
    + (x-tx)**2*(   c(57) + c(58)*(y-ty) + c(59)*(y-ty)**2 + c(60)*(y-ty)**3) &
    + (x-tx)**3*(   c(61) + c(62)*(y-ty) + c(63)*(y-ty)**2 + c(64)*(y-ty)**3))

dfx =             c(5)  + c(6)*(y-ty)  + c(7)*(y-ty)**2  + c(8)*(y-ty)**3 &
    + 2.0d0*(x-tx)*(   c(9)  + c(10)*(y-ty) + c(11)*(y-ty)**2 + c(12)*(y-ty)**3) &
    + 3.0d0*(x-tx)**2*( c(13) + c(14)*(y-ty) + c(15)*(y-ty)**2 + c(16)*(y-ty)**3) &
    + (z-tz)*((    c(21) + c(22)*(y-ty) + c(23)*(y-ty)**2 + c(24)*(y-ty)**3) &
    + 2.0d0*(x-tx)*(   c(25) + c(26)*(y-ty) + c(27)*(y-ty)**2 + c(28)*(y-ty)**3) &
    + 3.0d0*(x-tx)**2*( c(29) + c(30)*(y-ty) + c(31)*(y-ty)**2 + c(32)*(y-ty)**3)) &
    + (z-tz)**2*((  c(37) + c(38)*(y-ty) + c(39)*(y-ty)**2 + c(40)*(y-ty)**3) &
    + 2.0d0*(x-tx)*(   c(41) + c(42)*(y-ty) + c(43)*(y-ty)**2 + c(44)*(y-ty)**3) &
    + 3.0d0*(x-tx)**2*( c(45) + c(46)*(y-ty) + c(47)*(y-ty)**2 + c(48)*(y-ty)**3)) &
    + (z-tz)**3*((  c(53) + c(54)*(y-ty) + c(55)*(y-ty)**2 + c(56)*(y-ty)**3) &
    + 2.0d0*(x-tx)*(   c(57) + c(58)*(y-ty) + c(59)*(y-ty)**2 + c(60)*(y-ty)**3) &
    + 3.0d0*(x-tx)**2*( c(61) + c(62)*(y-ty) + c(63)*(y-ty)**2 + c(64)*(y-ty)**3))
 
dfy =             c(2)  + c(6)*(x-tx)  + c(10)*(x-tx)**2 + c(14)*(x-tx)**3 &
    + 2.0d0*(y-ty)*(   c(3)  + c(7)*(x-tx)  + c(11)*(x-tx)**2 + c(15)*(x-tx)**3) &
    + 3.0d0*(y-ty)**2*( c(4)  + c(8)*(x-tx)  + c(12)*(x-tx)**2 + c(16)*(x-tx)**3) &
    + (z-tz)*(     c(18) + c(22)*(x-tx) + c(26)*(x-tx)**2 + c(30)*(x-tx)**3 &
    + 2.0d0*(y-ty)*(   c(19) + c(23)*(x-tx) + c(27)*(x-tx)**2 + c(31)*(x-tx)**3) &
    + 3.0d0*(y-ty)**2*( c(20) + c(24)*(x-tx) + c(28)*(x-tx)**2 + c(32)*(x-tx)**3)) &
    + (z-tz)**2*(   c(34) + c(38)*(x-tx) + c(42)*(x-tx)**2 + c(46)*(x-tx)**3 &
    + 2.0d0*(y-ty)*(   c(35) + c(39)*(x-tx) + c(43)*(x-tx)**2 + c(47)*(x-tx)**3) &
    + 3.0d0*(y-ty)**2*( c(36) + c(40)*(x-tx) + c(44)*(x-tx)**2 + c(48)*(x-tx)**3)) &
    + (z-tz)**3*(   c(50) + c(54)*(x-tx) + c(58)*(x-tx)**2 + c(62)*(x-tx)**3 &
    + 2.0d0*(y-ty)*(   c(19) + c(55)*(x-tx) + c(59)*(x-tx)**2 + c(63)*(x-tx)**3) &
    + 3.0d0*(y-ty)**2*( c(20) + c(56)*(x-tx) + c(60)*(x-tx)**2 + c(64)*(x-tx)**3))

dfz =             c(17) + c(18)*(y-ty) + c(19)*(y-ty)**2 + c(20)*(y-ty)**3  &
    + (x-tx)*(     c(21) + c(22)*(y-ty) + c(23)*(y-ty)**2 + c(24)*(y-ty)**3) &
    + (x-tx)**2*(   c(25) + c(26)*(y-ty) + c(27)*(y-ty)**2 + c(28)*(y-ty)**3) &
    + (x-tx)**3*(   c(29) + c(30)*(y-ty) + c(31)*(y-ty)**2 + c(32)*(y-ty)**3) &
    + 2.0d0*(z-tz)*(   c(33) + c(34)*(y-ty) + c(35)*(y-ty)**2 + c(36)*(y-ty)**3  &
    + (x-tx)*(     c(37) + c(38)*(y-ty) + c(39)*(y-ty)**2 + c(40)*(y-ty)**3) &
    + (x-tx)**2*(   c(41) + c(42)*(y-ty) + c(43)*(y-ty)**2 + c(44)*(y-ty)**3) &
    + (x-tx)**3*(   c(45) + c(46)*(y-ty) + c(47)*(y-ty)**2 + c(48)*(y-ty)**3)) &
    + 3.0d0*(z-tz)**2*( c(49) + c(50)*(y-ty) + c(51)*(y-ty)**2 + c(52)*(y-ty)**3  &
    + (x-tx)*(     c(53) + c(54)*(y-ty) + c(55)*(y-ty)**2 + c(56)*(y-ty)**3) &
    + (x-tx)**2*(   c(57) + c(58)*(y-ty) + c(59)*(y-ty)**2 + c(60)*(y-ty)**3) &
    + (x-tx)**3*(   c(61) + c(62)*(y-ty) + c(63)*(y-ty)**2 + c(64)*(y-ty)**3))

end subroutine speva3


end module mod_spline
