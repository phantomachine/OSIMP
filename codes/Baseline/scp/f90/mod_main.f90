module mod_main


implicit none


contains


subroutine sus_iter(p1mat,p2mat,v1mat,v2mat,x1mat,x2mat,&
p1mat0,p2mat0,v1mat0,v2mat0,x1mat0,x2mat0,y1mat,y2mat,zmat,xmat,&
w1mat,w2mat,w11mat,w12mat,w21mat,w22mat,Gu,Pu,invT1,invT2,invTx)

use mod_foc
use mod_spline
use mod_parameters

implicit none

real(8), intent(out) :: p1mat(:,:,:,:,:), p2mat(:,:,:,:,:)
real(8), intent(out) :: v1mat(:,:,:,:,:), v2mat(:,:,:,:,:)
real(8), intent(out) :: x1mat(:,:,:,:,:), x2mat(:,:,:,:,:)
real(8), intent(in) :: p1mat0(:,:,:,:,:), p2mat0(:,:,:,:,:)
real(8), intent(in) :: v1mat0(:,:,:,:,:), v2mat0(:,:,:,:,:)
real(8), intent(in) :: x1mat0(:,:,:,:,:), x2mat0(:,:,:,:,:)
real(8), intent(inout) :: y1mat(:,:,:,:,:), y2mat(:,:,:,:,:), zmat(:,:,:,:,:), xmat(:,:,:,:,:)
real(8), intent(in) :: w1mat(:,:,:,:), w2mat(:,:,:,:), w11mat(:,:,:,:), w12mat(:,:,:,:), w21mat(:,:,:,:), w22mat(:,:,:,:)
real(8), intent(in) :: Gu(:,:), Pu(:,:), invT1(:,:), invT2(:,:), invTx(:,:)

integer flag
real(8) u1p, u2p, u1, u2, xp, y1p, y2p, w1, w11, w12, w2, w21, w22
real(8), allocatable :: p1cond(:,:,:), p2cond(:,:,:), v1cond(:,:,:), v2cond(:,:,:), x1cond(:,:,:), x2cond(:,:,:)
real(8), allocatable :: cp1(:,:,:,:), cp2(:,:,:,:)
real(8), allocatable :: cv1(:,:,:,:), cv2(:,:,:,:)
real(8), allocatable :: cx1(:,:,:,:), cx2(:,:,:,:)

real(8) xx0(2), xx1(3), f0(2), f1(3), y1, y2, z, x, p1, p2, xi1, xi2, &
xi1p, xi2p, v1, v2, v1f, v2f, dfx, dfy, dfz, p1f, p2f
integer iu, ju, ku, i1, i2, ix, i, rc, iter, j, j1, j2, jx

integer, external :: mexPrintf
integer k
character(len=80) line

allocate(p1cond(n1,n2,nx), p2cond(n1,n2,nx), v1cond(n1,n2,nx), v2cond(n1,n2,nx), x1cond(n1,n2,nx), x2cond(n1,n2,nx))
allocate(cp1(64,r1+1,r2+1,rx+1), cp2(64,r1+1,r2+1,rx+1))
allocate(cv1(64,r1+1,r2+1,rx+1), cv2(64,r1+1,r2+1,rx+1))
allocate(cx1(64,r1+1,r2+1,rx+1), cx2(64,r1+1,r2+1,rx+1))

!$omp parallel do private(ku,iu,ju,u1p,u2p,u1,u2, &
!$omp p1cond,p2cond,v1cond,v2cond,x1cond,x2cond,cp1,cp2,cv1,cv2,cx1,cx2)
do ku = 1,nu

    u1p = Gu(ku,1)
    u2p = Gu(ku,2)

    do iu = 1,nu

        u1 = Gu(iu,1)
        u2 = Gu(iu,2)

        p1cond = 0.0d0
        p2cond = 0.0d0
        v1cond = 0.0d0
        v2cond = 0.0d0
        x1cond = 0.0d0
        x2cond = 0.0d0

        do ju = 1,nu

            p1cond = p1cond + Pu(iu,ju)*p1mat0(iu,ju,:,:,:)
            p2cond = p2cond + Pu(iu,ju)*p2mat0(iu,ju,:,:,:)
            v1cond = v1cond + Pu(iu,ju)*v1mat0(iu,ju,:,:,:)
            v2cond = v2cond + Pu(iu,ju)*v2mat0(iu,ju,:,:,:)
            x1cond = x1cond + Pu(iu,ju)*x1mat0(iu,ju,:,:,:)
            x2cond = x2cond + Pu(iu,ju)*x2mat0(iu,ju,:,:,:)

        end do

        cp1 = spfit3(invT1,invT2,invTx,p1cond,r1,r2,rx,knots1,knots2,knotsx)
        cp2 = spfit3(invT1,invT2,invTx,p2cond,r1,r2,rx,knots1,knots2,knotsx)
        cv1 = spfit3(invT1,invT2,invTx,v1cond,r1,r2,rx,knots1,knots2,knotsx)
        cv2 = spfit3(invT1,invT2,invTx,v2cond,r1,r2,rx,knots1,knots2,knotsx)
        cx1 = spfit3(invT1,invT2,invTx,x1cond,r1,r2,rx,knots1,knots2,knotsx)
        cx2 = spfit3(invT1,invT2,invTx,x2cond,r1,r2,rx,knots1,knots2,knotsx)

        !$omp parallel do private(i1,i2,xx0,xx1,rc, &
        !$omp y1,y2,z,x,p1,p2,v1,v2,v1f,v2f,xi1,xi2,dfx,dfy,dfz, &
        !$omp xp,y1p,y2p,w1,w11,w12,w2,w21,w22) &
        !$omp firstprivate(u1p,u2p,u1,u2,cp1,cp2,cv1,cv2,cx1,cx2)
        do ix = 1,nx

            xp = knotsx(ix)

            do i1 = 1,n1

                y1p = knots1(i1)

                do i2 = 1,n2

                    y2p = knots2(i2)
                    w1 = w1mat(ku,iu,i1,i2)
                    w11 = w11mat(ku,iu,i1,i2)
                    w12 = w12mat(ku,iu,i1,i2)
                    w2 = w2mat(ku,iu,i1,i2)
                    w21 = w21mat(ku,iu,i1,i2)
                    w22 = w22mat(ku,iu,i1,i2)

                    ! first try
                    xx0 = (/y1mat(ku,iu,i1,i2,ix), y2mat(ku,iu,i1,i2,ix)/)
                    call csolve(xx0,rc,1d-6,50,0,u1p,u2p,u1,u2,xp,y1p,y2p,&
                    w1,w11,w12,w2,w21,w22,cp1,cp2,cv1,cv2,cx1,cx2)
                    y1 = xx0(1)
                    y2 = xx0(2)
                    z = 1.0d0
                    x = xp

                    call calcp(p1,p2,y1,y2,z,x,u1p,u2p,u1,u2,xp,y1p,y2p,0.0d0,0.0d0,cx1,cx2)
                    call speva3(cv1,y1,y2,x,r1,r2,rx,knots1,knots2,knotsx,v1f,dfx,dfy,dfz)
                    call speva3(cv2,y1,y2,x,r1,r2,rx,knots1,knots2,knotsx,v2f,dfx,dfy,dfz)
                    ! eq. (1)
                    v1 = -(gam*(eta+rho)*(y1 - (1.0d0-gam)/gam/(eta+rho)*u1)**2 &
                    + (1.0d0-gam)*(eta+rho)*(y2 + 1.0d0/(eta+rho)*u2)**2 &
                    + gam*(1.0d0-gam)*(1.0d0-rho)*(y1-y2)**2 &
                    + gam*sig/kap1*p1**2 + (1.0d0-gam)*sig/kap2*p2**2) + bet*v1f
                    ! eq. (2)
                    v2 = -((1.0d0-gam)*(eta+rho)*(y2 - gam/(1.0d0-gam)/(eta+rho)*u2)**2 &
                    + gam*(eta+rho)*(y1 + 1.0d0/(eta+rho)*u1)**2 &
                    + gam*(1.0d0-gam)*(1.0d0-rho)*(y2-y1)**2 &
                    + gam*sig/kap1*p1**2 + (1.0d0-gam)*sig/kap2*p2**2) + bet*v2f

                    xi1 = 0.0d0
                    xi2 = 0.0d0

                    ! check the sustainability constraints is binding.
                    ! note either of the constraints binds at a time.
                    ! if Home constraint is binding
                    if (v1<=w1) then

                        xx1 = (/y1mat(ku,iu,i1,i2,ix), y2mat(ku,iu,i1,i2,ix), zmat(ku,iu,i1,i2,ix)/)
                        call csolve(xx1,rc,1d-6,50,1,u1p,u2p,u1,u2,xp,y1p,y2p,&
                        w1,w11,w12,w2,w21,w22,cp1,cp2,cv1,cv2,cx1,cx2)
                        y1 = xx1(1)
                        y2 = xx1(2)
                        z = xx1(3)
                        x = 1.0d0-z*(1.0d0-xp)

                        call calcp(p1,p2,y1,y2,z,x,u1p,u2p,u1,u2,xp,y1p,y2p,w11,w12,cx1,cx2)
                        call speva3(cv1,y1,y2,x,r1,r2,rx,knots1,knots2,knotsx,v1f,dfx,dfy,dfz)
                        call speva3(cv2,y1,y2,x,r1,r2,rx,knots1,knots2,knotsx,v2f,dfx,dfy,dfz)
                        ! eq. (1)
                        v1 = -(gam*(eta+rho)*(y1 - (1.0d0-gam)/gam/(eta+rho)*u1)**2 &
                        + (1.0d0-gam)*(eta+rho)*(y2 + 1.0d0/(eta+rho)*u2)**2 &
                        + gam*(1.0d0-gam)*(1.0d0-rho)*(y1-y2)**2 &
                        + gam*sig/kap1*p1**2 + (1.0d0-gam)*sig/kap2*p2**2) + bet*v1f
                        ! eq. (2)
                        v2 = -((1.0d0-gam)*(eta+rho)*(y2 - gam/(1.0d0-gam)/(eta+rho)*u2)**2 &
                        + gam*(eta+rho)*(y1 + 1.0d0/(eta+rho)*u1)**2 &
                        + gam*(1.0d0-gam)*(1.0d0-rho)*(y2-y1)**2 &
                        + gam*sig/kap1*p1**2 + (1.0d0-gam)*sig/kap2*p2**2) + bet*v2f
                        !v1 = w1

                        xi1 = (1.0d0/z-1.0d0)*w11
                        xi2 = (1.0d0/z-1.0d0)*w12

                    ! if Foreign constraint is binding
                    elseif (v2<=w2) then

                        xx1 = (/y1mat(ku,iu,i1,i2,ix), y2mat(ku,iu,i1,i2,ix), zmat(ku,iu,i1,i2,ix)/)
                        call csolve(xx1,rc,1d-6,50,2,u1p,u2p,u1,u2,xp,y1p,y2p,&
                        w1,w11,w12,w2,w21,w22,cp1,cp2,cv1,cv2,cx1,cx2)
                        y1 = xx1(1)
                        y2 = xx1(2)
                        z = xx1(3)
                        x = z*xp

                        call calcp(p1,p2,y1,y2,z,x,u1p,u2p,u1,u2,xp,y1p,y2p,w21,w22,cx1,cx2)
                        call speva3(cv1,y1,y2,x,r1,r2,rx,knots1,knots2,knotsx,v1f,dfx,dfy,dfz)
                        call speva3(cv2,y1,y2,x,r1,r2,rx,knots1,knots2,knotsx,v2f,dfx,dfy,dfz)
                        ! eq. (1)
                        v1 = -(gam*(eta+rho)*(y1 - (1.0d0-gam)/gam/(eta+rho)*u1)**2 &
                        + (1.0d0-gam)*(eta+rho)*(y2 + 1.0d0/(eta+rho)*u2)**2 &
                        + gam*(1.0d0-gam)*(1.0d0-rho)*(y1-y2)**2 &
                        + gam*sig/kap1*p1**2 + (1.0d0-gam)*sig/kap2*p2**2) + bet*v1f
                        ! eq. (2)
                        v2 = -((1.0d0-gam)*(eta+rho)*(y2 - gam/(1.0d0-gam)/(eta+rho)*u2)**2 &
                        + gam*(eta+rho)*(y1 + 1.0d0/(eta+rho)*u1)**2 &
                        + gam*(1.0d0-gam)*(1.0d0-rho)*(y2-y1)**2 &
                        + gam*sig/kap1*p1**2 + (1.0d0-gam)*sig/kap2*p2**2) + bet*v2f
                        !v2 = w2

                        xi1 = (1.0d0/z-1.0d0)*w21
                        xi2 = (1.0d0/z-1.0d0)*w22

                    end if

                    if (z<1d-4 .or. z>1.0d0 .or. rc/=0) then

                        y1 = y1mat(ku,iu,i1,i2,ix)
                        y2 = y2mat(ku,iu,i1,i2,ix)
                        z = zmat(ku,iu,i1,i2,ix)
                        x = xmat(ku,iu,i1,i2,ix)
                        p1 = p1mat0(ku,iu,i1,i2,ix)
                        p2 = p2mat0(ku,iu,i1,i2,ix)
                        v1 = v1mat0(ku,iu,i1,i2,ix)
                        v2 = v2mat0(ku,iu,i1,i2,ix)
                        xi1 = x1mat0(ku,iu,i1,i2,ix)
                        xi2 = x2mat0(ku,iu,i1,i2,ix)

                    end if

                    y1mat(ku,iu,i1,i2,ix) = y1
                    y2mat(ku,iu,i1,i2,ix) = y2
                    zmat(ku,iu,i1,i2,ix)  = z
                    xmat(ku,iu,i1,i2,ix)  = x
                    p1mat(ku,iu,i1,i2,ix) = p1
                    p2mat(ku,iu,i1,i2,ix) = p2
                    v1mat(ku,iu,i1,i2,ix) = v1
                    v2mat(ku,iu,i1,i2,ix) = v2
                    x1mat(ku,iu,i1,i2,ix) = xi1
                    x2mat(ku,iu,i1,i2,ix) = xi2

                end do

            end do

        end do
        !$omp end parallel do

    end do

end do
!$omp end parallel do


end subroutine sus_iter


end module mod_main
