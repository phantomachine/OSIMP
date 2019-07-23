module mod_calcpar
    

implicit none


contains
    
    
subroutine calcpar(y1mat,y2mat,p1mat,p2mat,v1mat,v2mat,Gu,Pu,invT1,invT2)

use mod_focpar
use mod_spline
use mod_parameters
implicit none

real(8), intent(out) :: y1mat(nu,nu,n1,n2), y2mat(nu,nu,n1,n2), p1mat(nu,nu,n1,n2), p2mat(nu,nu,n1,n2), &
v1mat(nu,nu,n1,n2), v2mat(nu,nu,n1,n2) 
real(8), intent(in) :: Gu(nu,2), Pu(nu,nu), invT1(r1,r2), invT2(r1,r2)

real(8) p1mat0(nu,nu,n1,n2), p2mat0(nu,nu,n1,n2), v1mat0(nu,nu,n1,n2), v2mat0(nu,nu,n1,n2)

real(8) u1p, u2p, u1, u2, y1p, y2p, w1, w11, w12, w2, w21, w22
real(8), allocatable :: p1cond(:,:), p2cond(:,:), v1cond(:,:), v2cond(:,:)
real(8), allocatable :: cpn1(:,:,:), cpn2(:,:,:)
real(8), allocatable :: cvn1(:,:,:), cvn2(:,:,:)

real(8) xx0(2), f0(2), y1, y2, p1, p2, xi1, xi2, xi1p, xi2p, mu1, mu2, mu1p, mu2p, &
v1, v2, v1f, v2f, dfx, dfy
real(8) diff1, diff2, diff, crit
integer iu, ju, ku, i1, i2, i, rc

allocate(p1cond(n1,n2), p2cond(n1,n2), v1cond(n1,n2), v2cond(n1,n2))
allocate(cpn1(16,r1+1,r2+1), cpn2(16,r1+1,r2+1))
allocate(cvn1(16,r1+1,r2+1), cvn2(16,r1+1,r2+1))

! calc policy functions
! init values
p1mat0 = 0.0d0
p2mat0 = 0.0d0

crit = 1d-5
diff = 1d+4

do while (diff>crit)

    do ku = 1,nu
    
        u1p = Gu(ku,1)
        u2p = Gu(ku,2)

        do iu = 1,nu

            u1 = Gu(iu,1)
            u2 = Gu(iu,2)
        
            p1cond = 0.0d0
            p2cond = 0.0d0
            
            do ju = 1,nu

                p1cond = p1cond + Pu(iu,ju)*p1mat0(iu,ju,:,:)
                p2cond = p2cond + Pu(iu,ju)*p2mat0(iu,ju,:,:)
                
            end do
            
            cpn1 = spfit2(invT1,invT2,p1cond,r1,r2,knots1,knots2)
            cpn2 = spfit2(invT1,invT2,p2cond,r1,r2,knots1,knots2)
            
            do i1 = 1,n1

                y1p = knots1(i1)

                do i2 = 1,n2

                    y2p = knots2(i2)
                    
                    xx0 = 0.0d0
                    call csolve(xx0,rc,1d-6,50,u1p,u2p,u1,u2,y1p,y2p,cpn1,cpn2)
                    y1 = xx0(1)
                    y2 = xx0(2)

                    ! trade-off eqs.
                    mu1 = (lam-gam)*u1
                    mu2 = (lam-gam)*u2
                    mu1p = (lam-gam)*u1p
                    mu2p = (lam-gam)*u2p

                    xi1  = ((eta+rho+gam*(1.0d0-rho))*mu1 -gam*(1.0d0-rho)*mu2 )/(gam*(1.0d0+eta)*(rho+eta))
                    xi2  = -((eta+rho+(1.0d0-gam)*(1.0d0-rho))*mu2 -(1.0d0-gam)*(1.0d0-rho)*mu1 )/((1.0d0-gam)*(1.0d0+eta)*(rho+eta))
                    xi1p = ((eta+rho+gam*(1.0d0-rho))*mu1p-gam*(1.0d0-rho)*mu2p)/(gam*(1.0d0+eta)*(rho+eta))
                    xi2p = -((eta+rho+(1.0d0-gam)*(1.0d0-rho))*mu2p-(1.0d0-gam)*(1.0d0-rho)*mu1p)/((1.0d0-gam)*(1.0d0+eta)*(rho+eta))
                    p1 = -1.0d0/sig*(y1 - y1p - xi1 + xi1p)
                    p2 = -1.0d0/sig*(y2 - y2p - xi2 + xi2p)
                                                                
                    y1mat(ku,iu,i1,i2) = y1
                    y2mat(ku,iu,i1,i2) = y2
                    p1mat(ku,iu,i1,i2) = p1
                    p2mat(ku,iu,i1,i2) = p2
                                                                      
                end do

            end do
            
        end do
        
    end do

    diff1 = maxval(maxval(maxval(maxval(abs(p1mat-p1mat0),1),1),1),1)
    diff2 = maxval(maxval(maxval(maxval(abs(p2mat-p2mat0),1),1),1),1)
    diff = max(diff1,diff2)
    print *, diff

    p1mat0 = p1mat
    p2mat0 = p2mat
    
end do


! calc welfare values
v1mat0 = 0.0d0
v2mat0 = 0.0d0

crit = 1d-5
diff = 1d+4

do while (diff>crit)
    
    do ku = 1,nu
        
        u1p = Gu(ku,1)
        u2p = Gu(ku,2)
    
        do iu = 1,nu
            
            u1 = Gu(iu,1)
            u2 = Gu(iu,2)
            
            v1cond = 0.0d0
            v2cond = 0.0d0

            do ju = 1,nu

                v1cond = v1cond + Pu(iu,ju)*v1mat0(iu,ju,:,:)
                v2cond = v2cond + Pu(iu,ju)*v2mat0(iu,ju,:,:)

            end do

            cvn1 = spfit2(invT1,invT2,v1cond,r1,r2,knots1,knots2)
            cvn2 = spfit2(invT1,invT2,v2cond,r1,r2,knots1,knots2)

            do i1 = 1,n1

                do i2 = 1,n2

                    y1 = y1mat(ku,iu,i1,i2)
                    y2 = y2mat(ku,iu,i1,i2)
                    p1 = p1mat(ku,iu,i1,i2)
                    p2 = p2mat(ku,iu,i1,i2)

                    call speva2(cvn1,y1,y2,r1,r2,knots1,knots2,v1f,dfx,dfy)
                    call speva2(cvn2,y1,y2,r1,r2,knots1,knots2,v2f,dfx,dfy)                    
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

                    v1mat(ku,iu,i1,i2) = v1
                    v2mat(ku,iu,i1,i2) = v2

                end do

            end do

        end do
        
    end do
    
    diff1 = maxval(maxval(maxval(maxval(abs(v1mat-v1mat0),1),1),1),1)
    diff2 = maxval(maxval(maxval(maxval(abs(v2mat-v2mat0),1),1),1),1)
    diff = max(diff1,diff2)
    print *, diff

    v1mat0 = v1mat
    v2mat0 = v2mat
    
end do


end subroutine calcpar


end module mod_calcpar