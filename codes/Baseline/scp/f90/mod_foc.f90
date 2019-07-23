module mod_foc


implicit none


contains


function foc(xx0,flag,u1p,u2p,u1,u2,xp,y1p,y2p,w1,w11,w12,w2,w21,w22,cp1,cp2,cv1,cv2,cx1,cx2) result(f0)


    use mod_parameters
    use mod_spline

    implicit none
    real(8), intent(in) :: xx0(:)

    integer, intent(in) :: flag
    real(8), intent(in) :: u1p, u2p, u1, u2, xp, y1p, y2p, w1, w11, w12, w2, w21, w22, &
cp1(:,:,:,:), cp2(:,:,:,:), cv1(:,:,:,:), cv2(:,:,:,:), cx1(:,:,:,:), cx2(:,:,:,:)

    real(8) f0(size(xx0))
    real(8) y1, y2, z, x, p1, p2, p1f, p2f, v1f, v2f, xi1, xi2, xi1p, xi2p, dfx, dfy, dfz


    if (flag==0) then

        y1 = xx0(1)
        y2 = xx0(2)
        z = 1.0d0
        x = xp

        call calcp(p1,p2,y1,y2,z,x,u1p,u2p,u1,u2,xp,y1p,y2p,0.0d0,0.0d0,cx1,cx2)

        call speva3(cp1,y1,y2,x,r1,r2,rx,knots1,knots2,knotsx,p1f,dfx,dfy,dfz)
        call speva3(cp2,y1,y2,x,r1,r2,rx,knots1,knots2,knotsx,p2f,dfx,dfy,dfz)

        ! eqs. (3) and (4)
        f0(1) = -p1 + bet*p1f + kap1*u1 + kap1*((rho+eta)*y1 + (1.0d0-gam)*(1.0d0-rho)*(y1-y2))
        f0(2) = -p2 + bet*p2f + kap2*u2 + kap2*((rho+eta)*y2 - gam*(1.0d0-rho)*(y1-y2))

    elseif (flag==1) then

        y1 = xx0(1)
        y2 = xx0(2)
        z = xx0(3)
        x = 1.0d0-z*(1.0d0-xp)

        call calcp(p1,p2,y1,y2,z,x,u1p,u2p,u1,u2,xp,y1p,y2p,w11,w12,cx1,cx2)

        call speva3(cp1,y1,y2,x,r1,r2,rx,knots1,knots2,knotsx,p1f,dfx,dfy,dfz)
        call speva3(cp2,y1,y2,x,r1,r2,rx,knots1,knots2,knotsx,p2f,dfx,dfy,dfz)
        call speva3(cv1,y1,y2,x,r1,r2,rx,knots1,knots2,knotsx,v1f,dfx,dfy,dfz)

        ! eqs. (3) and (4)
        f0(1) = -p1 + bet*p1f + kap1*u1 + kap1*((rho+eta)*y1 + (1.0d0-gam)*(1.0d0-rho)*(y1-y2))
        f0(2) = -p2 + bet*p2f + kap2*u2 + kap2*((rho+eta)*y2 - gam*(1.0d0-rho)*(y1-y2))
        ! eq. (1)
        f0(3) = -w1 -(gam*(eta+rho)*(y1 - (1.0d0-gam)/gam/(eta+rho)*u1)**2 &
        + (1.0d0-gam)*(eta+rho)*(y2 + 1.0d0/(eta+rho)*u2)**2 &
        + gam*(1.0d0-gam)*(1.0d0-rho)*(y1-y2)**2 &
        + gam*sig/kap1*p1**2 + (1.0d0-gam)*sig/kap2*p2**2) + bet*v1f

    elseif (flag==2) then

        y1 = xx0(1)
        y2 = xx0(2)
        z = xx0(3)
        x = z*xp

        call calcp(p1,p2,y1,y2,z,x,u1p,u2p,u1,u2,xp,y1p,y2p,w21,w22,cx1,cx2)

        call speva3(cp1,y1,y2,x,r1,r2,rx,knots1,knots2,knotsx,p1f,dfx,dfy,dfz)
        call speva3(cp2,y1,y2,x,r1,r2,rx,knots1,knots2,knotsx,p2f,dfx,dfy,dfz)
        call speva3(cv2,y1,y2,x,r1,r2,rx,knots1,knots2,knotsx,v2f,dfx,dfy,dfz)

        ! eqs. (3) and (4)
        f0(1) = -p1 + bet*p1f + kap1*u1 + kap1*((rho+eta)*y1 + (1.0d0-gam)*(1.0d0-rho)*(y1-y2))
        f0(2) = -p2 + bet*p2f + kap2*u2 + kap2*((rho+eta)*y2 - gam*(1.0d0-rho)*(y1-y2))
        ! eq. (2)
        f0(3) = -w2 -((1.0d0-gam)*(eta+rho)*(y2 - gam/(1.0d0-gam)/(eta+rho)*u2)**2 &
        + gam*(eta+rho)*(y1 + 1.0d0/(eta+rho)*u1)**2 &
        + gam*(1.0d0-gam)*(1.0d0-rho)*(y2-y1)**2 &
        + gam*sig/kap1*p1**2 + (1.0d0-gam)*sig/kap2*p2**2) + bet*v2f

    end if


end function foc


subroutine calcp(p1,p2,y1,y2,z,x,u1p,u2p,u1,u2,xp,y1p,y2p,w11,w12,cx1,cx2)


    use mod_parameters
    use mod_spline

    implicit none
    real(8), intent(out) :: p1, p2
    real(8), intent(in) :: y1, y2, z, x
    real(8), intent(in) :: u1p, u2p, u1, u2, xp, y1p, y2p, w11, w12, cx1(:,:,:,:), cx2(:,:,:,:)

    real(8) xi1, xi2, eta1p, eta2p, zet1p, zet2p, xi1f, xi2f, zet1, zet2, eta1, eta2, dfx, dfy, dfz

    ! trade-off eqs.
    xi1 = (1.0d0/z-1.0d0)*w11
    xi2 = (1.0d0/z-1.0d0)*w12
    eta1p = (xp-gam)*u1p - bet/2.0d0*xi1
    eta2p = (xp-gam)*u2p + bet/2.0d0*xi2
    zet1p = ((eta+rho+gam*(1.0d0-rho))*eta1p-gam*(1.0d0-rho)*eta2p)/(gam*(1.0d0+eta)*(rho+eta))
    zet2p = -((eta+rho+(1.0d0-gam)*(1.0d0-rho))*eta2p-(1.0d0-gam)*(1.0d0-rho)*eta1p)/((1.0d0-gam)*(1.0d0+eta)*(rho+eta))

    call speva3(cx1,y1,y2,x,r1,r2,rx,knots1,knots2,knotsx,xi1f,dfx,dfy,dfz)
    call speva3(cx2,y1,y2,x,r1,r2,rx,knots1,knots2,knotsx,xi2f,dfx,dfy,dfz)
    eta1 = (x-gam)*u1 - bet/2.0d0*xi1f
    eta2 = (x-gam)*u2 + bet/2.0d0*xi2f
    zet1 = ((eta+rho+gam*(1.0d0-rho))*eta1-gam*(1.0d0-rho)*eta2)/(gam*(1.0d0+eta)*(rho+eta))
    zet2 = -((eta+rho+(1.0d0-gam)*(1.0d0-rho))*eta2-(1.0d0-gam)*(1.0d0-rho)*eta1)/((1.0d0-gam)*(1.0d0+eta)*(rho+eta))

    p1 = -1.0d0/sig*(y1-zet1-z*(y1p-zet1p))
    p2 = -1.0d0/sig*(y2-zet2-z*(y2p-zet2p))


end subroutine calcp


subroutine csolve(x,rc,crit,itmax,flag,u1p,u2p,u1,u2,xp,y1p,y2p,&
w1,w11,w12,w2,w21,w22,cp1,cp2,cv1,cv2,cx1,cx2)

    implicit none

    real(8), intent(inout) :: x(:)
    integer, intent(out) :: rc
    real(8), intent(in) :: crit
    integer, intent(in) :: itmax

    integer, intent(in) :: flag
    real(8), intent(in) :: u1p, u2p, u1, u2, xp, y1p, y2p, w1, w11, w12, w2, w21, w22, &
cp1(:,:,:,:), cp2(:,:,:,:), cv1(:,:,:,:), cv2(:,:,:,:), cx1(:,:,:,:), cx2(:,:,:,:)

    integer iv, nv
    real(8), allocatable :: x0(:), xmin(:), f(:), f0(:), f1(:), fmin(:), dx(:), dx0(:)
    real(8), allocatable :: grad(:,:), grad_try(:,:), tvec(:,:)
    integer done, itct, verbose, shrink, subDone
    real(8) af0, af00, lambda, lambdamin, af, afmin, delta, alpha, dxSize, factor, gnorm, rcond

    nv = size(x,1)
    allocate(x0(nv), xmin(nv), f(nv), f0(nv), f1(nv), fmin(nv), dx(nv), dx0(nv))
    allocate(grad(nv,nv), grad_try(nv,nv), tvec(nv,nv))

    verbose = 0
    delta = 1d-6
    alpha = 1d-3

    done = 0

!    f0 = func(x,varargin)
!    f0 = func(x)
    f0 = foc(x,flag,u1p,u2p,u1,u2,xp,y1p,y2p,w1,w11,w12,w2,w21,w22,cp1,cp2,cv1,cv2,cx1,cx2)
    af0 = sum(abs(f0))
    af00 = af0
    itct = 0
    rc = 0

    do while (done == 0)

    !    if itct>3 & af00-af0<crit*max(1,af0) & rem(itct,2)==1
    !        randomize=1
    !    else
    !
    !    if ~analyticg

        do iv = 1,nv

            x0 = x
            x0(iv) = x(iv) + delta
            f1 = foc(x0,flag,u1p,u2p,u1,u2,xp,y1p,y2p,w1,w11,w12,w2,w21,w22,cp1,cp2,cv1,cv2,cx1,cx2)
            grad(:,iv) = (f1-f0)/delta

        end do

    !    else % use analytic gradient
    !        grad=feval(gradfun,x,varargin{:});
    !    end
    !
    !    if isreal(grad)
    !         grad
        ! condition number of grad: rcond(grad) in matlab
        grad_try = grad
!        call dgetrf(nv,nv,grad_try,nv,IPIV,INFO)
!        gnorm = maxval(sum(abs(grad_try),2),1) ! 1-norm
!        call dgecon('1',nv,grad_try,nv,gnorm,rcond,WORK,IWORK,INFO)

!        if (rcond<1d-12) then

!            do iv = 1,nv
!                grad(iv,iv)=grad(iv,iv)+delta
!            end do

!            call dgetrf(nv,nv,grad,nv,IPIV,INFO)

!        else

!            grad = grad_try

!        end if

!        call dgetri(nv,grad,nv,IPIV,WORK,nv,INFO)
		grad = myinv(grad_try,nv)
        dx0 = -matmul(grad,f0)
        ! for 1-dim
        !dx0 = -f0(1)/grad(1,1)

    !    randomize = 0
    !
    !        else
    !            if(verbose),disp('gradient imaginary'),end
    !            randomize=1;
    !        end
    !    end
    !
    !    if (randomize==1) then
    !
    !        if(verbose),fprintf(1,'\n Random Search'),end
    !        dx0=norm(x)./randn(size(x));
    !
    !    end if

        lambda = 1.0d0
        lambdamin = 1.0d0
        fmin = f0
        xmin = x
        afmin = af0
        dxSize = mynorm2(dx0,nv)

        factor = .6d0
        shrink = 1
        subDone = 0

        do while (subDone == 0)

            dx = lambda*dx0
            f = foc(x+dx,flag,u1p,u2p,u1,u2,xp,y1p,y2p,w1,w11,w12,w2,w21,w22,cp1,cp2,cv1,cv2,cx1,cx2)
            af = sum(abs(f))

            if (af<afmin) then

                afmin = af
                fmin = f
                lambdamin = lambda
                xmin = x+dx

            end if

            if ( ( (lambda > 0.0d0) .and. (af0-af < alpha*lambda*af0) ) .or. ( (lambda < 0.0d0) .and. (af0-af < 0.0d0) ) ) then

                if (shrink == 0) then

                    factor = factor**.6d0
                    shrink = 1

                end if

                if (abs(lambda*(1.0d0-factor))*dxSize > .1d0*delta) then

                    lambda = factor*lambda

                elseif ((lambda > 0.0d0) .and. (factor==.6d0)) then ! %i.e., we've only been shrinking

                    lambda = -.3d0

                else

                    subDone = 1

                    if (lambda > 0.0d0) then
                        if (factor == .6d0) then
                            rc = 2
                        else
                            rc = 1
                        end if
                    else
                        rc = 3
                    end if

                end if

            elseif ( (lambda > 0.0d0) .and. (af-af0 > (1-alpha)*lambda*af0) ) then

                if (shrink == 1) then
                    factor = factor**0.6d0
                    shrink = 0
                end if
                lambda=lambda/factor

            else ! good value found

                subDone = 1
                rc = 0

            end if

        end do ! while ~subDone

        itct = itct + 1

        x = xmin
        f0 = fmin
        af00 = af0
        af0 = afmin
        !write(*,"('  itct ', I3, ', af ', F8.5, ', lambda ', F8.5, ', rc ', F8.5)") itct, afmin, lambdamin, rc
        !write(*,"('  x ', F8.5, ' ', F8.5, ' ', F8.5, ' ', F8.5)") xmin(1), xmin(2), xmin(3), xmin(4)
        !write(*,"('  f ', F8.5, ' ', F8.5, ' ', F8.5, ' ', F8.5)") fmin(1), fmin(2), fmin(3), fmin(4)

        if (itct >= itmax) then

            done = 1
            rc = 4

        elseif (af0 < crit) then

            done = 1
            rc = 0

        end if

    end do

    !deallocate(x0, xmin, f, f0, f1, fmin, dx, dx0)
    !deallocate(grad, grad_try, tvec)
    !deallocate(IPIV, WORK)


end subroutine csolve


function mynorm2(x,d) result(n)

real(8), intent(in) :: x(:)
integer, intent(in) :: d
real(8) n
integer i

n = 0.0d0

do i = 1,d

	n = n + x(i)**2

end do

n = sqrt(n)

end function mynorm2


function myinv(A,d) result(invA)

real(8), intent(in) :: A(:,:)
integer, intent(in) :: d
real(8) invA(d,d), detA, a11, a12, a13, a21, a22, a23, a31, a32, a33


!if (n>4)

!    error('the dimension of the matrix should be equal or less than 4');

!end

if (d==1) then

    invA(1,1) = 1.0d0/A(1,1)

elseif (d==2) then

    a11 = A(1,1)
    a12 = A(1,2)
    a21 = A(2,1)
    a22 = A(2,2)

    detA = a11*a22 - a12*a21

    invA(1,1) = a22
    invA(1,2) = -a12
    invA(2,1) = -a21
    invA(2,2) = a11

    invA = invA/detA

elseif (d==3) then

    a11 = A(1,1)
    a12 = A(1,2)
    a13 = A(1,3)
    a21 = A(2,1)
    a22 = A(2,2)
    a23 = A(2,3)
    a31 = A(3,1)
    a32 = A(3,2)
    a33 = A(3,3)

    detA = a11*a22*a33+a21*a32*a13+a31*a12*a23 &
        - (a11*a32*a23+a31*a22*a13+a21*a12*a33)

    invA(1,1) = a22*a33-a23*a32
    invA(1,2) = a13*a32-a12*a33
    invA(1,3) = a12*a23-a13*a22
    invA(2,1) = a23*a31-a21*a33
    invA(2,2) = a11*a33-a13*a31
    invA(2,3) = a13*a21-a11*a23
    invA(3,1) = a21*a32-a22*a31
    invA(3,2) = a12*a31-a11*a32
    invA(3,3) = a11*a22-a12*a21

    invA = invA/detA

end if


end function myinv


end module mod_foc
