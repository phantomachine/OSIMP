module mod_foccp

implicit none

contains


function foccp(xx0,u1,u2,y1p,y2p,cpn1,cpn2) result(f0)

    use mod_spline
	use mod_parameters

    implicit none
    real(8), intent(in) :: xx0(:)

    real(8), intent(in) :: u1, u2, y1p, y2p, cpn1(:,:,:), cpn2(:,:,:)

    real(8) f0(size(xx0))
    real(8) y1, y2, p1, p2, p1f, p2f, xi1, xi2, xi1p, xi2p, dfx, dfy


    y1 = xx0(1)
    y2 = xx0(2)

    call speva2(cpn1,y1,y2,r1,r2,knots1,knots2,p1f,dfx,dfy)
    call speva2(cpn2,y1,y2,r1,r2,knots1,knots2,p2f,dfx,dfy)

    ! trade-off eqs.
    p1 = -1.0d0/sig*(y1 - y1p)
    p2 = -1.0d0/sig*(y2 - y2p)

    ! eqs. (3) and (4)
    f0(1) = -p1 + bet*p1f + kap1*u1 + kap1*((rho+eta)*y1 + (1.0d0-gam)*(1.0d0-rho)*(y1-y2))
    f0(2) = -p2 + bet*p2f + kap2*u2 + kap2*((rho+eta)*y2 - gam*(1.0d0-rho)*(y1-y2))


end function foccp


! csolve in fortran - needs to be externalized?
! see https://scicomp.stackexchange.com/questions/11701/varargin-in-fortran
subroutine csolve(x,rc,crit,itmax,u1,u2,y1p,y2p,cpn1,cpn2)

    implicit none

    real(8), intent(inout) :: x(:)
    integer, intent(out) :: rc
    real(8), intent(in) :: crit
    integer, intent(in) :: itmax

    real(8), intent(in) :: u1, u2, y1p, y2p, cpn1(:,:,:), cpn2(:,:,:)

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
    f0 = foccp(x,u1,u2,y1p,y2p,cpn1,cpn2)
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
            f1 = foccp(x0,u1,u2,y1p,y2p,cpn1,cpn2)
            grad(:,iv) = (f1-f0)/delta

        end do

       ! else % use analytic gradient
       !     grad=feval(gradfun,x,varargin{:});
       ! end
	   !
       ! if isreal(grad)
       !      grad
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
            f = foccp(x+dx,u1,u2,y1p,y2p,cpn1,cpn2)
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


end module mod_foccp
