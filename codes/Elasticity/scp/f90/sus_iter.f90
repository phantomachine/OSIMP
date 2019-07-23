#include "fintrf.h"

subroutine mexFunction(nlhs, plhs, nrhs, prhs)

    use mod_main
	use mod_parameters
    implicit none

    mwPointer plhs(*), prhs(*)
    integer nlhs, nrhs
    mwPointer mxGetM, mxGetN

    mwPointer mxCreateNumericArray, mxGetPr
    integer mxClassIDFromClassName
    mwSize ndim, dims(5) ! integer(8) with -largeArrayDims option

    real(8), allocatable :: Gu(:,:), Pu(:,:), invT1(:,:), invT2(:,:), invTx(:,:)
    real(8), allocatable :: y1mat(:,:,:,:,:), y2mat(:,:,:,:,:), zmat(:,:,:,:,:), xmat(:,:,:,:,:)
    real(8), allocatable :: p1mat(:,:,:,:,:), p2mat(:,:,:,:,:), p1mat0(:,:,:,:,:), p2mat0(:,:,:,:,:)
    real(8), allocatable :: v1mat(:,:,:,:,:), v2mat(:,:,:,:,:), v1mat0(:,:,:,:,:), v2mat0(:,:,:,:,:)
    real(8), allocatable :: x1mat(:,:,:,:,:), x2mat(:,:,:,:,:), x1mat0(:,:,:,:,:), x2mat0(:,:,:,:,:)
    real(8), allocatable :: w1mat(:,:,:,:), w2mat(:,:,:,:), w11mat(:,:,:,:), w12mat(:,:,:,:), w21mat(:,:,:,:), w22mat(:,:,:,:)

    integer, external :: mexPrintf
    integer k, s
    character(len=80) line

! call omp_set_dynamic(.false.)
! call omp_set_num_threads(24)
! call kmp_set_stacksize_s(10240000)

! get the stack size in each thread
!omp parallel private(s)
!call kmp_get_stacksize_s()
!write(line,*) s
!k = mexPrintf(line//achar(10))
!omp end parallel

    nu = mxGetM(prhs(18))
    r1 = mxGetM(prhs(19))
    r2 = mxGetM(prhs(20))
    rx = mxGetM(prhs(21))
    n1 = mxGetM(prhs(29))
    n2 = mxGetM(prhs(30))
    nx = mxGetM(prhs(31))

	allocate(knots1(n1),knots2(n2),knotsx(nx))
	allocate(Gu(nu,2),Pu(nu,nu),invT1(r1,r1),invT2(r2,r2),invTx(rx,rx))
    allocate(y1mat(nu,nu,n1,n2,nx),y2mat(nu,nu,n1,n2,nx),zmat(nu,nu,n1,n2,nx),xmat(nu,nu,n1,n2,nx))
    allocate(p1mat(nu,nu,n1,n2,nx),p2mat(nu,nu,n1,n2,nx),p1mat0(nu,nu,n1,n2,nx),p2mat0(nu,nu,n1,n2,nx))
    allocate(v1mat(nu,nu,n1,n2,nx),v2mat(nu,nu,n1,n2,nx),v1mat0(nu,nu,n1,n2,nx),v2mat0(nu,nu,n1,n2,nx))
    allocate(x1mat(nu,nu,n1,n2,nx),x2mat(nu,nu,n1,n2,nx),x1mat0(nu,nu,n1,n2,nx),x2mat0(nu,nu,n1,n2,nx))
	allocate(w1mat(nu,nu,n1,n2),w2mat(nu,nu,n1,n2),w11mat(nu,nu,n1,n2),w12mat(nu,nu,n1,n2),w21mat(nu,nu,n1,n2),w22mat(nu,nu,n1,n2))

    call mxCopyPtrToReal8(mxGetPr(prhs(1)),p1mat0,nu*nu*n1*n2*nx)
    call mxCopyPtrToReal8(mxGetPr(prhs(2)),p2mat0,nu*nu*n1*n2*nx)
    call mxCopyPtrToReal8(mxGetPr(prhs(3)),v1mat0,nu*nu*n1*n2*nx)
    call mxCopyPtrToReal8(mxGetPr(prhs(4)),v2mat0,nu*nu*n1*n2*nx)
    call mxCopyPtrToReal8(mxGetPr(prhs(5)),x1mat0,nu*nu*n1*n2*nx)
    call mxCopyPtrToReal8(mxGetPr(prhs(6)),x2mat0,nu*nu*n1*n2*nx)
    call mxCopyPtrToReal8(mxGetPr(prhs(7)),y1mat,nu*nu*n1*n2*nx)
    call mxCopyPtrToReal8(mxGetPr(prhs(8)),y2mat,nu*nu*n1*n2*nx)
    call mxCopyPtrToReal8(mxGetPr(prhs(9)),zmat,nu*nu*n1*n2*nx)
    call mxCopyPtrToReal8(mxGetPr(prhs(10)),xmat,nu*nu*n1*n2*nx)
    call mxCopyPtrToReal8(mxGetPr(prhs(11)),w1mat,nu*nu*n1*n2)
    call mxCopyPtrToReal8(mxGetPr(prhs(12)),w2mat,nu*nu*n1*n2)
    call mxCopyPtrToReal8(mxGetPr(prhs(13)),w11mat,nu*nu*n1*n2)
    call mxCopyPtrToReal8(mxGetPr(prhs(14)),w12mat,nu*nu*n1*n2)
    call mxCopyPtrToReal8(mxGetPr(prhs(15)),w21mat,nu*nu*n1*n2)
    call mxCopyPtrToReal8(mxGetPr(prhs(16)),w22mat,nu*nu*n1*n2)
    call mxCopyPtrToReal8(mxGetPr(prhs(17)),Gu,nu*2)
    call mxCopyPtrToReal8(mxGetPr(prhs(18)),Pu,nu*nu)
    call mxCopyPtrToReal8(mxGetPr(prhs(19)),invT1,r1*r1)
    call mxCopyPtrToReal8(mxGetPr(prhs(20)),invT2,r2*r2)
    call mxCopyPtrToReal8(mxGetPr(prhs(21)),invTx,rx*rx)

    call mxCopyPtrToReal8(mxGetPr(prhs(22)),bet,1)
    call mxCopyPtrToReal8(mxGetPr(prhs(23)),eta,1)
    call mxCopyPtrToReal8(mxGetPr(prhs(24)),rho,1)
    call mxCopyPtrToReal8(mxGetPr(prhs(25)),sig,1)
    call mxCopyPtrToReal8(mxGetPr(prhs(26)),kap1,1)
    call mxCopyPtrToReal8(mxGetPr(prhs(27)),kap2,1)
    call mxCopyPtrToReal8(mxGetPr(prhs(28)),iota,1)
    call mxCopyPtrToReal8(mxGetPr(prhs(29)),knots1,n1)
    call mxCopyPtrToReal8(mxGetPr(prhs(30)),knots2,n2)
    call mxCopyPtrToReal8(mxGetPr(prhs(31)),knotsx,nx)

    call sus_iter(p1mat,p2mat,v1mat,v2mat,x1mat,x2mat,&
    p1mat0,p2mat0,v1mat0,v2mat0,x1mat0,x2mat0,y1mat,y2mat,zmat,xmat,&
    w1mat,w2mat,w11mat,w12mat,w21mat,w22mat,Gu,Pu,invT1,invT2,invTx)

	dims(1) = nu
	dims(2) = nu
	dims(3) = n1
	dims(4) = n2
    dims(5) = nx

    plhs(1) = mxCreateNumericArray(5, dims, mxClassIDFromClassName('double'), 0)
    plhs(2) = mxCreateNumericArray(5, dims, mxClassIDFromClassName('double'), 0)
    plhs(3) = mxCreateNumericArray(5, dims, mxClassIDFromClassName('double'), 0)
    plhs(4) = mxCreateNumericArray(5, dims, mxClassIDFromClassName('double'), 0)
    plhs(5) = mxCreateNumericArray(5, dims, mxClassIDFromClassName('double'), 0)
    plhs(6) = mxCreateNumericArray(5, dims, mxClassIDFromClassName('double'), 0)
    plhs(7) = mxCreateNumericArray(5, dims, mxClassIDFromClassName('double'), 0)
    plhs(8) = mxCreateNumericArray(5, dims, mxClassIDFromClassName('double'), 0)
    plhs(9) = mxCreateNumericArray(5, dims, mxClassIDFromClassName('double'), 0)
    plhs(10) = mxCreateNumericArray(5, dims, mxClassIDFromClassName('double'), 0)

    call mxCopyReal8ToPtr(p1mat,mxGetPr(plhs(1)),nu*nu*n1*n2*nx)
    call mxCopyReal8ToPtr(p2mat,mxGetPr(plhs(2)),nu*nu*n1*n2*nx)
    call mxCopyReal8ToPtr(v1mat,mxGetPr(plhs(3)),nu*nu*n1*n2*nx)
    call mxCopyReal8ToPtr(v2mat,mxGetPr(plhs(4)),nu*nu*n1*n2*nx)
    call mxCopyReal8ToPtr(x1mat,mxGetPr(plhs(5)),nu*nu*n1*n2*nx)
    call mxCopyReal8ToPtr(x2mat,mxGetPr(plhs(6)),nu*nu*n1*n2*nx)
    call mxCopyReal8ToPtr(y1mat,mxGetPr(plhs(7)),nu*nu*n1*n2*nx)
    call mxCopyReal8ToPtr(y2mat,mxGetPr(plhs(8)),nu*nu*n1*n2*nx)
    call mxCopyReal8ToPtr(zmat, mxGetPr(plhs(9)),nu*nu*n1*n2*nx)
    call mxCopyReal8ToPtr(xmat, mxGetPr(plhs(10)),nu*nu*n1*n2*nx)

	deallocate(knots1,knots2,knotsx)

    return

end subroutine mexFunction
