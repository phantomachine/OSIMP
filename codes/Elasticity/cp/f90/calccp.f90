#include "fintrf.h"
! NOTE: all the integers passed to lapack routines must be double precision with -largeArrayDims option

subroutine mexFunction(nlhs, plhs, nrhs, prhs)

    use mod_calccp
	use mod_parameters
    implicit none

    mwPointer plhs(*), prhs(*)
    integer nlhs, nrhs
    mwPointer mxGetM, mxGetN

    mwPointer mxCreateNumericArray, mxGetPr
    integer mxClassIDFromClassName
    mwSize ndim, dims(3) ! integer(8) with -largeArrayDims option

    real(8), allocatable :: Gu(:,:), Pu(:,:), invT1(:,:), invT2(:,:)
    real(8), allocatable :: y1mat(:,:,:), y2mat(:,:,:), p1mat(:,:,:), p2mat(:,:,:), v1mat(:,:,:), v2mat(:,:,:)

    integer, external :: mexPrintf
    integer k
    character(len=80) line

    nu = mxGetM(prhs(2))
    r1 = mxGetM(prhs(3))
    r2 = mxGetM(prhs(4))
    n1 = mxGetM(prhs(12))
    n2 = mxGetM(prhs(13))

	allocate(Gu(nu,2),Pu(nu,nu),invT1(r1,r1),invT2(r2,r2))
	allocate(y1mat(nu,n1,n2),y2mat(nu,n1,n2),p1mat(nu,n1,n2),p2mat(nu,n1,n2),v1mat(nu,n1,n2),v2mat(nu,n1,n2))
    allocate(knots1(n1),knots2(n2))

    call mxCopyPtrToReal8(mxGetPr(prhs(1)),Gu,nu*2)
    call mxCopyPtrToReal8(mxGetPr(prhs(2)),Pu,nu*nu)
    call mxCopyPtrToReal8(mxGetPr(prhs(3)),invT1,r1*r1)
    call mxCopyPtrToReal8(mxGetPr(prhs(4)),invT2,r2*r2)
    call mxCopyPtrToReal8(mxGetPr(prhs(12)),knots1,n1)
    call mxCopyPtrToReal8(mxGetPr(prhs(13)),knots2,n2)

    call mxCopyPtrToReal8(mxGetPr(prhs(5)),bet,1)
    call mxCopyPtrToReal8(mxGetPr(prhs(6)),eta,1)
    call mxCopyPtrToReal8(mxGetPr(prhs(7)),rho,1)
    call mxCopyPtrToReal8(mxGetPr(prhs(8)),sig,1)
    call mxCopyPtrToReal8(mxGetPr(prhs(9)),kap1,1)
    call mxCopyPtrToReal8(mxGetPr(prhs(10)),kap2,1)
    call mxCopyPtrToReal8(mxGetPr(prhs(11)),iota,1)

!    write(line,*) bet
!    k = mexPrintf(line//achar(10))

	call calccp(y1mat,y2mat,p1mat,p2mat,v1mat,v2mat,Gu,Pu,invT1,invT2)

	dims(1) = nu
	dims(2) = n1
	dims(3) = n2

    plhs(1) = mxCreateNumericArray(3, dims, mxClassIDFromClassName('double'), 0)
    plhs(2) = mxCreateNumericArray(3, dims, mxClassIDFromClassName('double'), 0)
    plhs(3) = mxCreateNumericArray(3, dims, mxClassIDFromClassName('double'), 0)
    plhs(4) = mxCreateNumericArray(3, dims, mxClassIDFromClassName('double'), 0)
    plhs(5) = mxCreateNumericArray(3, dims, mxClassIDFromClassName('double'), 0)
    plhs(6) = mxCreateNumericArray(3, dims, mxClassIDFromClassName('double'), 0)

    call mxCopyReal8ToPtr(y1mat,mxGetPr(plhs(1)),nu*n1*n2)
    call mxCopyReal8ToPtr(y2mat,mxGetPr(plhs(2)),nu*n1*n2)
    call mxCopyReal8ToPtr(p1mat,mxGetPr(plhs(3)),nu*n1*n2)
    call mxCopyReal8ToPtr(p2mat,mxGetPr(plhs(4)),nu*n1*n2)
    call mxCopyReal8ToPtr(v1mat,mxGetPr(plhs(5)),nu*n1*n2)
    call mxCopyReal8ToPtr(v2mat,mxGetPr(plhs(6)),nu*n1*n2)

    deallocate(knots1,knots2)

    return

end subroutine mexFunction
