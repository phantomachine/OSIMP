module mod_parameters

implicit none

real(8), save :: bet, eta, rho, sig, kap1, kap2, neu, sig1, sig2
integer, save :: nu, n1, n2, nx, r1, r2, rx
real(8), save, allocatable :: knots1(:), knots2(:), knotsx(:)

end module mod_parameters
