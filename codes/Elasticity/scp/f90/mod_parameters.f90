module mod_parameters

implicit none

real(8), save :: bet, eta, rho, sig, kap1, kap2, iota, sig1, sig2
integer, save :: nu, n1, n2, nx, r1, r2, rx
real(8), save, allocatable :: knots1(:), knots2(:), knotsx(:)
! chi = 2*aH*(1-aH)*(1-rho*iota) = 0.5*(1-rho*iota);
! varpi = 1/(4*aH*(1-aH)*rho*iota + (2*aH-1)^2) = 1/rho/iota;
! real(8), parameter :: chi = 0.5d0*(1.0d0-rho*iota)
! real(8), parameter :: varpi = 1.0d0/rho/iota

end module mod_parameters
