! =================================================
!   Module to set global variables.
! =================================================

module params

    ! Pysical parameters
    real*8, parameter :: Lx  = 3.
    real*8, parameter :: Ly  = 5.
    real*8, parameter :: rho = 7820.
    real*8, parameter :: c   = 460.
    real*8, parameter :: lbd = 15.
    real*8, parameter :: Ts  = 500.
    real*8, parameter :: Tn  = 300.
    real*8, parameter :: qw  = 800.
    real*8, parameter :: qe  = 800.

    ! Mesh settings 
    integer*4, parameter :: Nx = 30
    integer*4, parameter :: Ny = 50

    ! Time, space and temperature
    real*8                   :: dx, dy
    real*8, dimension(Ny,Nx) :: x, y
    real*8, dimension(Ny,Nx) :: T

    ! Dimensionless scalars
    real*8                   :: dx_, dy_
    real*8, dimension(Ny,Nx) :: x_, y_
    real*8, dimension(Ny,Nx) :: T_
    real*8                   :: qw_, qe_

    ! Time-marching settings
    real*8            :: alpha    = 0.25
    real*8            :: tau_
    real*8, parameter :: tau_step = 1.
    real*8, parameter :: tau_max  = 1000.

    ! Convegence settings
    real*8            :: omega     = 1.
    real*8, parameter :: steady_cr = 1e-4
    real*8, parameter :: conv_cr   = 1e-4
    
end module