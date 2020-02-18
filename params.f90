! =================================================
!   Module to set global variables.
! =================================================

module params

    ! Pysical parameters
    real*8, parameter :: Lx  = 3.0
    real*8, parameter :: Ly  = 4.5
    real*8, parameter :: rho = 7820.
    real*8, parameter :: c   = 460.
    real*8, parameter :: lbd = 15.
    real*8, parameter :: Ts  = 500.
    real*8, parameter :: Tn  = 300.
    real*8, parameter :: qw  = 800.
    real*8, parameter :: qe  = 800.

    ! Mesh settings 
    integer*4, parameter :: Nx = 10
    integer*4, parameter :: Ny = 15

    ! Time, space and temperature
    real*8                   :: dx, dy
    real*8, dimension(Nx,Ny) :: x, y
    real*8, dimension(Nx,Ny) :: T

    ! Dimensionless scalars
    real*8                   :: dx_, dy_
    real*8, dimension(Nx,Ny) :: x_, y_
    real*8, dimension(Nx,Ny) :: T_
    real*8                   :: qw_, qe_

    ! Time-marching settings
    real*8            :: tau_
    real*8, parameter :: tau_step = 0.001
    real*8, parameter :: tau_max  = 1000.

    ! Convegence settings
    real*8               :: omega     = 1.
    real*8,    parameter :: steady_cr = 1e-4
    real*8,    parameter :: conv_cr   = 1e-4
    integer*4, parameter :: conv_max  = 100
    
end module