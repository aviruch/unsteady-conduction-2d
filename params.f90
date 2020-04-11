! =================================================
!   Define global variables.
!   Note: All variables with suffix "_" represent 
!     non-dimensional scalars.
! =================================================

! Pysical parameters
module physical_def

    real*8, parameter :: Lx  = 3.0
    real*8, parameter :: Ly  = 4.5
    real*8, parameter :: rho = 7820.
    real*8, parameter :: c   = 460.
    real*8, parameter :: lbd = 15.
    real*8, parameter :: Ts  = 500.
    real*8, parameter :: Tn  = 300.
    real*8, parameter :: qw  = 800.
    real*8, parameter :: qe  = 800.

    real*8 :: Lx_, Ly_
    real*8 :: qw_, qe_
    
end module physical_def


! Mesh settings 
module coordinates
    
    integer*4, parameter :: Nx = 10
    integer*4, parameter :: Ny = 15
    
    real*8                   :: dx,  dy
    real*8                   :: dx_, dy_
    real*8, dimension(Nx,Ny) :: x,  y
    real*8, dimension(Nx,Ny) :: x_, y_

end module coordinates


! Time-marching settings
module time_config

    real*8, parameter :: tau_step = 0.001
    real*8, parameter :: tau_max  = 1000.

end module time_config


! Convegence settings
module conv_config

    real*8               :: omega     = 1.
    real*8,    parameter :: steady_cr = 1e-4
    real*8,    parameter :: conv_cr   = 1e-4
    integer*4, parameter :: conv_max  = 100

end module conv_config


! Output settings
module log_config
    
    integer*4, parameter :: log_freq = 20
    logical,   parameter :: log_conv_show = .false.

end module log_config
