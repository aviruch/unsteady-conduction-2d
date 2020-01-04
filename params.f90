module params

    ! Pysical Parameters
    real(kind=8), parameter :: H = 3.0
    real(kind=8), parameter :: rho = 7820
    real(kind=8), parameter :: c = 460
    real(kind=8), parameter :: lambda = 15
    real(kind=8), parameter :: T0 = 250.0
    real(kind=8), parameter :: T1 = 400.0
    real(kind=8), parameter :: qw = 750.0
    real(kind=8), parameter :: qe = 750.0

    ! Mesh Settings
    integer(kind=4), parameter :: N = 10

    ! Dimensionless Scalars
    real(kind=8) :: X_s, theta_s, tau_s, q_s

    ! Iteration Settings
    real(kind=8) :: alpha, omega
    real(kind=8) :: dX
    real(kind=8) :: taumax
    real(kind=8) :: dtau
    real(kind=8), parameter :: steady_cr = 1e-4
    real(kind=8), parameter :: conv_cr = 1e-4

    ! Iterative Variables
    real(kind=8), allocatable :: x(:), y(:)
    real(kind=8), allocatable :: T(:)
    real(kind=8), allocatable :: theta(:), theta0(:)
    
end module