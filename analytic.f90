! =================================================
!   Analytical solution.
! =================================================

subroutine analytic

    use params 
    implicit none

    character(len=:), allocatable :: dir
    real*8,           parameter   :: PI = 4*atan(1.0_8)
    real*8                        :: L
    integer*4                     :: S, Smax
    real*8                        :: A2, B2, A3, B3
    real*8                        :: T_1(Ny,Nx)
    real*8                        :: T_2(Ny,Nx)
    real*8                        :: T_3(Ny,Nx)
    real*8                        :: x_tmp2(Ny,Nx)
    real*8                        :: y_tmp2(Ny,Nx)
    real*8                        :: x_tmp3(Ny,Nx)
    real*8                        :: y_tmp3(Ny,Nx)

    integer*4 :: I,J

    ! Initialization
    call init()

    ! Seperate variable method
    L    = .5*(Lx + Ly)
    T_1 = L/Ly*y_

    Smax = 10
    T_2 = 0.0_8
    T_3 = 0.0_8
    do S = 1, Smax, 2
        x_tmp2 = S*PI*L/Ly*x_
        y_tmp2 = S*PI*L/Ly*y_
        A2  = -2*Ly**2*qw_/(S*PI*L)**2
        B2  = -A2/tanh(S*PI*Lx/Ly)
        T_2 = T_2 + sin(y_tmp2)* & 
            & (A2*sinh(x_tmp2) + B2*cosh(x_tmp2))
        x_tmp3 = S*PI*L/Ly*(Lx/L - x_)
        y_tmp3 = S*PI*L/Ly*y_
        A3  = -2*Ly**2*qe_/(S*PI*L)**2
        B3  = -A3/tanh(S*PI*Lx/Ly)
        T_3 = T_3 + sin(y_tmp3)* & 
            & (A2*sinh(x_tmp3) + B2*cosh(x_tmp3))
    end do

    T_ = T_1 + T_2 + T_3
    T  = T_*(Tn - Ts) + Ts

    ! Visualization
    dir  = ".\\out\\analytic\\"
    call output(dir, 0, "T", T)

end subroutine