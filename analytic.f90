! =================================================
!   Analytical solution.
! =================================================

module analytic_solver

    implicit none
    public :: calc_analytic
    public :: comp_analytic
    
contains

    subroutine calc_analytic(T, T_)

        use physical_def, only: Lx, Ly, Tn, Ts, Lx_, Ly_, qw_, qe_
        use coordinates,  only: Nx, Ny, x_, y_
        use conv_config,  only: conv_cr, conv_max

        real*8, intent(out) :: T(Nx,Ny), T_(Nx,Ny)

        character(len=:), allocatable :: dir
        real*8,           parameter   :: PI = 4*atan(1.0_8)
        integer*4                     :: n, INF = conv_max
        real*8                        :: T0_(Nx,Ny), err
        
        ! Series summation
        T0_ = y_/Ly_
        do n = 1, INF
            T_ = T0_ - 2*(-1 + (-1)**n)*Ly_/(n*PI)**2                    &
            &   *(qe_*cosh(n*PI/Ly_*x_) + qw_*cosh(n*PI/Ly_*(Lx_ - x_))) &
            &   /sinh(n*PI*Lx_/Ly_)*sin(n*PI/Ly_*y_)
            call max_deviation(T_, T0_, err)
            if (err > 0 .and. err < conv_cr) exit
            T0_ = T_
        end do
        T = T_*(Tn - Ts) + Ts

        ! Visualization
        dir = ".\\out\\analytic\\"
        call check_dir(dir)
        call save_tecplot_file(dir, 0, "Temperature", T)
        
    end subroutine calc_analytic
    
    subroutine comp_analytic(mat_n, mat_a, avg, std)

        use coordinates, only: Nx, Ny
        
        real*8, intent(in)  :: mat_n(Nx,Ny)
        real*8, intent(in)  :: mat_a(Nx,Ny)
        real*8, intent(out) :: avg, std
        
        real*8, allocatable :: vec_n(:)
        real*8, allocatable :: vec_a(:)
        real*8, allocatable :: err(:)
        
        vec_n = reshape(mat_n, (/Nx*Ny/))
        vec_a = reshape(mat_a, (/Nx*Ny/))
        err = (vec_n - vec_a)/vec_a
        avg = sum(err)/(Nx*Ny)
        std = sqrt(sum((err - avg)**2))/(Nx*Ny)
    
    end subroutine comp_analytic
    
    
end module 
