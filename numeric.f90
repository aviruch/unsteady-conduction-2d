! =================================================
!   Analytical solution, including explicit method
!     and implicit method
! =================================================

module numeric_solver

    use coordinates, only: Nx, Ny
    implicit none
    
    private
    
    real*8 :: T(Nx,Ny)
    real*8 :: T_(Nx,Ny), T0_(Nx,Ny)

    character(len=:), allocatable :: dir
    real*8                        :: tic, toc
    integer*4                     :: iter = 0
    
    real*8 :: Alpx, Alpy
    real*8 :: tau_ = 0
    real*8 :: steady_err
    
    public :: iterate

contains

    subroutine init()
    
        use physical_def, only: Ly, Tn, Ts
        use coordinates,  only: y

        ! Initalize temperature: linear distributed
        T  = (Tn - Ts)*y/Ly + Ts
        T_ = (T - Ts)/(Tn - Ts)
        
        ! Check if output directory exists
        call check_dir(dir)
        call save_tecplot_file(dir, iter, "Temperature", T)
    
    end subroutine init
    
    subroutine update_ex()
    
        use physical_def, only: qw_, qe_
        use coordinates,  only: Nx, Ny, dx_, dy_
        use time_config,  only: tau_step
        integer*4 :: I, J
        
        iter = iter + 1
        tau_ = tau_ + tau_step
        
        ! Interior
        do I = 2, Nx - 1
            do J = 2, Ny - 1
                T_(I,J) = (1 - 2*Alpx - 2*Alpy)*T0_(I,J)     &
                    &   + Alpx*T0_(I-1,J) + Alpx*T0_(I+1,J)  &
                    &   + Alpy*T0_(I,J-1) + Alpy*T0_(I,J+1)
            end do
        end do
        
        ! South
        do I = 2, Nx - 1
            T_(I,1)  = (1 - 2*Alpx - 3*Alpy)*T0_(I,1)        &
                &    + Alpx*T0_(I-1,1) + Alpx*T0_(I+1,1)     &
                &    + Alpy*T0_(I,2)
        end do

        ! North
        do I = 2, Nx - 1
            T_(I,Ny) = (1 - 2*Alpx - 3*Alpy)*T0_(I,Ny)       &
                &    + Alpx*T0_(I-1,Ny) + Alpx*T0_(I+1,Ny)   &
                &    + Alpy*T0_(I,Ny-1) + 2*Alpy
        end do

        ! West
        do J = 2, Ny - 1
            T_(1,J)  = (1 - Alpx - 2*Alpy)*T0_(1,J)          &
                &    + Alpx*dx_*qw_ + Alpx*T0_(2,J)          &
                &    + Alpy*T0_(1,J-1) + Alpy*T0_(1,J+1)
        end do

        ! East
        do J = 2, Ny - 1
            T_(Nx,J) = (1 - Alpx - 2*Alpy)*T0_(Nx,J)         &
                &    + Alpx*T0_(Nx-1,J) + Alpx*dx_*qe_       &
                &    + Alpy*T0_(Nx,J-1) + Alpy*T0_(Nx,J+1)
        end do
    
        ! Southwest
        T_(1,1)   = (1 - Alpx - 3*Alpy)*T0_(1,1)             &
            &     + Alpx*dx_*qw_ + Alpx*T0_(2,1)             &
            &     + Alpy*T0_(1,2)
        
        ! Northwest
        T_(1,Ny)  = (1 - Alpx - 3*Alpy)*T0_(1,Ny)            &
            &     + Alpx*dx_*qw_ + Alpx*T0_(2,Ny)            &
            &     + Alpy*T0_(1,Ny-1) + 2*Alpy
        
        ! Southeast
        T_(Nx,1)  = (1 - Alpx - 3*Alpy)*T0_(Nx,1)            &
            &     + Alpx*T0_(Nx-1,1) + Alpx*dx_*qe_          &
            &     + Alpy*T0_(Nx,2)
            
        ! Northeast
        T_(Nx,Ny) = (1 - Alpx - 3*Alpy)*T0_(Nx,Ny)           &
            &     + Alpx*T0_(Nx-1,Ny) + Alpx*dx_*qe_         &
            &     + Alpy*T0_(Nx,Ny-1) + 2*Alpy

    end subroutine update_ex
    
    subroutine update_im(mode)
    
        use physical_def, only: qw_, qe_
        use coordinates,  only: Nx, Ny, dx_, dy_
        use time_config,  only: tau_step
        use conv_config,  only: conv_cr, conv_max
        
        character(len=*), intent(in) :: mode
        
        real*8    :: A(Nx*Ny,Nx*Ny)
        real*8    :: b(Nx*Ny)
        real*8    :: T_vec(Nx*Ny)
        integer*4 :: I, J, P
        
        iter = iter + 1
        tau_ = tau_ + tau_step
        
        A = 0.
        b = 0.
        
        ! Interior
        do I = 2, Nx - 1
            do J = 2, Ny - 1
                P         = (I - 1)*Ny + J
                A(P,P)    = 1 + 2*Alpx + 2*Alpy
                A(P,P-Ny) = -Alpx
                A(P,P+Ny) = -Alpx
                A(P,P-1)  = -Alpy
                A(P,P+1)  = -Alpy
                b(P)      = T0_(I,J)
            end do
        end do
        
        ! South
        do I = 2, Nx - 1
            P         = (I - 1)*Ny + 1
            A(P,P)    = 1 + 2*Alpx + 3*Alpy
            A(P,P-Ny) = -Alpx
            A(P,P+Ny) = -Alpx
            A(P,P+1)  = -Alpy
            b(P)      = T0_(I,1)
        end do

        ! North
        do I = 2, Nx - 1
            P         = I*Ny
            A(P,P)    = 1 + 2*Alpx + 3*Alpy
            A(P,P-Ny) = -Alpx
            A(P,P+Ny) = -Alpx
            A(P,P-1)  = -Alpy
            b(P)      = T0_(I,Ny) + 2*Alpy
        end do

        ! West
        do J = 2, Ny - 1
            P         = J
            A(P,P)    = 1 + Alpx + 2*Alpy
            A(P,P+Ny) = -Alpx
            A(P,P-1)  = -Alpy
            A(P,P+1)  = -Alpy
            b(P)      = T0_(1,J) + Alpx*dx_*qw_
        end do

        ! East
        do J = 2, Ny - 1
            P         = (Nx - 1)*Ny + J
            A(P,P)    = 1 + Alpx + 2*Alpy
            A(P,P-Ny) = -Alpx
            A(P,P-1)  = -Alpy
            A(P,P+1)  = -Alpy
            b(P)      = T0_(Nx,J) + Alpx*dx_*qe_
        end do

        ! Southwest
        P         = 1
        A(P,P)    = 1 + Alpx + 3*Alpy
        A(P,P+Ny) = -Alpx
        A(P,P+1)  = -Alpy
        b(P)      = T0_(1,1) + Alpx*dx_*qw_

        ! Northwest
        P         = Ny
        A(P,P)    = 1 + Alpx + 3*Alpy
        A(P,P+Ny) = -Alpx
        A(P,P-1)  = -Alpy
        b(P)      = T0_(1,Ny) + Alpx*dx_*qw_ + 2*Alpy

        ! Southeast
        P         = (Nx - 1)*Ny + 1
        A(P,P)    = 1 + Alpx + 3*Alpy
        A(P,P-Ny) = -Alpx
        A(P,P+1)  = -Alpy
        b(P)      = T0_(Nx,1) + Alpx*dx_*qe_

        ! Northeast
        P         = Nx*Ny
        A(P,P)    = 1 + Alpx + 3*Alpy
        A(P,P-Ny) = -Alpx
        A(P,P-1)  = -Alpy
        b(P)      = T0_(Nx,Ny) + Alpx*dx_*qe_ + 2*Alpy
        
        ! Solute linear equation
        select case (mode)
            case ("iter_ja_p")
                call iter_ja_p(A, b, T_vec, conv_cr, conv_max)
            case ("iter_gs_p")
                call iter_gs_p(A, b, T_vec, conv_cr, conv_max)
            case ("iter_ja_l")
                call iter_ja_l(A, b, T_vec, conv_cr, conv_max)
        end select
        
        ! Vector -> Matrix
        P = 0
        do I = 1, Nx
            do J = 1, Ny
                P = P + 1
                T_(I,J) = T_vec(P)
            end do
        end do
        
    end subroutine update_im
    
    subroutine check_steady()
    
        use conv_config, only: steady_cr
        use log_config,  only: log_freq
        
        call max_deviation(T_, T0_, steady_err)
        if (mod(iter, log_freq) == 0 .or. steady_err < steady_cr) then
            print *, "Iter:", iter, "Err:", steady_err
        end if
    
    end subroutine check_steady

    subroutine iterate_ex()
    
        use physical_def, only: Tn, Ts
        use time_config,  only: tau_max
        use conv_config,  only: steady_cr
        
        if (Alpx + Alpy > 0.5) then 
            print *, "Error: diffusion number exceeds:"
            print *, "Alpha=", Alpx + Alpy
            return
        end if
        
        call init()
        T0_ = T_
        
        do while (tau_ < tau_max)
            call update_ex()
            T = T_*(Tn - Ts) + Ts
            call save_tecplot_file(dir, iter, "Temperature", T)
            call check_steady()
            if (steady_err < steady_cr) exit
            T0_ = T_
        end do
        
    end subroutine iterate_ex
    
    subroutine iterate_im(mode)
    
        use physical_def, only: Tn, Ts
        use time_config,  only: tau_max
        use conv_config,  only: steady_cr
    
        character(len=*), intent(in) :: mode
        
        call init()
        T0_ = T_
        
        do while (tau_ < tau_max)
            call update_im(mode)
            T = T_*(Tn - Ts) + Ts
            call save_tecplot_file(dir, iter, "Temperature", T)
            call check_steady()
            if (steady_err < steady_cr) exit
            T0_ = T_
        end do
    
    end subroutine iterate_im
    
    subroutine iterate(res, res_, method)
    
        use coordinates, only: Nx, Ny, dx_, dy_
        use time_config, only: tau_step
        
        integer*4, intent(in)  :: method
        real*8,    intent(out) :: res(Nx,Ny), res_(Nx,Ny)
        
        call cpu_time(tic)
        
        Alpx = tau_step/dx_**2
        Alpy = tau_step/dy_**2
        
        select case (method)
            case (0)
                print *, "Method: explicit."
                dir = ".\\out\\explicit\\"
                call iterate_ex()
            case (1)
                print *, "Method: implicit, Jacobi."
                dir = ".\\out\\implicit_ja_p\\"
                call iterate_im("iter_ja_p")
            case (2)
                print *, "Method: implicit, Gauss-Seidel."
                dir = ".\\out\\implicit_gs_p\\"
                call iterate_im("iter_gs_p")
            case (3)
                print *, "Method: implicit, Jacobi (block) + TDMA."
                dir = ".\\out\\implicit_ja_l\\"
                call iterate_im("iter_ja_l")
            case default
                return
        end select
        
        call cpu_time(toc)
        print *, "Converged!"
        print *, "Time consumed:", toc - tic, "s."
        
        res  = T
        res_ = T_
    
    end subroutine iterate
    
end module
