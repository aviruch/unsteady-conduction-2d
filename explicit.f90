! =================================================
!   Explicit method.
! =================================================

subroutine explicit

    use params 
    implicit none

    character(len=:), allocatable :: dir
    real*8                        :: tic, toc
    integer*4                     :: iter
    real*8                        :: T0(Nx,Ny)
    real*8                        :: T0_(Nx,Ny)
    real*8                        :: a_x, a_y
    real*8                        :: steady_err
    real*8                        :: avg, std
    integer*4                     :: I, J

    write(*, *) "Explicit Solver Starts!"
    call cpu_time(tic)

    ! Initialization
    call init()

    ! Check if output directory exists
    iter = 0
    dir  = ".\\out\\explicit\\"
    call checkdir(dir)
    call output(dir, iter, "T", T)

    ! Diffusion number judgement
    a_x = tau_step/dx_**2
    a_y = tau_step/dy_**2
    if (a_x + a_y > 0.5) then
        write(*, *) "Error: diffusion number exceeds."
        write(*, *) "alpha=", a_x + a_y
        return
    end if

    ! Time Marching
    tau_ = 0
    T0_  = T_
    do while (tau_ < tau_max)

        iter = iter + 1
        tau_ = tau_ + tau_step
        
        ! Boundary conditions
        T_(1,1)   = (1 - 2*a_x - 2*a_y)*T0_(1,1)   &
            &     + a_x*dx_*qw_ + a_x*T0_(2,1)     &
            &     + a_y*T0_(1,2)
        T_(1,Ny)  = (1 - 2*a_x - 2*a_y)*T0_(1,Ny)  &
            &     + a_x*dx_*qw_ + a_x*T0_(2,Ny)    &
            &     + a_y*T0_(1,Ny-1) + a_y*2
        T_(Nx,1)  = (1 - 2*a_x - 2*a_y)*T0_(Nx,1)  &
            &     + a_x*T0_(Nx-1,1) + a_x*dx_*qe_  &
            &     + a_y*T0_(Nx,2)
        T_(Nx,Ny) = (1 - 2*a_x - 2*a_y)*T0_(Nx,Ny) &
            &     + a_x*T0_(Nx-1,Ny) + a_x*dx_*qe_ &
            &     + a_y*T0_(Nx,Ny-1) + a_y*2

        ! do I = 2, Nx - 1
        !     T_(I,1)  = T_(I,2)/3
        !     T_(I,Ny) = T_(I,Ny-1)/3 + 2/3
        ! end do

        ! write(*,*) T_(:,Ny)
        ! write(*,*) T_(:,Ny-1)
        ! write(*,*) ""

        do I = 2, Nx - 1
            T_(I,1)  = 0
            T_(I,Ny) = 1
        end do

        do J = 2, Ny - 1
            T_(1,J)  = (1 - a_x - 2*a_y)*T0_(1,J)      &
                &    + a_x*(0.5*dx_*qw_ + T0_(2,J))    &
                &    + a_y*(T0_(1,J-1) + T0_(1,J+1))
            T_(Nx,J) = (1 - a_x - 2*a_y)*T0_(Nx,J)     &
                &    + a_x*(T0_(Nx-1,J) + 0.5*dx_*qe_) &
                &    + a_y*(T0_(Nx,J-1) + T0_(Nx,J+1))
        end do

        ! Interior
        do I = 2, Nx - 1
            do J = 2, Ny - 1
                T_(I,J) = (1 - 2*a_x - 2*a_y)*T0_(I,J)  &
                    &   + a_x*(T0_(I-1,J) + T0_(I+1,J)) &
                    &   + a_y*(T0_(I,J-1) + T0_(I,J+1))
            end do
        end do

        ! Dimensionalize & output
        T = T_*(Tn - Ts) + Ts
        call output(dir, iter, "T", T)

        ! Steady criteria
        call maxerr(Nx, Ny, T_, T0_, steady_err)
        write(*, '(A8, I4, A8, E12.4E2)') &
            & "Iter:", iter,  "Err:", steady_err
        if (steady_err < steady_cr) exit

        ! Update
        T0_ = T_

    end do

    call cpu_time(toc)
    write(*, *) "Explicit Solver Ends!"
    write(*, *) "Time Consumed:     ", toc - tic, "s."

    call analytic(T0)
    call cmperr(Nx, Ny, T, T0, avg, std)
    write(*, *) "Average error:     ", avg
    write(*, *) "Standard deviation:", std

end subroutine