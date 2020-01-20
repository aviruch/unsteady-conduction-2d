! =================================================
!   Explicit discrete method.
! =================================================

subroutine explicit

    use params 
    implicit none

    character(len=:), allocatable :: dir
    real*8                        :: tic, toc
    integer*4                     :: iter
    real*8                        :: T0_(Ny,Nx)
    real*8                        :: steady_err
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

    stop

    ! Time Marching
    tau_ = 0
    T0_  = T_
    do while (tau_ < tau_max)

        iter = iter + 1
        tau_ = tau_ + tau_step
        
        

        ! Dimensionalize & output
        T = T_*(Tn - Ts) + Ts
        call output(dir, iter, "T", T)
        
        ! Steady criteria
        call maxerr(T_, T0_, steady_err)
        write(*, '(A8, I4, A8, E12.4E2)') &
            & "Iter:", iter,  "Err:", steady_err
        if (steady_err < steady_cr) exit

        ! Update
        T0_ = T_

    end do

    call cpu_time(toc)
    write(*, *) "Explicit Solver Ends!"
    write(*, *) "Time Consumed:", toc - tic, "s."

end subroutine