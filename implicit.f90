! =================================================
!   Implicit method.
! =================================================

subroutine implicit

    use params 
    implicit none

    character(len=:), allocatable :: dir
    real*8                        :: tic, toc
    integer*4                     :: iter
    real*8                        :: T0(Nx,Ny)
    real*8                        :: T0_(Nx,Ny)
    real*8                        :: a_x, a_y
    real*8                        :: A(Nx*Ny,Nx*Ny)
    real*8                        :: b(Nx*Ny)
    real*8                        :: T_1d(Nx*Ny)
    real*8                        :: steady_err
    real*8                        :: avg, std
    integer*4                     :: I, J

    write(*, *) "Implicit Solver Starts!"
    call cpu_time(tic)

    ! Initialization
    call init()

    ! Check if output directory exists
    iter = 0
    dir  = ".\\out\\implicit\\"
    call checkdir(dir)
    call output(dir, iter, "T", T)

    ! Time Marching
    tau_ = 0
    T0_  = T_
    do while (tau_ < tau_max)

        iter = iter + 1
        tau_ = tau_ + tau_step

        ! Iteration
        T_1d = reshape(T0_, (/Nx*Ny/))
        write(*,*) T_1d
        call getAb(A, b, T_1d, a_x, a_y)
        call iter_ja_p(A, b, T_1d, conv_cr, conv_max)
        T_ = reshape(T_1d, (/Nx,Ny/))
        
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
    write(*, *) "Implicit Solver Ends!"
    write(*, *) "Time Consumed:     ", toc - tic, "s."

    call analytic(T0)
    call cmperr(Nx, Ny, T, T0, avg, std)
    write(*, *) "Average error:     ", avg
    write(*, *) "Standard deviation:", std

end subroutine


! =================================================
!   Initialize matrix 'A' and vector 'b' in linear 
!       equation system 'Ax = b'.
! =================================================

subroutine getAb(A, b, x0, a_x, a_y)

    use params, only: Nx, Ny, dx_, dy_, qw_, qe_
    implicit none

    real*8, intent(out) :: A(Nx*Ny,Nx*Ny)
    real*8, intent(out) :: b(Nx*Ny)
    real*8, intent(in)  :: x0(Nx,Ny)
    real*8, intent(in)  :: a_x, a_y
    integer*4           :: I, J, P

    ! Global initialization
    A = 0
    b = 0

    ! Boundary conditions
    P = 1
    A(P,P)    = 1 + 2*a_x + 2*a_y
    A(P,P+Ny) = -a_x
    A(P,P+1)  = -a_y
    b(P)      = x0(1,1) + a_x*dx_*qw_

    P = Ny
    A(P,P)    = 1 + 2*a_x + 2*a_y
    A(P,P+Ny) = -a_x
    A(P,P-1)  = -a_y
    b(P)      = x0(1,Ny) + a_x*dx_*qw_ + 2*a_y

    P = (Nx - 1)*Ny + 1
    A(P,P)    = 1 + 2*a_x + 2*a_y
    A(P,P-Ny) = -a_x
    A(P,P+1)  = -a_y
    b(P)      = x0(Nx,1) + a_x*dx_*qe_

    P = Nx*Ny
    A(P,P)    = 1 + 2*a_x + 2*a_y
    A(P,P-Ny) = -a_x
    A(P,P-1)  = -a_y
    b(P)      = x0(Nx,Ny) + a_x*dx_*qe_ + 2*a_y

    do I = 2, Nx - 1
        P = (I - 1)*Ny + 1
        A(P,P)   = 1.
        A(P,P+1) = -1./3
        b(P)     = 0.
        P = I*Ny
        A(P,P)   = 1.
        A(P,P-1) = -1./3
        b(P)     = 2./3
    end do

    do J = 2, Ny - 1
        P = J
        A(P,P)    = 1 + a_x + 2*a_y
        A(P,P+Ny) = -a_x
        A(P,P-1)  = -a_y
        A(P,P+1)  = -a_y
        b(P)      = x0(1, J) + 0.5*a_x*dx_*qw_
        P = (Nx - 1)*Ny + J
        A(P,P)    = 1 + a_x + 2*a_y
        A(P,P-Ny) = -a_x
        A(P,P-1)  = -a_y
        A(P,P+1)  = -a_y
        b(P)      = x0(Nx, J) + 0.5*a_x*dx_*qe_
    end do
    
    ! Interior
    do I = 2, Nx - 1
        do J = 2, Ny - 1
            P = (I - 1)*Ny + J
            A(P,P)    = 1 + 2*a_x + 2*a_y
            A(P,P-Ny) = -a_x
            A(P,P+Ny) = -a_x
            A(P,P-1)  = -a_y
            A(P,P+1)  = -a_y
            b(P)      = x0(I, J)
        end do
    end do

end subroutine


! =================================================
!   Jocobi iteration (point) method to solve linear 
!       equation system 'Ax = b'.
! =================================================

subroutine iter_ja_p(A, b, x, iter_cr, iter_max)

    use params, only: Nx, Ny
    implicit none

    real*8,    intent(in)    :: A(Nx*Ny,Nx*Ny)
    real*8,    intent(in)    :: b(Nx*Ny)
    real*8,    intent(inout) :: x(Nx*Ny)
    real*8,    intent(in)    :: iter_cr
    integer*4, intent(in)    :: iter_max

    real*8    :: x0(Nx*Ny)
    real*8    :: iter_err, tmp
    integer*4 :: iter_cnt
    integer*4 :: P, Q

    x0 = x
    iter_cnt = 0
    do while (iter_cnt < iter_max)

        iter_cnt = iter_cnt + 1

        do P = 1, Nx*Ny
            tmp = 0
            do Q = 1, Nx*Ny
                if (Q /= P) then
                    tmp = tmp + A(P,Q)*x0(Q)
                end if
            end do
            x(P) = (b(P) - tmp)/A(P,P)
        end do

        iter_err = maxval(abs(x - x0)/x0)
        write(*,*) iter_cnt, iter_err
        if (iter_err < iter_cr) exit

        x0 = x

    end do

end subroutine