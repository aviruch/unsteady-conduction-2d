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

    ! Diffusion number
    a_x = tau_step/dx_**2
    a_y = tau_step/dy_**2

    ! Time Marching
    tau_ = 0
    T0_  = T_
    do while (tau_ < tau_max)

        iter = iter + 1
        tau_ = tau_ + tau_step

        ! Iteration
        call vec1d(T0_, T_1d)
        call getAb(A, b, T_1d, a_x, a_y)
        call iter_ja_p(Nx*Ny, A, b, T_1d, &
            & conv_cr, conv_max)
        call mat2d(T_1d, T_)

        ! do J = 1, Nx*Ny
        !     write(*,*) A(:,J)
        ! end do
        ! write(*,*) ""

        ! do J = 1, Ny
        !     write(*,*) T_(:,J)
        ! end do

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

        ! stop

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
!   Vectorize 2d array to 1d array.
! =================================================

subroutine vec1d(x2d, x1d)

    use params, only: Nx, Ny
    implicit none

    real*8, intent(in)  :: x2d(Nx,Ny)
    real*8, intent(out) :: x1d(Nx*Ny)
    integer*4           :: I, J, P

    P = 0
    do I = 1, Nx
        do J = 1, Ny
            P = P + 1
            x1d(P) = x2d(I,J) 
        end do
    end do

end subroutine


! =================================================
!   Turn vectorized 1d array back to 2d array.
! =================================================

subroutine mat2d(x1d, x2d)

    use params, only: Nx, Ny
    implicit none

    real*8, intent(in)  :: x1d(Nx*Ny)
    real*8, intent(out) :: x2d(Nx,Ny)
    integer*4           :: I, J, P

    P = 0
    do I = 1, Nx
        do J = 1, Ny
            P = P + 1
            x2d(I,J) = x1d(P)
        end do
    end do

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
    A = 0.
    b = 0.

    ! Boundary conditions
    P = 1
    A(P,P)    = 1 + a_x + 3*a_y
    A(P,P+Ny) = -a_x
    A(P,P+1)  = -a_y
    b(P)      = x0(1,1) + a_x*dx_*qw_

    P = Ny
    A(P,P)    = 1 + a_x + 3*a_y
    A(P,P+Ny) = -a_x
    A(P,P-1)  = -a_y
    b(P)      = x0(1,Ny) + a_x*dx_*qw_ + 2*a_y

    P = (Nx - 1)*Ny + 1
    A(P,P)    = 1 + a_x + 3*a_y
    A(P,P-Ny) = -a_x
    A(P,P+1)  = -a_y
    b(P)      = x0(Nx,1) + a_x*dx_*qe_

    P = Nx*Ny
    A(P,P)    = 1 + a_x + 3*a_y
    A(P,P-Ny) = -a_x
    A(P,P-1)  = -a_y
    b(P)      = x0(Nx,Ny) + a_x*dx_*qe_ + 2*a_y

    J = 1
    do I = 2, Nx - 1
        P = (I - 1)*Ny + J
        A(P,P)    = 1 + 2*a_x + 3*a_y
        A(P,P-Ny) = -a_x
        A(P,P+Ny) = -a_x
        A(P,P-1)  = -a_y
        b(P)      = x0(Nx,Ny)
    end do

    J = Ny
    do I = 2, Nx - 1
        P = (I - 1)*Ny + J
        A(P,P)    = 1 + 2*a_x + 3*a_y
        A(P,P-Ny) = -a_x
        A(P,P+Ny) = -a_x
        A(P,P+1)  = -a_y
        b(P)      = x0(Nx,Ny) + 2*a_y
    end do

    I = 1
    do J = 2, Ny - 1
        P = (I - 1)*Ny + J
        A(P,P)    = 1 + a_x + 2*a_y
        A(P,P+Ny) = -a_x
        A(P,P-1)  = -a_y
        A(P,P+1)  = -a_y
        b(P)      = x0(1,J) + a_x*dx_*qw_
    end do

    I = Nx
    do J = 2, Ny - 1
        P = (I - 1)*Ny + J
        A(P,P)    = 1 + a_x + 2*a_y
        A(P,P-Ny) = -a_x
        A(P,P-1)  = -a_y
        A(P,P+1)  = -a_y
        b(P)      = x0(Nx,J) + a_x*dx_*qe_
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
            b(P)      = x0(I,J)
        end do
    end do

end subroutine


! =================================================
!   Jocobi iteration (point) method to solve linear 
!       equation system 'Ax = b'.
! =================================================

subroutine iter_ja_p(n, A, b, x, iter_cr, iter_max)

    implicit none

    integer*4, intent(in)    :: n
    real*8,    intent(in)    :: A(n,n)
    real*8,    intent(in)    :: b(n)
    real*8,    intent(inout) :: x(n)
    real*8,    intent(in)    :: iter_cr
    integer*4, intent(in)    :: iter_max

    real*8    :: x0(n)
    real*8    :: iter_err, tmp
    integer*4 :: iter_cnt
    integer*4 :: I, J

    x0 = x
    iter_cnt = 1
    do while (iter_cnt < iter_max)

        do I = 1, n
            tmp = 0.
            do J = 1, n
                if (J /= I) then
                    tmp = tmp + A(I,J)*x0(J)
                end if
            end do
            x(I) = (b(I) - tmp)/A(I,I)
        end do

        iter_err = maxval(abs(x - x0)/x0)
        write(*,*) iter_cnt, iter_err
        if (iter_err < iter_cr) exit

        x0 = x
        iter_cnt = iter_cnt + 1

    end do

end subroutine