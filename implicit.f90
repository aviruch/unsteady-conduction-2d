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
    integer*4                     :: I, J, P

    ! Initialization
    call cpu_time(tic)
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
        call getAb(A, b, T0_, a_x, a_y)
        ! call iter_ja_p(Nx*Ny, A, b, T_1d, &
        !     & conv_cr, conv_max)
        call iter_gs_p(Nx*Ny, A, b, T_1d, &
            & conv_cr, conv_max)

        ! Vector -> Matrix
        P = 0
        do I = 1, Nx
            do J = 1, Ny
                P = P + 1
                T_(I,J) = T_1d(P)
            end do
        end do

        ! Dimensionalize & output
        T = T_*(Tn - Ts) + Ts
        call output(dir, iter, "T", T)

        ! Steady criteria
        call maxerr(Nx, Ny, T_, T0_, steady_err)
        write(*, "(A5, I4, A8, E12.4E2)") &
            & "Iter:", iter,  "Err:", steady_err
        if (steady_err < steady_cr) exit

        ! Update
        T0_ = T_

    end do

    ! Print time consumption
    write(*, "(A20)") ""
    write(*, "(A20)") "Converged!          "
    call cpu_time(toc)
    write(*, "(A20, F8.3, A2)") &
        & "Time consumed:      ", toc - tic, "s."
    write(*, "(A20)") ""

    ! Print numerical solution
    write(*, "(A20)") "Numerical solution: "
    call print(T_)
    write(*, "(A20)") ""

    ! Print analytical solution
    write(*, "(A20)") "Analytical solution:"
    call analytic(T0)
    call cmperr(Nx, Ny, T, T0, avg, std)
    call print((T0 - Ts)/(Tn - Ts))
    write(*, "(A20)") ""
    write(*, "(A20, E12.4E2)") "Average error:      ", avg
    write(*, "(A20, E12.4E2)") "Standard deviation: ", std

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
        b(P)      = x0(I,1)
    end do

    J = Ny
    do I = 2, Nx - 1
        P = (I - 1)*Ny + J
        A(P,P)    = 1 + 2*a_x + 3*a_y
        A(P,P-Ny) = -a_x
        A(P,P+Ny) = -a_x
        A(P,P+1)  = -a_y
        b(P)      = x0(I,Ny) + 2*a_y
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
