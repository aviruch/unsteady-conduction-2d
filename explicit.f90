subroutine explicit

    use params
    implicit none

    character(len=:), allocatable :: dir
    real(kind=8) :: tic, toc
    integer(kind=4) :: iter
    real(kind=8) :: time, tau
    real(kind=8) :: steady_err
    integer(kind=4) :: I, J, P

    write(*, '(A)') "Explicit Solver Starts!"
    call cpu_time(tic)

    ! Temperature Initialization
    call init()
    dir = ".\\out\\explicit\\"
    call checkdir(dir)
    call visualize(dir, 0)

    ! Time Marching
    tau = 0
    iter = 0
    dtau = alpha*dX**2
    taumax = dtau*1000
    do while (tau < taumax)

        iter = iter + 1
        tau = tau + dtau
        
        ! Top Boundary
        J = 1
        do I = 2, N - 1
            P = (I - 1)*N + J
            theta(P) = theta(P+1)/3
        end do

        ! Bottom Boundary
        J = N
        do I = 2, N - 1
            P = (I - 1)*N + J
            theta(P) = (theta(P-1) + 2)/3
        end do

        ! Left Boundary
        I = 1
        do J = 2, N - 1
            P = (I - 1)*N + J
            theta(P) = alpha*(0.5*q_s*qw*dX + theta0(P+N) + theta0(P-1) + theta0(P+1)) + (1 - 3*alpha)*theta0(P)
        end do

        ! Right Boundary
        I = N
        do J = 2, N - 1
            P = (I - 1)*N + J
            theta(P) = alpha*(0.5*q_s*qe*dX + theta0(P-N) + theta0(P-1) + theta0(P+1)) + (1 - 3*alpha)*theta0(P)
        end do
        
        ! Corners Boundary
        I = 1; J = 1
        P = (I - 1)*N + J
        theta(P) = alpha*(q_s*qw*dX + theta0(P+N) + theta0(P+1)) + (1 - 4*alpha)*theta0(P)
        I = 1; J = N
        P = (I - 1)*N + J
        theta(P) = alpha*(q_s*qw*dX + theta0(P+N) + theta0(P-1) + 2) + (1 - 4*alpha)*theta0(P)
        I = N; J = 1
        P = (I - 1)*N + J
        theta(P) = alpha*(q_s*qe*dX + theta0(P-N) + theta0(P+1)) + (1 - 4*alpha)*theta0(P)
        I = N; J = N
        P = (I - 1)*N + J
        theta(P) = alpha*(q_s*qe*dX + theta0(P-N) + theta0(P-1) + 2) + (1 - 4*alpha)*theta0(P)

        ! Interior
        do I = 2, N - 1
            do J = 2, N - 1
                P = (I - 1)*N + J
                theta(P) = alpha*(theta0(P-N) + theta0(P+N) + theta0(P-1) + theta0(P+1)) + (1 - 4*alpha)*theta0(P)
            end do
        end do

        ! Dimensionless -> Dimension
        time = tau*tau_s
        T = theta*theta_s + T0

        ! Temperature Visualization
        call visualize(dir, iter)
        
        ! Steady Criteria
        steady_err = maxval(abs(theta - theta0)/theta0)
        write(*, '(A8, I4, A8, E12.4E2, A8, E12.4E2)') "Iter:", iter, "Time:", time, "Err:", steady_err
        if (steady_err < steady_cr) exit
        theta0 = theta

    end do

    call cpu_time(toc)
    write(*, '(A)') "Explicit Solver Ends!"
    write(*, '(A, F12.8, A, /)') "Time Consumed:", toc - tic, "s."

    deallocate(x, y, t, theta)

end subroutine