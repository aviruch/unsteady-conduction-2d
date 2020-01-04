subroutine init

    use params
    implicit none
    integer :: I, J, K

    ! Dimensionless Scalars
    X_s = H
    theta_s = T1 - T0
    tau_s = rho*c*H**2/lambda
    q_s = H/(lambda*theta_s)

    ! Nodes Coordinates
    dX = 1./N
    allocate(x(N*N), y(N*N))
    K = 0
    do I = 1, N
        x(K+1:K+N) = H/N*(I - 0.5)
        do J = 1, N
            K = K + 1
            y(K) = H/N*(N - J + 0.5)
        end do
    end do

    ! Temperature Initialization
    allocate(T(N*N), theta(N*N))
    T = T0 + (T1 - T0)/H*x
    T(1:N) = (2*T0 + T(N+1:2*N))/3
    T((N-1)*N+1:N*N) = (2*T1 + T((N-2)*N+1:(N-1)*N))/3
    theta0 = (T - T0)/(T1 - T0)
    theta = theta0

end subroutine