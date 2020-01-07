subroutine analytic

    use params
    implicit none

    real*8, parameter      :: PI = 4*atan(1.0_8)
    integer*4              :: s, smax
    real*8, allocatable    :: A2s(:), B2s(:)
    real*8, allocatable    :: A3s(:), B3s(:)
    real*8, dimension(N,N) :: X2d, Y2d
    real*8, dimension(N,N) :: theta1, theta2, theta3
    
    character(len=:), allocatable :: dir
    integer*4                     :: I, J, P

    call init()

    dir  = ".\\out\\analytic\\"
    call checkdir(dir)

    smax = 100
    allocate(A2s(smax), B2s(smax), A3s(smax), B3s(smax))
    do s = 1, smax
        if (mod(s, 2) == 1) then
            A2s(s) = -2*H*qw/lambda/(T1- T0)/(s*PI)**2
            B2s(s) = -A2s(s)/tanh(s*PI)
            A3s(s) = -2*H*qe/lambda/(T1- T0)/(s*PI)**2
            B3s(s) = -A3s(s)/tanh(s*PI)
        else
            A2s(s) = 0
            B2s(s) = 0
            A3s(s) = 0
            B3s(s) = 0
        end if
    end do

    X2d = reshape(x/H, (/N, N/))
    Y2d = reshape(y/H, (/N, N/))
    theta1 = 1 - Y2d
    theta2 = 0
    theta3 = 0
    do s = 1, smax
        theta2 = theta2 + (A2s(s)*sinh(s*PI*X2d) + B2s(s)*cosh(s*PI*X2d))*sin(s*PI*Y2d)
        theta3 = theta3 + (A3s(s)*sinh(s*PI*(1 - X2d)) + B3s(s)*cosh(s*PI*(1 - X2d)))*sin(s*PI*Y2d)
    end do

    theta = reshape(theta1 + theta2 + theta3, (/N*N/))
    T = theta*(T1 - T0) + T0
    call visualize(dir, 0)

    ! call output(dir, 0, "Temperature", T)

    deallocate(x, y, t, theta)

end subroutine