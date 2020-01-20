! =================================================
!   Subroutine to initialize and nondimensionalize 
!   time, space and temperature.
! =================================================

subroutine init

    use params
    implicit none

    real*8    :: L
    integer*4 :: I, J

    ! Mesh coordinates
    dx = Lx/Nx
    dy = Ly/Ny
    do I = 1, Ny
        do J = 1, Nx
            x(I,J) = (J - 0.5)*dx
            y(I,J) = (I - 0.5)*dy
        end do
    end do

    ! Inital temperature, linear distributed
    T = (Tn - Ts)*y/Ly + Ts

    ! Dimensionless scalars
    L   = .5*(Lx + Ly)
    dx_ = dx/L
    dy_ = dy/L
    x_  = x/L
    y_  = y/L
    T_  = (T - Ts)/(Tn - Ts)
    qw_ = qw*L/lbd/(Tn - Ts)
    qe_ = qe*L/lbd/(Tn - Ts)

end subroutine