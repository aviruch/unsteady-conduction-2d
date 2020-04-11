! =================================================
!   Initialize global variables.
! =================================================

subroutine init_global_vars()

    use physical_def
    use coordinates
    implicit none

    real*8    :: L
    integer*4 :: I, J

    ! Init coordinates
    dx = Lx/Nx
    dy = Ly/Ny
    do I = 1, Nx
        do J = 1, Ny
            x(I,J) = (I - 0.5)*dx
            y(I,J) = (J - 0.5)*dy
        end do
    end do

    ! Non-dimensionalize scalars
    L   = .5*(Lx + Ly)
    dx_ = dx/L
    dy_ = dy/L
    x_  = x/L
    y_  = y/L
    Lx_ = Lx/L
    Ly_ = Ly/L
    qw_ = qw*L/lbd/(Tn - Ts)
    qe_ = qe*L/lbd/(Tn - Ts)

end subroutine