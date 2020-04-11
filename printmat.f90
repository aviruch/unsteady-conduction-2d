! =================================================
!   Print matrix to console.
! =================================================

subroutine print_matrix(mat)

    use coordinates, only: Nx, Ny
    implicit none

    real*8,   intent(in) :: mat(Nx,Ny)
    integer*4            :: I, J

    do J = Ny, 1, -1
        do I = 1, Nx
            write(*, "(F8.3)", advance="no") mat(I,J)
        end do
        print *, ""
    end do

end subroutine