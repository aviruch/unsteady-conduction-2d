! =================================================
!   Print matrix to console.
! =================================================

subroutine print(mat)

    use params, only: Nx, Ny
    implicit none

    real*8,   intent(in) :: mat(Nx,Ny)
    integer*4            :: I, J

    do J = Ny, 1, -1
        do I = 1, Nx
            write(*, "(F8.3)", advance="no") mat(I,J)
        end do
        write(*,*) ""
    end do

end subroutine