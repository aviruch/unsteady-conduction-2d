! =================================================
!   Subroutine to calculate max deviation between
!   two matrices.
! =================================================

subroutine maxerr(Nx, Ny, val, val0, err)

    implicit none

    integer*4, intent(in)  :: Nx, Ny
    real*8,    intent(in)  :: val(Ny,Nx)
    real*8,    intent(in)  :: val0(Ny,Nx)
    real*8,    intent(out) :: err

    real*8, allocatable :: tmp(:)
    real*8, allocatable :: tmp0(:)

    tmp  = reshape(val, (/Ny*Nx/))
    tmp0 = reshape(val0, (/Ny*Nx/))
    err  = maxval(abs(tmp - tmp0)/tmp0)

end subroutine