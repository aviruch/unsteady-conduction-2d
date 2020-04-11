! =================================================
!   Calculate max deviation between two matrices.
! =================================================

subroutine max_deviation(mat, mat0, err)

    use coordinates, only: Nx, Ny
    use conv_config, only: conv_cr
    implicit none

    real*8, intent(in)  :: mat(Nx,Ny)
    real*8, intent(in)  :: mat0(Nx,Ny)
    real*8, intent(out) :: err

    real*8, dimension(Nx*Ny) :: vec, vec0, rand

    vec  = reshape(mat,  (/Nx*Ny/))
    vec0 = reshape(mat0, (/Nx*Ny/))
    call random_seed()
    call random_number(rand)
    err  = maxval(abs(vec - vec0)/(vec0 + 0.001*conv_cr*rand))

end subroutine