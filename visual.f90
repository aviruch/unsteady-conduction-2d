! =================================================
!   Output data files for post-processing.
! =================================================

subroutine check_dir(dir)

    implicit none

    character(len=*) :: dir

    call system("if exist "//dir//" (rmdir /s/q "//dir//")")
    call system("md "//dir)

end subroutine check_dir

subroutine save_tecplot_file(dir, iter, var, val)

    use coordinates, only: Nx, Ny, x, y
    implicit none

    character(len=*), intent(in)  :: dir, var
    integer*4,        intent(in)  :: iter
    real*8,           intent(in)  :: val(Nx,Ny)

    character(len=8)              :: iter_str
    character(len=:), allocatable :: file
    real*8                        :: x1d(Nx*Ny)
    real*8                        :: y1d(Nx*Ny)
    real*8                        :: f1d(Nx*Ny)
    integer*4                     :: I, J

    write(iter_str, '(I8.8)') iter
    file = dir//var//"."//iter_str//".dat"
    open(10, file=file, action="write", &
        & status="replace")
    write(10, *) "variables=x,y,"//var
    write(10, *) "zone i=", Ny, "j=", Nx, "f=point"

    x1d = reshape(transpose(x), (/Nx*Ny/))
    y1d = reshape(transpose(y), (/Nx*Ny/))
    f1d = reshape(transpose(val), (/Nx*Ny/))
    do I = 1, Nx*Ny
        write(10, *) x1d(I), y1d(I), f1d(I)
    end do

    close(10)

end subroutine save_tecplot_file