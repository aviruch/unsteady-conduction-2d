! =================================================
!   Subroutines to output data files, for post- 
!   processing.
! 
!   Usage: 
!   1. Check if output directory exists. Add the 
!      following lines outside the time-marching 
!      loop.
! 
!      character(len=:), allocatable :: dir
!      integer*4                     :: iter
!      iter = 0
!      dir  = ".\\out\\"
!      call checkdir(dir)
! 
!   2. Call subroutine to write data files. Add 
!      the following lines in the time-marching 
!      loop, note that "var" and "val" should 
!      be replaced with specific name and value
!      of your variale.
!
!      iter = iter + 1
!      call output(dir, iter, "var", val)
! 
! =================================================


subroutine checkdir(dir)

    implicit none

    character(len=*) :: dir
    logical          :: dir_exist

    inquire(file=dir, exist=dir_exist)
    if (dir_exist) then
        call system("rmdir /s/q "//dir)
    end if
    call system("md "//dir)

end subroutine


subroutine output(dir, iter, var, val)

    use params, only: Nx, Ny, x, y
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

end subroutine