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
!      loop, note that "name" and "value" should 
!      be replaced with specific name and value
!      of your variale.
!
!      iter = iter + 1
!      call output(dir, iter, "name", vaule)
! 
! =================================================


subroutine checkdir(dir)

    implicit none

    character(len=*) :: dir
    logical :: dir_exist

    inquire(file=dir, exist=dir_exist)
    if (dir_exist) then
        call system("rmdir /s/q "//dir)
    end if
    call system("md "//dir)

end subroutine


subroutine output(dir, iter, name, value)

    use params, only: N
    implicit none

    character(len=*), intent(in)  :: dir, name
    integer*4,        intent(in)  :: iter
    real*8,           intent(in)  :: value(N,N)
    character(len=8)              :: iter_str
    character(len=:), allocatable :: file_name
    integer*4                     :: I, J

    write(iter_str, '(I8.8)') iter
    file_name = dir//name//"."//iter_str//".dat"
    open(10, file=file_name, action="write", &
        & status="replace")
    write(10, *) "variables=x,y,"//name
    write(10, *) "zone i=", N, "j=", N, "f=point"
    do J = 1, N
        do I = 1, N
            write(10, *) I, J, value(I,J)
        end do
    end do
    close(10)

end subroutine