include './lib/VTKWRITER.f90'

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

subroutine visualize(dir, iter)

    use params
    use vtkwriter
    implicit none

    type(writer) :: vtk_writer
    character(len=*) :: dir
    integer(kind=4) :: iter
    character(len=8) :: iter_str

    write(iter_str, '(I8.8)') iter
    call vtk_writer%init(dir // "data." // iter_str // ".vtk")
    call vtk_writer%write2d(x, y, (/N, N/), "Temperature", T)

end subroutine