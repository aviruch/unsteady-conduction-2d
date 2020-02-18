! =================================================
!   Jocobi iteration (point) method to solve linear 
!       equation system 'Ax = b'.
! =================================================

subroutine iter_ja_p(n, A, b, x, iter_cr, iter_max)

    implicit none

    integer*4, intent(in)    :: n
    real*8,    intent(in)    :: A(n,n)
    real*8,    intent(in)    :: b(n)
    real*8,    intent(inout) :: x(n)
    real*8,    intent(in)    :: iter_cr
    integer*4, intent(in)    :: iter_max

    real*8    :: x0(n)
    real*8    :: iter_err, tmp
    integer*4 :: iter_cnt
    integer*4 :: I, J

    x0 = 0
    iter_cnt = 1
    do while (iter_cnt < iter_max)

        do I = 1, n
            tmp = 0.
            do J = 1, n
                if (J /= I) then
                    tmp = tmp + A(I,J)*x0(J)
                end if
            end do
            x(I) = (b(I) - tmp)/A(I,I)
        end do

        iter_err = maxval(abs(x - x0)/x0)
        write(*, "(A5, I4, A8, E12.4E2)") &
            & "", iter_cnt,  "", iter_err
        if (iter_err < iter_cr) exit

        x0 = x
        iter_cnt = iter_cnt + 1

    end do

end subroutine