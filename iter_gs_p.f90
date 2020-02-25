! =================================================
!   Gauss-Seidel iteration (point) method to solve  
!       linear equation system 'Ax = b'.
! =================================================

subroutine iter_gs_p(n, A, b, x, iter_cr, iter_max)

    implicit none

    integer*4, intent(in)    :: n
    real*8,    intent(in)    :: A(n,n)
    real*8,    intent(in)    :: b(n)
    real*8,    intent(inout) :: x(n)
    real*8,    intent(in)    :: iter_cr
    integer*4, intent(in)    :: iter_max

    real*8    :: x0(n)
    real*8    :: iter_err, tmp1, tmp2
    integer*4 :: iter_cnt
    integer*4 :: I, J

    x0 = 0
    iter_cnt = 1
    do while (iter_cnt < iter_max)

        do I = 1, n
            tmp1 = 0.
            do J = 1, I - 1
                tmp1 = tmp1 + A(I,J)*x(J)
            end do
            tmp2 = 0.
            do J = I + 1, n
                tmp2 = tmp2 + A(I,J)*x0(J)
            end do
            x(I) = (b(I) - tmp1 - tmp2)/A(I,I)
        end do

        iter_err = maxval(abs(x - x0)/x0)
        write(*, "(A5, I4, A8, E12.4E2)") &
            & "", iter_cnt,  "", iter_err
        if (iter_err < iter_cr) exit

        x0 = x
        iter_cnt = iter_cnt + 1

    end do

end subroutine