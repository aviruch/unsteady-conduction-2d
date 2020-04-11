! =================================================
!   Jocobi iteration (line) method to solve linear 
!       equation system 'Ax = b'.
! =================================================

subroutine iter_ja_l(A, b, x, iter_cr, iter_max)

    use coordinates, only: Nx, Ny
    use log_config,  only: log_conv_show
    implicit none

    integer*4, parameter     :: n = Nx*Ny
    real*8,    intent(in)    :: A(n,n)
    real*8,    intent(in)    :: b(n)
    real*8,    intent(inout) :: x(n)
    real*8,    intent(in)    :: iter_cr
    integer*4, intent(in)    :: iter_max

    real*8    :: x0(n)
    real*8    :: iter_err
    integer*4 :: iter_cnt
    integer*4 :: I, J, P

    real*8    :: aa(Ny), bb(Ny), cc(Ny), dd(Ny)
    real*8    :: xx(Ny)

    call random_seed()
    call random_number(x0)
    iter_cnt = 1
    do while (iter_cnt < iter_max)

        do I = 1, Nx

            aa(1) = 0.
            do J = 2, Ny
                P = (I - 1)*Ny + J
                aa(J) = A(P,P-1)
            end do

            do J = 1, Ny - 1
                P = (I - 1)*Ny + J
                cc(J) = A(P,P+1)
            end do
            cc(Ny) = 0.

            do J = 1, Ny
                P = (I - 1)*Ny + J
                bb(J) = A(P,P)
                if (I == 1)      then
                    dd(J) = A(P,P+Nx)*x0(P+Nx) + b(P)
                elseif (I == Nx) then
                    dd(J) = A(P,P-Nx)*x0(P-Nx) + b(P)
                else
                    dd(J) = A(P,P-Nx)*x0(P-Nx) &
                        & + A(P,P+Nx)*x0(P+Nx) + b(P)
                end if
            end do
        
            call TDMA(Ny, aa, bb, cc, dd, xx)
            do J = 1, Ny
                P = (I - 1)*Ny + J
                x(P) = xx(J)
            end do

        end do
        
        ! do I = 1, n
        !     write(*,*) x(I)
        ! end do
        ! stop

        ! Update: x -> x0
        iter_err = maxval(abs(x - x0)/x0)
        if (log_conv_show) print *, "", iter_cnt,  iter_err
        if (iter_err < iter_cr) exit

        x0 = x
        iter_cnt = iter_cnt + 1

    end do

end subroutine