! =================================================
!   TDMA algorithm
! =================================================

subroutine TDMA(n, a, b, c, d, x)

    implicit none

    integer*4, intent(in)  :: n
    real*8,    intent(in)  :: a(n)
    real*8,    intent(in)  :: b(n)
    real*8,    intent(in)  :: c(n)
    real*8,    intent(in)  :: d(n)
    real*8,    intent(out) :: x(n)

    real*8    :: P(n) 
    real*8    :: Q(n) 
    real*8    :: R(n) 
    integer*4 :: I

    P(1) = -c(1)/b(1)
    Q(1) =  d(1)/b(1)
    do I = 2, n
        R(I) =  a(I)*P(I-1) + b(I)
        P(I) = -c(I)/R(I)
        Q(I) = (d(I) - a(I)*Q(I-1))/R(I)
    end do

    x(N) = Q(N)
    do I = N - 1, 1, -1
        x(I) = P(I)*x(I+1) + Q(I)
    end do

end subroutine