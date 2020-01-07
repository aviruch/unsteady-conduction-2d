program main

    use params, only: alpha, omega
    implicit none

    alpha = 0.25
    omega = 1.20

    ! Analytic solution
    call analytic()

    ! Explicit format
    call explicit()

    ! Implicit format
    call implicit_ja_p()

end program