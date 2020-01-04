include 'params.f90'
include 'init.f90'
include 'visual.f90'
include 'analytic.f90'
include 'explicit.f90'
include 'implicit_ja_p.f90'

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