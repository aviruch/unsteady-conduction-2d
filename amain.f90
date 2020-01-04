include 'params.f90'
include 'init.f90'
include 'visual.f90'
include 'analytic.f90'
include 'explicit.f90'

program main

    use params
    implicit none

    alpha = 0.25

    ! Analytic solution
    call analytic()

    ! Explicit format
    call explicit()

end program