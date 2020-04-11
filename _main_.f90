! =================================================
!   Main program.
! =================================================

program main

    use coordinates,    only: Nx, Ny
    use numeric_solver, only: calc_numeric => iterate
    use analytic_solver
    implicit none
    
    real*8 :: Tn(Nx,Ny), Tn_(Nx,Ny)
    real*8 :: Ta(Nx,Ny), Ta_(Nx,Ny)
    real*8 :: avg, std
    
    call init_global_vars()
    
    ! Numerical solution
    call calc_numeric(Tn, Tn_, method=3)
    print *, "Numerical solution:"
    call print_matrix(Tn_)
    
    ! Analytical solution
    call calc_analytic(Ta, Ta_)
    print *, "Analytical solution:"
    call print_matrix(Ta_)
    
    ! Compare numerical solution with analytical solution
    call comp_analytic(Tn_, Ta_, avg, std)
    print *, "Average error:     ", avg
    print *, "Standard deviation:", std
    
    pause
    
end program