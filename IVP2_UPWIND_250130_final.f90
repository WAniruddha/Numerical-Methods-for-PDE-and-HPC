! IVP2 Implementation 
program IVP2_Upwind
    implicit none

    ! Declare variables
    integer :: nx, nx_analytical, i, iterations
    real :: dx, dt, CFL, L, t_final, c, x_min, x_max, dx_analytical, dt_analytical, x_shifted
    real, allocatable :: x(:), f(:,:), f_analytical(:)
    
    ! Constants
    nx = 2000             
    nx_analytical = 4000
    t_final = 01000.0       
    c = 1.0
    x_min = -1.0
    x_max = 1.0
    L = x_max - x_min
    CFL = 01.0  ! Set CFL to 1.0 for validation

    ! Grid and time step
    dx = L / (nx - 1)
    dx_analytical = L / (nx_analytical - 1)
    dt = CFL * dx / c
    dt_analytical = CFL * dx_analytical / c
    iterations = int(t_final / dt)

    ! Allocate arrays
    allocate(x(nx_analytical), f(nx, 2), f_analytical(nx_analytical))

    ! Initialize x(i) properly
    do i = 1, nx_analytical
        x(i) = x_min + (i - 1) * dx_analytical
    end do

    ! Compute Analytical Solution with Periodic Shift
    do i = 1, nx_analytical
        x_shifted = x(i) - c * t_final

        ! Correct periodic boundary conditions
        do while (x_shifted < x_min)
            x_shifted = x_shifted + L
        end do
        do while (x_shifted > x_max)
            x_shifted = x_shifted - L
        end do

        ! Assign Initial Condition
        if (x_shifted >= -0.8d0 .and. x_shifted <= -0.6d0) then
            f_analytical(i) = exp(log(2.0d0) * ((x_shifted + 0.7d0)**2) / 0.0009d0)
        else if (x_shifted >= -0.4 .and. x_shifted <= -0.2) then
            f_analytical(i) = 1.0
        else if (x_shifted >= 0.0 .and. x_shifted <= 0.2) then
            f_analytical(i) = 1.0 - abs(10.0 * x_shifted - 1.0)
        else if (x_shifted >= 0.4 .and. x_shifted <= 0.6) then
            f_analytical(i) = sqrt(max(0.0, 1.0 - 100.0 * (x_shifted - 0.5)**2))
        else
            f_analytical(i) = 0.0
        end if
    end do

    ! Write Analytical Solution
    open(10, file='Analytical_Solution.dat')
    do i = 1, nx_analytical
        write(10, *) x(i), f_analytical(i)
    end do
    close(10)

    ! Deallocate and reallocate for numerical scheme
    deallocate(x, f, f_analytical)
    allocate(x(nx), f(nx, 2))

    ! Initialize x(i) for numerical scheme
    do i = 1, nx
        x(i) = x_min + (i - 1) * dx
    end do

    ! Initialize function values
    do i = 1, nx
        if (x(i) >= -0.8d0 .and. x(i) <= -0.6d0) then
            f(i,1) = exp(log(2.0d0) * ((x(i) + 0.7d0)**2) / 0.0009d0)
        else if (x(i) >= -0.4 .and. x(i) <= -0.2) then
            f(i,1) = 1.0
        else if (x(i) >= 0.0 .and. x(i) <= 0.2) then
            f(i,1) = 1.0 - abs(10.0 * x(i) - 1.0)
        else if (x(i) >= 0.4 .and. x(i) <= 0.6) then
            f(i,1) = sqrt(max(0.0, 1.0 - 100.0 * (x(i) - 0.5)**2))
        else
            f(i,1) = 0.0
        end if
    end do



    ! Time-stepping loop for Upwind Scheme
    do iterations = 1, int(t_final / dt)
        ! Apply Periodic Boundary Condition
        f(1, 1) = f(nx, 1)

        do i = 2, nx
            ! Apply Upwind Scheme Formula
            f(i, 2) = f(i, 1) - CFL * (f(i, 1) - f(i - 1, 1))
        end do

        ! Apply Periodic Boundary Condition
        f(1, 2) = f(nx, 1)

        ! Update values
        f(:, 1) = f(:, 2)
    end do

    ! Write Upwind Scheme Solution
    open(11, file='Upwind_Solution.dat')
    do i = 1, nx
        write(11, *) x(i), f(i,2)
    end do
    close(11)

    deallocate(x,f)

end program IVP2_Upwind