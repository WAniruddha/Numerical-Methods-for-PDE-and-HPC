! IVP2 Lax-Wendroff Implementation
program ivp2_lax_wendroff
    implicit none

    ! Declare variables
    integer :: nx, nx_analytical, i, iterations
    real :: dx, dt, CFL, L, t_final, c, x_min, x_max, dx_analytical, dt_analytical, x_shifted, t_c
    real, allocatable :: x(:), f(:,:), f_analytical(:)

    ! Constants
    nx = 2000               ! Number of grid points
    nx_analytical = 4000   ! Higher resolution for analytical solution
    t_final = 1000.0       ! Final simulation time
    c = 1.0                ! Wave speed
    x_min = -1.0
    x_max = 1.0
    L = x_max - x_min      ! Domain length
    CFL = 1.0              ! CFL fixed at 1.0

    ! Grid and time step
    dx = L / (nx - 1)
    dx_analytical = L / (nx_analytical - 1)
    dt = CFL * dx / c      ! Ensures CFL = 1 always
    dt_analytical = CFL * dx_analytical / c
    iterations = int(t_final / dt)

    ! Allocate arrays
    allocate(x(nx_analytical), f_analytical(nx_analytical))

    ! Initialize analytical solution grid
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

    deallocate(x, f_analytical) 
    
    ! Allocate arrays for numerical scheme
    allocate(x(nx), f(nx, 2))

    ! Initialize numerical grid
    do i = 1, nx
        x(i) = x_min + (i - 1) * dx
    end do

    ! Initialize function values for Lax-Wendroff
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

    ! Time-stepping loop for Lax-Wendroff Scheme
    t_c = 0.0  ! Initialize time counter

    do iterations = 1, int(t_final / dt)
        ! Apply Lax-Wendroff Scheme
        do i = 2, nx-1
            f(i,2) = f(i,1) - 0.5 * CFL * (f(i+1,1) - f(i-1,1)) + 0.5 * CFL**2 * (f(i+1,1) - 2.0*f(i,1) + f(i-1,1))
        end do

        ! Apply periodic boundary conditions
        f(1,2) = f(nx-1,1)
        f(nx,2) = f(2,1)

        ! Update time
        t_c = t_c + dt

        ! Update values for the next step
        f(:,1) = f(:,2)

        ! Exit loop if we reached t_final exactly - add for debuging
        !if (abs(t_c - t_final) < 1.0e-6) exit
    end do

    ! Write Lax-Wendroff Scheme Solution
    open(11, file='Lax_Wendroff_Solution.dat')
    do i = 1, nx
        write(11, *) x(i), f(i,2)
    end do
    close(11)

    ! Deallocate arrays
    deallocate(x, f)

end program ivp2_lax_wendroff
