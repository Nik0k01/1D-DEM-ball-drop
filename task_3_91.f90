PROGRAM task_3_91
! 
!   Purpose:
!       Calculate the minimum of a function describing the distance between a sphere and an elypse.
!       The distance function: F(theta, phi) = (x_c + R * cos(phi) - a * cos(theta)) ** 2 +
!                                            + (y_c + R * sin(phi) - b * sin(theta)) ** 2
!       Algorithm implemented - gradient descent  
! 
    IMPLICIT NONE
    ! Constants
    REAL, PARAMETER :: pi = 3.14159265359       ! Pi value
    REAL, PARAMETER :: tol = 1.e-6              ! Solution tolerance
    INTEGER, PARAMETER :: iter_max = 1000       ! Maximum number of iteration
    REAL, PARAMETER :: damping_coeff = 0.2      ! Relaxation parameter
    ! Data dictionary: declare calling parameters, types & units
    REAL :: x_c, y_c                            ! Sphere coordinates
    REAL :: a, b                                ! Elipse parameters
    REAL :: R                                   ! Sphere radius
    REAL :: theta, phi, theta_0, phi_0          ! Angles used in parametric formulation, radians
    REAL :: f_prime_phi_val                     ! Value of the function gradient
    REAL :: f_prime_theta_val                   ! Value of the function gradient
    REAL :: initial_f_val, final_f_val          ! Function values for comaprision 
    INTEGER :: iter_count = 0                   ! Iteration counter

    ! Example values to test the program
    WRITE (*,*) "Provide values for x_c, y_c, a, b, R:"
    READ (*,*) x_c, y_c, a, b, R
    ! Initial guess for the values of the angles
    WRITE (*,*) "Provide initial guess for phi and theta (degrees):"
    READ (*,*) phi, theta

    ! Convert degrees to radians
    phi = phi / 180 * pi 
    theta = theta / 180 * pi

    ! Initial value of the function
    initial_f_val = (x_c + R * cos(phi) - a * cos(theta)) ** 2 &
                  + (y_c + R * sin(phi) - b * sin(theta)) ** 2

    min_loop: DO
        iter_count = iter_count + 1

        ! Calcualte the derivative with respect to theta, and then update the value
        f_prime_theta_val = a * sin(theta) * 2. * (x_c + R * cos(phi) - a * cos(theta)) &
                          - b * cos(theta) * 2. * (y_c + R * sin(phi) - b * sin(theta))
        
        ! Calculate the derivative with respect to phi, and then update the value
        f_prime_phi_val = R * cos(phi) * 2. * (y_c + R * sin(phi) - b * sin(theta)) & 
                        - R * sin(phi) * 2. * (x_c + R * cos(phi) - a * cos(theta)) 

        theta_0 = theta
        theta = theta - damping_coeff * f_prime_theta_val
        ! Check if the value of the angle is within bounds (0, 2pi)
        IF (theta > 2. * pi) THEN
            theta = 2. * pi
        ELSEIF (theta < 0.) THEN
            theta = 0.
        END IF

        phi_0 = phi
        phi = phi - damping_coeff * f_prime_phi_val
        ! Check if the value of the angle is within bounds (0, 2pi)
        IF (phi > 2 * pi) THEN
            phi = 2 * pi
        ELSEIF (phi < 0.) THEN
            phi = 0.
        END IF

        ! Loop exit conditions
        IF (abs(phi - phi_0) <= tol .AND. abs(theta - theta_0) <= tol .OR. iter_count == iter_max) EXIT
    END DO min_loop

    ! Final value of the function
    final_f_val = (x_c + R * cos(phi) - a * cos(theta)) ** 2 &
                  + (y_c + R * sin(phi) - b * sin(theta)) ** 2

    ! Display the results
    WRITE (*, '("Iteration count ", I3)') iter_count
    WRITE (*, '("phi = ", F8.2, " theta = ", F8.2)') phi, theta
    WRITE (*, '("Initial value of the function = ", ES10.2)') initial_f_val
    WRITE (*, '("Final value of the function = ", ES10.2)') final_f_val

END PROGRAM task_3_91