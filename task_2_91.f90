program task_2_91
! 
!   Purpose:
!       Calculate the equation roots of equation a*x^2 + b*x + c = 0 using the Newton-Rapson method

implicit none

! Data dictionary, define variable types & definitons
REAL :: a, b, c                                     ! Coefficients of the quadratic equation
REAL :: x_0                                         ! Initial guess 
REAL :: r                                           ! Root of the quadratic equation
REAL :: f_val, f_prime_val                          ! Values of the function and its derivative at x0
INTEGER :: iter_count = 0                           ! Iteration counter
REAL, PARAMETER :: TOL = 1.E-5, COUNT_MAX = 100     ! Tolerance of the solution convergance and maximum No of iterations

! Read the coefficients of the quadratic equation 
WRITE (*,*) "Please provide the coefficients of the quadratic equation " &
            //"a, b, c in that order."
READ (*,*) a, b, c
! Ask the user for the initial guess
WRITE (*,*) "Please provide the initial guess for x0!"
READ (*,*) x_0

DO 
    iter_count = iter_count + 1         ! Update the iteration counter
    f_val = x_0**2. * a + x_0 * b + c   ! Calcualte f(x0)
    f_prime_val = 2 * x_0 * a + b       ! Calculate f'(x0)
    IF (abs(f_prime_val) < 1.E-8) THEN
        WRITE(*,*) "The derivative is too small; Newton-Raphson method cannot proceed."
        EXIT
    END IF

    r = x_0  - f_val / f_prime_val      ! Calcualte the root

    ! Check if the soultion has converged or the iteration counter reached maximu
    if1: IF (abs(f_val) <= TOL .OR. iter_count > COUNT_MAX) THEN
        WRITE(*, '("The estimated root is ", F5.2)') r
        WRITE(*, '("The value of the function at root is ", ES13.5)') f_val
        WRITE(*, '("Total iterations: ", I5)') iter_count
        EXIT    
    END IF if1

    ! Update the x value for the second guess
    x_0 = r
END DO
end program task_2_91