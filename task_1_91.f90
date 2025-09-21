program task_1_91
! 
!   Purpose:
!       Calculate the equation roots of equation a*x^2 + b*x + c = 0 algebraically 

implicit none

! Data dictionary, define variable types & definitons
REAL :: a, b, c             ! Coefficients of the quadratic equation
REAL :: delta               ! Delta of quadratic equation
REAL :: r1, r2              ! Roots of the quadratic equation
REAL :: r1_imag, r2_imag    ! Imaginary parts of the roots
REAL :: r1_real, r2_real    ! Real parts of the roots
COMPLEX :: r1_c, r2_c       ! Complex roots

WRITE (*,*) "Please provide the coefficients of the quadratic equation " &
            //"a, b, c in that order."
READ (*,*) a, b, c

delta = b**2. - 4. * a * c
IF (delta > 0.) THEN
    WRITE (*,*) 'There are two real roots!'
    r1 = (-b - sqrt(delta)) / (2. * a)
    r2 = (-b + sqrt(delta)) / (2. * a)
    WRITE(*, '("The roots of the equation are " F5.2 " and " F5.2)') r1, r2
ELSE IF (abs(delta) < 1.E-5) THEN
    WRITE (*,*) 'There is one repeated root!'
    r1 = -b / (2. * a)
    WRITE(*, '("The repeated root is ", F5.2)') r1
ELSE
    WRITE (*,*) 'There are two complex roots!'
    r1_real = -b / (2. * a)
    r1_imag = sqrt(abs(delta)) / (2. * a)
    r2_real = r1_real
    r2_imag = -r1_imag
    r1_c = CMPLX(r1_real, r1_imag)
    r2_C = CMPLX(r2_real, r2_imag)
    WRITE(*, '("The roots of the equation are ", F5.2, " + ", F5.2, "i and ", F5.2, " - ", F5.2, "i")') &
            r1_real, r1_imag, r2_real, abs(r2_imag)
END IF
end program task_1_91