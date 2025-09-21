module distributions
    implicit none
    REAL(8), PARAMETER :: pi = 3.14159265359d0
    
contains

    FUNCTION distribution_array(dist_type, n, params, max_val, min_val) 
        ! Purpose:
        !   Returns a 1000-element 2D matrix with arguments and values for specified distributions.

        REAL(8), DIMENSION(1000, 2) :: distribution_array       ! Return array
        CHARACTER(*) :: dist_type                               !  Distribution type: 'NORM', 'LOG-NORM', 'UNI'
        INTEGER :: n                                            ! Number of parameters
        REAL(8), DIMENSION(n) :: params                         ! Distribution parameters
        REAL(8) :: max_val, min_val                             ! Argument range
        REAL(8) :: incr                                         ! Argument increment
        INTEGER :: i                                            ! Loop variable

        ! Validate input range
        IF (max_val <= min_val) THEN
            WRITE(*,*) "Error: max_val must be greater than min_val!"
            distribution_array = 0.0
            RETURN
        END IF

        ! Initialize increment
        incr = (max_val - min_val) / 999.0

        ! Distribution selection
        SELECT CASE (dist_type)
        CASE ('NORM')
            IF (n < 2) THEN
                WRITE(*,*) "Error: Normal distribution requires at least 2 parameters (mean, std)."
                distribution_array = 0.0
                RETURN
            END IF
            DO i = 1, 1000
                distribution_array(i, 1) = min_val + (i - 1) * incr
                distribution_array(i, 2) = 1.0 / (params(2) * sqrt(2 * pi)) * &
                                            EXP(-0.5 * ((distribution_array(i, 1) - params(1)) / params(2))**2)
            END DO

        CASE ('LOG-NORM')
            IF (n < 2) THEN
                WRITE(*,*) "Error: Log-normal distribution requires at least 2 parameters (mean, std)."
                distribution_array = 0.0
                RETURN
            END IF
            DO i = 1, 1000
                distribution_array(i, 1) = min_val + (i - 1) * incr
                IF (distribution_array(i, 1) <= 0) THEN
                    distribution_array(i, 2) = 0.0
                ELSE
                    distribution_array(i, 2) = 1.0 / (distribution_array(i, 1) * params(2) * sqrt(2 * pi)) * &
                                            EXP(-(LOG(distribution_array(i, 1)) - params(1))**2 / (2 * params(2)**2))
                END IF
            END DO

        CASE ('UNI')
            DO i = 1, 1000
                distribution_array(i, 1) = min_val + (i - 1) * incr
                distribution_array(i, 2) = 1.0 / (max_val - min_val)
            END DO

        CASE DEFAULT
            WRITE(*,*) "Error: Distribution type not recognized!"
            distribution_array = 0.0
        END SELECT
    END FUNCTION

    SUBROUTINE data_file_generator(dist_type, min_value, max_value, mean, std)
        ! Purpose:
        !   Generate a data file with values of a given distribution in a specified 
        !   argument range. 

        ! Data dictionary
        CHARACTER(*), INTENT(IN) :: dist_type                   ! Distribution type
        REAL(8), INTENT(IN) :: min_value                        ! Minimum value
        REAL(8), INTENT(IN) :: max_value                        ! Maximum value
        REAL(8), INTENT(IN) :: mean                             ! Distribution mean - can be whatever for uniform
        REAL(8), INTENT(IN) :: std                              ! Standard deviation - can be whatever for uniform
        REAL(8), DIMENSION(1000, 2) :: value_matrix             ! Matrix with arguments and values
        REAL(8), DIMENSION(2) :: params                         ! Array for distribution parameters
        INTEGER :: file_id = 9                                  ! File id
        CHARACTER(10) :: date                                   ! Date
        CHARACTER(12) :: time                                   ! Time
        CHARACTER(50) :: file_name                              ! File name
        INTEGER :: ierror                                       ! File status
        CHARACTER(80) :: err_string                             ! Error message

        ! Use the distribution_array function to generate matrix with arguments and value
        params = [mean, std]
        value_matrix = distribution_array(dist_type, 2, params, max_value, min_value)

        ! Generate a file name
        CALL date_and_time(date, time)
        file_name = date // '_' // time // '_' // dist_type // '.dat'
        ! Open a file 
        OPEN (UNIT=file_id, FILE=file_name, STATUS='NEW', ACTION='WRITE', &
              IOSTAT=ierror, IOMSG=err_string)
        

    END SUBROUTINE
    
end module distributions