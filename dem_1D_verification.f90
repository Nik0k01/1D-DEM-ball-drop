PROGRAM dem_1D_verification
    ! Purpose:
    !   To verify if the results of the 1D DEM simulation are in agreement with 
    !   the analytical solution. Only part of the particle motion is going to be simulated.
    
    IMPLICIT NONE

    ! Physical constants, properties
    REAL, PARAMETER :: g = 9.8067                           ! Gravitational constant m/s^2
    REAL, PARAMETER :: k = 800.0                            ! Stiffness N/m 
    REAL, PARAMETER :: r = 0.1                              ! Sphere radius m
    REAL, PARAMETER :: m = 0.1                              ! Sphere mass kg
    REAL, PARAMETER :: c = 0.5                              ! Viscous damping force

    ! Simulation parameters
    INTEGER(8), PARAMETER :: steps_no = 1000                  ! Number of iterations
    REAL(8), PARAMETER :: total_time = 0.8d0                  ! Total simulation time
    REAL(8), PARAMETER :: time_step = total_time / REAL(steps_no) ! Time step 

    ! Data dictionary: types, variables, units
    CHARACTER(50) :: data_wo_damping                        ! DEM simulation result without damping 
    CHARACTER(50) :: data_w_damping                         ! DEM simulation result with damping 
    CHARACTER(50) :: results_formula                        ! Name of file with the analytical results 
    REAL(8), DIMENSION(steps_no) :: position                ! Position array
    REAL(8), DIMENSION(steps_no) :: time                    ! Time array
    REAL(8) :: time_free_fall                               ! Time of free fall 
    REAL(8) :: init_height                                  ! Initial height
    INTEGER(8) :: i                                         ! Loop counter
    REAL(8) :: c_1, c_2                                     ! Integration constants
    REAL(8) :: omega                                        ! Temporal frequency of the spring
    REAL(8) :: rest_factor_analy                            ! Analytially calculated restitution factor
    REAL(8) :: rest_factor_num                              ! Numerically calculated restitution factor

    ! Initial conditions
    init_height = 1.0d0    
    position(1) = init_height  
    time(1) = 0d0

    ! Calculate the time in which the particle is falling
    time_free_fall = sqrt(2 * (init_height - r) / g)
    ! Calculate the spring temporal frequency
    omega = sqrt(k / m)

    ! Main calculation loop - analytical position vs time
    main_loop: DO i = 2, steps_no
            ! Update the time
            time(i) = time(i-1) + time_step
            ! Function of position in free fall is different than in contact with the ground
            time_check: IF (time(i) < time_free_fall) THEN
                ! Update the position
                position(i) = -g * time(i) ** 2 / 2.d0 + init_height
            ELSE
                c_2 = (-sqrt(2 * g * (init_height - r)) + m * g / k * tan(omega * time_free_fall)) &
                * (omega * cos(omega * time_free_fall) + tan(omega * time_free_fall) * omega   &
                * sin(omega * time_free_fall)) ** (-1) 
                ! 1st integration constant
                c_1 = m * g / k * cos(omega * time_free_fall) ** (-1) - c_2 * tan(omega * time_free_fall)
                ! Position function is a solution of nonhomogenous 2nd order ordinary differential equation
                position(i) = c_1 * cos(omega * time(i)) + c_2 * sin(omega * time(i)) + r - m * g / k
            END IF time_check
        END DO main_loop

    ! Select the data file with simulation results
    data_wo_damping = '20241205_092036.406_results.dat'
    data_w_damping = '20241205_123409.795_results.dat'

    ! Write analytical results to a file
    results_formula = 'analytical_data_file.dat'        ! Name of the file with analytical results
    CALL file_writer(time, position, steps_no, results_formula)

    ! Calculate the restitution coefficients
    rest_factor_analy = resitution_factor(k, m, c)                          ! Analytical calculation
    CALL resitution_factor_numerical(data_w_damping, rest_factor_num, r)    ! Numerical calculation
    ! Display the results to the users
    WRITE(*, '("Numerical restitution factor:", F10.5, " Analytical restitution factor:", F10.5)') rest_factor_num, rest_factor_num

    ! Create comparision plots - case without damping coefficient
    CALL plot_compare(results_formula, data_wo_damping)

    CONTAINS

    FUNCTION resitution_factor(k, m, c)
        ! Purpose:
        !   Calculate the restitution factor analytically 
        ! Function return
        REAL(8) :: resitution_factor
        ! Constants
        REAL(8), PARAMETER :: pi = 3.14159265359d0
        ! Variables
        REAL :: k, m, c                              ! Stifness coefficient, mass, damping factor
        resitution_factor = exp(pi * c / sqrt(4 * k * m - c ** 2))
    END FUNCTION

    SUBROUTINE resitution_factor_numerical(data_file_name, rest_factor, r)
        CHARACTER(50), INTENT(IN) :: data_file_name                     ! Name of the file with numerical data
        REAL(8), INTENT(OUT) :: rest_factor                             ! Numerical rest factor
        CHARACTER(80) :: msg                                            ! Error message
        INTEGER :: nvals = 0                                            ! Number of values to read in
        INTEGER :: status                                               ! Error flag when reading the file
        REAL(8), ALLOCATABLE, DIMENSION(:, :) :: result_matrix          ! Matrix result
        CHARACTER(len=256) :: line                                      ! Line in the dat file
        INTEGER :: i, j                                                 ! Control loop integers
        REAL(8) :: vel_before, vel_after                                ! Velocity before and after the collision
        REAL :: r                                                       ! Radius of the sphere - needed to detect collison

        ! Open the file
        OPEN (UNIT=3, FILE=data_file_name, STATUS='OLD', ACTION='READ', &
              IOSTAT=status, IOMSG=msg)
        ! Read the contents of the file 
        openif: IF (status == 0) THEN
            
            ! Open was ok, read the value
            countloop: DO
                READ(3, *, IOSTAT=status) line      ! Read the line
                IF (status /= 0) EXIT               ! EXIT if not valid
                nvals = nvals + 1
            END DO countloop

        ! Rewind the file to read again
            rewind(3)
        ! Allocate the memory for the result matrix
            ALLOCATE(result_matrix(nvals, 5))
        ! Read the data into the matrix
            DO i = 1, nvals
                READ(3, *) (result_matrix(i, j), j = 1, 5)
            END DO
            CLOSE (3)
        END IF openif

        ! Detect the moment just before the collision
        DO i = 1, nvals
            ! Check the position of the sphere
            IF (result_matrix(i, 2) <= r) THEN
                vel_before = result_matrix(i-1, 3)            ! Velocity before collision is one position above in the result matrix
                EXIT                                          ! Velocity has been found - exit the loop
            END IF
        END DO
        ! Detect the moment just after the collision
        DO j = i, nvals                                       ! Loop through matrix from the position of the first collision
            ! Check if the sphere is no longer touching the ground
            IF (result_matrix(j, 2) > r) THEN
                vel_after = result_matrix(j, 3)               ! Velocity after the collision is when the sphere is no longer touching the ground
                EXIT                                          ! Velocity has been found - exit the loop
            END IF
        END DO
        DEALLOCATE(result_matrix)                             ! Free up the memory
        ! Calculate the restitution factor
        rest_factor = abs(vel_after / vel_before)
    END SUBROUTINE

    SUBROUTINE file_writer(time, position, steps_no, file_name) 
        ! Writes calculation results to a file
        REAL(8), DIMENSION(steps_no), INTENT(IN) :: time
        REAL(8), DIMENSION(steps_no), INTENT(IN) :: position
        CHARACTER(*), INTENT(IN) :: file_name   ! Name of the file with results
        INTEGER(8), INTENT(IN) :: steps_no      ! Number of time steps
        INTEGER :: ierror                       ! File status
        INTEGER :: file_id                      ! File id
        CHARACTER(80) :: err_string             ! Error message
        INTEGER(8) :: i                         ! Counter 

        ! Id of the file
        file_id = 20
        ! Open a file
        OPEN (UNIT=file_id, FILE=file_name, STATUS='REPLACE', ACTION='WRITE', &
              IOSTAT=ierror, IOMSG=err_string)
        
        ! Check if the file was opened correctly
        openif: IF (ierror == 0) THEN
            1010 FORMAT (F10.5, ",", F10.5)
            ! Open was ok. Write values
            writeloop: DO i = 1, steps_no
                WRITE(file_id, 1010) time(i), position(i)
            END DO writeloop
        ! If openining of the file was not correct
        ELSE openif
            WRITE(*, 1040) ierror
            1040 FORMAT ('Error openin file: IOSTAT = ', I6) 
            WRITE(*, 1050) TRIM(err_string)
            1050 FORMAT (A)  
        END IF openif 
        
        ! Close file
        CLOSE (UNIT=file_id)
    END SUBROUTINE

    SUBROUTINE plot_compare(analytical_file, numerical_file)
        ! Generate a Gnuplot script and plot the results
        CHARACTER(50), INTENT(IN) :: numerical_file    ! File with numerical results
        CHARACTER(50), INTENT(IN) :: analytical_file    ! File with analytical results
        CHARACTER(9) :: script_filename                 ! Gnuplot script name
        CHARACTER(50) :: plot_name                      ! Array with plot names
        CHARACTER(LEN=:), ALLOCATABLE :: command        ! Command made by Gnuplot
        CHARACTER(1), DIMENSION(3) :: column_pairs      ! Pairs of columns for plotting
        CHARACTER(LEN=230) :: plot_command_1
        CHARACTER(LEN=230) :: plot_command_2

    
        ! Name of the Gnuplot script
        script_filename = 'script.gp'
        
        ! Set up the output plot names
        plot_name = 'Plots/position_comparision_plot.png' 
        
        ! Specify the columns to use for each plot
        column_pairs = ['4', '2', '3'] ! (time:force, time:position, time:velocity)
        
        ! Write Gnuplot commands to the script
        OPEN (UNIT=20, FILE=script_filename, STATUS='REPLACE')
        WRITE (20, '(A)') 'set terminal pngcairo size 800,600 enhanced font "Arial,14"'
        WRITE (20, '(A)') 'set grid'
        WRITE (20, '(A)') 'set style line 1 lt 1 lw 2 lc rgb "blue"'
        WRITE (20, '(A)') 'set style line 2 lt 2 lw 2 lc rgb "red"'
        WRITE (20, '(A)') 'set style line 3 lt 3 lw 2 lc rgb "green"'
        WRITE (20, '(A)') 'set xrange [0:0.6]'
        WRITE (20, '(A)') 'set output "' // plot_name // '"'
        WRITE (20, '(A)') 'set xlabel "Time (s)"'
        WRITE (20, '(A)') 'set title "Position vs Time"'
        WRITE (20, '(A)') 'set ylabel "Position (m)"'
        plot_command_1 = 'plot "' // trim(analytical_file) // '" using 1:2 with lines ls 2 title "Position Analytical",' 
        plot_command_2 = '"' // trim(numerical_file) // '" using 1:2 with lines ls 1 dt 2 title "Position Numerical"'
        WRITE (20, '(A, A)') plot_command_1,  plot_command_2
  

        CLOSE(UNIT=20)
    
        ! Execute the Gnuplot script
        command = 'gnuplot ' // script_filename
        CALL execute_command_line(command)
    END SUBROUTINE

END PROGRAM dem_1D_verification