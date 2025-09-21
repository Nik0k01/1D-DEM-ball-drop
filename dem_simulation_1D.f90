PROGRAM dem_simulation_1D
    ! 
    ! Purpose: 
    !   Simulate the 1d movement of a ball falling from a distance h.

    IMPLICIT NONE

    ! Physical constants, properties
    REAL, PARAMETER :: g = 9.8067                           ! Gravitational constant m/s^2
    REAL, PARAMETER :: k = 1200.0                           ! Stiffness N/m 
    REAL, PARAMETER :: r = 0.1                              ! Sphere radius m
    REAL, PARAMETER :: m = 0.5                              ! Sphere mass kg
    REAL, PARAMETER :: c = 0.1                              ! Viscous damping force

    ! Simulation parameters
    INTEGER(8), PARAMETER :: steps_no = 10000000            ! Number of iterations
    REAL(8), PARAMETER :: total_time = 5d0                  ! Total simulation time
    REAL(8), PARAMETER :: time_step = total_time / REAL(steps_no) ! Time step           

    ! Data dictionary: types, variables, units
    REAL(8) :: init_height                                  ! Initial height
    REAL(8), DIMENSION(steps_no) :: velocity                ! Velocity vector
    REAL(8), DIMENSION(steps_no) :: position                ! Position vector
    REAL(8), DIMENSION(steps_no) :: time                    ! Time vector
    REAL(8), DIMENSION(steps_no) :: net_force               ! Net force vector
    REAL(8), DIMENSION(steps_no) :: overlap                 ! Overlap
    REAL(8) :: contact_force                                ! Contact force
    REAL(8) :: damping_force                                ! Viscous drag force, opposed to the motion of the particle 
    INTEGER(8) :: i                                         ! Iteration counter                 
    REAL(8), DIMENSION(steps_no, 5) :: result_matrix        ! Result matrix
    CHARACTER(50) :: data_file_name                         ! Name of the file with data for plotting

    ! Initial conditions
    init_height = 1.0d0                                     ! Ball is at specific height
    velocity(1) = 0d0                                       ! Ball is at rest before dropping
    position(1) = init_height                               ! First position is its initial height
    time(1) = 0d0
    net_force(1) = -1.d0 * m * g                            ! Negative gravity force (Jump will be against the gravity)
    overlap(1) = 0
    
    ! Main calculation loop
    main_loop: DO i = 2, steps_no

        ! Update the velocity
        velocity(i) = velocity(i-1) + net_force(i-1) * time_step / m
        ! Update the position
        position(i) = position(i-1) + velocity(i-1) * time_step
        ! Update the time
        time(i) = time(i-1) + time_step

        ! Check if the contact force is active
        contact: IF (r > position(i)) THEN
            contact_force = k * (r - position(i))           ! Contact is detected if the radius is below the z axis
            ! Implement the damping force 
            damping_force = - c * velocity(i)               ! Spring contstant - overlap (Radius and the distance of displaced ball to the center of wall)
            ! Calculate the overlap
            overlap(i) = r - position(i)                    ! Damping coefficient must be greater than 0
        ELSE
            overlap(i) = 0
            contact_force = 0 
            damping_force = 0
        END IF contact

        ! Calculate the net force 
        net_force(i) = -1.d0 * m * g + contact_force + damping_force        ! Add up all of the forces
    END DO main_loop

    ! Write the results in a matrix
    result_matrix(:, 1) = time; result_matrix(:, 2) = position              ! Creating matrix to see and plot the results 
    result_matrix(:, 3) = velocity; result_matrix(:, 4) = net_force
    result_matrix(:, 5) = overlap
    CALL file_writer(result_matrix, steps_no, data_file_name)               ! Creating a file for result

    ! Plot the results 
    CALL plot_results(data_file_name)
    
    ! Plot every energy type in the system
    CALL energy_conservation(result_matrix, steps_no, m, g, c, k, r)

    CONTAINS

    SUBROUTINE file_writer(result_matrix, steps_no, file_name) 
        ! Writes calculation results to a file
        REAL(8), DIMENSION(steps_no, 5), INTENT(IN) :: result_matrix
        INTEGER(8), INTENT(IN) :: steps_no
        INTEGER :: ierror                       ! File status
        INTEGER :: file_id                      ! File id
        CHARACTER(80) :: err_string             ! Error message
        CHARACTER(10) :: date                   ! Date
        CHARACTER(12) :: time                   ! Time
        CHARACTER(50), INTENT(OUT) :: file_name ! Name of the file
        INTEGER(8) :: i                         ! Counter 

        ! Generate a file name
        CALL date_and_time(date, time)          ! Writing a dynamic file name by concatenating the current date and time with a fixed suffix 'results.dat'
        file_name = trim(date) // '_' // trim(time) // '_' // 'results.dat'

        ! Open a file
        OPEN (UNIT=file_id, FILE=file_name, STATUS='NEW', ACTION='WRITE', &
              IOSTAT=ierror, IOMSG=err_string)      ! If the operation is successful, `ierror` is set to 0; otherwise, it holds a non-zero error code
        
        ! Check if the file was opened correctly
        openif: IF (ierror == 0) THEN
            1010 FORMAT (F10.5, ",", F10.5, ",", F10.5, ",", F10.5, "," F10.5)
            ! Open was ok. Write values
            writeloop: DO i = 1, steps_no, 50           ! Write to a file every 50th result
                WRITE(file_id, 1010) result_matrix(i, 1), result_matrix(i, 2), &
                                     result_matrix(i, 3), result_matrix(i, 4), &
                                     result_matrix(i, 5)
            END DO writeloop
        ! If openining of the file was not correct
        ELSE openif
            WRITE(*, 1040) ierror                       ! Writes ierror if there is error in file opening
            1040 FORMAT ('Error openin file: IOSTAT = ', I6) 
            WRITE(*, 1050) TRIM(err_string)
            1050 FORMAT (A)  
        END IF openif 
        
        ! Close file
        CLOSE (UNIT=file_id)
    END SUBROUTINE

    SUBROUTINE plot_results(data_file_name)
        ! Generate a Gnuplot script and plot the results
        CHARACTER(50), INTENT(IN) :: data_file_name
        CHARACTER(10) :: date                       ! Date
        CHARACTER(12) :: time                       ! Time
        CHARACTER(9) :: script_filename             ! Gnuplot script name
        CHARACTER(50) :: base_plot_name             ! Common base name for each plot
        CHARACTER(50), DIMENSION(5) :: plot_names   ! Array with plot names
        CHARACTER(LEN=:), ALLOCATABLE :: command    ! Command made by Gnuplot
        INTEGER :: i                                ! Loop counter
    
        ! Name of the Gnuplot script
        script_filename = 'script.gp'
        
        ! Generate names for the plots of important variables
        CALL date_and_time(date, time)
        base_plot_name = 'Plots/' // trim(date) // '_' // trim(time) 
        
        ! Set up the output plot names
        plot_names(1) = trim(base_plot_name) // '_force_plot.png'           ! Generating png files for graphs
        plot_names(2) = trim(base_plot_name) // '_position_plot.png'
        plot_names(3) = trim(base_plot_name) // '_velocity_plot.png'
        plot_names(4) = trim(base_plot_name) // '_overlap_plot.png'
        plot_names(5) = trim(base_plot_name) // '_overlap_force_plot.png'
         
        ! Write Gnuplot commands to the script
        OPEN (UNIT=20, FILE=script_filename, STATUS='REPLACE')
        WRITE (20, '(A)') 'set terminal pngcairo size 800,600 enhanced font "Arial,14"'
        WRITE (20, '(A)') 'set grid'
        WRITE (20, '(A)') 'set style line 1 lt 1 lw 2 lc rgb "blue"'
        WRITE (20, '(A)') 'set style line 2 lt 2 lw 2 lc rgb "red"'
        WRITE (20, '(A)') 'set style line 3 lt 3 lw 2 lc rgb "green"'
    
        plot: DO i = 1, 5
            WRITE (20, '(A)') 'set output "' // plot_names(i) // '"'
            WRITE (20, '(A)') 'set xlabel "Time (s)"'
            SELECT CASE (i)
                CASE (1)
                    WRITE (20, '(A)') 'set title "Force vs Time"'           ! Plot force vs time
                    WRITE (20, '(A)') 'set ylabel "Force (N)"'
                    WRITE (20, '(A)') 'plot "' // data_file_name // '" using 1:4 with lines ls 1 title "Force"'
                CASE (2)
                    WRITE (20, '(A)') 'set title "Position vs Time"'        ! Plot position vs time
                    WRITE (20, '(A)') 'set ylabel "Position (m)"'
                    WRITE (20, '(A)') 'plot "' // data_file_name // '" using 1:2 with lines ls 2 title "Position"'
                CASE (3)
                    WRITE (20, '(A)') 'set title "Velocity vs Time"'        ! Plot velocity vs time
                    WRITE (20, '(A)') 'set ylabel "Velocity (m/s)"'
                    WRITE (20, '(A)') 'plot "' // data_file_name // '" using 1:3 with lines ls 3 title "Velocity"'
                CASE (4)
                    WRITE (20, '(A)') 'set title "Overlap vs Time"'         ! Plot overlap vs time
                    WRITE (20, '(A)') 'set ylabel "Overlap (m)"'
                    WRITE (20, '(A)') 'plot "' // data_file_name // '" using 1:5 with lines ls 3 title "Overlap"'
                CASE (5)
                    WRITE (20, '(A)') 'set title "Force vs Overlap"'
                    WRITE (20, '(A)') 'set xlabel "Overlap (m)"'
                    WRITE (20, '(A)') 'set ylabel "Force (N)"'
                    WRITE (20, '(A)') 'plot "' // data_file_name // '" using 5:4 with lines ls 3 title "Overlap"'
            END SELECT
        END DO plot
        CLOSE(UNIT=20)
    
        ! Execute the Gnuplot script
        command = 'gnuplot ' // script_filename
        CALL execute_command_line(command)
    END SUBROUTINE
        
    SUBROUTINE energy_conservation(result_matrix, steps_no, m, g, c, k, r)
        ! Purpose:
        !   Plot the total energy of the system vs time to check the energy conservation law
        INTEGER(8), INTENT(IN) :: steps_no                                  ! Number of steps 
        REAL(8), DIMENSION(steps_no, 5), INTENT(IN) :: result_matrix        ! Result matrix

        ! Physical constants, properties
        REAL, INTENT(IN) :: g                                               ! Gravitational constant m/s^2
        REAL, INTENT(IN) :: k                                               ! Stiffness N/m 
        REAL, INTENT(IN) :: r                                               ! Sphere radius m
        REAL, INTENT(IN) :: m                                               ! Sphere mass kg
        REAL, INTENT(IN) :: c                                               ! Viscous damping force

        ! Data dictionary
        REAL(8), DIMENSION(steps_no) :: kinetic_energy                      ! Kinetic energy
        REAL(8), DIMENSION(steps_no) :: potential_energy                    ! Potential energy
        REAL(8), DIMENSION(steps_no) :: elastic_energy                      ! Elastic energy
        REAL(8), DIMENSION(steps_no) :: dissipated_energy                   ! Dissipated energy
        REAL(8), DIMENSION(steps_no) :: total_energy                        ! Total energy
        REAL(8), DIMENSION(steps_no, 6) :: energy_matrix                    ! Total energy
        INTEGER(8) :: i, j                                                  ! Loop counter
        CHARACTER(30) :: command                                            ! Gnuplot command

        ! Initialize the values of energies
        potential_energy(1) = position(1) * m * g
        kinetic_energy(1) = velocity(1) ** 2 * m / 2
        elastic_energy(1) = 0d0
        dissipated_energy(1) = 0d0
 
        ! For each time step calculate the value of each type of energy, dissipated energy can only increase - work done
        DO i = 2, steps_no
            potential_energy(i) = position(i) * m * g                                                      ! Potential energy
            kinetic_energy(i) = velocity(i) ** 2 * m / 2                                                ! Kinetic energy
            IF (r > result_matrix(i, 2)) THEN
                elastic_energy(i) = 0.5 * k * (r - result_matrix(i, 2)) ** 2                            ! Energy stored in the spring
                dissipated_energy(i) = dissipated_energy(i-1) + c * result_matrix(i, 3) ** 2 &          ! Work done by theviscous forces
                                     * (result_matrix(i, 1) - result_matrix(i-1, 1))
            ELSE
                elastic_energy(i) = 0
                dissipated_energy(i) = dissipated_energy(i-1)
            END IF
        END DO
        ! Calculate the total energy of the system
        total_energy = potential_energy + kinetic_energy + elastic_energy + dissipated_energy   
        energy_matrix(:, 1) = result_matrix(:, 1); energy_matrix(:, 2) = potential_energy
        energy_matrix(:, 3) = kinetic_energy; energy_matrix(:, 4) = elastic_energy
        energy_matrix(:, 5) = dissipated_energy; energy_matrix(:, 6) = total_energy
        
        ! Create a new data for the results
        OPEN (3, FILE='Energy_data.dat', STATUS='REPLACE', ACTION='WRITE')
        ! Write the data to the file
        1015 FORMAT (F10.5, ",", F10.5, ",", F10.5, ",", F10.5, "," F10.5, "," F10.5)
        DO i = 1, steps_no, 50
            WRITE(3, 1015) (energy_matrix(i, j), j = 1, 6)
        END DO
        ! Close the file 
        CLOSE(3)

        ! Plot the results
        OPEN (20, FILE='gnuplot.gp', STATUS='REPLACE')
        WRITE (20, '(A)') 'set output "Plots/Energy_plot.png"'
        WRITE (20, '(A)') 'set terminal pngcairo size 800,600 enhanced font "Arial,14"'
        WRITE (20, '(A)') 'set grid'
        WRITE (20, '(A)') 'set style line 1 lt 1 lw 2 lc rgb "blue"'
        WRITE (20, '(A)') 'set style line 2 lt 2 lw 2 lc rgb "red"'
        WRITE (20, '(A)') 'set style line 3 lt 3 lw 2 lc rgb "green"'
        WRITE (20, '(A)') 'set style line 4 lt 3 lw 2 lc rgb "black"'
        WRITE (20, '(A)') 'set style line 5 lt 3 lw 2 lc rgb "pink"'
        WRITE (20, '(A)') 'set xlabel "Time (s)"'
        WRITE (20, '(A)') 'set ylabel "Energy (J)"'
        WRITE (20, '(A)') 'set title "Energy vs Time"'
        WRITE (20, '(A)') 'plot "Energy_data.dat" using 1:2 with lines ls 1 title "Potential",' // & 
                          ' "Energy_data.dat" using 1:3 with lines ls 2 title "Kinetic",' // &
                          ' "Energy_data.dat" using 1:4 with lines ls 3 title "Elastic",' // &
                          ' "Energy_data.dat" using 1:5 with lines ls 4 title "Dispersed",' // &
                          ' "Energy_data.dat" using 1:6 with lines ls 5 title "Total",'
        CLOSE(20)

        command = 'gnuplot gnuplot.gp'
        CALL execute_command_line(command)
    END SUBROUTINE
END PROGRAM dem_simulation_1D