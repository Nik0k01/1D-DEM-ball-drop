PROGRAM dem_simulation_2D
    ! 
    ! Purpose: 
    !   Simulate the 2d movement of a ball falling from a distance h.

    IMPLICIT NONE

    ! Physical constants, properties
    REAL, PARAMETER :: g = 9.8067                               ! Gravitational constant m/s^2
    REAL, PARAMETER :: k = 5000.0                                ! Stiffness N/m 
    REAL, PARAMETER :: r = 0.2                                  ! Sphere radius m
    REAL, PARAMETER :: m = 0.2                                  ! Sphere mass 
    REAL(8), PARAMETER :: mu = 0.5d0                            ! Friction factor
    REAL(8), PARAMETER :: c = 0.5                               ! Translation damping factor

    ! Simulation parameters
    INTEGER(8), PARAMETER :: steps_no = 10000000                    ! Number of iterations
    REAL(8), PARAMETER :: total_time = 5d0                          ! Total simulation time
    REAL(8), PARAMETER :: time_step = total_time / REAL(steps_no)   ! Time step           

    ! Data dictionary: types, variables, units
    REAL(8) :: init_height                                      ! Initial height
    REAL(8), DIMENSION(steps_no, 2) :: velocity                 ! Velocity vector in x and y direction
    REAL(8), DIMENSION(steps_no, 2) :: position                 ! Position vector in x and y direction
    REAL(8), DIMENSION(steps_no) :: time                        ! Time vector   
    REAL(8), DIMENSION(steps_no, 2) :: net_force                ! Net force vector in x and y direction  
    REAL(8), DIMENSION(steps_no) :: torque                      ! Torque resulting from the friction force  
    REAL(8), DIMENSION(steps_no) :: angular_vel                 ! Angular velocity of the sphere
    REAL(8) :: inertia                                          ! Intertia  
    REAL(8) :: contact_force                                    ! Contact force in x direction
    REAL(8) :: damping_force_x                                  ! Damping force in the x direction
    REAL(8) :: damping_force_y                                  ! Damping force in the y drection
    REAL(8) :: damping_rotation                                 ! Damping torque 
    REAL(8) :: friction_force                                   ! Friction force in y direction
    INTEGER(8) :: i                                             ! Iteration counter                 
    REAL(8), DIMENSION(steps_no, 7) :: result_matrix            ! Result matrix
    CHARACTER(50) :: data_file_name                             ! Name of the file with data for plotting
    REAL(8) :: c_r                                              ! Rotation damping coefficient
    REAL(8), DIMENSION(steps_no, 3) :: marker_position          ! Marker position (x, y) on the disk 

    ! Initial conditions
    init_height = 2.0d0                                         ! Height of the ball
    velocity(1, :) = [0d0 , 10d0]                               ! Velocity in x is zero, and 1 m/s for y direction
    position(1, :) = [init_height, 0d0]                         ! Initial position with respect to x and y axis
    time(1) = 0d0
    net_force(1, 1) = -1.d0 * m * g
    net_force(1, 2) = 0.d0                                      ! No force in y direction
    angular_vel(1) = 0.d0                                       ! Sphere is not rotating in the beginning
    torque(1) = 0.d0                                            ! No tangential force acting on the sphere
    inertia = 2. / 5. * m * r ** 2                              ! Inertia of a sphere
    c_r = 10 * c * r ** 2                                       ! Definition of the damping coefficient for angular movement
    marker_position(1, 1) = time(1)
    marker_position(1, 2) = init_height + r
    marker_position(1, 3) = 0

    ! Main calculation loop
    main_loop: DO i = 2, steps_no

        ! Update the velocity in the x direction
        velocity(i, 1) = velocity(i-1, 1) + net_force(i-1, 1) * time_step / m
        ! Update the velocity in the y direction
        velocity(i, 2) = velocity(i-1, 2) + net_force(i-1, 2) * time_step / m
        IF (velocity(i, 2) < 0) THEN ! In the terminal case the force may be large enough to produce negatve velocity
            velocity(i, 2) = 0.d0
        END IF
        ! Update the position in the x direction
        position(i, 1) = position(i-1, 1) + velocity(i-1, 1) * time_step
        ! Update the position in the y direction 
        position(i, 2) = position(i-1, 2) + velocity(i-1, 2) * time_step
        ! Update the angular velocity 
        angular_vel(i) = angular_vel(i-1) - torque(i-1) / inertia * time_step 
        ! Update the time
        time(i) = time(i-1) + time_step
        ! Calculate the marker position
        marker_position(i, 1) = time(i)             ! Time vector
        marker_position(i, 2) = position(i-1, 1) + (r + 0.5) * COS(angular_vel(i-1) * time(i-1))
        marker_position(i, 3) = position(i-1, 2) + (r + 0.5) * SIN(angular_vel(i-1) * time(i-1))

        ! Check if the contact force is active
        contact: IF (r > position(i, 1)) THEN
            contact_force = k * (r - position(i, 1))
            damping_force_x = - c * velocity(i, 1)                      ! Damping force in the x direction
            ! Friction force cannot excess the maximum value, the second part is the force required to 
            ! stop the movement of the particle in the y direction
            vel_check: IF (velocity(i, 2) > 0.d0) THEN ! The friction force acts only if there is a movement in y velocity
                friction_force = -1.0d0 * min(mu * contact_force, m * velocity(i, 2) / time_step)
                damping_force_y = - c * velocity(i, 2)                  ! Damping force in the y direction
                damping_rotation = - c_r * angular_vel(i)               ! Damping of the particle rotation
                torque(i) = friction_force * r + damping_rotation       ! Net torque acting on the particle
            END IF vel_check
        ELSE
            torque(i) = 0.d0
            contact_force = 0 
            friction_force = 0
            damping_force_x = 0
            damping_force_y = 0
        END IF contact

        ! Calculate the net force in the x and y direction 
        net_force(i, :) = [-1.d0 * m * g + contact_force + damping_force_x, friction_force + damping_force_y]
    END DO main_loop

    ! Write the results in a matrix
    result_matrix(:, 1) = time; result_matrix(:, 2:3) = position
    result_matrix(:, 4:5) = velocity; result_matrix(:, 6:7) = net_force
    CALL file_writer(result_matrix, steps_no, data_file_name, angular_vel)

    ! ! Plot the results 
    ! CALL plot_results(data_file_name)

    ! Generate the gif
    CALL gif_maker(position, marker_position, steps_no)
    
    CONTAINS

    SUBROUTINE file_writer(result_matrix, steps_no, file_name, angular_vel) 
        ! Writes calculation results to a file
        REAL(8), DIMENSION(steps_no, 7), INTENT(IN) :: result_matrix
        REAL(8), DIMENSION(steps_no), INTENT(IN) :: angular_vel
        INTEGER(8), INTENT(IN) :: steps_no
        INTEGER :: ierror                       ! File status
        INTEGER :: file_id                      ! File id
        CHARACTER(80) :: err_string             ! Error message
        CHARACTER(10) :: date                   ! Date
        CHARACTER(12) :: time                   ! Time
        CHARACTER(50), INTENT(OUT) :: file_name ! Name of the file
        INTEGER(8) :: i                         ! Counter 

        ! Generate a file name
        CALL date_and_time(date, time)
        file_name = trim(date) // '_' // trim(time) // '_' // 'results_2D.dat'

        ! Open a file
        OPEN (UNIT=file_id, FILE=file_name, STATUS='NEW', ACTION='WRITE', &
              IOSTAT=ierror, IOMSG=err_string)

        openif: IF (ierror == 0) THEN
            WRITE(file_id, '(A)') "Time, Pos(x), Pos(y), Vel(x), Vel(y), Force(x), Force(y), Angular"
            1010 FORMAT (F10.5, ",", F10.5, ",", F10.5, "," F10.5, "," F10.5, "," F10.5, "," F10.5, F10.5)
            ! Open was ok. Write values
            writeloop: DO i = 1, steps_no, 10
                WRITE(file_id, 1010) result_matrix(i, 1), result_matrix(i, 2), &
                                     result_matrix(i, 3), result_matrix(i, 4), &
                                     result_matrix(i, 5), result_matrix(i, 6), &
                                     result_matrix(i, 7), angular_vel(i)
            END DO writeloop
        
        ELSE openif
            WRITE(*, 1040) ierror
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
        CHARACTER(50), DIMENSION(6) :: plot_names   ! Array with plot names
        CHARACTER(LEN=:), ALLOCATABLE :: command    ! Command made by gnuplot
        CHARACTER(1), DIMENSION(6) :: column_pairs  ! Pairs of columns for plotting
        INTEGER(8) :: i                             ! Loop counter
    
        ! Name of the gnuplot script
        script_filename = 'script.gp'
        
        ! Generate names for the plots of important variables
        CALL date_and_time(date, time)
        base_plot_name = 'Plots/' // trim(date) // '_' // trim(time) 
    
        ! Set up the output plot names
        plot_names(1) = trim(base_plot_name) // '_pos_x_2d_plot.png'
        plot_names(2) = trim(base_plot_name) // '_pos_y_2d_plot.png'
        plot_names(3) = trim(base_plot_name) // '_vel_x_2d_plot.png'
        plot_names(4) = trim(base_plot_name) // '_vel_y_2d_plot.png'
        plot_names(5) = trim(base_plot_name) // '_for_x_2d_plot.png'
        plot_names(6) = trim(base_plot_name) // '_for_y_2d_plot.png'
    
        ! Specify the columns to use for each plot
        column_pairs = ['2', '3', '4', '5', '6', '7'] ! (time:force, time:position, time:velocity)
    
        ! Write Gnuplot commands to the script
        OPEN (UNIT=20, FILE=script_filename, STATUS='REPLACE')
        plot: DO i = 1, 6
            WRITE (20, '(A)') 'set terminal png'
            WRITE (20, '(A)') 'set output "' // trim(plot_names(i)) // '"'
            WRITE (20, '(A)') 'set xlabel "Time (s)"'
            SELECT CASE (i)
                CASE (1)
                    WRITE (20, '(A)') 'set title "Position x vs Time"'
                    WRITE (20, '(A)') 'set ylabel "Position x (m)"'
                CASE (2)
                    WRITE (20, '(A)') 'set title "Position y vs Time"'
                    WRITE (20, '(A)') 'set ylabel "Position y (m)"'
                CASE (3)
                    WRITE (20, '(A)') 'set title "Velocity x vs Time"'
                    WRITE (20, '(A)') 'set ylabel "Velocity x (m/s)"'
                CASE (4)
                    WRITE (20, '(A)') 'set title "Velocity y vs Time"'
                    WRITE (20, '(A)') 'set ylabel "Velocity y (m/s)"'
                CASE (5)
                    WRITE (20, '(A)') 'set title "Force x vs Time"'
                    WRITE (20, '(A)') 'set ylabel "Force x (N)"'
                CASE (6)
                    WRITE (20, '(A)') 'set title "Force y vs Time"'
                    WRITE (20, '(A)') 'set ylabel "Force y (N)"'
            END SELECT
            WRITE (20, '(A, I1, A, I1, A)') 'plot "' // data_file_name // '" using 1:' // column_pairs(i) // ' with lines'
        END DO plot
        CLOSE(UNIT=20)
    
        ! Execute the Gnuplot script
        command = 'gnuplot ' // script_filename
        CALL execute_command_line(command)
    END SUBROUTINE

    SUBROUTINE gif_maker(position, marker, steps_no)
        IMPLICIT NONE
        REAL(8), DIMENSION(steps_no, 2), INTENT(IN) :: position     ! xy position of the centre of mass
        REAL(8), DIMENSION(steps_no, 3), INTENT(IN) :: marker       ! xy position of the marker
        INTEGER(8), INTENT(IN) :: steps_no                          ! Number of steps
        INTEGER(8) :: step                                          ! Loop counter
        CHARACTER(100) :: command                                   ! Gnuplot command
        CHARACTER(50) :: frame_file                                 ! File with the data for the frames
        CHARACTER(50) :: frame_gp                                   ! Temporary gnuplot script file
        INTEGER :: frame_counter                                    ! Frame counter
        REAL :: radius = 0.7
    
        frame_counter = 0
        ! Write frame data into a file
        DO step = 1, steps_no, 10000
            frame_counter = frame_counter + 1
            WRITE(frame_file, '("frame_data_", I4.4, ".dat")') frame_counter
            OPEN (UNIT=90, FILE=frame_file, STATUS='REPLACE', ACTION='WRITE')
            WRITE(90, '(F10.5, F10.5, F10.5, F10.5, F10.5)') position(step, 1), position(step, 2), marker(step, 2), marker(step, 3), radius
            CLOSE(90)
    
            ! Generate gnuplot script
            WRITE(frame_gp, '("frame_", I4.4, ".gp")') frame_counter
            OPEN(UNIT=20, FILE=frame_gp, STATUS='REPLACE')
            WRITE(20, '(A)') 'set terminal pngcairo size 800, 800'
            WRITE(20, '(A, A)') 'set output "', 'frame_' // TRIM(ADJUSTL(frame_gp)) // '.png"'
            WRITE(20, '(A)') 'unset key'
            WRITE(20, '(A)') 'set xrange [0:12]'
            WRITE(20, '(A)') 'set yrange [0:5]'
            WRITE(20, '(A)') 'set size ratio -1'
            WRITE(20, '(A)') 'plot "' // TRIM(ADJUSTL(frame_file)) // '" using 2:1:5 with circles lc rgb "blue" title "Disk",' // &
                              ' "' // TRIM(ADJUSTL(frame_file)) // '" using 4:3 with points pt 7 lc rgb "red" title "Marker"'
            CLOSE(20)
    
            ! Call gnuplot to generate the frame
            command = 'gnuplot ' // TRIM(ADJUSTL(frame_gp))
            CALL execute_command_line(command)
        END DO
    
        ! Generate GIF from PNG frames
        CALL execute_command_line('magick -delay 2 -loop 0 frame_*.png simulation.gif')
    
        ! Clean up temporary files
        CALL execute_command_line('del frame_*.png frame_*.gp frame_data_*.dat')
    
    END SUBROUTINE
    
END PROGRAM dem_simulation_2D