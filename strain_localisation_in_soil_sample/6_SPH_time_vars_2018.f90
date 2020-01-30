
!--------------------------------------------------------


                 MODULE SPH_time_vars_2018 
                 

!--------------------------------------------------------
                      
use variable_types

implicit none

real (kind=DP) :: time, dt		! time increment, time, 
real (kind=DP) :: time_sph, dt_sph	! local time and dt for sph

integer:: current_ts		! absolute number of time step. Accounts several restarts
integer:: nstart		! number of timesteps in previous cycle
integer:: itimestep		! local time step in present restart cycle
integer:: itimestep_sph	! same for sph (there is an inner loop)

real (kind=DP) :: dt_limit_inf	! minimum allowed dt, set to 1.e-7*dt
real (kind=DP) :: time_end		! max time for computation
integer :: maxtimestep	! maximum number of timesteps in computation

integer :: print_step     	!  Print every print_step    (On Screen)
integer :: save_step      	!  Save  every  save_step    (To Disk File)
integer :: plot_step      	!  Plot  every  plot_step    (For GID v.7)
integer :: moni_part      	!  Number of control particle (only one)

real   :: time_plot		! we plot every time_plot time
real   :: t_plot_reset	! time counter reset to zero after plotting
real   :: time_save		! we save every time_plot time
real   :: t_save_reset	! time counter reset to zero after saving
real   :: time_print		! we print every time_print time
real   :: t_print_reset	! time counter reset to zero after printing
						
integer:: ntcurves 				!     number of time curves
integer:: mptstcurves			!     max number of time pts in all t curves
integer, allocatable:: nptstcurves(:)	!     nr. of time pts in each curve
real (kind=DP), allocatable::   ttcurves(:,:)	!     time values for all time curves
real, allocatable::   ftcurves(:,:)	!     factors at all times. 
						!     Multiply a0 in prer variable
END MODULE SPH_time_vars_2018 
