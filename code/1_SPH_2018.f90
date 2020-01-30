             
!-------------------------------------------------------------------------------
       
                      PROGRAM SPH_2018                

!-------------------------------------------------------------------------------
!contains simulation evolution and output controls   
!-------------------------------------------------------------------------------

use variable_types       

use SPH_time_vars_2018      

use SPH_global_vars_2018 
                                
                                
use SPH_main_2018,  only:  Init_SPH,              	&
                           time_integration,  	& 
                           Out_plot_sph,          	&
                           Out_print_sph,         	&  
                           clean_up_sph                                         

use SPH_material_2018



implicit none

input_file = 997
chk_file = 998

!open simulation input and output files
!-------------------------------------------------------------------------------

open(input_file,file='input.txt')
open(output_file,file='output.txt')

!initialise SPH
!-------------------------------------------------------------------------------

call Init_sph                 
                
!initialise time stepping
!-------------------------------------------------------------------------------

nstart   = 0.0; current_ts = 0        
time     = 0.0; dt = 0.001 
time_sph = 0.0; itimestep_sph = 0 

DTGT0: do while (dt > 0.0)

        read (input_file,*)  text
        read (input_file,*)  dt, time_end, maxtimestep 
        write(chk_file,*)  'dt, time_end, maxtimestep'  
        write(chk_file,*)    dt, time_end, maxtimestep

        if (dt.le.0) EXIT DTGT0
        read(input_file,*) text
        read (input_file,*) print_step, save_step, plot_step 
        write(chk_file,*)  'print_step, save_step, plot_step'
        write(chk_file,*)   print_step, save_step, plot_step 

        dt_sph = dt
   
        call OutputRes !outputs initial results (at time zero)
    
        time_print    = print_step*dt; time_save  = save_step*dt; time_plot  = plot_step*dt
        t_print_reset = 0.; t_save_reset = 0.; t_plot_reset = 0. 
   
        dt_limit_inf  = 1.e-7*dt  ! minimum allowed time step
   
        do itimestep = nstart+1, nstart + maxtimestep   
             
                current_ts    = current_ts + 1      
                t_print_reset = t_print_reset + dt
                t_plot_reset  = t_plot_reset  + dt
                t_save_reset  = t_save_reset  + dt
      
                if (dt.lt.0.0.or.dt.lt.dt_limit_inf) EXIT                  
	
                itimestep_sph = itimestep_sph +1                    
                call time_integration
                time_sph = time_sph + dt_sph   
       
                time = time + dt
                
                if (time > time_end) EXIT 
   
                !output files
   
                if (t_plot_reset >= time_plot) then
                        call out_plot                      
                        t_plot_reset = 0.0 
                end if 
   
                if (t_print_reset >= time_print) then
                        call out_prn                       
                        t_print_reset = 0.0 
                end if
            
   
        end do  

        nstart = current_ts

end do DTGT0

call clean_up_sph   
   
CONTAINS

!-------------------------------------------------------------------

subroutine out_plot  

!Output sph results files     
!-------------------------------------------------------------------

implicit none

call OutputRes
 
end subroutine out_plot 

!-------------------------------------------------------------------

subroutine out_prn  
 
!Output sph data for debugging/monitoring          
!-------------------------------------------------------------------

implicit none

call out_print_sph

end subroutine out_prn

END PROGRAM SPH_2018

