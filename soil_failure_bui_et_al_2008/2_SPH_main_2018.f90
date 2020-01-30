!-------------------------------------------------------------------------------             
!-------------------------------------------------------------------------------
      
              MODULE SPH_main_2018

! ------------------------------------------------------------------------------ 
! ------------------------------------------------------------------------------
!contains main SPH routine, including time integration, SPH discretisation and
!node--stress-point interpolation
! ------------------------------------------------------------------------------               
! ------------------------------------------------------------------------------ 

use variable_types        !double precision defined here.
use SPH_time_vars_2018      !contains time variables 
use SPH_global_Vars_2018     !contains main SPH variables 

use SPH_material_2018    !Problem initialisation and problem-specific routines (as well as other miscellaneous routines)

PRIVATE       !No variable is accesible from outside the module, except for the public routines
             
public:: Init_SPH            
public:: time_integration
public:: Out_plot_sph
public:: Out_print_sph
public:: clean_up_sph 

CONTAINS

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

        SUBROUTINE Init_sph 

!initialises SPH simulation
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

implicit none
       
call problem_input_data  !set up problem (located in 3_SPH_material_2018.f90)

x00 = x0    !------  Store initial particle position 
       
call OutputMesh  !mesh for GiD visualisation purposes
       
        END SUBROUTINE Init_sph    

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

       SUBROUTINE Out_print_sph

!prints simulation info to output file for monitoring/debugging purposes
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

implicit none

integer :: ipoin
real (kind=DP) :: total_mass

total_mass=0.0
do ipoin = 1,ntotal
        total_mass = total_mass + mass(ipoin)
end do

write(output_file,*)'.****  SPH ******........................................'
write(output_file,*)' time step_sph is ', itimestep_sph, ' time_sph =', time_sph, '  dt_sph=', dt_sph , ' MASS = ', total_mass   
write(output_file,*) 'x','velocity'      
write(output_file,*)  x(1:ndimn,nnode),vel(1:ndimn,nnode) 
write(output_file,*)'.........................................'
         
        END SUBROUTINE Out_print_sph     

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

        SUBROUTINE time_integration

!updates problem variables in time
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

implicit none

real (kind=DP) :: t0 

integer :: i, j, idimn, iamaTG, k, ncrit
real :: vel_half(ndimn,ntotal)

ncrit = props(1,2)

call Check_Out_Domain 
 
call grid_find_new (itimestep_sph)

call Pint_Update   
  
if (SPH_shift) then
       if (itimestep_sph > 1) then
               if (mod(itimestep_sph-1,shift_update) == 0) then

                       call stress_point_update 
      
                       if (ncrit == 12) call adapt_stress2 
                
               end if
      end if
end if

!store variables at beginning of time step
!-------------------------------------------------------------------------------

t0 = time_sph; vx0(:,1:ntotal) = vel(:,1:ntotal); x0(:,1:ntotal) = x(:,1:ntotal)
vel_half = 0.0

!RK4 discretisation loop             
!-------------------------------------------------------------------------------

call RK4

!get deviatoric plastic strain 
!------------------------------------------------------------------------------- 
 
call update_strain(dt)
  
!update velocity on stress-points and stress on nodes
!------------------------------------------------------------------------------- 

call stress_point_update 
      
if (SP_SPH) then !only necessary to adapt stress and bcs if stress-particle SPH is used
        if (ncrit == 12) call adapt_stress2 !adapt stress in elastoplastic model, if stress-particle SPH is used
        if (no_bcs > 0) call BCs !update fixed boundary conditions
end if

!updates particle positions
!-------------------------------------------------------------------------------
	
if (update_x) then

        if (XSPH) then        
                call XSPH_update    !XSPH particle position update            
        else        
                vel_half = 0.5*(vx0(:,1:ntotal) + vel(:,1:ntotal))
                x(:,1:ntotal) = x0(:,1:ntotal)+ (vel_half(:,1:ntotal))*dt_sph  !standard particle position update             
        end if
                            
        !locate free surface nodes
        !-----------------------------------------------------------------------
        
        if (ndimn == 2) then
                call get_nodes_on_free_surface
        end if
        	
        !update stress-point positions in the outside approach
        !-----------------------------------------------------------------------	
        	
        if ((SP_SPH) .and. (.not. inside_approach)) then 	        		                
                call shift_stress_points !shift stress_points after a specified number of time steps                                      
        end if
        
        !ensure stress-points = nodes for Standard SPH (should be true anyway)
        !-----------------------------------------------------------------------
        
        if (.not. SP_SPH) x(:,nnode+1:ntotal) = x(:,1:nnode)
        
        !total displacement so far
        !-----------------------------------------------------------------------
                                                
	displ(:,1:nnode) = x(:,1:nnode) - x00(:,1:nnode)          
       
        !update boundary conditions
        !-----------------------------------------------------------------------
        !apply stress free 
        !if (ncrit == 12) call adapt_stress2			
                          		
else 	
		
	displ(:,1:nnode) = displ(:,1:nnode) + 0.5*(vx0(:,1:nnode)+vel(:,1:nnode))*dt_sph 
                        
end if
			
END SUBROUTINE  time_integration

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

        SUBROUTINE XSPH_update

!XSPH routine to update particle velocites 
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

implicit none

real (kind=DP) :: eps
real (kind=DP) :: vel_sum(ndimn,ntotal), x0(ndimn,ntotal)
integer:: i, j

 vel_sum = 0.0; x0 = x(:,1:ntotal)
 eps = 0.5
 n_int = 0.0
 conc = 0.0

 current => last
do while (associated(current))
               i = current%pair_i !stress point
               j = current%pair_j !node 
                
               if (current%pint_type == 2) then !stress-point-stress-point
               
                       vel_sum(:,i) = vel_sum(:,i) + (mass(j)/rho(j))*(vel(:,j) - vel(:,i))*current%w
                       vel_sum(:,j) = vel_sum(:,j) + (mass(i)/rho(i))*(vel(:,i) - vel(:,j))*current%w                        
                       
               else if (current%pint_type == 3) then !node-node
               
                       vel_sum(:,j) = vel_sum(:,j) + (mass(i)/rho(i))*(vel(:,i) - vel(:,j))*current%w
                       vel_sum(:,i) = vel_sum(:,i) + (mass(j)/rho(j))*(vel(:,j) - vel(:,i))*current%w
                        
                       n_int(i) = n_int(i) + 1 !number of node-node interactions
                       n_int(j) = n_int(j) + 1
                        
        		if (bc_or_not(i) == 2 .and. bc_or_not(j) /= 2) then
        			bc_or_not(j) = 3
        		end if

			if (bc_or_not(j) == 2 .and. bc_or_not(i) /= 2) then
				bc_or_not(i) = 3
			end if
   
               end if                
               current => current%next
end do

x(:,1:nnode) = x0(:,1:nnode) + dt_sph*(vel(:,1:nnode)+eps*vel_sum(:,1:nnode))
x(:,nnode+1:ntotal) = x0(:,nnode+1:ntotal) + dt_sph*(vel(:,nnode+1:ntotal)+eps*vel_sum(:,nnode+1:ntotal))
 
        END SUBROUTINE XSPH_update

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

        SUBROUTINE shift_stress_points

!updates the position of the stress-points with respect to the updated node positions
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

implicit none

integer :: i,j,k
integer :: clumping(nnode), tensile(nnode), varI1(nnode)
real :: rij, disp_0, cos_theta(nnode), sin_theta(nnode) 
real :: r_x2(ndimn,nnode)
logical :: switch

if (mod(itimestep_sph,shift_update) == 0) then !shift stress-points every specified number of time steps 

!if velocity vector method is used to orientate stress-point positions in the outside approach
!-------------------------------------------------------------------------------
if (vel_vector) then
 
        disp_x(1,:) = x(1,1:nnode) - x_10(1,1:nnode)
        disp_x(2,:) = x(2,1:nnode) - x_10(2,1:nnode)
        disp_10 = sqrt((x(1,1:nnode) - x_10(1,1:nnode))**2 + (x(2,1:nnode) - x_10(2,1:nnode))**2)
        disp_10 = disp_10/dx

        x_10 = x(:,1:nnode)
  
        cos_theta = vel(1,1:nnode)/sqrt((vel(1,1:nnode))**2 + (vel(2,1:nnode))**2)
        sin_theta = vel(2,1:nnode)/sqrt((vel(1,1:nnode))**2 + (vel(2,1:nnode))**2)
 
        r_x2(1,:) = (dx/2)*cos_theta
        r_x2(2,:) = (dx/2)*sin_theta
        
        k = nnode+1
        do i = 1,nnode
                if (r_x2(1,i) > 0 .and. r_x2(1,i) < dx/5) r_x2(1,i) = r_x2(1,i) + dx/3.
                if (r_x2(2,i) > 0 .and. r_x2(2,i) < dx/5) r_x2(2,i) = r_x2(2,i) + dx/3.
                if (r_x2(1,i) < 0 .and. r_x2(1,i) > -dx/5) r_x2(1,i) = r_x2(1,i) - dx/3.
                if (r_x2(2,i) < 0 .and. r_x2(2,i) > -dx/5) r_x2(2,i) = r_x2(2,i) - dx/3.
           
                if (disp_10(i) > disp_tol) then 
                        x(1,k) = x(1,i) + r_x2(1,i)
                        x(2,k) = x(2,i) + r_x2(2,i)
                        x(1,k+1) = x(1,i) - r_x2(1,i)
                        x(2,k+1) = x(2,i) - r_x2(2,i)
                end if             
                k = k+2           
        end do

!standard outside approach *** take care with stress-point orientation ***
!-------------------------------------------------------------------------------         
else
 
        if (npoints == 1) then
                    
                x(1,nnode+1:ntotal) = x(1,1:nnode) + r_x
                x(2,nnode+1:ntotal) = x(2,1:nnode) + r_y
        
        else if (npoints == 2) then
    
                k = nnode+1
                do i = 1,nnode                     
                        x(1,k) = x(1,i) + r_x
                        x(2,k) = x(2,i) + r_y
                        x(1,k+1) = x(1,i) - r_x
                        x(2,k+1) = x(2,i) - r_x
                        
                        k = k+2        
                end do         
       
        else if (npoints == 3) then
    
                k = nnode+1
                do i = 1,nnode        
                        x(1,k) = x(1,i) 
                        x(2,k) = x(2,i) + r_y
                        x(1,k+1) = x(1,i) - r_x
                        x(2,k+1) = x(2,i) - r_y
                        x(1,k+2) = x(1,i) + r_x
                        x(2,k+2) = x(2,i) - r_y
           
                        k = k+3           
                end do       
       
   end if  
   
end if

end if

!-------------------------------------------------------------------------------
!adapt the stress-points so that they align with their associated nodes (essentially turning them off)
!-------------------------------------------------------------------------------

!adapts them in 2 cases: 
!1) near the boundary (for nodes travelling faster than a certain velocity)
!2) for isolated nodes (those with neighbours less than 2)
                
!isolated nodes are calculated in artificial viscosity and XSPH subroutines (as node-node interactions are already considered there)
!if these aren't used, call isolated nodes:
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

if ( (alpha == 0 .and. beta == 0) .and. (.not. XSPH) ) call isolated_nodes 
                
abs_vel = sqrt(vel(1,1:nnode)**2 + vel(2,1:nnode)**2)  
                 
k = nnode+1
do i = 1,nnode
        if (bc_int(i) == 1) then 
                if (abs_vel(i) > 0.4) then               
                        x(:,k) = x(:,i) 
                        if (npoints > 1) x(:,k+1) = x(:,i)
                        if (npoints > 2) x(:,k+2) = x(:,i)
                end if
        end if
        if (n_int(i) < 2) then               
                x(:,k) = x(:,i) 
                if (npoints > 1) x(:,k+1) = x(:,i)
                if (npoints > 2) x(:,k+2) = x(:,i)
        end if
        k=k+npoints   
end do
             
        END SUBROUTINE

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

        SUBROUTINE isolated_nodes

!calculates number of node-node interactions to determine any isolated nodes
!(not necessary if artificial viscosity or XSPH is used)
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

implicit none

integer:: i, j

n_int = 0.0
 
 current => last
do while (associated(current))
        i = current%pair_i !stress point
        j = current%pair_j !node 

        if (current%pint_type == 3) then !node-node
                n_int(i) = n_int(i) + 1
                n_int(j) = n_int(j) + 1
        end if                
        current => current%next
 end do
 
        END SUBROUTINE isolated_nodes

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

        SUBROUTINE stress_point_update

!interpolates the stresses on the nodes from the stress-points, and the velocities on
!the stress-points from the nodes
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

implicit none

integer :: i, j, idimn, istre
real (kind=DP) :: vel_temp(ndimn,ntotal), stress_temp(nstre,ntotal), internal_vars_temp(1,ntotal)  !temporary storage for interpolation calculations
real (kind=DP) :: rho_temp(ntotal), hsml_temp(ntotal), cspm_norm(ntotal), h1, h2

!interpolates variables via SPH-CSPM method
!-------------------------------------------------------------------------------

if (SP_SPH) then
        !initialise temporary variables   
        stress_temp = 0; rho_temp = 0; hsml_temp = 0; vel_temp = 0
        internal_vars_temp = 0; cspm_norm = 0
                       
        current => last
        do while (associated(current))
                j = current%pair_i !stress point
                i = current%pair_j !node 
                
                if (current%pint_type == 1) then                
                        h1 = (mass(j)/rho(j))*current%w
                        h2 = (mass(i)/rho(i))*current%w                               
                        
                        vel_temp(:,j) = vel_temp(:,j) + vel(:,i)*h2 !SPH velocity approximation                        
                        stress_temp(:,i) = stress_temp(:,i) + stress(:,j)*h1 !SPH stress approximation
                        internal_vars_temp(1,i) = internal_vars_temp(1,i) + internal_vars(1,j)*h1 !SPH dev. strain approximation
                         
                        if (cont_density) rho_temp(i) = rho_temp(i) + rho(j)*h1 !SPH density approximation   
                                
                        cspm_norm(j) = cspm_norm(j) + (current%w*mass(i))/rho(i) !cspm normalisation
                        cspm_norm(i) = cspm_norm(i) + (current%w*mass(j))/rho(j)
                                                                                                                                   
                end if
                current => current%next
        end do

        do i = nnode+1,ntotal                
                if (cspm_norm(i) == 0) cycle !don't do the interpolation with normalisation if there aren't enough interpolation points
                vel(1,i) = vel_temp(1,i)/cspm_norm(i)
                vel(2,i) = vel_temp(2,i)/cspm_norm(i)
        end do
        
        do i = 1,nnode
                if (cspm_norm(i) /= 0) then
                        stress(1,i) = stress_temp(1,i)/cspm_norm(i)
                        stress(2,i) = stress_temp(2,i)/cspm_norm(i)
                        stress(3,i) = stress_temp(3,i)/cspm_norm(i)
                        stress(4,i) = stress_temp(4,i)/cspm_norm(i)
                        internal_vars(1,i) = internal_vars_temp(1,i)/cspm_norm(i)
                else
                        vel(:,i) = 0 !if there aren't any interacting stress-points, make the node stop moving (only necessary for cohesive soil simulations with the inside appraoch)
                end if
        end do
                
        if (cont_density) then                        
                rho(1:nnode) = rho_temp(1:nnode)/cspm_norm(1:nnode)
        !        if (sle == 2) then
        !                hsml(k) = hsml_temp(k)/no_int_stress(k)
        !        end if
        end if
         
         
else if (.not. SP_SPH) then !Standard SPH: nodes = sps
                                
                vel(:,nnode+1:ntotal) = vel(:,1:nnode)
                stress(:,1:nnode) = stress(:,nnode+1:ntotal)
                internal_vars(1,1:nnode) = internal_vars(1,nnode+1:ntotal)   
                
                if (cont_density) rho(1:nnode) = rho(nnode+1:ntotal)
                          
end if

        END SUBROUTINE stress_point_update

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

        SUBROUTINE get_derivatives(fi1,div1,fi2,div2)

!calculates the divergence of f_u (f1) and f_sigma (f2) 
! first calculates the gradients of velocity and stress
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

implicit none

real (kind=DP), intent(in) :: fi1(:,:), fi2(:,:) !input: variables (stress and velocity)
real (kind=DP), intent(out) :: div1(:,:), div2(:,:) !output: divergence of f_u and f_sigma (functions of stress and velocity)
integer :: i, j, iamaTG, idimn, istre
real (kind=DP) :: grad1_tmp(ndimn,ndimn,ntotal), grad2_tmp(nstre-1,ndimn,ntotal), h1, h2, g1, g2 !gradient terms 
real (kind=DP) :: AE(5,ntotal), h1b, h2b  !CSPM terms
real (kind=DP) :: da, db, beta, beta_max, vel_wall, wall, dummy_vel(ndimn,ntotal2) !dummy node parameters

grad1_tmp = 0.0; grad2_tmp = 0.0; div1  = 0.0; div2 = 0.0 !gradient and divergence terms
AE = 0.0 !CSPM matrix
dummy_vel = 0.0 !dummy node velocity
bc_int = 0 !node--dummy node interactions

!loops through all interacting pairs of particles
!-------------------------------------------------------------------------------			         			           			        
 current => last 
do while (associated(current))          			      
        j = current%pair_i !stress-point
        i = current%pair_j !node                                       
        h1 = current%dwdx*mass(i)/rho(i)        
        if (current%pint_type == 1) then !node--stress-point
         
                !calculate velocity gradient on stress-points
                !---------------------------------------------------------------      
                do idimn = 1,ndimn                                
                        grad1_tmp(idimn,1,j) =  grad1_tmp(idimn,1,j) + (fi1(idimn,i)- fi1(idimn,j))*h1                                        						              	
                        if (ndimn == 2) then                                        						
                                h2 = current%dwdy*mass(i)/rho(i)
                                grad1_tmp(idimn,2,j) =  grad1_tmp(idimn,2,j) + (fi1(idimn,i)- fi1(idimn,j))*h2						                
                        end if                                        
                end do
        
                !calculate stress gradient on nodes
                !---------------------------------------------------------------                        
                do istre = 1,nstre-1                                                                           
                        g1 = current%dwdx*(fi2(istre,i)/rho(i)**2 + fi2(istre,j)/rho(j)**2)                                        
                        grad2_tmp(istre,1,i) =  grad2_tmp(istre,1,i) - mass(j)*g1                                                                                 
                        if (ndimn == 2) then                                                                                
                                g2 = current%dwdy*(fi2(istre,i)/rho(i)**2 + fi2(istre,j)/rho(j)**2)                                                        
                                grad2_tmp(istre,2,i) =  grad2_tmp(istre,2,i) - mass(j)*g2                                                        
                        end if                                
                end do  
                
                !CSPM gradient normalisation (according to Chen et al. 2000)
                if (CSPM) then
                        h1b = -current%dwdx*mass(j)/rho(j)
                        h2b = -current%dwdy*mass(j)/rho(j)
                        AE(1,j) = AE(1,j) + (x(1,i) - x(1,j))*h1
                        AE(2,j) = AE(2,j) + (x(2,i) - x(2,j))*h1
                        AE(3,j) = AE(3,j) + (x(1,i) - x(1,j))*h2
                        AE(4,j) = AE(4,j) + (x(2,i) - x(2,j))*h2                         
                        AE(1,i) = AE(1,i) + (x(1,j) - x(1,i))*h1b
                        AE(2,i) = AE(2,i) + (x(2,j) - x(2,i))*h1b
                        AE(3,i) = AE(3,i) + (x(1,j) - x(1,i))*h2b
                        AE(4,i) = AE(4,i) + (x(2,j) - x(2,i))*h2b
                end if
                                                                               				
        else if (current%pint_type == 9) then ! stress-point--dummy node                        
                ! i  is stress point
                ! j is dummy node
                beta_max = 1.5 !dummy node parameters
                vel_wall = 0.0  
                wall = wall_position(j)                              
                if (horizontal_or_not(j) == 1) then !horizontal wall
                        da = abs(x(2,i)-wall)
                        db = abs(x(2,j)-wall)
                else                                !vertical wall (left or right)
                        da = abs(x(1,i)-wall)
                        db = abs(x(1,j)-wall)
                end if  
                                                                                             
                beta = min(beta_max, 1+(db/da))                                
                dummy_vel(:,j) = vel(:,i)*(1-beta) + beta*vel_wall
                h1 = current%dwdx*mass(j)/rho(j)                                             
                do idimn = 1,ndimn
                        grad1_tmp(idimn,1,i) =  grad1_tmp(idimn,1,i) + (fi1(idimn,i)- dummy_vel(idimn,j))*h1		                
                        if (ndimn == 2) then
                                h2 = current%dwdy*mass(j)/rho(j)
                                grad1_tmp(idimn,2,i) =  grad1_tmp(idimn,2,i) + (fi1(idimn,i)- dummy_vel(idimn,j)) * h2                                                
                        end if  
                end do 
                                
        else if (current%pint_type == 6) then !node--dummy node
                        
                        bc_int(i) = 1 !stores node--dummy node interactions, used when shifting stress-points near the boundary in the outside approach
                        stress(:,j) = stress(:,i)                                                      
                        do istre = 1,nstre-1                         
                                h1 = current% dwdx*(fi2(istre,i)/rho(i)**2 + stress(istre,j)/rho(j)**2)                                        
                                grad2_tmp(istre,1,i) =  grad2_tmp(istre,1,i) - mass(j)*h1
                                if (ndimn == 2) then 
                                        h2 = current% dwdy*(fi2(istre,i)/rho(i)**2 + stress(istre,j)/rho(j)**2)                                                        
                                        grad2_tmp(istre,2,i) =  grad2_tmp(istre,2,i) - mass(j)*h2
                                end if                                                                                               
                        end do                                                                                                    		                       
       end if                                                                                     
       current=>current%next         
end do

!matrix for CSPM normalisation
!-------------------------------------------------------------------------------

if (CSPM) then
        do i = 1,ntotal
                AE(5,i) = AE(1,i)*AE(4,i) - AE(2,i)*AE(3,i)
                if (abs(AE(5,i)) < 1e-07) then
                        AE(5,i) = 1; AE(1,i) = 1; AE(2,i) = 0; AE(3,i) = 0; AE(4,i) = 1 !reduces to no CSPM if AE(5,i) is too small
                else
                        AE(5,i) = 1/AE(5,i)
                end if
        end do
end if	

!get divergence terms for the momentum and constitutive equations
!-------------------------------------------------------------------------------				        		                        		                        
		                                		                        			         			          
if (ndimn == 1) then 
                                                     
        div1(1,nnode+1:ntotal) = -young*grad1_tmp(1,1,nnode+1:ntotal)                   
        div2(1,1:nnode) = -grad2_tmp(1,1,1:nnode) 
                                                     
else if (ndimn == 2) then

        if (CSPM) then
                do i = 1,ndimn
                        grad1_tmp(i,1,nnode+1:ntotal) = (AE(5,nnode+1:ntotal))*(AE(1,nnode+1:ntotal)*grad1_tmp&
                                        (i,1,nnode+1:ntotal)+ AE(2,nnode+1:ntotal)*grad1_tmp(i,2,nnode+1:ntotal))
                        grad1_tmp(i,2,nnode+1:ntotal) = (AE(5,nnode+1:ntotal))*(AE(3,nnode+1:ntotal)*grad1_tmp&
                                        (i,1,nnode+1:ntotal)+ AE(4,nnode+1:ntotal)*grad1_tmp(i,2,nnode+1:ntotal)) 
                        grad2_tmp(i,1,1:nnode) = (AE(5,1:nnode))*(AE(1,1:nnode)*grad2_tmp(i,1,1:nnode) &
                                                + AE(2,1:nnode)*grad2_tmp(i,2,1:nnode))
                        grad2_tmp(i,2,1:nnode) = (AE(5,1:nnode))*(AE(3,1:nnode)*grad2_tmp(i,1,1:nnode) &
                                                + AE(4,1:nnode)*grad2_tmp(i,2,1:nnode))
                end do
        end if
        
        !calculate divergence of f_u on stress-points
        !-----------------------------------------------------------------------
                  
        div1(1,nnode+1:ntotal) = -(D11*grad1_tmp(1,1,nnode+1:ntotal) + D12*grad1_tmp(2,2,nnode+1:ntotal))
        div1(2,nnode+1:ntotal) = -(D12*grad1_tmp(1,1,nnode+1:ntotal) + D22*grad1_tmp(2,2,nnode+1:ntotal))
        div1(3,nnode+1:ntotal) = -(D33*grad1_tmp(2,1,nnode+1:ntotal) + D33*grad1_tmp(1,2,nnode+1:ntotal))
        div1(4,nnode+1:ntotal) = -(D41*grad1_tmp(1,1,nnode+1:ntotal) + D42*grad1_tmp(2,2,nnode+1:ntotal))
        
        !calculate divergence of f_sigma on nodes
        !-----------------------------------------------------------------------
                                            
        div2(1,1:nnode) = -(grad2_tmp(1,1,1:nnode) + grad2_tmp(3,2,1:nnode))
        div2(2,1:nnode) = -(grad2_tmp(3,1,1:nnode) + grad2_tmp(2,2,1:nnode))
        
        grad_u = grad1_tmp !store velocity gradient for use in other subroutines
        
end if

        END SUBROUTINE get_derivatives

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

        SUBROUTINE RK4

!RK4 time integration loop, to update problem variables in time
!velocity is updated on the SPH nodes, stress is updated on the SPH stress-points
!(in Standard SPH nodes = stress-points)  
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


implicit none

integer :: i, j, idimn, istre, ncrit
real (kind=DP) :: f1rk(4), f2rk(4), divf1(nstre,ntotal), divf2(ndimn,ntotal)
real (kind=DP) :: RK_dev_strain(ntotal), RK_vel(ndimn,ntotal), RK_stress(nstre,ntotal), RK_rho(ntotal), RK_hsml(ntotal)
real (kind=DP) :: art_force(ndimn,nnode), RHS_x(ndimn,ntotal), G_local(nstre)
real (kind=DP) :: der_intvars_local(nint_Vars), RHS_1(nstre,ntotal), RHS_2(ndimn,ntotal), spin(nstre,ntotal)
real (kind=DP) :: variab(nstre), vel0(ndimn,ntotal), stress0(nstre,ntotal), source_stress(nstre,ntotal), source_grav(ndimn,nnode) 
real (kind=DP) :: omega(ndimn,ntotal), G_grav(ndimn,nnode)
real (kind=DP) :: rho0(nstress), RHS_rho(ntotal), hsml0(nstress), RHS_hsml(ntotal)


!RK4 integrator coefficients
!-------------------------------------------------------------------------------

data f1rk/0., 0.5, 0.5, 1.0/				
data f2rk/1., 2., 2., 1.0/			


vel0 = vel(:,1:ntotal); stress0 = stress(:,1:ntotal); rho0 = rho(nnode+1:ntotal); hsml0 = hsml(nnode+1:ntotal) !initial conditions for the start of each RK4 routine

RK_stress=0.0; RK_vel = 0.0; RK_rho = 0; RK_hsml = 0; RK_dev_strain = 0 !initialise variable RK4 terms

RHS_1 = 0.0; RHS_2 = 0.0; RHS_rho = 0; RHS_hsml = 0 !initialise RHS RK4 terms

source_stress = 0.0; source_grav = 0.0; divf1 = 0.0; divf2 = 0.0 !initialise terms within RHS 
art_visc = 0.0; art_force = 0.0; G_grav = 0.0; G_local = 0.0; spin = 0.0; omega = 0.0

vel = 0.0; stress = 0.0 !initialise variables

ncrit = props(1,2) !constitutive model

!main RK4 loop
!-------------------------------------------------------------------------------
do i = 1,4       
        !update velocity and stress within each RK4 integration step
        !-----------------------------------------------------------------------
                        
	vel(:,1:nnode) = vel0(:,1:nnode) + f1rk(i)*(dt_sph)*RHS_2(:,1:nnode)
	stress(:,nnode+1:ntotal) = stress0(:,nnode+1:ntotal) + f1rk(i)*(dt_sph)*RHS_1(:,nnode+1:ntotal)
		
	!updates density and smoothing length if turned on
	!-----------------------------------------------------------------------
		
	if (cont_density) then		
	        call density_update(RHS_rho)
	        rho(nnode+1:ntotal) = rho0 + f1rk(i)*(dt_sph)*RHS_rho(nnode+1:ntotal)
	        if (sle == 2) then
	                RHS_hsml(nnode+1:ntotal) = -(hsml0(1:nstress)/(rho(nnode+1:ntotal)*ndimn))*RHS_rho(nnode+1:ntotal)
	                hsml(nnode+1:ntotal) = hsml0(1:nstress) + f1rk(i)*(dt_sph)*RHS_hsml(nnode+1:ntotal)
	        end if		        
	end if	

        if (ncrit == 12) call adapt_stress2 !adapt the stress state for elastoplastic constitutive law (as described in Bui et al. 2008)
        if (no_bcs > 0) call BCs !update fixed boundary conditions
                
        call stress_point_update !transfer information from nodes to stress-points (and vice versa)
      
        if (SP_SPH) then !only necessary to adapt stress and bcs if stress-particle SPH is used
                if (ncrit == 12) call adapt_stress2 !adapt stress in elastoplastic model, if stress-particle SPH is used
                if (no_bcs > 0) call BCs !update fixed boundary conditions
        end if
                     				
	call get_derivatives(vel,divf1,stress,divf2) !div(vel) at stress-points

        !get plastic terms for constitutive equation
        !----------------------------------------------------------------------- 
                
        do j = nnode+1,ntotal
                call plastic_terms(j,ndimn,stress(:,j),G_local,der_intvars_local)
                source_stress(:,j) = G_local(:)
                RK_dev_strain(j) = RK_dev_strain(j) + der_intvars_local(1)*f2rk(i)
        end do
                
        !get gravity force for equation of momentum
        !-----------------------------------------------------------------------
                
        call gravity_force(ndimn,vel(:,1:nnode),G_grav)
        source_grav = G_grav
                 
        if (inside_approach) call boundary_forces !boundary repulsive force - only required for inside approach
                
        if (alpha > 0 .or. beta > 0) call artificial_viscosity(vel) !artificial viscosity (if used)
                                                
        if (art_stress) call artificial_force(vel,art_force) !artificial stress (treatment for the tensile instability)
                
        !gets spin rate tensor to replace stress with Jaumann stress for large deformation problems
        !-----------------------------------------------------------------------
                                
        if (update_x .and. ndimn == 2) then                  
                call get_spin_rate_tensor(omega)                
                spin(1,nnode+1:ntotal) = 2*omega(1,nnode+1:ntotal)*stress(3,nnode+1:ntotal)
                spin(2,nnode+1:ntotal) = 2*omega(2,nnode+1:ntotal)*stress(3,nnode+1:ntotal)
                spin(3,nnode+1:ntotal) = omega(2,nnode+1:ntotal)*stress(1,nnode+1:ntotal) &
                                        + omega(1,nnode+1:ntotal)*stress(2,nnode+1:ntotal)                  
        end if  
                
        !right hand side terms of equations for RK4 integration
        !-----------------------------------------------------------------------             
                
        RHS_1(:,nnode+1:ntotal) = -divf1(:,nnode+1:ntotal) + spin(:,nnode+1:ntotal) + source_stress(:,nnode+1:ntotal)
        RHS_2(:,1:nnode) = -divf2(:,1:nnode) + source_grav(:,1:nnode) &
                          + art_visc(:,1:nnode) + f_bound(:,1:nnode) + art_force(:,:)

        !RK4 terms required to update variables after main RK4 loop
        !---------------------------------------------------------------
                
	RK_vel(:,1:nnode) = RK_vel(:,1:nnode) + f2rk(i)*RHS_2(:,1:nnode)
	RK_stress(:,nnode+1:ntotal) = RK_stress(:,nnode+1:ntotal) + f2rk(i)*RHS_1(:,nnode+1:ntotal)
		
	if (cont_density) then
		RK_rho = RK_rho + f2rk(i)*RHS_rho		
		if (sle == 2) then
		        RK_hsml = RK_hsml + f2rk(i)*RHS_hsml
                end if
	end if
		   
end do

!update stress and velocity at time t = n+1
!-------------------------------------------------------------------------------

vel(:,1:nnode) = vel0(:,1:nnode) + (dt_sph/6)*RK_vel(:,1:nnode)
stress(:,nnode+1:ntotal) = stress0(:,nnode+1:ntotal) + (dt_sph/6)*RK_stress(:,nnode+1:ntotal)
	
if (ncrit == 12) call adapt_stress2 !adapts stress state for elastoplastic material
if (no_bcs > 0) call BCs !update fixed boundary conditions
!updates density and smoothing length if turned on
!-------------------------------------------------------------------------------
	
if (cont_density) then
	rho(nnode+1:ntotal) = rho0 + (dt_sph/6.)*RK_rho(nnode+1:ntotal)	
	if (sle == 2) then
	        hsml(nnode+1:ntotal) = hsml0 + (dt_sph/6.)*RK_hsml(nnode+1:ntotal)
	end if
end if
	
Ddev_strn(:) = RK_dev_strain(:)/6 !to obtain deviatoric plastic strain


        END SUBROUTINE RK4  

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

        SUBROUTINE density_update(rho0)

!calculates right hand side of continuity equation if updating the density (in the continuity approach)
!calculated on stress-points as velocity gradient is already stored there
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

implicit none

real (kind=DP), intent(out) :: rho0(:)

rho0(nnode+1:ntotal) = -rho(nnode+1:ntotal)*(grad_u(1,1,nnode+1:ntotal) + grad_u(2,2,nnode+1:ntotal))


        END SUBROUTINE

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

        SUBROUTINE artificial_viscosity(vel0)

!calculates artificial viscosity (via node-node interactions)
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

implicit none

integer :: i, j, iamaTG, idimn, ipoin, switch, i0, j0, i2, k

real ::  hf1, dwdx, w, gradw, sigma, h1, h2, wall
real ::  h, theta, div_u, c, rho2, young, visc, xij, yij
real :: div_u_out_temp(nnode), art_visc_temp(ndimn,nnode)
real :: art_visc_temp_bound(ndimn,nnode)
real(kind=DP) :: vel0(ndimn,nnode), dummy_vel(ndimn,ntotal2)
real :: da, db, beta2, y0, beta_max

n_int = 0.0

art_visc_temp= 0.0; visc = 0.0; dummy_vel = 0.0; art_visc_temp_bound = 0.0
xij = 0.0; div_u = 0.0
		
young = props(1,3) !E = young's modulus
			             
 ! approximate div(art_viscosity) on the nodes, using information from the stress-points
 
 current => last
do while (associated(current)) 
        i = current%pair_i !i = node
        j = current%pair_j !j = stress point
        xij = x(1,i) - x(1,j)                                                
        yij = 0.0 
                                                                                       
        if (ndimn == 2) then                                                
                yij = x(2,i) - x(2,j)                                                        
        end if 
                                                                                                                                                           
        h = 0.5*(hsml(i) + hsml(j))                                                
        rho2 = 0.5*(rho(i) + rho(j))                                                
        !c = 0.5*(sqrt(props(1,3)/rho(i)) + sqrt(props(1,3)/rho(j))) 
        c = 600                                       
                                                                                                                                                                        
        if (current%pint_type == 3) then !node - node
        
                n_int(i) = n_int(i) + 1 !calculates number of node-node interactions for each node (for use in the outside approach, useful to have anyway)
                n_int(j) = n_int(j) + 1
        
                div_u = xij*(vel0(1,i) - vel0(1,j))
                                        
                if (ndimn == 2) then
                        div_u = div_u + yij*(vel0(2,i) - vel0(2,j))
                end if  
                
                theta = (h*div_u)/((sqrt(xij**2 +yij**2))**2 + 0.01*h**2)
                                        
                if (div_u < 0) then                                                                               
                        visc = (-alpha*c*theta + beta*theta**2)/rho2                                                
                else if (div_u >= 0) then                                                
                        visc = 0                                                     
                end if  
                                                
                art_visc_temp(1,i) = art_visc_temp(1,i) + visc*current%dwdx*mass(j)
                art_visc_temp(1,j) = art_visc_temp(1,j) - visc*current%dwdx*mass(i)
                                                
                if (ndimn == 2) then
                        art_visc_temp(2,i) = art_visc_temp(2,i) + visc*current%dwdy*mass(j)
                        art_visc_temp(2,j) = art_visc_temp(2,j) - visc*current%dwdy*mass(i)
                end if                
                
        end if
		                       
        current => current%next
		                       
end do
		                       
art_visc = -art_visc_temp                
        

        END SUBROUTINE artificial_viscosity

!---------------------------------------------------------------------

        SUBROUTINE artificial_force(vel0,art_force)

!---------------------------------------------------------------------
!applies artificial force method (Monaghan 2000), that repels particles that are too close together

real (kind =DP) :: vel0(:,:), art_force(:,:)
real (kind =DP) :: s12i, s12j, theta, theta2, eps, n, f, w2, gradw2(ndimn), hsml0, h1, h2, dx2(ndimn)
real (kind =DP) :: sigma2(ndimn,nnode), R2(ndimn,nnode), R(nstre-1,nnode), art_force_temp(nstre-1,ndimn,nnode)
integer :: i, j, i2, istre

sigma2 = 0; R = 0; R2 = 0
art_force_temp = 0
eps = 0.1
n = 2.55

w2 = 0.0
gradw2 = 0.0
hsml0 = 1.2*dx

dx2(1) = dx
dx2(2) = dy

 call kernel(dx,dx2,hsml0,w2,gradw2)

 current => last
 
 do while (associated(current))
 
        if (current%pint_type == 3) then !node-node
        
                i = current%pair_i
                j = current%pair_j
                
                s12i = (stress(1,i) - stress(2,i))
                s12j = (stress(1,j) - stress(2,j))
                
                if (s12i >= 1e-08) then
                        theta = 0.5*atan(2*stress(3,i)/s12i)
                else
                        theta = 0.0
                end if
                
                if (s12j >= 1e-08) then 
                        theta2 = 0.5*atan(2*stress(3,j)/s12j)
                else
                        theta2 = 0.0
                end if
                
                sigma2(1,i) = (cos(theta)**2)*stress(1,i) + 2*cos(theta)*sin(theta)*stress(3,i) + (sin(theta)**2)*stress(2,i)
                
                sigma2(1,j) = (cos(theta2)**2)*stress(1,j) + 2*cos(theta2)*sin(theta2)*stress(3,j) + (sin(theta2)**2)*stress(2,j)
                
                sigma2(2,i) = (sin(theta)**2)*stress(1,i) - 2*cos(theta)*sin(theta)*stress(3,i) + (cos(theta)**2)*stress(2,i)
                
                sigma2(2,j) = (sin(theta2)**2)*stress(1,j) - 2*cos(theta2)*sin(theta2)*stress(3,j) + (cos(theta2)**2)*stress(2,j)
                
                f = (current%w)/w2
                
                do i2 = 1,2
                                if (sigma2(i2,i) > 0) then
                                        R2(i2,i) = -eps*(sigma2(i2,i)/(rho(i)**2))
                                else
                                        R2(i2,i) = 0
                                end if
                                if (sigma2(i2,j) > 0) then
                                        R2(i2,j) = -eps*(sigma2(i2,j)/(rho(j)**2))
                                else
                                        R2(i2,j) = 0
                                end if
                end do
                
                R(1,i) = R2(1,i)*(cos(theta)**2) + R2(2,i)*(sin(theta)**2)
                R(2,i) = R2(1,i)*((cos(theta)**2) + (sin(theta)**2))
                R(3,i) = (R2(1,i)-R2(2,i))*(cos(theta)*sin(theta))
                        
                R(1,j) = R2(1,j)*(cos(theta2)**2) + R2(2,j)*(sin(theta2)**2)
                R(2,j) = R2(1,j)*((cos(theta2)**2) + (sin(theta2)**2))
                R(3,j) = (R2(1,j)-R2(2,j))*(cos(theta2)*sin(theta2))
                
                do istre = 1,nstre-1
                
                        h1 = current% dwdx*(f**n)*(R(istre,i) + R(istre,j))                        
                        art_force_temp(istre,1,i) = art_force_temp(istre,1,i) + mass(j)*h1 
                        art_force_temp(istre,1,j) = art_force_temp(istre,1,j) - mass(i)*h1
                        
                        if (ndimn == 2) then
                                
                                h2 = current% dwdy*(f**n)*(R(istre,i) + R(istre,j))
                                art_force_temp(istre,2,i) = art_force_temp(istre,2,i) + mass(j)*h2 
                                art_force_temp(istre,2,j) = art_force_temp(istre,2,j) - mass(i)*h2 
                            
                        end if
                        
                end do
                
          end if
          
          current => current%next
          
  end do
  
  do i = 1,nnode
        
        art_force(1,i) = art_force_temp(1,1,i) + art_force_temp(3,2,i)
        art_force(2,i) = art_force_temp(3,1,i) + art_force_temp(2,2,i)
        
  end do
                
        END SUBROUTINE

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

        SUBROUTINE get_spin_rate_tensor(omega2)
        
!calculates the spin rate tensor for large deformation problems
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

implicit none

real (kind=DP) :: omega2(:,:)

omega2(1,:) = 0.5*(grad_u(1,2,:) - grad_u(2,1,:))
omega2(2,:) = -0.5*(grad_u(1,2,:) - grad_u(2,1,:))

        END SUBROUTINE
        
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

        SUBROUTINE boundary_forces

!calculates a repulsive boundary force to prevent particles penetrating the boundary
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

implicit none

real :: beta, k, w, r(ndimn), r2, r0, r_bound, d(ndimn), d0, c, f, f2, h
integer :: i, j, idimn, p1, p2, test, start, finish


r2 = 0.0

f_bound = 0.0
r = 0.0
d0 = 0.0
f = 0; f2 = 0

 !c = 0.5*(sqrt(props(1,3)/rho(i)) + sqrt(props(1,3)/rho(j))) 
 c = 20 !Speed of sound, 600 is max in a soil

test = 2

if (test == 1) then

 current => last
 do while (associated(current))
        if (current%pint_type == 5) then !node - boundary node
                i = current%pair_i !node 
                j = current%pair_j !boundary node
                do idimn = 1,ndimn
                        r(idimn) = x(idimn,i) - x(idimn,j)
                        d(idimn) = x00(idimn,i) - x00(idimn,j)                                               
                end do
                
                
                r2 = sqrt(r(1)**2 + r(2)**2) !current distance between particles
                d0 = sqrt(d(1)**2 + d(2)**2) !initial distance
                f2 = 1 - (r2/d0)
                h = 0.5*(hsml(i) + hsml(j))
                r_bound = r2/(0.75*h)
                
                if (r2 < d0) then
                
                        if (0 < r_bound .and. r_bound <= 2./3.) then
                
                                f = 2./3.                                        
                        
                        else if (2./3. < r_bound .and. r_bound <= 1) then
                        
                                f = 2*r_bound - 1.5*r_bound**2
                        
                        else if (1 < r_bound .and. r_bound < 2) then
                
                                f = 0.5*(2-r_bound)**2
                        
                        else 
                
                                f = 0.
                                
                        end if
                        
                        f_bound(:,i) = f_bound(:,i) + (0.01*c**2)*f2*f*(r(:)/(r2**2))
                        
                else
                
                
                        f_bound(:,i) = 0.0
                        
                                                
                        
               end if
                
        end if
        current => current%next
 end do 
 
 else if (test == 2) then
 
 start = ntotal + 1
 finish = ntotal + ndummy2
 do j = start,finish !boundary
        do i = 1,nnode !node
                do idimn = 1,ndimn
                        r(idimn) = x(idimn,i) - x(idimn,j)
                        d(idimn) = x00(idimn,i) - x00(idimn,j)                                               
                end do                                
                
                r2 = sqrt(r(1)**2 + r(2)**2) !current distance between particles
                !d0 = sqrt(d(1)**2 + d(2)**2) !initial distance
                d0 = dx/2.
                f2 = 1 - (r2/(1.5*d0))
                h = 0.5*(hsml(i) + hsml(j))
                r_bound = r2/(0.75*h)
                
                if (r2 > 0 .and. r2 < 1.5*d0) then                
                        f2 = 1 - (r2/(1.5*d0))                        
                else                
                        f2 = 0.0                        
                end if
                
                        if (0 < r_bound .and. r_bound <= 2./3.) then
                
                                f = 2./3.                                        
                        
                        else if (2./3. < r_bound .and. r_bound <= 1) then
                        
                                f = 2*r_bound - 1.5*r_bound**2
                        
                        else if (1 < r_bound .and. r_bound < 2) then
                
                                f = 0.5*(2-r_bound)**2
                        
                        else 
                
                                f = 0.
                                
                        end if
                        
                 f_bound(:,i) = f_bound(:,i) + (0.01*c**2)*f2*f*(r(:)/(r2**2))
        end do
 end do
 
 end if
                
 end subroutine                               

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
 
       SUBROUTINE Check_Out_Domain  

!checks if each particle is outside the computational domain
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

implicit none

integer :: ipoin, idimn
real (kind=DP) ::  dxx 

!      ------  Initialization 

! If_Out_Domain = 0               CHANGED 19 July 2007

!      ------  

DO ipoin = 1,ntotal2
   Do idimn = 1, ndimn
      dxx   = (x(idimn,ipoin)-Xmin_Domain(idimn))*(x(idimn,ipoin)-Xmax_Domain(idimn))
      if ( dxx.gt.0.0 ) If_Out_Domain(ipoin) = 1
   Enddo
Enddo

        END SUBROUTINE Check_Out_Domain

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
 
       SUBROUTINE grid_find_NEW  (itimestep_sph)  

! Routine which finds the particles that are interacting with one another
! The particles are divided into a smaller number of grids, and each grid is searched for interacting particles 
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
 
implicit none

type (pairs), SAVE, pointer::keepit0, keepit1, current 
 
integer, SAVE:: mdivt, ndivt, m_pairs, n_pairs 
integer, SAVE, allocatable:: which_cell(:), info_picell(:,:), list_picell(:)

real (kind=DP), SAVE, allocatable:: xmin(:), xmax(:), deltx(:)
real (kind=DP), SAVE, allocatable:: dxiac(:), tdwdx(:)
real (kind=DP) ::  length, length_new 
real (kind=DP) ::    driac, r, mhsml, w 
integer :: ndivx(3), idivx(3), jdivx(3), itimestep_sph, idimn, itotal, jtotal 
integer :: idivt, jdivt, idest, ipost, i0, idelt, i, j, k, scale_k, countiac_real(ntotal2) 
integer :: idx, idy, idz, jdx, jdy, jdz, idx0, idy0, idz0, idx1, idy1, idz1
integer :: sumiac, maxiac, miniac, noiac, maxp, minp, npicsi, npicsj, jj, jdest, j0 
 
countiac_real = 0
!      ------  The first time we have not allocated anything. Do it

if (itimestep_sph.eq.1) then
   allocate ( xmin(ndimn),  xmax(ndimn),  deltx(ndimn) )
   allocate ( dxiac(ndimn), tdwdx(ndimn) )
endif

xmin = 1.e+10; xmax = -xmin; deltx = 0.0; ndivx = 1; idivx = 1; jdivx = 1
 
if (skf.eq.1) then 
   scale_k = 2 
else if (skf.eq.2) then 
   scale_k = 3 
else if (skf.eq.3) then 
   scale_k = 3 
endif 

niac =0; countiac = 0

!      ------  Task 1: Setup Grid parameters

  
Do idimn = 1,ndimn
   Do itotal = 1, ntotal2
      if ( If_Out_Domain(itotal).eq.1 ) CYCLE
      if (x(idimn,itotal).lt.xmin(idimn)) xmin(idimn)  = x(idimn,itotal)
      if (x(idimn,itotal).gt.xmax(idimn)) xmax(idimn)  = x(idimn,itotal)
      if (deltx(idimn).lt.hsml(itotal))   deltx(idimn) = hsml(itotal)
   enddo
   deltx(idimn) = deltx(idimn)*2       
   length       = xmax(idimn) - xmin(idimn)
   ndivx(idimn) = (length/deltx(idimn)) + 1
   length_new   = ndivx(idimn)*deltx(idimn)
   xmin(idimn)  = xmin(idimn) - (length_new-length)/2 - 0.001*length
   xmax(idimn)  = xmax(idimn) + (length_new-length)/2 + 0.001*length
Enddo

ndivt = ndivx(1)*ndivx(2)*ndivx(3) 

!      ------  Now we can allocate or reallocate some control arrays

if (itimestep_sph.eq.1) then
    mdivt = ndivt
    allocate ( which_cell(ntotal2), info_picell(3,mdivt),  list_picell(ntotal2) )
else if (ndivt.gt.mdivt) then
    deallocate(info_picell)
    mdivt = max (2*mdivt, ndivt)
    allocate  (info_picell(3,mdivt) )
endif

which_cell = 0; info_picell=0; list_picell=0

!      ------  Task 2: Fill Grid structure with particles

DO itotal = 1, ntotal2
   if ( If_Out_Domain(itotal).eq.1 ) CYCLE
   DO idimn = 1, ndimn
      idivx(idimn) = (x(idimn,itotal)-xmin(idimn))/deltx(idimn) + 1
      if (idivx(idimn).gt.ndivx(idimn)) idivx(idimn)=ndivx(idimn)
   enddo
   idivt =  ndivx(1)*ndivx(2)*(idivx(3)-1) + ndivx(1)*(idivx(2)-1) + idivx(1)
   which_cell(itotal)   = idivt
   info_picell(1,idivt) = info_picell(1,idivt) +1
Enddo

!      ------  in 2nd row we store, for each cell, the position in array list_picell(ntotal)
!              where we start storing particles in this cell

ipost = 1
do idivt = 1,ndivt
   info_picell(2,idivt) = ipost
   ipost = ipost + info_picell(1,idivt)
enddo

info_picell(3,:) = -1

do itotal = 1, ntotal2
   if ( If_Out_Domain(itotal).eq.1 ) CYCLE
   idivt  = which_cell  (itotal)
   info_picell(3,idivt) = info_picell(3,idivt) + 1 
   idest  = info_picell(2,idivt) + info_picell(3,idivt)
   list_picell(idest) = itotal
enddo

!      ------  Task 3: Get list of interactions

if (itimestep_sph.eq.1) then
   m_pairs = 0
   n_pairs = 0
   nullify (last)
   allocate (keepit0, keepit1)
else
   if (n_pairs.lt.m_pairs) then
       keepit0%next => keepit1
       current => last
   endif
endif   

n_pairs = 0
DO idz = 1, ndivx(3)
DO idy = 1, ndivx(2)
DO idx = 1, ndivx(1)
   idivt  =  ndivx(1)*ndivx(2)*(idz-1) + ndivx(1)*(idy-1) + idx 
   npicsi = info_picell(1,idivt)
   DO i = 1, npicsi
      i0     = info_picell(2,idivt)
      idest  = i0 + i -1 
      itotal = list_picell(idest)
      idx0   = max (       1,idx-1); idy0   = max (       1,idy-1); idz0   = max (       1,idz-1)
      idx1   = min (ndivx(1),idx+1); idy1   = min (ndivx(2),idy+1); idz1   = min (ndivx(3),idz+1)
      DO jdz = idz, idz1        !  before it was idz0, idz1
      DO jdy = idy, idy1        !  before it was idy0, idy1
      DO jdx = idx0,idx1 
         jdivt = ndivx(1)*ndivx(2)*(jdz-1) + ndivx(1)*(jdy-1) + jdx
         if (jdivt.lt.idivt) CYCLE
         jj = 1
         if (idivt.eq.jdivt) jj = i + 1
         npicsj = info_picell(1,jdivt)
         DO j = jj,npicsj
            j0     = info_picell(2,jdivt)
            jdest  = j0 + j -1
            jtotal = list_picell(jdest)
            dxiac(1) = x(1,itotal)-x(1,jtotal)          !          note jdivt is always greater than idivt
            driac    = dxiac(1)*dxiac(1)
        	do idimn = 2, ndimn
               dxiac(idimn) = x(idimn,itotal)-x(idimn,jtotal)
               driac = driac + dxiac(idimn)*dxiac(idimn)
            enddo
            mhsml =  (hsml(itotal)+hsml(jtotal))/2.
            r = sqrt(driac)
            if (r.lt.scale_k*mhsml) then
               n_pairs = n_pairs + 1
               countiac(itotal) = countiac(itotal) + 1
               countiac(jtotal) = countiac(jtotal) + 1
               If (itype(itotal).ne.itype(jtotal)) then
			   countiac_real(itotal) = countiac_real(itotal) +1
			   countiac_real(jtotal) = countiac_real(jtotal) +1 
		Endif
		if (n_pairs > m_pairs) then
                   allocate(current)                    !      create next pair structure
                   !allocate (current%dwdx(ndimn))       !      allocate dwdx
                   !current% dwdx = 0.0
                   !current% dwdy = 0.0
                   current%next => last                 !      
                   last         => current              !      
                   m_pairs      =  m_pairs + 1          ! 
               end if                                    !    
               
               call kernel(r, dxiac, mhsml, w, tdwdx)  !      using r, dxiac and mshml obtains w and its gradient
               current%niac      = n_pairs             !      this is the pair number 
               current%pair_i    = itotal              !      first particle
               current%pair_j    = jtotal              !      and second particle
               current%Pint_type = 0                   !      Sets interaction type to default (same material)
               current%w         = w                   !      Weight
               current%dwdx   = tdwdx(1)            !      Gradient of weight
               current%dwdy   = tdwdx(2) 
               
               if (n_pairs.lt.m_pairs) then            !
                   keepit0      => current             !               
                   current      => current%next        ! 
               endif
               
            endif
         ENDDO
         
      ENDDO
      ENDDO
      ENDDO
      
   ENDDO
   
ENDDO
ENDDO
ENDDO

If (n_pairs.lt.m_pairs) then
   keepit1 => current
   nullify (keepit0%next)
endif

current=> last
 
!      ------  Statistics for the interaction

sumiac = 0
maxiac = 0
miniac = 1000
noiac  = 0
DO i = 1, ntotal2
   sumiac = sumiac + countiac(i)  
   if (countiac(i).gt.maxiac) then
      maxiac = countiac(i)
      maxp = i
   endif
   if (countiac(i).lt.miniac) then 
      miniac = countiac(i)
      minp = i
   endif
   if (countiac(i).eq.0)      noiac  = noiac + 1
enddo
 
if (mod(itimestep_sph,print_step).eq.0) then
       write(output_file,*) ' >> Statistics: interactions per particle:', time, dt_sph
       write(output_file,*)'**** Particle:',maxp, ' maximal interactions:',maxiac
       write(output_file,*)'**** Particle:',minp, ' minimal interactions:',miniac
       write(output_file,*)'**** Average :',real(sumiac)/real(ntotal2)
       write(output_file,*)'**** Total pairs : ',niac
       write(output_file,*)'**** Particles with no interactions:',noiac
endif    



        END SUBROUTINE  grid_find_NEW 

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

        SUBROUTINE kernel(r,dx2,hsml,w,dwdx)   

!calculates SPH kernel and kernel gradient values
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!   Subroutine to calculate the smoothing kernel wij and its 
!              derivatives dwdxij.
!   if skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
!          = 2, Gauss kernel   (Gingold and Monaghan 1981) 
!          = 3, Quintic kernel (Morris 1997)

!   r    : Distance between particles i and j                     [in]
!   dx   : x-, y- and z-distance between i and j                  [in]  
!   hsml : Smoothing length                                       [in]
!   w    : Kernel for all interaction pairs                      [out]
!   dwdx : Derivative of kernel with respect to x, y and z       [out]

implicit none
      
real (kind=DP)::   r, dx2(ndimn), hsml, w, dwdx(ndimn)    
real  (kind=DP)::  q, factor
integer :: idimn  

q = r/hsml 
w       = 0. 
dwdx    = 0.

if (skf.eq.1) then 
    if (ndimn.eq.1) then
       factor = 1.e0/hsml
    elseif (ndimn.eq.2) then
       factor = 15.e0/(7.e0*pi*hsml*hsml)
    elseif (ndimn.eq.3) then
       factor = 3.e0/(2.e0*pi*hsml*hsml*hsml)
    else
       print *,' >>> Error <<< : Wrong dimension: Dim =',ndimn
       stop
    endif                                           
    if (q.ge.0.and.q.le.1.e0) then          
       w = factor * (2./3. - q*q + q**3 / 2.)
       do idimn = 1, ndimn
          dwdx(idimn) = factor * (-2.+3./2.*q)/hsml**2 * dx2(idimn)       
       enddo   
    else if (q.gt.1.e0.and.q.le.2) then          
       w = factor * 1.e0/6.e0 * (2.-q)**3 
       do idimn = 1, ndimn
          dwdx(idimn) =-factor * 1.e0/6.e0 * 3.*(2.-q)**2/hsml * (dx2(idimn)/r)        
       enddo              
    else
       w    = 0.0
       dwdx = 0.0           
    endif  
    
else if (skf.eq.2) then      
    factor = 1.e0 / (hsml**ndimn * pi**(ndimn/2.))      
    if(q.ge.0.and.q.le.3) then
       w = factor * exp(-q*q)
       do idimn = 1, ndimn
          dwdx(idimn) = w * ( -2.* dx2(idimn)/hsml/hsml)
       enddo 
    else
    w    = 0.
    dwdx = 0.0     
endif 

else if (skf.eq.3) then      
    if (ndimn.eq.1) then
        factor = 1.e0 / (120.e0*hsml)
    elseif (ndimn.eq.2) then
        factor = 7.e0 / (478.e0*pi*hsml*hsml)
    elseif (ndimn.eq.3) then
        factor = 1.e0 / (120.e0*pi*hsml*hsml*hsml)
    else
        print *,' >>> Error <<< : Wrong dimension: Dim =',ndimn
        stop
    endif              
    if(q.ge.0.and.q.le.1) then
       w = factor * ( (3-q)**5 - 6*(2-q)**5 + 15*(1-q)**5 )
       do idimn= 1, ndimn
         dwdx(idimn) = factor * ( (-120 + 120*q - 50*q**2) / hsml**2 * dx2(idimn) ) 
       enddo 
    else if(q.gt.1.and.q.le.2) then
       w = factor * ( (3-q)**5 - 6*(2-q)**5 )
       do idimn= 1, ndimn
          dwdx(idimn) = factor * (-5*(3-q)**4 + 30*(2-q)**4) / hsml * (dx2(idimn)/r) 
       enddo 
    else if(q.gt.2.and.q.le.3) then
       w = factor * (3-q)**5 
       do idimn= 1, ndimn
          dwdx(idimn) = factor * (-5*(3-q)**4) / hsml * (dx2(idimn)/r) 
       enddo 
    else   
       w    = 0.
       dwdx = 0.0
    endif       
endif 
        
        END SUBROUTINE kernel

!-------------------------------------------------------------------------------
!------------------------------------------------------------------------------- 

          SUBROUTINE clean_up_sph

!closes open files      
!-------------------------------------------------------------------------------  
!-------------------------------------------------------------------------------

implicit none   

 close(dat_file)
 close(gid_res)  
 close(gid_msh)

        END SUBROUTINE clean_up_sph 

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

END MODULE SPH_main_2018


 
