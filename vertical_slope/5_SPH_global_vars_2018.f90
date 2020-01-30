
!--------------------------------------------------------


                 MODULE SPH_global_vars_2018  


!--------------------------------------------------------
                      
use variable_types

implicit none

TYPE pairs
     integer :: niac                             ! Used in a linked list describing part interactions
     integer :: pint_type			   
     integer :: pint_type2 !for ghost particles
     integer :: pair_i, pair_j                   ! two linked particles
     real    :: w                                ! Weight function and
     real    :: dwdx    !         derivatives
     real    :: dwdy 
     type    (pairs), pointer :: next   ! Pointer to next interacting pair of particles 
END TYPE pairs

type (pairs), pointer :: current, last              ! List described by pointers to current and last of list

real (kind=DP), allocatable:: x0       (:,:)   !  initial coordinates of particles. Updated at every tn 
real (kind=DP), allocatable:: x00      (:,:)   !  initial coordinates of particles. Used for relative displacements 
real (kind=DP), allocatable:: x        (:,:)   !  coordinates of particles 
real (kind=DP), allocatable:: vx0      (:,:)   !  velocities at time n 
real (kind=DP), allocatable:: vx       (:,:)   !  velocities of particles 
real (kind=DP), allocatable:: dvx      (:,:)   !  dvx = dvx/dt, force per unit mass 

real (kind=DP), allocatable:: mass     (:)     !  mass of particles 
real (kind=DP), allocatable:: mass0    (:)     !  mass of particles at time 0 
real (kind=DP), allocatable:: rho0     (:)     !  density at time tn   
real (kind=DP), allocatable:: rho      (:)     !  densities of particles
real (kind=DP), allocatable:: drho     (:)     !  drho =  drho/dt 

real (kind=DP), allocatable:: hsml     (:)     !  smoothing lengths of particles
real (kind=DP), allocatable:: c        (:)     !  sound velocity of particles

real (kind=DP), allocatable:: Dx_RK (:,:)	     ! These variables only for RK4 integrator
real (kind=DP), allocatable:: Dvx_RK(:,:) 
real (kind=DP), allocatable:: Du_RK   (:)
real (kind=DP), allocatable:: Drho_RK (:) 
real(kind=DP), allocatable :: grad_u(:,:,:), art_visc_out(:), disp_10(:),x_10(:,:)
real(kind=DP), allocatable :: disp_x(:,:),vec_x(:,:)
real(kind=DP), allocatable :: art_visc(:,:), adapt_stress(:), subset(:)
integer, allocatable :: keep_or_not(:), BC_or_not(:), bc_int(:)
real (kind=DP), allocatable :: f_drucker(:)
integer :: npoints
real (kind=DP) :: r_x, r_y
logical :: inside_approach, SP_SPH, art_stress, particle_shift
logical :: SPH_update, SPH_shift
integer :: n_update, shift_update
integer ::  ifsigman		

real (kind=DP) :: sml, alpha, beta, dx, dy
real (kind=DP) :: x_left, x_right, y_bottom, y_top
real (kind=DP) :: t0, t1, t2, t3, ft0, ft1, ft2, ft3
real (kind=DP) :: m1, m2, m3, c1, c2, c3

integer, allocatable:: itype    (:)     !  type of particles 
integer, allocatable:: countiac (:)     !  List of particles connected to each one.

integer :: vel_file, stress_file, time_file, node_file, node_bed_file, stress_point_file, geom_file, free_surf_file  
integer :: input_file, output_file
    
integer :: type_of_problem!  1 for SW, 2 for NS, 3 for DF, 11 for TGSolid

integer :: pa_sph         !  SPH algorithm for particle approximation (pa_sph)
                            !      pa_sph = 1 : (e.g. (p(i)+p(j))/(rho(i)*rho(j))
                            !               2 : (e.g. (p(i)/rho(i)**2+p(j)/rho(j)**2)
                            !               3 : uses 2 + alpha correction factor 		** NEW 15th Oct 2006
integer :: nnps           !  Nearest neighbor particle searching (nnps) method
                            !  nnps = 1 : Simplest and direct searching
                            !         2 : Sorting grid linked list
integer :: sle            !  Smoothing length evolution (sle) algorithm
                            !  sle = 0 : Keep unchanged,
                            !        1 : h = fac * (m/rho)^(1/dim)
                            !        2 : dh/dt = (-1/dim)*(h/rho)*(drho/dt)
                            !        3 : Other approaches (e.g. h = h_0 * (rho_0/rho)**(1/dim) )
integer :: skf            !  Smoothing kernel function 
                            !  skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
                            !        2, Gauss kernel   (Gingold and Monaghan 1981) 
                            !        3, Quintic kernel (Morris 1997)
logical :: summation_density, cont_density   ! . Use density summation model or density continuity model (or don't update density) 
                               
logical :: ghost_nodes, dummy_nodes   

logical :: MLS                      

integer :: ntotal, ntotal2, npoin, nvirt         ! ntotal  : real + virtual particles
                                          ! npoin   : number of real particles                                          ! nvirt   : virtual particles for bcs
integer :: ndimn                        ! dimension of the problem 1,2 3
integer :: niac                         ! Maximum number of interacting pairs.  
 
integer :: dat_file, chk_file, gid_msh, gid_res, pts_file    ! files for input/output
integer :: res_file1,res_file2, res_file3          ! files for input/output
integer :: init_dat, init_chk , node_n,stress_n, node_stressa_n, node_stressb_n, stressa_n, stressb_n
integer :: ndivx, ndivy, ndivxa, ndivya, ndivxb, ndivyb


integer, allocatable:: If_Out_Domain(:)                ! value is 1 if particle is outside control domain 
integer, allocatable:: If_Out_Domain_LAST(:)           ! las known value, to check wheteher part is just exited              
real (kind=DP), allocatable:: Xmin_Domain(:), Xmax_Domain(:)  ! Xmin_Domain(ndimn) and Xmax_Domain(ndimn)

real (kind=DP) :: xmin_geom, xmax_geom        ! max and min. X coordinates of domain
real (kind=DP) :: ymin_geom, ymax_geom        ! max and min. Y coordinates of domain
real (kind=DP) :: zmin_geom, zmax_geom        ! max and min. Z coordinates of domain

real (kind=DP) :: pi

character(15) :: problem_name            ! name of problem being solved.
character(60) :: text                    ! general purpose char string

integer :: n_avail_nodes                                !  we store here the list of nodes which we can inject. 
integer, allocatable :: list_avail_nodes(:)  	  !  set in set global vars, updated in injection
integer :: npoin0					  !  sometimes nodes 1 to npoin 0 are special (like behind gates)


integer :: ntype_solid			!  Type of solid problem 
						!       0 ... 1D
						!       1 ... plane stress
						!       2 ... plane strain
						!       3 ... axial symm.
						!       4 ... 3D

integer :: nprop, nmats                      !  number of mat props and number of mats
real (kind=DP), allocatable:: props       (:,:)   !  Material properties nmats * nprop
integer, allocatable:: matno         (:)   !  Material types of all nodes


integer :: nint_Vars                         !  number of internal vars
real (kind=DP), allocatable:: Internal_Vars (:,:)!  Internal Vars at nodes               nint_Vars * npoin1
real (kind=DP), allocatable:: Ddev_strn (:), Dvol_strn(:),  DPc(:)

integer :: IC_norm_TG				!  0 nothing 1 fi normalized 2 Fi and gradient normalized

integer  namaTG, nstre, nstr1, nstress, nnode, nmass, nelem, npoin1
!integer :: nelem, nbound, nbounds, nmiddle, nstress_mid
integer :: ndummy, ndummy_s, nrow, ndummy_total, nghost, ndummy2
integer, allocatable :: element_coord(:,:)

real (kind=DP), allocatable:: stress (:,:), vel(:,:), node(:,:)  !  vector of unks at main nodes						
real (kind=DP), allocatable::  RHS(:,:), stress_pt(:,:), conc(:,:), conc0(:)   !  vector of unks at intermediate nodes				
real (kind=DP)  , allocatable:: f1(:,:,:), f2(:,:,:) 													   
real (kind=DP), allocatable:: divef1(:,:), divef2(:,:)   !  div of fluxes a intr nodes	

real (kind=DP), allocatable :: f_bound(:,:), normal(:,:) 				

real, allocatable :: horizontal_or_not(:)
real, allocatable :: wall_position(:), n_int(:), abs_vel(:)

real, allocatable::  x_ghost(:,:), rho_ghost(:), mass_ghost(:), hsml_ghost(:), vel_ghost(:,:), stress_ghost(:,:)
integer, allocatable::  itype_ghost(:)
real, allocatable :: nghost_wall(:),nghost_sym(:)
real :: x_wall, x_sym
real, allocatable :: x_wall2(:), x_sym2(:)

integer, allocatable:: BNneighb_nod(:,:), BN_count_nod(:)
integer, allocatable:: BNnei_vct_nod(:,:)
integer  nmlspts_nod
integer, allocatable:: BNneighb_ele(:,:), BN_count_ele(:)
integer, allocatable:: BNnei_vct_ele(:,:)
integer :: nmlspts_ele

logical :: update_x, cspm, XSPH, vel_vector

real (kind=DP), allocatable:: displ(:,:)
real (kind=DP) :: disp_tol
integer :: no_bcs             !number of variables used to define prescribed values

integer, allocatable:: bc_info(:,:)   !  relation between node and the BCs applied
integer, allocatable:: intmaTG(:,:)    ! we store interm. nodes at elems intmaTG(nnodeTG,nelemTG)
integer, allocatable :: vel_out(:), stress_out(:)
integer :: rho_out, sml_out, disp_out, strain_out

!constitutive parameters
!-------------------------------------------------------------------------------
real (kind=DP) :: poiss, young, K_mod, G_mod, D11, D22, D12, D33, D41, D42

END MODULE SPH_global_vars_2018 
