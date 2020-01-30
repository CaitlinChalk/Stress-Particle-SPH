!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


                MODULE  SPH_material_2018

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!contains material and problem-specific information, including problem initialisation & boundary conditions
!Also contains particle-particle interaction definitions, and results output routines
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

use variable_types
use SPH_time_vars_2018
use SPH_global_vars_2018 
use SPH_material_vars_2018

implicit none

PRIVATE       !No variables are accesible from outside the module unless declared public 

public:: problem_input_data
public:: pint_update 
public:: OutputMesh, OutputRes
public:: plastic_terms
public:: BCs
public:: update_strain
public:: gravity_force 
public:: get_nodes_on_free_surface 
public:: apply_stress_free
public:: normal_bcs 
public:: adapt_stress2                  

CONTAINS

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

       SUBROUTINE problem_input_data

!reads problem information from data (material information) and input (main system controls) files      
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

implicit none

integer :: len1, nline, icunknoS, lmat, imats, iprop, ncrit 
integer :: iline, i, j, itotv, ipoin, idimn, iiabs, iivn0, k

real (kind=DP) :: twopi 

!Assign units to files, read problem name, open files
!-------------------------------------------------------------------------------

dat_file = 10
chk_file = 21
gid_msh = 41
gid_res = 42
pts_file = 73

read(input_file,*) text
read(input_file,*) problem_name
len1 = len_trim(problem_name)

open(dat_file,  file=problem_name(1:len1)//'.dat')
open(chk_file,  file=problem_name(1:len1)//'.chk')
open(gid_msh,  file=problem_name(1:len1)//'.post.msh')
open(gid_res,  file=problem_name(1:len1)//'.post.res')

!read lines of text describing the problem 
!-------------------------------------------------------------------------------

read (dat_file,*)  nline
write(chk_file,*)  nline    
do iline = 1,nline
   read (dat_file,*) text
   write(chk_file,*) text
enddo

!problem dimensions (the current code is only applicable for 2D problems...
!... but this keeps the option open for future developments)
!-------------------------------------------------------------------------------

read (dat_file,*) text
write(chk_file,*) text
read (dat_file,*) ndimn   
write(chk_file,*) ndimn  

!open .pts file (containing geometry information)
!-------------------------------------------------------------------------------

open(pts_file,  file=problem_name(1:len1)//'.pts')
  
!type of problem --- plane stress or plane strain (code currently only implemented in plaine strain)
!nstre = number of non-zero stress componenets
!-------------------------------------------------------------------------------
        
read (dat_file,*) text
write(chk_file,*) text
read (dat_file,*) nstre, nmats  ! nstre = 4 -- plane strain (only current functioning option)
write(chk_file,*) nstre, nmats  ! nstre = 3 = plane stress

nprop = 20  

!           For solid problems: 
!           props...Material properties (nmats * nprop)

!             1..ntype eco 1..elastic 2..EP  3 VP
!             2..ncrit (2 Von-Mises 5 Cam Clay 12 Drucker-Prager)
!             3..Young Modulus
!             4..Poisson
!             5..
!             6..density
!             7..Yield stress    !
!             8..Hardening modulus!                     
!             9..Phi Mohr (degrees)
!             10..gamma (Perzyna's viscoplasticity)
!             11..delta (Perzyna)
!             12..n     (Perzyna)
!             13..Cam Clay  M  (related to phi) / tan(phi)
!             14..          lambda / cohesion
!             15            kappa
!             16            Pc0 
!             17            e0 
! 

allocate ( props(nmats,nprop) )       
read (dat_file,*) text
write(chk_file,*) text
DO imats=1,nmats
   read (dat_file,*) lmat, (props(lmat,iprop),iprop=1,12)
   write(chk_file,*) lmat, (props(lmat,iprop),iprop=1,12)
   ncrit = props(lmat,2)
   if (ncrit == 5 .or. ncrit == 12 ) then  !  This is for Cam Clay or Drucker-Prager: read the extra model parameters.   
      read (dat_file,*) text    !  Provides compatibility 
      write(chk_file,*) text    !    with previous versions
      read (dat_file,*) (props(lmat,iprop),iprop=13,nprop)
      write(chk_file,*) (props(lmat,iprop),iprop=13,nprop)
   endif
enddo


 pi = 4*atan(1.0)
   
 Call setup_particles    !  ------  Setup nodes and stress-points
    
 Call Setup_Global_Arrays       !  ------  Allocate global arrays, pass local to global, deallocate, initialize

	countiac = 0        

 
!     prescribed values : bc_info : node nr. & variable
!                         bc_list : (a0+a1*sin(wt-fi))*(1-exp(-t/tf))
!                                  a0,a1,w,fi,tf
!     if time curves are used: 
!        v3 is given a negative value  (frequency never is)
!        v1 is the time curve being used
!        v2 is the factor by which we multiply all values f(t)
!        time will be read later, in this section
!        f    is a f(t). Value of unknown is v2*f(t)

read (dat_file,*) text
write(chk_file,*) text
read (dat_file,*) no_bcs, ifsigman   
write(chk_file,*)  no_bcs, ifsigman

if (no_bcs.gt.0) then
	
        no_bc_vars = 8
        allocate ( bc_list(no_bc_vars,no_bcs) )
        twopi = acos(0.0)*4.0/360.
        read (dat_file,*) text
        write(chk_file,*) text

        do i = 1,no_bcs
                read (dat_file,*) (bc_list(j,i),j=1,8)
                write(chk_file,*) (bc_list(j,i),j=1,8)
                bc_list(7,i)=bc_list(7,i)*twopi
        enddo

endif

read (dat_file,*) text
write(chk_file,*) text
read (dat_file,*) no_segments_bc
write(chk_file,*) no_segments_bc
If (no_segments_bc.gt.0) then
	read (dat_file,*) text
	write(chk_file,*) text
	allocate (Segment_BCs(5,no_segments_bc))
	do i = 1,no_segments_bc
	        read (dat_file,*) (Segment_BCs(j,i),j=1,5)
	        write(chk_file,*) (Segment_BCs(j,i),j=1,5)
	Enddo
Endif

read (dat_file,*) text
write(chk_file,*) text
read (dat_file,*) no_nodes_bc
write(chk_file,*) no_nodes_bc
If (no_nodes_bc.gt.0) then
	read (dat_file,*) text
	write(chk_file,*) text
	allocate (Nodal_BCs(2,no_nodes_bc)) 
	do i = 1,no_nodes_bc
	        read (dat_file,*) (Nodal_BCs(j,i),j=1,2)
	        write(chk_file,*) (Nodal_BCs(j,i),j=1,2)
	Enddo
Endif

Call Get_BCs_on_node


read (dat_file,*) text                 !ic_tcurve=0 ==> no time curve ; ic_tcurve=1 ==> time curve
write(chk_file,*) text
read (dat_file,*) ic_tcurve   
write(chk_file,*) ic_tcurve   


    if (ic_tcurve==1) then                    !    read time curves
                                        !    variables belong to sph_gfl_time_vars
       read (dat_file,*) text
       write(chk_file,*) text
       read (dat_file,*) ntcurves, mptstcurves ! number of time curves
       write(chk_file,*) ntcurves, mptstcurves ! and max nr of pts in all
       
       allocate ( nptstcurves(ntcurves) )               !
       allocate ( ttcurves(ntcurves,mptstcurves) )      !
       allocate ( ftcurves(ntcurves,mptstcurves) )      !
       
       do i = 1, ntcurves
          read (dat_file,*) text
          write(chk_file,*) text
          read (dat_file,*) nptstcurves(i)     !    nr of pts in each time curve
          write(chk_file,*) nptstcurves(i)
          read (dat_file,*) text
          write(chk_file,*) text
          read (dat_file,*) (ttcurves(i,j), j=1,nptstcurves(i))        !  t values
          read (dat_file,*) (ftcurves(i,j), j=1,nptstcurves(i))        !  f values
          write(chk_file,*) (ttcurves(i,j), j=1,nptstcurves(i))  
          write(chk_file,*) (ftcurves(i,j), j=1,nptstcurves(i)) 
       enddo
          
     endif
     
     
     t0 = ttcurves(1,1)
     t1 = ttcurves(1,2)
     t2 = ttcurves(1,3)
     t3 = ttcurves(1,4)
     
     ft0 = ftcurves(1,1)
     ft1 = ftcurves(1,2)
     ft2 = ftcurves(1,3)
     ft3 = ftcurves(1,4)
     
     m1 = (ft1 - ft0)/(t1 - t0)
     c1 = ft1 - m1*t1
     
     m2 = (ft2 - ft1)/(t2 - t1)
     c2 = ft2 - m2*t2
     
     m3 = (ft3 - ft2)/(t3 - t2)
     c3 = ft3 - m3*t3
     
!      ------  Read the domain limits. If nodes get outside, they will be deactivated

read (dat_file,*) text
write(chk_file,*) text  
read (dat_file,*)  (Xmin_Domain(idimn), idimn=1,ndimn),(Xmax_Domain(idimn), idimn=1,ndimn)
write(chk_file,*)  (Xmin_Domain(idimn), idimn=1,ndimn),(Xmax_Domain(idimn), idimn=1,ndimn)  
         
!  ****  Initial CONDITIONS: unknoS 
nint_Vars = 10     ! 1..Pc  2..ep_vol  3..ep_shear  Keep track at nodes & elems
allocate (internal_Vars(nint_Vars,ntotal) )     ;internal_Vars = 0.0 


 Call Initial_conditions
           
!      ------  read control parameters

read (dat_file,*) text
write(chk_file,*) text
read (dat_file,*)  pa_sph, nnps, sle, skf, CSPM, update_x, XSPH
write(chk_file,*)  pa_sph, nnps, sle, skf, CSPM, update_x, XSPH
                         
!      ------  read constants
        
read (dat_file,*) text
write(chk_file,*) text
read (dat_file,*) summation_density, cont_density
write(chk_file,*) summation_density, cont_density

read(dat_file,*) text
read(dat_file,*) DampingTG

read(dat_file,*) text
read(dat_file,*) alpha, beta

!      ------  read gravity force
read (dat_file,*) text
write(chk_file,*) text
read (dat_file,*) ic_grav
write(chk_file,*) ic_grav

If (ic_grav==1) then
    read (dat_file,*) text          
    write(chk_file,*) text          
    read (dat_file,*)    (cgrav(i),i=1,ndimn), tcurve_grav, ft_grav
    write(chk_file,*)    (cgrav(i),i=1,ndimn), tcurve_grav, ft_grav
Endif  

read(dat_file,*) text
read(dat_file,*) text
read(dat_file,*) stress_out(1), stress_out(2), stress_out(3), stress_out(4), & 
                vel_out(1), vel_out(2), strain_out, rho_out, sml_out, disp_out  


  
END SUBROUTINE problem_input_data

! ----------------------------------------------------------------------   

subroutine setup_particles

! ----------------------------------------------------------------------   

implicit none      

integer :: icunk_s, icc 

read (dat_file,*) text
write(chk_file,*) text
read (dat_file,*) icunk_s   
write(chk_file,*) icunk_s 

call Read_2DMesh

 
end subroutine  setup_particles

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

        subroutine Read_2DMesh
        
!set up all particles, initialise everything, read simulation information 
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


implicit none

real (kind=DP) :: x1, x2, x3, x4, y1, y2, y3, y4
real (kind=DP) :: lx, ly, tol, rx_factor, ry_factor
integer :: i, j, k, k1, l, idimn, geom_type
real (kind=DP) , allocatable :: area_1(:), area_s(:)
integer, allocatable :: no_int_node(:), no_int_stress(:)

!read simulation information (i.e. Standard SPH or S-P SPH)
!-------------------------------------------------------------------------------

read(input_file,*) text
read(input_file,*) SP_SPH, art_stress, particle_shift
if (SP_SPH) then
        read(input_file,*) text
        read(input_file,*) inside_approach !the inside or outside approach? read from input file
        read(input_file,*) text
        read(input_file,*) npoints !number of stress-points: per node for the outside approach, inside each virtual quadrilateral for the inside approach
        if (.not. (inside_approach)) then        
                read(input_file,*) text
                read(input_file,*) SPH_shift, vel_vector, shift_update, rx_factor, ry_factor, disp_tol
        end if
end if

if (vel_vector) npoints = 2 !2 stress-points per node required for the velocity vector S-P SPH approach
if (.not. SP_SPH) npoints = 1 !stress-points = nodes in Standard SPH

!read geometry 
!-------------------------------------------------------------------------------

read(pts_file,*) text
write(chk_file,*) text
read(pts_file,*) geom_type
write(chk_file,*) geom_type

!geom_type - defines initial geometry of material
!currently only a rectangle and trapezium have been implemented
!geom_type = 1 : rectangle, geom_type = 2 : trapezium
 
read(pts_file,*) text
write(chk_file,*) text
read(pts_file,*) x1, x2, x3, x4, dx
write(chk_file,*) x1, x2, x3, x4, dx
read(pts_file,*) text
write(chk_file,*) text
read(pts_file,*) y1, y2, y3, y4, dy
write(chk_file,*) y1, y2, y3, y4, dy
!coordinate order : (x1,y1) = bottom left corner, (x2,y2) = bottom right corner, (x3,y3) = top left corner, (x4,y4) = top right corner
!for rectangle: y1 = bottom, y4 = top, x1 = left, x4 = right 
read(pts_file,*) text
write(chk_file,*) text
read(pts_file,*) dummy_nodes
write(chk_file,*) dummy_nodes

x_left = x1; x_right = x4
y_bottom = y1; y_top = y4

!width and height of material and number of divisions in each
!-------------------------------------------------------------------------------        
lx = x4 - x1
ly = y4 - y1
        
ndivx = lx/dx + 1
ndivy = ly/dy + 1
        
!total number of nodes and stress-points
!-------------------------------------------------------------------------------

nnode = ndivx*ndivy
nelem = (ndivx-1)*(ndivy-1) !number of elements (required for gid output)

if ((inside_approach)) then

        if (npoints == 1) then
                nstress = (ndivx-1)*(ndivy-1)
        else if (npoints == 2) then
                nstress = 2*(ndivx-1)*(ndivy-1) 
        else if (npoints == 3) then
                nstress = 3*(ndivx-1)*(ndivy-1)
        else if (npoints == 4) then
                nstress = 4*(ndivx-1)*(ndivy-1)
        end if
        
else        

        nstress = npoints*(nnode)
        
end if

ntotal = nnode + nstress

!read smoothing length from input file 
!-------------------------------------------------------------------------------

read(input_file,*) text
read(input_file,*) sml 

!set-up dummy nodes for no-slip wall boundaries (if included)
!-------------------------------------------------------------------------------

if (dummy_nodes) then        
        call set_up_dummy_nodes
else        
    ndummy = 0   
end if

!allocate and initialise nodes
!-------------------------------------------------------------------------------
!no_int_node/no_int_stress defines whether the particle is interior (8), lies on the boundary (4) or lies on a corner (2), in order to distribute the mass appropriately  

allocate (x_1(ndimn,nnode), rho_1(nnode), mass_1(nnode), hsml_1(nnode))
allocate (itype_1(nnode), stress_1(nstre,nnode))
allocate(vel_1(ndimn,nnode),no_int_node(nnode),area_1(nnode))

itype_1 = 2 !marker for nodes
no_int_node = 0
vel_1= 0.0; stress_1 = 0.0; area_1 = 0.0 

!allocate and initialise elements for gid output
!-------------------------------------------------------------------

allocate(element_coord(4,nelem),x_out(ndimn,nnode))

element_coord = 0

!allocate and initialize stress points
!-------------------------------------------------------------------- 

allocate ( x_s(ndimn,nstress),  rho_s(nstress), mass_s(nstress), hsml_s(nstress) )
allocate ( itype_s(nstress), stress_s(nstre,nstress) )
allocate(vel_s(ndimn,nstress),no_int_stress(nstress),area_s(nstress))

itype_s = 1 !marker for stress-points
no_int_stress = 8;
vel_s = 0.0; stress_s = 0.0; area_s = 0.0


!-------------------------------------------
!define initial particle positions
!--------------------------------------------

!node positions
!--------------


k=1
do i = 1,ndivx
        do j = 1,ndivy
                x_1(1,k) =  x1 + (i-1)*dx
                x_1(2,k) =  y1 + (j-1)*dy                
                if (x_1(1,k) == x1 .or. x_1(1,k) == x4 .or. x_1(2,k) == y4) then
                        no_int_node(k) = 4
                else
                        no_int_node(k) = 8
                end if
                if ((x_1(1,k) == x1 .and. (x_1(2,k) == y1 .or. x_1(2,k) == y4)) .or. &
                        (x_1(1,k) == x4 .and. (x_1(2,k) == y1 .or. x_1(2,k) == y4))) then
                        no_int_node(k) = 2
                end if                                   
                k=k+1
        end do
end do


!elements for gid output
!---------------------------------------


i = 2
k = ndivy+ 2
j=1
do while (i <= ndivx)
        do while (k <= i*ndivy)
                element_coord(1,j) = k 
                element_coord(2,j) = k - ndivy
                element_coord(3,j) = k - ndivy - 1
                element_coord(4,j) = k - 1
                k=k+1
                j = j+1
        end do
        i = i+1
        k=k+1
end do


!stress-point positions
!--------------------------------------------------------

if (inside_approach) then

        if (npoints == 1) then
                
                k = 1
                do i = 2,ndivx
                        do j = 2,ndivy
                                x_s(1,k) =  x1 - dx/2. + (i-1)*dx
                                x_s(2,k) =  y1 - dx/2. + (j-1)*dy                   
                                k=k+1
                        end do
                end do
                
        else if (npoints == 2) then

                k=1
                do j = 1,ndivx-1
                        do i = 1,ndivy-1
                                x_s(1,k) = (x_1(1,i+(j-1)*ndivy) + x_1(1,i+1+(j-1)*ndivy) + x_1(1,i+ndivy+(j-1)*ndivy))/3.
                                x_s(2,k) = (x_1(2,i) + x_1(2,i+1) + x_1(2,i+ndivy))/3.
                                
                                x_s(1,k+1) = (x_1(1,i+1+(j-1)*ndivy) + x_1(1,i+ndivy+(j-1)*ndivy) + x_1(1,i+ndivy+(j-1)*ndivy))/3.
                                x_s(2,k+1) = (x_1(2,i+1) + x_1(2,i+ndivy) + x_1(2,i+ndivy+1))/3.
                                k=k+2
                        end do
                end do
                
         else if (npoints == 3) then
         
                k=1
                do j = 1,ndivx-1
                        do i = 1,ndivy-1
                                x_s(1,k) = (x_1(1,i+(j-1)*ndivy) + x_1(1,i+1+(j-1)*ndivy) + x_1(1,i+ndivy+(j-1)*ndivy))/3.
                                x_s(2,k) = (x_1(2,i) + x_1(2,i+1) + x_1(2,i+ndivy))/3.              
                
                                x_s(1,k+1) = (x_1(1,i+(j-1)*ndivy) + x_1(1,i+ndivy+(j-1)*ndivy) + x_1(1,i+1+ndivy+(j-1)*ndivy))/3.
                                x_s(2,k+1) = (x_1(2,i) + x_1(2,i+ndivy) + x_1(2,i+ndivy+1))/3.
                
                                x_s(1,k+2) = (x_1(1,i+1+(j-1)*ndivy) + x_1(1,i+1+ndivy+(j-1)*ndivy) + (x_1(1,i+ndivy+(j-1)*ndivy)&
                                + x_1(1,i+(j-1)*ndivy))/2.)/3.
                                x_s(2,k+2) = (x_1(2,i+1) + x_1(2,i+1+ndivy) + (x_1(2,i)+x_1(2,i+1))/2.)/3.                
                                k=k+3
                        end do
                end do
                
         else if (npoints == 4) then
                         
                k=1
                do j = 1,ndivx-1
                        do i = 1,ndivy-1
                                x_s(1,k) = (x_1(1,i+(j-1)*ndivy) + x_1(1,i+(j-1)*ndivy) + x_1(1,i+1+ndivy+(j-1)*ndivy))/3.
                                x_s(2,k) = (x_1(2,i) + x_1(2,i+1) + x_1(2,i+ndivy))/3.
                                x_s(1,k+1) = (x_1(1,i+1+(j-1)*ndivy) + x_1(1,i+ndivy+(j-1)*ndivy) + x_1(1,i+ndivy+1+(j-1)*ndivy))/3.
                                x_s(2,k+1) = (x_1(2,i+1) + x_1(2,i+ndivy) + x_1(2,i+ndivy+1))/3.
                
                                x_s(1,k+2) = (x_1(1,i+(j-1)*ndivy) + x_1(1,i+1+(j-1)*ndivy) + x_1(1,i+1+ndivy+(j-1)*ndivy))/3.
                                x_s(2,k+2) = (x_1(2,i) + x_1(2,i+1) + x_1(2,i+1+ndivy))/3.
                                x_s(1,k+3) = (x_1(1,i+(j-1)*ndivy) + x_1(1,i+ndivy+(j-1)*ndivy) + x_1(1,i+1+ndivy+(j-1)*ndivy))/3.
                                x_s(2,k+3) = (x_1(2,i) + x_1(2,i+ndivy) + x_1(2,i+ndivy+1))/3.
                
                                k=k+4
                        end do
                end do
                
          end if
                
else 

        r_x = dx*rx_factor; r_y = dx*ry_factor !define shifting distances for the outside approach, in terms of dx
        
        if (npoints == 1) then
        
                x_s(1,:) = x_1(1,:) + r_x
                x_s(2,:) = x_1(2,:) + r_y
                
        else if (npoints == 2) then
        
                k=1
                do i = 1,nnode
                        x_s(1,k) = x_1(1,i) + r_x
                        x_s(2,k) = x_1(2,i) + r_y
                        x_s(1,k+1) = x_1(1,i) - r_x
                        x_s(2,k+1) = x_1(2,i) - r_y
                        k=k+2
                end do 
                                             
        else if (npoints == 3) then
        
                k=1
                do i = 1,nnode
                        x_s(1,k) = x_1(1,i)
                        x_s(2,k) = x_1(2,i) + r_y
                        x_s(1,k+1) = x_1(1,i) - r_x
                        x_s(2,k+1) = x_1(2,i) - r_y
                        x_s(1,k+2) = x_1(1,i) + r_x
                        x_s(2,k+2) = x_1(2,i) - r_y
                        k = k+3
                end do  
                
         else if (npoints == 4) then
         
                k = 1
                do i = 1,nnode
                        x_s(1,k) = x_1(1,i) - r_x
                        x_s(2,k) = x_1(2,i) - r_y
                        x_s(1,k+1) = x_1(1,i) - r_x
                        x_s(2,k+1) = x_1(2,i) + r_y
                        x_s(1,k+2) = x_1(1,i) + r_x
                        x_s(2,k+2) = x_1(2,i) + r_y
                        x_s(1,k+3) = x_1(1,i) + r_x
                        x_s(2,k+3) = x_1(2,i) - r_y
                        k = k+4
                end do
                      
         end if
         
end if         


!Initial density, mass, smoothing length
!-------------------------------------------------------
        
do i = 1,nnode
        rho_1 (i) = props (1,6)        
        hsml_1(i) = sml*dx
        area_1(i) = (rho_1(i)*dx*dy)/8 
        mass_1(i) = (no_int_node(i)*area_1(i))
end do


do i = 1,nstress
        rho_s(i) = props(1,6)
        hsml_s(i) = sml*dx
        area_s(i) = (rho_s(i)*dx*dy)/8      
end do


!calculate the mass of the stress-points (sum of stress-point mass = sum of nodes = total mass (approx))

if (inside_approach) then

        mass_s = no_int_stress*area_s
        mass_s = mass_s/npoints
     
else
        if (npoints /= 2) then
        
                k=1
                do i = 1,nnode
                        do j = 1,npoints
                                mass_s(k+(j-1)) = mass_1(i)/npoints
                        end do
                        k=k+npoints
                end do               
        
        
        else if (npoints == 2) then !"turning off" outer stress-points is currently only applied for 2 points per node
        
        
        k=1
        do i = 1,nnode
                do j = 1,npoints
                        mass_s(k+(j-1)) = mass_1(i)/npoints
                end do
                k=k+npoints
        end do 
        
         end if
        
end if        

!ntotal2 is the total number of particles - nodes, stress-points and dummy nodes

ntotal2 = ntotal + ndummy

!read material constants 
!-------------------------------------------------------------------------------

young = props(1,3) !young's modulus
poiss  = props(1,4)  !poisson's ratio                                    
K_mod = young/(3*(1-2*Poiss)) !bulk modulus
G_mod = young/(2*(1+Poiss)) !shear modulus
!elastic constitutive matrix constants
!-------------------------------------
D11 = (4./3.)*G_mod + K_mod 
D22 = D11
D12 = -(2./3.)*G_mod + K_mod
D33 = G_mod
D41 = D12
D42 = D12
    
END subroutine Read_2DMesh

!---------------------------------------------------------

subroutine set_up_dummy_nodes

!----------------------------------------------------------
!-------------------------------------
!dummy particles for walls
!-------------------------------------

implicit none

!real :: wallh_position, wallv_position
real :: wallh_x1, wallh_x2, wallv_y1, wallv_y2
real :: l_wall, dx2, dy2
integer :: n_walls
integer, allocatable :: wall_id(:), ndiv_wall(:)
real, allocatable :: wallx(:,:), wall_position2(:)
integer :: i,j,k,l,m

!define 
!-------------------------------------

nrow = 3 !number of rows/columns of dummy particles

dx2 = dx !particle spacing for dummy particles
dy2 = dy

!read boundary info from file
!-------------------------------------
!wall_id defines the wall orientation:
!1 - vertical LHS wall
!2 - horizontal bottom wall
!3 - vertical RHS wall

read(pts_file,*) text
write(chk_file,*) text
read(pts_file,*) n_walls
write(chk_file,*) n_walls

allocate(wall_id(n_walls), wall_position2(n_walls), wallx(2,n_walls), ndiv_wall(n_walls))

ndummy2 = 0
do i = 1,n_walls
        read(pts_file,*) text
        read(pts_file,*) text
        read(pts_file,*) wall_id(i), wall_position2(i), wallx(1,i), wallx(2,i)
        l_wall = wallx(2,i) - wallx(1,i)
        ndiv_wall(i) = l_wall/(dx2) + 1
        ndummy2 = ndummy2 + ndiv_wall(i)
end do

ndummy = nrow*ndummy2 !total number of dummy nodes

!allocate and initialise dummy node arrays
!------------------------------------------

allocate(x_dummy(ndimn,ndummy),itype_dummy(ndummy),rho_dummy(ndummy) &
,mass_dummy(ndummy),hsml_dummy(ndummy),vel_dummy(ndimn,ndummy), stress_dummy(nstre,ndummy))
allocate(horizontal_or_not(ntotal+ndummy), wall_position(ntotal+ndummy))

itype_dummy = 25 !marker for dummy node
vel_dummy = 0.0
stress_dummy = 0.0
horizontal_or_not = 0.0

!define the inner most layer first, to make it easier to apply the boundary repulsive force if needed
k=1
l = ntotal+1
do i = 1,n_walls
        do j = 1,ndiv_wall(i)
                if (wall_id(i) == 1 .or. wall_id(i) == 3) then !vertical walls
                        x_dummy(1,k) = wall_position2(i)
                        x_dummy(2,k) = wallx(1,i) + (j-1)*dy2
                        if (wall_id(i) == 1) horizontal_or_not(l) = 2 !LHS vertical wall
                        if (wall_id(i) == 3) horizontal_or_not(l) = 22 !RHS vertical wall
                else if (wall_id(i) == 2) then !horizontal wall
                        x_dummy(1,k) = wallx(1,i) + (j-1)*dx2
                        x_dummy(2,k) = wall_position2(i)
                        horizontal_or_not(l) = 1
                end if
                wall_position(l) = wall_position2(i)
                k=k+1
                l=l+1
         end do
end do
                        
!define outer dummy layers 

do i = 1,n_walls
        do m = 2,nrow
                do j = 1,ndiv_wall(i)
                        if (wall_id(i) == 1 .or. wall_id(i) == 3) then !vertical walls
                                x_dummy(1,k) = wall_position2(i) - (m-1)*dx2
                                x_dummy(2,k) = wallx(1,i) + (j-1)*dy2
                                if (wall_id(i) == 1) horizontal_or_not(l) = 2 !LHS vertical wall
                                if (wall_id(i) == 3) horizontal_or_not(l) = 22 !RHS vertical wall
                        else if (wall_id(i) == 2) then !horizontal wall
                                x_dummy(1,k) = wallx(1,i) + (j-1)*dx2
                                x_dummy(2,k) = wall_position2(i) - (m-1)*dy2
                                horizontal_or_not(l) = 1
                        end if
                        wall_position(l) = wall_position2(i)
                        k=k+1
                        l=l+1
                end do
         end do
end do

!define dummy node information
!------------------------------

do k = 1,ndummy

        rho_dummy(k) = props(1,6)
        mass_dummy(k) = ((dx2)**2)*props(1,6) !mass of interior nodes
       ! do i = 1,n_walls
                !if (wall_id(i) == 1 .or. wall_id(i) == 3) then
                !        if (x_dummy(1,k) == wall_position(i)) then
                !                mass_dummy(k) = mass_dummy(k)/2. !reduce mass on 
                !else if (x_dummy(1,k) == wallv_position1 .or. x_dummy(1,k) == wallv_position2) then
                !mass_dummy(k) = mass_dummy(k)/2.
        !end if
        hsml_dummy(k) = sml*dx !0
        vel_dummy(:,k) = 0.0
        stress_dummy(:,k) = 0.0

end do
               

!write geometry to output file
!-----------------------------

geom_file = 57

open(unit=geom_file,file='geom.csv')

write(geom_file,*) 'x-coord,', 'y-coord,', 'z-coord'

do i = 1,ndummy
        write(geom_file,*) x_dummy(1,i), ',', x_dummy(2,i), ',', 0.0
end do

 close(unit=geom_file)

end subroutine


!----------------------------------------------------------------------------------

        Subroutine Setup_Global_Arrays 

!----------------------------------------------------------------------------------

!      ------  Allocate global arrays, pass local to global, deallocate, initialize

implicit none

integer :: itotv, ipoin, idimn,iiabs, iivn0
integer :: it1, it2, it3, i, istre

allocate ( x00(ndimn,ntotal2),  x0 (ndimn,ntotal2),  x(ndimn,ntotal2), vx0(ndimn,ntotal2), vx(ndimn,ntotal2) )
allocate ( mass(ntotal2), rho(ntotal2),  hsml(ntotal2), c(ntotal) )
allocate ( mass0  (ntotal2) )
allocate ( rho0(ntotal2), drho(ntotal2) ) 
allocate ( itype(ntotal2), countiac(ntotal2) )
allocate (If_Out_Domain(ntotal2))
allocate (Xmin_Domain(ndimn))
allocate (Xmax_Domain(ndimn))
allocate(art_visc(ndimn,nnode))

allocate ( matno(ntotal2) )
allocate( stress(nstre,ntotal2), vel(ndimn,ntotal2))
allocate( RHS (ndimn,ntotal))
allocate ( f1 (ndimn,nstre,ntotal), f2 (ndimn,nstre,ntotal) )
allocate ( divef1(nstre,ntotal), divef2(ndimn,ntotal) )
allocate (Ddev_strn(ntotal), Dvol_strn(ntotal),  DPc(ntotal))

allocate (displ(ndimn,nnode))
allocate (cgrav(ndimn))
allocate (f_bound(ndimn,nnode))
allocate(adapt_stress(ntotal))
allocate(keep_or_not(ntotal),disp_10(nnode),x_10(ndimn,nnode), disp_x(ndimn,nnode),vec_x(ndimn,nnode))

allocate (subset(ntotal))
allocate (normal(ndimn,ntotal))
allocate(grad_u(ndimn,ndimn,ntotal))
allocate(art_visc_out(nnode),n_int(nnode),abs_vel(nnode))
allocate(f_drucker(ntotal),bc_int(nnode),conc(ndimn,ntotal),conc0(nnode))


allocate(stress_out(nstre), vel_out(ndimn))

f_drucker = 0.0

mass0 = 0.0; disp_10 = 0.0; vec_x = 0.0; disp_x = 0.0

grad_u = 0.0
art_visc_out = 0.0
art_visc = 0.0
adapt_stress = 0.0
keep_or_not = 0
n_int = 0.0; abs_vel = 0.0
if_out_domain = 0

npoin0 = 0
        
!      ------  From local main nodes to global array, then deallocate

itotv = 0                           
  

!nodes
!--------------------------------- 
    do i = 1, nnode 
      itotv = itotv + 1
      do idimn = 1, ndimn
                x  (idimn,itotv) = x_1 (idimn,i)		
		vel(idimn,itotv) = vel_1(idimn,i)
      enddo
      do istre = 1,nstre
        stress(istre,itotv) = stress_1(idimn,i)
      end do
      itype (itotv) = itype_1(i)
      rho   (itotv) = rho_1  (i)
      mass  (itotv) = mass_1 (i)
      mass0 (itotv) = mass_1 (i)
      hsml  (itotv) = hsml_1 (i)
   enddo
   
  
   deallocate (x_1, itype_1, rho_1, mass_1, hsml_1)
  

!stress points
!------------------------------------- 

    do ipoin    = 1, nstress
      itotv    = itotv + 1
      do idimn = 1, ndimn
                x  (idimn,itotv) = x_s (idimn,ipoin)
		vel(idimn,itotv) = vel_s(idimn,ipoin)
      enddo
      do istre = 1,nstre
        stress(istre,itotv) = stress_s(istre,ipoin)
      end do
      itype (itotv) = itype_s(ipoin)
      rho   (itotv) = rho_s(ipoin)
      mass  (itotv) = mass_s(ipoin)
      hsml  (itotv) = hsml_s(ipoin)
   enddo
   
   deallocate (x_s, itype_s, rho_s, mass_s, hsml_s)
   
!ghost nodes
!-------------------------------------------
if (dummy_nodes) then

    do ipoin    = 1, ndummy
      itotv    = itotv + 1
      do idimn = 1, ndimn
                x  (idimn,itotv) = x_dummy (idimn,ipoin)
		vel(idimn,itotv) = vel_dummy(idimn,ipoin)
      enddo
      do istre = 1,nstre
        stress(istre,itotv) = stress_dummy(istre,ipoin)
      end do
      itype (itotv) = itype_dummy(ipoin)
      rho   (itotv) = rho_dummy(ipoin)
      mass  (itotv) = mass_dummy(ipoin)
      hsml  (itotv) = hsml_dummy(ipoin)
   enddo
   
   deallocate (x_dummy, itype_dummy, rho_dummy, mass_dummy, hsml_dummy)
end if  
        
!      ------  Initialize (changed, because we have to initialize x,vx,itype,rho,mass,hsml,u and p)
    
  x0       = x  ;  x00      =  x0;  vx0     = vx  
  rho0     = rho; drho     = 0.0  
  x_10 = x(:,1:nnode)
  f1    = 0.0;  f2   = 0.0
  divef1 = 0.0; divef2 = 0.0
  Ddev_strn = 0.0; Dvol_strn  = 0.0;   DPc = 0.0
  stress = 0.0; stress_n = 0.0; vel = 0.0

  displ = 0.0
  normal = 0.0
  subset = 0.0
  !nnode = nnode + ndummy
  !nstress = nstress+ndummy
 

END SUBROUTINE Setup_Global_Arrays

!------------------------------------------

        subroutine Get_BCs_on_node

!------------------------------------------

implicit none

integer :: ipoin, isegment, inodal, ic_node
real (kind=DP) :: xx1, yy1, xx2, yy2, xx, yy
real (kind=DP) :: xx1s, yy1s, xx2s, yy2s
real (kind=DP) :: ymax, ymin, xmax, xmin
real (kind=DP) :: z, z2, zs
integer :: k, no_bc, bc_segment
real (kind=DP) :: bc_var, tol

allocate(bc_info(8,ntotal))
allocate(BC_or_not(ntotal))

BC_or_not = 0
bc_info = 0.0
tol = dx/10000.

Do ipoin = 1, ntotal
bc_info(1,ipoin) = ipoin
Enddo

If (no_Segments_bc.gt.0) then
	Do isegment=1,no_Segments_bc
	        xx1	= Segment_BCs(1,isegment)
	        yy1	= Segment_BCs(2,isegment)
	        xx2	= Segment_BCs(3,isegment)
	        yy2	= Segment_BCs(4,isegment)
	        xmax = max(xx1,xx2)
	        xmin = min(xx1,xx2)
	        ymax = max(yy1,yy2)
	        ymin = min(yy1,yy2)
		Do ipoin = 1, ntotal
			no_bc = Segment_BCs(5,isegment)		! which BC has to be applied ?
			bc_var = bc_list(2,no_bc)		! On which variable is the BC applied ?
			xx = x(1,ipoin)
			If (ndimn==2) then
			        yy = x(2,ipoin)
			Endif
			z=(yy1-yy2)*(xx-xx1)+(xx2-xx1)*(yy-yy1) 
			If (z == 0.)  then     !ipoin is on the boundary
				If ((xx >= xmin) .AND. (xx <= xmax)) then
					If ((yy >= ymin) .AND. (yy <= ymax)) then
						bc_info(2,ipoin)= bc_info(2,ipoin)+1  !how many boundary conditions on this node?
						k=bc_info(2,ipoin)
						bc_info(k+2,ipoin)=no_bc
						if (BC_or_not(ipoin) == 0) then
						        BC_or_not(ipoin) = 1
						end if						
					Endif
			        Endif
			Endif
		Enddo
				
	Enddo
Endif


If  (no_nodes_bc.gt.0) then
	Do inodal = 1,no_nodes_bc
		ic_node = Nodal_BCs(1,inodal)
		no_bc   = Nodal_BCs(2,inodal)		! which BC has to be applied ?
		bc_var  = bc_list(2,no_bc)			! On which variable is the BC applied ?
		bc_info(2,ic_node)=bc_info(2,ic_node)+1
		k=bc_var+2
			If (bc_info(k,ic_node).ne.0 .AND. bc_info(k,ic_node).ne.no_bc) then
			        write(*,*) ic_node, 'has BC number', bc_info(k,ic_node), 'and', no_bc, 'applied for the same variable', bc_var
			        write(*,*) 'Which BC do you want to apply to this node?', bc_info(k,ic_node), 'or', no_bc
			        read (*,*) no_bc
			        bc_info(2,ic_node)=bc_info(2,ic_node)-1
			Endif
		bc_info(k,ic_node)=no_bc
		if (BC_or_not(ic_node) == 0) then
		        BC_or_not(ic_node) = 1
		end if								
	Enddo
Endif


END Subroutine Get_BCs_on_node

!-----------------------------------------

subroutine get_nodes_on_free_surface

!-----------------------------------------

!detects the nodes on the free surface and approximates their normal vectors

implicit none 

integer :: i, j, idimn, ipoin, i1,i11,i12,i2,i21,i22
real :: hf
real :: A_norm (5, ntotal)
real :: hf1i, hf1j, hf2i, hf2j, gradW
real :: lambda(2,npoin1)
real :: aa, bb, cc, dd
real :: lambda_min
real :: vv(3,ntotal), ff(2,ntotal)
real :: tt(2,ntotal)
real :: tau(2,ntotal)
real :: xji, xjt_norm, prod_scal
real :: xjt(2)
real :: xij, xit_norm
real :: xit(2)
real :: limit
real :: m, m2

real :: x1, x2, y1, y2, x_vect, y_vect, norm_vect
real, allocatable :: neighbour(:)

real :: f_int(ntotal)
real :: grad_f_int(ndimn,ntotal)
real ::  p_scal

subset = 0.0
!subset(nnode+1:ntotal) = 2.
hf = 0.0
vv = 0.0
ff = 0.0
A_norm = 0.0
tt = 0.0
p_scal = 0.0

allocate (neighbour(ntotal))
f_int = 0.0
grad_f_int = 0.0

!Step 1: calculate the renormalisation matrix

 current => last
 do while (associated(current))
        i = current%pair_i !stress point
        j = current%pair_j !node 
       if (current%pint_type == 2 .or. current%pint_type == 3) then !node-node or stress pt - stress pt
      !  if (current%pint_type == 1) then       
                gradW = current%dwdx
                hf1i =  mass(j)*gradW/rho(j)
	        hf1j = - mass(i)*gradW/rho(i)
                A_norm(1,i)=A_norm(1,i)+(x(1,j)-x(1,i))*hf1i		! A11 for node i
                A_norm(1,j)=A_norm(1,j)+(x(1,i)-x(1,j))*hf1j            ! A11 for node j
                ff(1,i)    =ff(1,i) + hf1i !+ hf1i		   
                ff(1,j)    =ff(1,j) + hf1j !+ 1.e9*hf1j
                            
                if (ndimn==2) then
      	                hf2i = mass(j)*current%dwdy/rho(j)
                        hf2j = -mass(i)*current%dwdy/rho(i)
                        A_norm(2,i)=A_norm(2,i)+(x(2,j)-x(2,i))*hf1i		    ! A12 for node i
                        A_norm(3,i)=A_norm(3,i)+(x(1,j)-x(1,i))*hf2i		    ! A21 for node i
                        A_norm(4,i)=A_norm(4,i)+(x(2,j)-x(2,i))*hf2i		    ! A22 for node i
                        A_norm(2,j)=A_norm(2,j)+(x(2,i)-x(2,j))*hf1j		    ! A12 for node j
                        A_norm(3,j)=A_norm(3,j)+(x(1,i)-x(1,j))*hf2j		    ! A21 for node j
                        A_norm(4,j)=A_norm(4,j)+(x(2,i)-x(2,j))*hf2j		    ! A22 for node j
                        ff(2,i)=ff(2,i) + hf2i !+ 1.e9*hf2i		   
                        ff(2,j)=ff(2,j) + hf2j !+ 1.e9*hf2j            
                end if
                
        else if (current%pint_type == 1) then !node-stress point
                gradW = current%dwdx
                hf1i =   mass(j)*gradW/rho(j)
	        hf1j =  - mass(i)*gradW/rho(i)
                A_norm(1,i)=A_norm(1,i)+(x(1,j)-x(1,i))*hf1i		! A11 for node i
                A_norm(1,j)=A_norm(1,j)+(x(1,i)-x(1,j))*hf1j            ! A11 for node j
                ff(1,i)    =ff(1,i) + hf1i		   
                ff(1,j)    =ff(1,j) + hf1j 
                            
                if (ndimn==2) then
      	                hf2i =  mass(j)*current%dwdy/rho(j)
                        hf2j =  - mass(i)*current%dwdy/rho(i)
                        A_norm(2,i)=A_norm(2,i)+(x(2,j)-x(2,i))*hf1i		    ! A12 for node i
                        A_norm(3,i)=A_norm(3,i)+(x(1,j)-x(1,i))*hf2i		    ! A21 for node i
                        A_norm(4,i)=A_norm(4,i)+(x(2,j)-x(2,i))*hf2i		    ! A22 for node i
                        A_norm(2,j)=A_norm(2,j)+(x(2,i)-x(2,j))*hf1j		    ! A12 for node j
                        A_norm(3,j)=A_norm(3,j)+(x(1,i)-x(1,j))*hf2j		    ! A21 for node j
                        A_norm(4,j)=A_norm(4,j)+(x(2,i)-x(2,j))*hf2j		    ! A22 for node j
                        ff(2,i)=ff(2,i) + hf2i !+ 1.e9*hf2i		   
                        ff(2,j)=ff(2,j) + hf2j !+ 1.e9*hf2j            
                end if
         end if
        current=>current%next
 end do
 
!Step 2: Approximate the normal vector 
 
do i = 1,ntotal

     !ff(:,i) = ff(:,i)*1.e-9
     
     do i2=1,4
         if (abs(A_norm(i2,i)) <= 1.e-8) then   
             A_norm(i2,i) = 0.0
         end if
     end do
     
     vv(1,i) = -(A_norm(1,i)*ff(1,i) + A_norm(2,i)*ff(2,i))
     vv(2,i) = -(A_norm(3,i)*ff(1,i) + A_norm(4,i)*ff(2,i))
     vv(3,i) = (vv(1,i)**2 + vv(2,i)**2)**0.5
     normal(1,i) = vv(1,i)/vv(3,i)
     normal(2,i) = vv(2,i)/vv(3,i) 
     
end do

!Step 3: Define scan region and check if particles are in the region

! Calculate the orthogonal vector to the normal	
do i = 1,ntotal

    tt(1,i) = x(1,i) + hsml(i)*normal(1,i)
    tt(2,i) = x(2,i) + hsml(i)*normal(2,i)
    tau(1,i) = -normal(2,i)
    tau(2,i) = normal(1,i)
    
end do

! Check if there is a SPH particle inside the control area

 current => last
    DO WHILE (associated(current))
	i = current%pair_i !stress-pt
	j = current%pair_j !node
	
	    if (current%pint_type == 2 .or. current%pint_type == 3 ) then    !node-node and stress-point stress-point or node-stress pt
	    
                    if (subset(i) == 0) then
                    
  		        xji = ( (x(1,j)-x(1,i))**2 + (x(2,j)-x(2,i))**2 )**0.5
  		        xjt(1) =x(1,j)-tt(1,i)
  		        xjt(2) =x(2,j)-tt(2,i)
  		        xjt_norm = ((xjt(1))**2 + (xjt(2))**2 )**0.5
  		        prod_scal = abs(normal(1,i)*xjt(1) + normal(2,i)*xjt(2)) + abs( tau(1,i)*xjt(1)+ tau(2,i)*xjt(2))
  		        limit = (2**0.5)*hsml(i)  
  		          
  		        if (xji >= limit .and. xjt_norm < hsml(i)) then  !belongs to set 1
  		                subset(i) = 2
  		                f_int(i) = -1
  		        else if (xji < limit .and. prod_scal < hsml(i)) then !belongs to set 2
   		                subset(i) = 2
   		                f_int(i) = -1
 		        end if
 		        
  		    end if

 		 if (subset(j) == 0) then
  		    
  		        xij = ( (x(1,i)-x(1,j))**2 + (x(2,i)-x(2,j))**2 )**0.5
  		        xit(1) =x(1,i)-tt(1,j)
  		        xit(2) =x(2,i)-tt(2,j)
  		        xit_norm = ((xit(1))**2 + (xit(2))**2 )**0.5
  		        prod_scal = abs(normal(1,j)*xit(1) + normal(2,j)*xit(2)) + abs( tau(1,j)*xit(1)+ tau(2,j)*xit(2))
  		        limit = (2**0.5)*hsml(j)  
  		          
  		        if (xij >= limit .and. xit_norm < hsml(j)) then !belongs to set 1
  		                subset(j) = 2
  		                f_int(j) = -1
  		        else if(xij < limit .and. prod_scal < hsml(j)) then !belongs to set 2
   		                subset(j) = 2
   		                f_int(j) = -1
 		        end if
 		        
  		    end if 
  		    
  	     end if
                
             if (current%pint_type == 6 .or. current%pint_type == 9 ) then 
          		    
  		    if (subset(j) == 0) then
  		    
  		        xij = ( (x(1,i)-x(1,j))**2 + (x(2,i)-x(2,j))**2 )**0.5
  		        xit(1) =x(1,i)-tt(1,j)
  		        xit(2) =x(2,i)-tt(2,j)
  		        xit_norm = ((xit(1))**2 + (xit(2))**2 )**0.5
  		        prod_scal = abs(normal(1,j)*xit(1) + normal(2,j)*xit(2)) + abs( tau(1,j)*xit(1)+ tau(2,j)*xit(2))
  		        limit = (2**0.5)*hsml(j)  
  		          
  		        if (xij >= limit .and. xit_norm < hsml(j)) then !belongs to set 1
  		                subset(j) = 2
  		                f_int(j) = -1
  		        else if(xij < limit .and. prod_scal < hsml(j)) then !belongs to set 2
   		                subset(j) = 2
   		                f_int(j) = -1
 		        end if
 		        
  		    end if  		     		    
  	end if
  		
    current=>current%next
    
    end do  
    
! We put nodes on the free surface in subset 1
do i = 1,ntotal
        if (subset(i) == 0) then
                subset(i) = 1
        end if
end do        

!we now have, if subset(i) = 1, i is on the free surface. If subset(i) = 2, it is not.

! find nodes which are on the surface, but have other boundary conditions

do i = 1,ntotal
        if (bc_or_not(i) /= 1) then
                bc_or_not(i) = 0
        end if
end do

 do i = 1,ntotal
        if (subset(i) == 2 .and. bc_or_not(i) /=1) then
                bc_or_not(i) = 0
        else if (subset(i) == 1 .and. bc_or_not(i) /= 1) then
                bc_or_not(i) = 2
        else if (subset(i) == 3 .and. bc_or_not(i) /= 1) then
                !bc_or_not(i) = 3
        else if (subset(i) == 1 .and. bc_or_not(i) == 1) then
                subset(i) = 2
        end if
 end do
 
 !Step 4: Get a better approximation to the normal
 
normal = 0.0
neighbour = 0
 current => last
DO WHILE (associated(current))
        i = current%pair_i 
	j = current%pair_j 
	if (current%pint_type == 2 .or. current%pint_type == 3 .or. current%pint_type == 1) then     !  i..node j..node .or. current%pint_type == 1
        if (subset(i) == 1 .and. subset(j) == 1) then           ! i and j on the boundary
		    x1 = x(1,i)
		    y1 = x(2,i)
		    x2 = x(1,j)
		    y2 = x(2,j)
		    x_vect = x2 - x1
		    y_vect = y2 - y1
		    normal(1,i) = normal(1,i) - y_vect
		    normal(2,i) = normal(2,i) + x_vect
		    normal(1,j) = normal(1,j) - y_vect
		    normal(2,j) = normal(2,j) + x_vect
		    neighbour(i) = neighbour(i) + 1
		    neighbour(j) = neighbour(j) + 1
	end if
	
	! Calculate the gradient of the function f_int
	gradW = current%dwdx
	hf1i =   mass(j)*gradW/rho(j)
	hf1j =  - mass(i)*gradW/rho(i)
        hf2i =   mass(j)*current%dwdy/rho(j)
        hf2j =  - mass(i)*current%dwdy/rho(i)
        grad_f_int(1,i) = grad_f_int(1,i) + (f_int(j) - f_int(i)) * hf1i
        grad_f_int(1,j) = grad_f_int(1,j) + (f_int(i) - f_int(j)) * hf1j
        grad_f_int(2,i) = grad_f_int(2,i) + (f_int(j) - f_int(i)) * hf2i
        grad_f_int(2,j) = grad_f_int(2,j) + (f_int(i) - f_int(j)) * hf2j
        
	end if
	
    current=>current%next
    
end do

do i = 1,ntotal
        
    if (subset(i)==1) then
    
        normal(:,i) = normal(:,i)/neighbour(i)
        norm_vect = ( normal(1,i)**2 + normal(2,i)**2 ) ** 0.5
        normal(:,i) = normal(:,i)/norm_vect 
        p_scal = grad_f_int(1,i)* normal(1,i) + grad_f_int(2,i)* normal(2,i)
    
    if (p_scal < 0) then
        normal(:,i) = -normal(:,i)
    end if
 
    end if
    
enddo       
 
 

end subroutine
!------------------------------------------

        subroutine Initial_conditions

!------------------------------------------

!      ------  Get initial conditions for all components of UNKNO
!              read ic_unkno.....1 for 1D problems
!                                2 read from file
implicit none
 
integer :: icunknoTG  
   
!      ------  Read  

read (dat_file,*) text
write(chk_file,*) text
read (dat_file,*) icunknoTG
write(chk_file,*) icunknoTG

if (icunknoTG==2) then
   call Get_Init2D
elseif (icunknoTG==30) then
   call Get_Init_sigma0
else   
   write(*,*) ' the type of initial condition has not been implemented yet '
   !pause
endif


END subroutine Initial_conditions

!------------------------------------------

        subroutine Get_Init2D

!------------------------------------------

!      ------  Get initial conditions for all components of UNKNO
!              read ic_unkno.....1 for 1D problems
!

implicit none

integer :: ictype 
integer :: ipoin, iamat , test, istre, idimn 
real (kind=DP) :: x0, y0, R0, value, xx, yy, RR, factR, a ,b ,c, xmin, xmax, pi

pi=4*atan(1.0)

DO istre = 1, nstre

   read (dat_file,*) text
   write(chk_file,*) text
   read (dat_file,*)    ICtype			!  ICs for var iamat
   write(chk_file,*)    ICtype 

   if (ICtype==0) then
      stress(istre,:) = 0.0
   end if
   
enddo  

DO idimn = 1, ndimn

   read (dat_file,*) text
   write(chk_file,*) text
   read (dat_file,*)    ICtype			!  ICs for var iamat
   write(chk_file,*)    ICtype 

   if (ICtype==0) then
      vel(idimn,:) = 0.0
   end if
   
enddo 

END Subroutine Get_Init2D

!------------------------------------------

        subroutine Get_Init_sigma0

!------------------------------------------

!      ------   We will read: 
!               sigma0...       nstre components at nodes
!               Internal vars...        Internal_vars(nint_vars, nelemT) 
!               velocities...... Assumed to be zero
!               Internal Vars...
!                               1.... evpstn accumulated deviatoric plastic strain
!                               2.... plastic volumetric strain 
!                               3.... destn(nstre)  plastic increment of strain
!                               4.... size of yield surface, Pc for Cam Clay model

implicit none
 
integer :: type_sigma0, ipoin, ia, ielem, ivar, len1
character(60)  text
real (kind=DP) :: stress_ic(nstre), vars(nint_Vars)

stress        = 0.0
internal_Vars = 0.0
len1 = len_trim(problem_name)

   
!      ------  Read  text, sigma0, internal Vars 

read (dat_file,*) text
write(chk_file,*) text
read (dat_file,*) type_sigma0 
write(chk_file,*) type_sigma0 


if ( type_sigma0 ==1) then        ! constant stress and int vars state
   
   read  (dat_file,*) text
   write (chk_file,*) text
   read  (dat_file,*)    (stress_ic (ia), ia = 1,nstre)
   write (chk_file,*)    (stress_ic (ia), ia = 1,nstre)
   do ipoin = 1, ntotal
      do ia = 1,nstre
         stress(ia,ipoin) = stress_ic(ia)
      enddo
   enddo
      
   read  (dat_file,*) text  
   write (chk_file,*) text
   read  (dat_file,*)    (vars (ia), ia = 1,nint_Vars)
   write (chk_file,*)    (vars (ia), ia = 1,nint_Vars)
   
   do ipoin = 1,ntotal
      do ivar = 1,nint_Vars
         internal_Vars(ivar,ipoin)= vars(ivar)
      enddo
   enddo

end if 

END subroutine  Get_Init_sigma0


! -----------------------------------------------------------------
 
        SUBROUTINE Pint_Update
 
! -----------------------------------------------------------------

!   Subroutine to update interactions between the particles  
!   Stress-point - 1
!   Node - 2
!   Boundary node - 25
 
implicit none

integer :: i, j, ii, jj, isumm, idiff


current=>last

do while (associated(current))

      i = current%pair_i
      j = current%pair_j
      ii    = iabs(itype(i))
      jj    = iabs(itype(j))
      isumm = itype(i) + itype(j)
      idiff = itype(i) - itype(j)
      
      if (ii == 2  .and. jj == 1) then                        !      ------ First check. Valid combinations (1,2);(1,-3);(1,-12);(2,-3),(2,-12))
          current%pair_i = j  !so i is stress point
          current%pair_j = i  !j is boundary
          current%dwdx = -current%dwdx !because w'_ij = - w'_ji
          current%dwdy = -current%dwdy !because w'_ij = - w'_ji
      else if (ii == 2 .and. jj == 25) then !stress-point - no-slip dummy node
          current%pair_i = j  !i is dummy node
          current%pair_j = i  !j is node
          current%dwdx = -current%dwdx   
          current%dwdy = -current%dwdy  
      else if (ii == 1 .and. jj == 25) then !dummy node - node
           current%pair_i = j !i is node
           current%pair_j = i !j is dummy node
           current%dwdx= -current%dwdx  
           current%dwdy = -current%dwdy                  
      end if
     
        
      if (isumm == 3) then    !node - stress pt 
                current%pint_type = 1
      else if (isumm == 2) then  !node - node
                current%pint_type = 2
      else if (isumm == 4) then !stress point - stress point
                current%pint_type = 3
      else if (isumm == 27) then !dummy node - stress point
                current%pint_type = 6     
      else if (isumm == 26) then !dummy node - node
                current%pint_type = 9 
      end if
          
      current=>current%next
      
ENDDO


END subroutine Pint_Update




!-------------------------------------------------------------------

       Subroutine BCs
       
!-------------------------------------------------------------------
implicit none

real (kind=DP) :: ic_time

if     (no_bcs > 0) then   
  
  call Normal_BCs(ic_time)
  
  if (ifsigman == 1) then
        
        call apply_stress_free(stress)
        
  end if

end if

END SUBROUTINE BCs


!-------------------------------------------------------------------

       Subroutine Normal_BCs(ic_time)
       
!-------------------------------------------------------------------

implicit none


integer :: ipoin, iprer, it_curves, j, k, l, ipts_tcurves,start, bc_seg
integer :: it_var, it_prer, nber_BC, bcs_sigman,i,bc_type,ndivxs, ndivys
real (kind=DP) :: t_actual, ic_time, tt0, tt1, xi, unk_pres, bc_value, bc_var
real  (kind=DP) :: a0, a1, w, phi,tt, fact, argum



t_actual = time_sph+ic_time*dt_sph


k = nnode+1
do ipoin = 1,nnode
!k = k+1
nber_BC = bc_info(2,ipoin)				! Nber of BCs on this node
!bcs_sigman = bc_info(8,ipoin)			! =1 if sigmanormal imposed ==> BCs already applied

if (bc_or_not(ipoin) == 2) cycle
if (bc_or_not(ipoin) == 0) cycle

if (bc_or_not(ipoin) == 1) then  
If (nber_BC==0 ) cycle
	DO i=1,5								! First or second BCs
	bc_type   = bc_info(i+2,ipoin)     ! Which Bcs is applied ?
		If (bc_type==0)  cycle
			it_curves = bc_list(3,bc_type)		! Is it a time curve or sinusoidal time function
			If (it_curves==0) then			! BCs with a0 + a1 sin(omega*t+phi)
				a0    = bc_list(5,bc_type)     
				a1    = bc_list(4,bc_type)           
				w     = bc_list(6,bc_type)    
				phi   = bc_list(7,bc_type)
				tt    = bc_list(8,bc_type)
				fact = 1.0-exp(-t_actual/tt)
				argum = w*t_actual - phi
				bc_value = (a0+a1*sin(argum))*fact

			Elseif (it_curves.gt.0) then			! BCs with a1 * f(t)  where f is a time curve
				a1    = bc_list(4,bc_type)           
				do ipts_tcurves = 1, nptstcurves(it_curves)-1
				        tt0 = ttcurves(it_curves,ipts_tcurves)
				        tt1 = ttcurves(it_curves,ipts_tcurves+1)
				if ( t_actual.ge.tt0.AND.t_actual.le.tt1 ) EXIT
					tt0 = -1000.
				enddo
				if (tt0.ge.0.0) then
					xi = (t_actual - tt0) /(tt1-tt0)
					bc_value = (1.-xi) * ftcurves(it_curves , ipts_tcurves)  &
					+ xi  * ftcurves(it_curves , ipts_tcurves+1) 
				else
					bc_value = 0.0
				endif
				bc_value = bc_value * a1    

			Endif

			bc_var = bc_list(2,bc_type)
			
			if (ndimn == 1) then	       

			        if (bc_var == 2) then 
			                vel(1,ipoin) = bc_value 
			                vel(1,ipoin+nnode) = bc_value
			                vel(1,ipoin+2*nnode) = bc_value
			        else if (bc_var == 1) then
			                stress(1,ipoin) = bc_value
			                stress(1,ipoin+nnode) = bc_value
			                stress(1,ipoin+2*nnode) = bc_value
			        end if
			        
		        else if (ndimn == 2) then
		        
		                if (bc_var == 5) then !ux
		                        vel(1,ipoin) = bc_value
		                        do j = 1,3
                                               !vel(1,k+(j-1)) = bc_value
                                        end do		                         		                        
		                else if (bc_var == 6) then !uy
                                        vel(2,ipoin) = bc_value
                                        do j = 1,3
                                              ! vel(2,k+(j-1)) = bc_value
                                        end do	
                                else if (bc_var == 1) then !sigmaxx
                                        stress(1,ipoin) = bc_value
		                else if (bc_var == 3) then !sigmaxy
		                        stress(3,ipoin) = bc_value
		                else if (bc_var == 2) then !sigmayy
		                        stress(2,ipoin) = bc_value		                        
		                end if		                		     		                
		                
		        end if
		                
	ENDDO	
end if	
		                       
	
k = k+3	
ENDDO

l = 2

if (l ==1) then

stress(3,1:ndivy) = 0
vel(1,1:ndivy) = 0

stress(3,nnode-ndivy+1:nnode) = 0
stress(1,nnode-ndivy+1:nnode) = 0

l=1
k=ndivy
do i = 1,ndivx
        vel(:,l) = 0
        vel(1,k) = 0
        if (time_sph < 0.0005) then
                vel(2,k) = 2000*time_sph
        else
                vel(2,k) = 1
        end if
        k=k+ndivy
        l=l+ndivy
end do


end if


END SUBROUTINE Normal_BCs


!------------------------------------------

        subroutine apply_stress_free(fi)

!------------------------------------------

!      ------  Impose tractions on inclined boundaries 

!      ifsigmanT		... 1 if there are inclined boundaries
!                             	       with prescribed tractions
!      nsigmanT			... number of nodes
!      bsigmanT(nsigmanT)	... list of nodes
!      sigmanT (ndimnT, nsigmanT).. 1D  normal stress
!                                   2D  normal and tangential (aclockwise)
!				    3D  normal and 2 tangentials


implicit none
 

integer :: ipoin, iprer, kind_bc, k, j

real (kind=DP) :: sigmatt, sigmann, sigmant, sigman, sigmaxx0, sigmayy0, sigmaxy0
real(kind=DP) :: fi(:,:), tol, sigmaxx, sigmayy, sigmaxy,costh, sinth, sq, s2, c2, sc

tol = dx/10000.

if (ndimn==1) then
   write(*,*) ' we cannot apply yet normal stresses in 1D and 3D '
   write(*,*) ' in 1d is not necessary at all '
   !pause
   stop
elseif (ndimn==2) then
k=1
 DO ipoin=1,nnode
        !kind_bc=bc_or_not(ipoin)
        if (bc_or_not(ipoin) == 2 .and. bc_int(ipoin) /= 1) then       ! BC of SigmaN and TauN
                    if (isnan(normal(1,ipoin)) .or. isnan(normal(2,ipoin))) cycle ! then
                    !            normal(1,ipoin) = 0
                    !            normal(2,ipoin) = 1
                    !end if                                        ! Implemented only SigmaN =0 ==> Not implemented TauN =0
	            !sigman = 0.0
!	            sigman = (time_sph)*(-5.e5)
!	            Endif
!	        Else
!	        sigman=  0.0             ! unpre_TG(4,iprer)                     ! Which is the prescribed value for SigmaN (generally 0)
!	        Endif
		    costh =  normal(1,ipoin)   ! unpre_TG(5,iprer)
                    sinth =  normal(2,ipoin)   ! unpre_TG(6,iprer)
		    s2    = sinth*sinth
		    c2    = costh*costh
		    sc    = sinth*costh
          
	            sigmaxx0 = fi(1,ipoin)                               ! Value of the tension in main reference system
		    sigmayy0 = fi(2,ipoin)
		    sigmaxy0 = fi(3,ipoin)
		        		
		    sigmatt = s2*sigmaxx0 - 2*sc*sigmaxy0 + c2*sigmayy0    ! Value of the Tau in the new reference system
            
    	            sigmaxx =  s2 * sigmatt                            ! Corresponding value of the tension in the main reference system
		    sigmayy =  c2*sigmatt   
		    sigmaxy =  -sc * sigmatt
		    		      

            stress(1,ipoin) = sigmaxx                               ! Give the value to the vector of unknows
            stress(2,ipoin) = sigmayy
            stress(3,ipoin) = sigmaxy
            
            if (nstre == 4) then
		        stress(4,ipoin) = c2*sigmatt
            end if
            
            !do j = 1,3
              ! stress(:,k+nnode+(j-1)) = stress(:,ipoin)
              ! vel(:,k+nnode+(j-1)) = vel(:,ipoin)
              ! bc_or_not(k+nnode+(j-1)) = 2
            !end do
            
        endif
        k=k+3
   ENDDO
   
endif
 
 
END SUBROUTINE  apply_stress_free



! --------------------------------------------------------------------

        subroutine update_strain(dta)

! -------------------------------------------------------------------


implicit none

integer :: ipoin, idimn, i
real (kind=DP) :: dta

do ipoin = nnode+1,ntotal
	internal_vars(1,ipoin) = internal_vars(1,ipoin) + dta * Ddev_strn(ipoin)
	
end do


end subroutine update_strain

!-------------------------------------------------------------------------------

       Subroutine plastic_terms(ie,nd,fi,Gs, der_intvars) 
            
!------------------------------------------------------------------------------

 

!  input is: lmat    ... mat. number
!            ie      ... element (used to get accum.evp)
!            nd      ... ndimn
!            na      ... namat nr of components in Unkno and Unkne
!            Fi      ... Local (node or elem) vector of unknowns (namatT->na)
!			 if_integrate =1 if integrate the sources in this subroutine ortherwise integrate Gs in RK subroutines
!  output is Gs      ... Local Sources (na)
!            debar   ... variation of the evp strain, output for RK subroutines
 
implicit none

integer :: nd, na, ntype_eco,ie, istre, jstre, ip, if_integrate,nintvars, ncrit
real (kind=DP) :: fi(:), Gs(:), stress2(nstre), vivel(nstre), G2(nstre), der_intvars(:)
real (kind=DP) :: Lint_vars(nint_Vars), evpstn, devpstn, YieldH, sigma_out, dt_source
real (kind=DP) :: dens, Young, Poiss, Yield, Hards, gamma, delta, nvp, Dmatx(nstre,nstre)   
 

nintvars = nint_Vars
Gs = 0.0; G2 = 0.0
stress2 = 0.0
vivel = 0.0
ntype_eco = props (1,1)
ncrit = props(1,2)

IF (ntype_eco.gt.1) THEN
	Gs = 0.0
	Poiss = props(1,4)
    Lint_vars(:) = Internal_Vars (:,ie) !  Local Internal vars of the element
	
	Do istre = 1,nstre
		stress2(istre) = fi(istre)
	enddo
	
	
	if (ncrit <= 5) then

	        Call Get_Vivel (ie, stress2, nintvars, Lint_vars, vivel)         ! obtains viscoplastic strain rate
	        Call Get_Dmatx (dmatx)                         ! get De matrix
   
                Do istre = 1,nstre
		        Do jstre = 1,nstre
			        Gs (istre) = Gs(istre) - Dmatx(istre,jstre)*vivel(jstre)
		        Enddo
	        Enddo
	        
	else if (ncrit == 12) then
	
	        call drucker_prager(ie, stress2, vivel, G2)
	        	        
	        Gs = -G2
	        
	end if

	Call Get_derivative_intvars(ie, vivel, der_intvars)
ELSE
	Gs = 0.0
	der_intvars = 0.0
ENDIF


End Subroutine plastic_terms

!------------------------------------------------------------------------------------

subroutine drucker_prager(ie, stress2, vivel, G)

!------------------------------------------------------------------------------------

implicit none

integer :: ncrit, i, j, ie

real(kind=DP) :: G(:), stress2(:), vivel(:), stress0(nstre), f0, df
real :: devia(nstre), dsigma(nstre)
real :: root3, hards, fric, coh, tol
real :: smean, varj2, varj3, steff, sint3, theta, yield
real :: phira, snphi, vari1
real :: evpstn, fdatm0, fdatm, fact
real :: M0,M           ! M0 CSL slope in compression M=M(Th)
real :: p,q, Pc, Py    !  p =-smean  q= root(3*J2)  
                            ! Py size of YS for CamClay
integer :: nintvars
real :: K_mod, G_mod, young, poiss
real :: emean, eps11, eps12, eps22
real ::  s_eps, kc, f1, alpha2, tanfi
real :: G1(4), G2, lambda_1, lambda_2, lambda_3
integer :: delta(4)

adapt_stress = 0

do i = 1,4
        delta(i) = 1
end do

delta(3) = 0


!if (time_sph == 0) then
 !       f0 = 0
!end if

f0 = f_drucker(ie)

tol = 10e-4

!---------------------------
!bui rheology
!---------------------------

  tanfi  = props(1,13)
  coh = props(1,14)
  young = props(1,3)
  poiss = props(1,4)


  smean = (stress2(1)+stress2(2)+stress2(4))/3.0
  
  devia(1) = stress2(1) - smean
  devia(2) = stress2(2) - smean
  devia(3) = stress2(3)
  devia(4) = stress2(4) - smean

  varj2    = devia(3)*devia(3) + 0.5*(devia(1)*devia(1) +       &
                       devia(2)*devia(2) + devia(4)*devia(4) )
                       
  vari1 = 3*smean                     
  !f1 = sqrt(varj2)
  !-------------------------------------------------                     
  eps11 = grad_u(1,1,ie)
  eps12 = 0.5*(grad_u(1,2,ie) + grad_u(2,1,ie))
  eps22 = grad_u(2,2,ie)
  
  emean = eps11 + eps22
  
  alpha2 = tanfi/(sqrt(9 + 12*tanfi**2))
  kc = (3*coh)/(sqrt(9 + 12*(tanfi**2)))
  
  yield = -alpha2*vari1 + kc
  f_drucker(ie) = sqrt(varj2) - yield
  f1 = f_drucker(ie)
  df = f1 - f0
 
  !----------------------------------
  !things needed for plastic flow
  
  G_mod = young/(2.*(1.+Poiss))
  K_mod = young/(3.*(1.-2.*Poiss))
  
  s_eps = devia(1)*eps11 + 2*devia(3)*eps12 + devia(2)*eps22
  
  !if (abs(-yield + f1) < 10e-08 .and. abs(f1) > 10e-09) then .and. df >= 0 
   !if (yield <= f1 .and. abs(f1) > 10e-09 ) then
   !if (abs(-yield + f1) < 10e-08 .and. abs(f1) > 10e-09) then  .and. abs(df) <= 0
   if (time_sph .gt. 0 .and. f1 >= 0  .and. df >= 0 .and. sqrt(varj2) >= 10e-06) then 
        do i = 1,4
                !G1(i) = 3*alpha2*K_mod*delta(i) + (G_mod/f1)*devia(i)
                G1(i) = (G_mod/sqrt(varj2))*devia(i) !non-associated
        end do
        
        lambda_1 = 3*alpha2*K_mod*emean 
        lambda_2 = (G_mod/sqrt(varj2))*s_eps
       !lambda_3 = 9*alpha2*K_mod + G_mod
        lambda_3 = G_mod
        
        G2 = (lambda_1 + lambda_2)/lambda_3
        
        G = G1*G2
        
        !-------------------------
        !plastic strain output
        !-------------------------
        
        
        
        vivel(1) = (1./6.)*(1./sqrt(varj2))*(2*devia(1) - devia(2) - devia(4)) !+ alpha2
        vivel(2) = (1./6.)*(1./sqrt(varj2))*(2*devia(2) - devia(1) - devia(4)) !+ alpha2
        vivel(3) = (1./sqrt(varj2))*devia(3)
        vivel(4) = (1./6.)*(1/sqrt(varj2))*( 2*devia(4)- devia(1) - devia(2)) !+ alpha2
        
        vivel = vivel*G2
         
  else 
  
        G = 0.0
        vivel = 0.0
        
  end if
        

end subroutine 

!-----------------------------------------------

subroutine adapt_stress2
!-----------------------------------------------
!to adapt the stress state for the drucker prager model, if it leaves the yield surface

real(kind=DP) :: smean, kc, rn
real(kind=DP) :: devia(4), varj2, fric, coh
real(kind=DP) :: alpha2, tanfi, yield
integer :: i

  tanfi = props(1,13)
  coh = props(1,14)
         
  alpha2 = tanfi/(sqrt(9 + 12*(tanfi**2)))
  kc = (3*coh)/(sqrt(9 + 12*(tanfi**2)))


  do i = 1,ntotal
        smean = (stress(1,i)+stress(2,i)+stress(4,i))/3.0
  
        devia(1) = stress(1,i) - smean
        devia(2) = stress(2,i) - smean
        devia(3) = stress(3,i)
        devia(4) = stress(4,i) - smean

  varj2    = devia(3)*devia(3) + 0.5*(devia(1)*devia(1) +       &
                       devia(2)*devia(2) + devia(4)*devia(4) )
  
  yield = -alpha2*3*smean + kc                     
  

        if (yield < 0) then 

                stress(1,i) = stress(1,i) - smean + kc/(3*alpha2)

                stress(2,i) = stress(2,i) - smean + kc/(3*alpha2)
        
                stress(4,i) = stress(4,i) - smean + kc/(3*alpha2)
                
                smean = (stress(1,i)+stress(2,i)+stress(4,i))/3.0
  
                devia(1) = stress(1,i) - smean
                devia(2) = stress(2,i) - smean
                devia(3) = stress(3,i)
                devia(4) = stress(4,i) - smean

                varj2    = devia(3)*devia(3) + 0.5*(devia(1)*devia(1) +   &
                       devia(2)*devia(2) + devia(4)*devia(4) )
  
                yield = -alpha2*3*smean + kc 
                
        end if

        if (yield < sqrt(varj2)) then
        
                rn = (-3*alpha2*smean + kc)/(sqrt(varj2))
  
                if (sqrt(varj2) <= 10e-06) then
                        rn = 0
                end if

                stress(1,i) = rn*devia(1) + smean
        
                stress(2,i) = rn*devia(2) + smean
        
                stress(4,i) = rn*devia(4) + smean
        
                stress(3,i) = rn*devia(3)

        end if
        
  end do



end subroutine

! -----------------------------------------------------------------------

   subroutine Get_Vivel (ie, stress2,nintvars, Lintvars, vivel)
                    
!----------------------------------------------------------------------

!      Obtains  (i)   Invariants info: Devia, varj2, sint3, theta, steff and yield          
!               (ii)  Direction of plastic flow dg/ds
!               (iii) vector of evp rate (1..nstre)  (nstre+1) is accum dev evp strain 

!      Input  is lmat           (= mat type), 
!                evpstn         (= acc dev evp strain)  
!                stress(nstre)  (note dim is nstre)

!      Output is vivel (nstre)  (= accum evp. strain )
!				 debar    (= derivative of acc evp. strain) used in RK4 for sources

implicit none

integer ::  nintvars, ie
real (kind=DP) :: Lintvars(:), avect(nstre), smean, varj2, sint3, theta, steff
real (kind=DP) :: stress2 (:), vivel(:), devia(nstre), fdatm0, fdatm, hards 
real (kind=DP) :: frict, root3, yield, devpstn, Young, gamma, delta, nvp, evpstn, sigma_out   





       call invar09 (nintVars, LintVars, stress2, devia, smean, varj2, sint3, theta, steff, & 
                     yield, fdatm)
              
	   if (yield.gt.fdatm) then                         !   only if yielded for localization in vertical slope
                        !   only if yielded   for localization in shear band test
           call yieldf09  ( nintvars, LintVars,                   & 
                            devia, smean, varj2, sint3, theta, steff, avect ) 

           call flowvp09 ( nintvars,LintVars,yield, fdatm, avect,vivel )
    
       else

          vivel   = 0.0
       
	   endif
   

end subroutine Get_Vivel


!----------------------------------------------------------------------

   subroutine invar09 ( nintvars, Lintvars, stress2, devia, smean, varj2, sint3, theta, steff, yield, fdatm )
                    
!----------------------------------------------------------------------


!              Uses  : (i)  stress  stress at point
!                      (ii) ncrit  1..Tresca 2..VMises  3..MCoulomb 4..D.Prager
!                                  5..Mod.Cam Clay
!                      (iii) LintVars (nintvars)... Internal variables at this elem/nodes
!                                    1... accum deviatoric viscoplastic strai
!                                    2... volumetric plastic strain
!                                    3... size of yield surface (Pc for CamClay)

!              Obtains:(a)  Deviatoric stresses DEVIA
!                      (b)  Invariants:  smean, varj2, sint3,theta
!                                        steff 1D equiv. stress
!                      (c)  Yield: Size of F passing through Stress point
!                      (d)  fdatm: Size of current Yield Surface


implicit none 

integer :: ncrit, nintvars

real (kind=DP) :: devia(:), stress2(:), Lintvars(:), root3, hards, phira, snphi
real (kind=DP) :: smean, varj2, varj3, steff, sint3, theta, yield
real (kind=DP) :: evpstn, fdatm0, fdatm, fact, M0, M, p, q, Pc, Py  


root3 = sqrt(3.00)
ncrit = props (1,2)
 

!      ------ Invariants for plane stress conditions


  smean = (stress2(1)+stress2(2)+stress2(4))/3.0

  devia(1) = stress2(1) - smean
  devia(2) = stress2(2) - smean
  devia(3) = stress2(3)
  devia(4) = stress2(4) - smean

  varj2    = devia(3)*devia(3) + 0.5*(devia(1)*devia(1) +       &
                       devia(2)*devia(2) + devia(4)*devia(4) )
  varj3    = devia(4)*(devia(4)*devia(4) - varj2)
  


steff = sqrt(varj2)
if (steff.ne.0.0)  then
    sint3 = -3.0*root3*varj3/(2.0*varj2*steff)
    if(sint3.gt.1.0) sint3=1.0
else
    sint3=0.0
endif
if(sint3.lt.-1.0) sint3=-1.0
if(sint3.gt. 1.0) sint3= 1.0
theta = asin(sint3)/3.0

if     (ncrit==1) then        !      ------  tresca

   yield = 2.0*cos(theta)*steff
   
elseif (ncrit==2) then        !      ------  von mises
 
    yield = root3*steff
   
elseif (ncrit==3) then        !      ------  mohr-coulomb
 
    phira = props(1,9)*0.017453292
    snphi = sin(phira)
    yield = smean*snphi+steff*(cos(theta)-sin(theta)*snphi/root3)
   
elseif (ncrit==4) then        !      ------  drucker-prager
 
    phira = props(1,9)*0.017453292
    snphi = sin(phira)
    yield = 6.0*smean*snphi/(root3*(3.0-snphi))+steff
   
elseif (ncrit==5) then        !      ------  Modified Cam Clay

    M0    = props (1,13)
    snphi = 3.*M0/(6.+ M0)
    M     = 6.*snphi/(3-snphi*sint3)
    p  = - smean
    q  = root3*steff
    if (p.gt.1.e-6) then
        Py = p + (1./p)*(q/M)**2.
    else
        Py = p
    endif
    Yield  = Py
   
endif   


!      ------  Obtain size of yield surface. Before it was computed in flow

if (ncrit.ge.1.and.ncrit.le.4) then     ! Tresca, VM, MC and DP with linear Hard.

   evpstn = LintVars(1)         ! this is accum. deviat. vp strain
   fdatm0 = props(1,7)       ! Initial yield stress 
   hards  = props(1,8)       ! hardening

   if(fdatm0.gt.0.001) then             ! Yield surface should not be smaller than
      fdatm = fdatm0 + hards*evpstn     !     0.1 its initial size. Limit to soft
      fact = abs(fdatm)/abs(fdatm0)
      if(fact.lt.0.1) fact=0.1
      fdatm = fdatm0*fact
   else
      write(*,*) 'with initial yield surface size=',fdatm0
      write(*,*) ' Too small..Hit CR to stop the program '
      !pause
      Stop
   endif

elseif (ncrit==5) then                ! Cam Clay
   fdatm = LintVars(3)                  ! This is Pc. Yield is Py, already obtained
endif
 
End Subroutine invar09


! --------------------------------------------------------------------

  subroutine yieldf09 ( nintvars, Lintvars, devia, smean, varj2, sint3, theta, steff, avect )

!  --------------------------------------------------------------------

!  Obtain plastic vector  df/ds  avect(nstre)

!  ***MP Revision 18th November  for CamClay Model ncrit = 5
!     avoid entering the routine for Tresca and VMises with J2 = 0!!!!

!                props  Material properties nmats * nprop

!                       1..ntype eco 1..elastic 2..EP  3 VP
!                       2..ncrit (1 Tresca 2 VMises 3 MohrC 4 DPrager  5 Cam Clay)
!                       3..Young Modulus
!                       4..Poisson
!                       5..
!                       6..density
!                       7..Yield stress    !
!                       8..Hardening modulus!                     
!                       9..Phi Mohr (degrees)
!                       10..gamma (Perzyna's viscoplasticity)
!                       11..delta (Perzyna)
!                       12..n     (Perzyna)
!                       13..Cam Clay  M  (related to phi)
!                       14..          lambda
!                       15            kappa
!                       16            Pc0 
!                       17            e0 
!
!               LintVars Local Internal Variables
!                        
!                        1..accum. deviatoric plastic strain
!                        2. volum. plastic strain
!                        3. size of yield surface (Pc for instance)

implicit none

integer :: ncrit, istr1, nintvars
real (kind=DP) :: avect(:), smean, varj2, steff, veca1(6), veca2(6), veca3(6) 
real (kind=DP) ::  frict, tanth, tant3, sinth, costh, sint3, cost3, theta  
real (kind=DP) ::  root3, cons1, cons2, cons3, abthe, plumi, snphi
real (kind=DP) :: devia(:), LintVars(:), M0, dfdM, dMdth, dfdth, M, p, Pc

avect = 0.0
veca1 = 0.0
veca2 = 0.0
veca3 = 0.0
ncrit  = props (1,2)
frict  = props (1,9)        
tanth  = tan(theta)
tant3  = tan(3.0*theta)
sinth  = sin(theta)
costh  = cos(theta)
cost3  = cos(3.0*theta)
root3  = sqrt(3.0)

!      ------  2D situations


    veca1(1)=1.0             !      ------  Get vector a1=d(Sm)/ds
    veca1(2)=1.0
    veca1(3)=0.0
    veca1(4)=1.0

	If (steff.gt.0) then
	   do istr1=1,nstre         !      ------  Get vector a2=d(Steff)/ds
          veca2(istr1)=devia(istr1)/(2.0*steff)
       enddo
       veca2(3)=devia(3)/steff
	Else
	   veca2 = 0.0
	Endif       
                                !      ------  Get vector a3
    veca3(1)=devia(2)*devia(4)+varj2/3.0 
    veca3(2)=devia(1)*devia(4)+varj2/3.0
    veca3(3)=-2.0*devia(3)*devia(4)
    veca3(4)=devia(1)*devia(2)-devia(3)*devia(3)+varj2/3.0



IF (ncrit==1 ) THEN        !      ------  tresca

   cons1  =  0.0
   abthe  =  abs(theta*57.29577951308)
   if (abthe.ge.29.0) then
      cons2 = root3
      cons3 = 0.0
   else
      cons2 = 2.0*(costh+sinth*tant3)
      cons3 = root3*sinth/(varj2*cost3)
   endif

ELSEIF (ncrit==2) THEN    !      ------  von mises

    cons1 = 0.0
    cons2 = root3
    cons3 = 0.0

ELSEIF (ncrit==3) THEN     !      ------  mohr-coulomb

    cons1 = sin(frict*0.017453292)/3.0
    abthe = abs(theta*57.29577951308)
    if(abthe.ge.29.0) then
       cons3 = 0.0
       plumi = 1.0
       if(theta.gt.0.0) plumi=-1.0
       cons2 = 0.5*(root3+plumi*cons1*root3)
    else
       cons2 = costh*((1.0+tanth*tant3)+cons1*(tant3-tanth)*root3)
       cons3 = (root3*sinth+3.0*cons1*costh)/(2.0*varj2*cost3)
    endif

ELSEIF (ncrit==4) THEN     !      ------  drucker-prager

    snphi = sin(frict*0.017453292)
    cons1 = 2.0*snphi/(root3*(3.0-snphi))
    cons2 = 1.0
    cons3 = 0.0

ELSEIF (ncrit==5) THEN     !      ------  Cam Clay

    Pc    = LintVars(3)
    p     = - smean 
    M0    = props (1,13)
    snphi = 3.*M0/(6.+ M0)
    M     = 6.*snphi/(3-snphi*sint3)
    dMdTh = 18*snphi*snphi*cost3/((3-snphi*sint3)*(3-snphi*sint3))
    dfdM  = (2.*M/Pc)*p*(p-pc)
    dfdTh = dfdM * dmdth
   
    cons1 = -(1./3.)*(2.*p - Pc)*M*M/Pc
    cons2 = 6.*steff/Pc - (tant3*dfdTh/steff)
    cons3 = -(root3/(2.*cost3*steff**3.))* dfdth

ELSE
   write(*,*) ' ncrit wrong in yieldf routine !pauseD '
   !pause

ENDIF

DO istr1 = 1,nstre
   avect(istr1) = cons1*veca1(istr1) + cons2*veca2(istr1) + cons3*veca3(istr1)
ENDDO

end subroutine Yieldf09

! -------------------------------------------------------------

   subroutine flowvp09 ( nintvars, Lintvars, yield, fdatm, avect, vivel)

! -----------------------------------------------------------------
 
!      ------  this subroutine evaluates the viscoplastic strain rate
!              and stores it in vivel(nstre)
!
!              Update internal Vars done in another routine using vivel
 
implicit none


real (kind=DP) :: vivel(:), avect(:), LintVars(:), yield, evpstn, fdatm, fdatm0
real (kind=DP) :: allow, hards, frict, gamma, delta, nflow, root3, fact, fcurr, fnorm, cmult      
integer :: istr1, ncrit, nintvars
 
allow = 0.01
vivel = 0.0

ncrit  = props(1,2)
fdatm0 = props(1,7)      !   this is the reference yield stress
hards  = props(1,8)
frict  = props(1,9)
gamma  = props(1,10)
delta  = props(1,11)
nflow  = props(1,12)
root3  = 1.73205080757
frict  = frict*0.017453292

if(ncrit==3) fdatm0 = fdatm0*cos(frict)

if(ncrit==4) then
	fdatm0 = 6.0*fdatm0*cos(frict)/(root3*(3.0-sin(frict)))
endif

fcurr = yield - fdatm
fnorm = fcurr/fdatm

if(fnorm.ge.allow) then
	if (nflow.ne.1) then
		cmult = gamma*(exp(delta*fnorm)-1.0)
	else
		cmult = gamma*(fnorm**delta)
	endif
	do istr1 = 1, nstre
		vivel(istr1) = cmult*avect(istr1)      ! store in avect evp velocity
	enddo
endif

END subroutine  flowvp09


!----------------------------------------------------------------------

      subroutine Get_Dmatx (dmatx)

!-----------------------------------------------------------------------

!      ------  Obtain the constitutive elastic matrix (Isotropic)

implicit none

real (kind=DP) :: dmatx(nstre,nstre), Young, Poiss, const, conss, alfa, beta, G 

young = props(1,3)
poiss = props(1,4)

Dmatx = 0.0                                     !      ------  Initializes De


     const      = young*(1.0-poiss)/((1.0+poiss)*(1.0-2.0*poiss))
     dmatx(1,1) = const
     dmatx(2,2) = const
     dmatx(1,2) = const*poiss/(1.0-poiss)
     dmatx(2,1) = const*poiss/(1.0-poiss)
     dmatx(3,3) = (1.0-2.0*poiss)*const/(2.0*(1.0-poiss)) 
     dmatx(1,4) = const*poiss/(1.0-poiss)
     dmatx(2,4) = const*poiss/(1.0-poiss)
     dmatx(4,1) = const*poiss/(1.0-poiss)
     dmatx(4,2) = const*poiss/(1.0-poiss)
     dmatx(4,4) = const
        

 

END subroutine Get_Dmatx



! -----------------------------------------------------------------------

   subroutine Get_derivative_intvars(ie, vivel, der_intvars)
                    
!----------------------------------------------------------------------



!      Calculate derivative of Int_Variables (NintVars, ie)   
!
!                          1..acc. dev.vp strain
!                          2..volum. plastic strain
!                          3..Size of Yield Surface Pc   
!
!      Input  is 
!                lmat...          (= mat type)
!				 ie....    element/node 
!                vivel(nstre) ...  vp strain rate
!
!	   Output is
!				 der_intvars: derivative of acc. dev. vp strain
!							  derivative of volum. plastic strain
!							  derivative of the size of yield surface Pc 		
!       from Main vars:

!                props  Material properties nmats * nprop

!                       1..ntype eco 1..elastic 2..EP  3 VP
!                       2..ncrit (1 Tresca 2 VMises 3 MohrC 4 DPrager  5 Cam Clay)
!                       3..Young Modulus
!                       4..Poisson
!                       5..
!                       6..density
!                       7..Yield stress    !
!                       8..Hardening modulus!                     
!                       9..Phi Mohr (degrees)
!                       10..gamma (Perzyna's viscoplasticity)
!                       11..delta (Perzyna)
!                       12..n     (Perzyna)
!                       13..Cam Clay  M  (related to phi)
!                       14..          lambda
!                       15            kappa
!                       16            Pc0 
!                       17            e0  

implicit none

integer :: ie, ncrit
real (kind=DP) :: vivel(:), der_intvars(:), ddev_strn, dvol_strn
real (kind=DP) :: Pc, e0, lambda, kappa, dPc

der_intvars = 0.0

ncrit = props (1,2)


if (ncrit==5) then
   Pc     = internal_vars (3,ie)
   lambda = props (1,14)
   kappa  = props (1,15)
   e0     = props (1,17)
endif

!      ------  Obtain deviatoric vp strain


    ddev_strn = sqrt((2.0*(vivel(1)*vivel(1) + vivel(2)*vivel(2)       &
                         + vivel(4)*vivel(4))+ vivel(3)*vivel(3))/3.0)
    dvol_strn = vivel(1) + vivel(2) + vivel(4)
    
    der_intvars(1) = ddev_strn
    der_intvars(2) = dvol_strn

    if (ncrit==5) then       
       dPc = - ((1+e0)/(lambda-kappa))*Pc*dvol_strn
	   der_intvars(3) = dPc
    endif  
    

        
End subroutine Get_derivative_intvars  


!----------------------------------------------------------------------------------
 
    Subroutine gravity_force( nd, fi,Gs) 
             
!----------------------------------------------------------------------------------

!  INTEGRATE Gravity sources and damping

!  Input is: nd      ...        ndimnT
!            Fi      ...        Local (node or elem) vector of unknowns (namatT->na)

!  output is Gs      ...        Local Sources (na)
 
!      ------  Obtain fluxes Gs for TG algorithm in solid problems

implicit none

integer :: nd, ip, idimn, nd_aux, it_curves, ipts_tcurves  
real(kind=DP) :: fi(:,:), Gs(:,:), tt0, tt1, xi, factg            !  

Gs = 0.0
!      ------  Gravity Loading 

if (ic_grav.ne.0) then 

   it_curves = tcurve_grav

   do ipts_tcurves = 1, nptstcurves(it_curves)-1
      tt0 = ttcurves(it_curves,ipts_tcurves)
      tt1 = ttcurves(it_curves,ipts_tcurves+1)
      if ( time_sph.ge.tt0.AND.time_sph.le.tt1 ) EXIT
      tt0 = -1000.
   enddo

   if (tt0.ge.0.0) then
      xi    = (time_sph - tt0) /(tt1-tt0)
      factg = (1.-xi) * ftcurves(it_curves , ipts_tcurves)    &
                + xi  * ftcurves(it_curves , ipts_tcurves + 1) 
   else
      factg = 0.0
   endif  
   factg = factg * ft_grav

   if (ic_grav==1) then                      !  Apply gravity forced dv/dt = g + ...
   
          do idimn = 1,nd
!             Gs (nstre + idimn) = Gs (nstre + idimn) + factg * cgrav (idimn)
             Gs (idimn,:)  =  factg * cgrav (idimn)
          enddo
    endif

endif

!      ------ Damping

nd_aux = nd
!if (ntype_solid==12) nd_aux = 2


do idimn = 1,nd_aux
  Gs (idimn,:) = Gs (idimn,:) - DampingTG* fi(idimn,:)
enddo 


End Subroutine gravity_force


!-------------------------------------------------------------------

       Subroutine OutputMesh
       
!-------------------------------------------------------------------

implicit none

integer :: idimn, ielem, ipoin, inode, ndummy3 


write (gid_msh,*) 'MESH    dimension 3 ElemType Quadrilateral  Nnode 4 '

write (gid_msh,*) ' coordinates'


do ipoin = 1,nnode!nnode+1, ntotal - ndummy2
   write(gid_msh,*) ipoin,(x(idimn,ipoin),idimn=1,ndimn), 0.0
enddo

write(gid_msh,*) ' End coordinates' 

write (gid_msh,*) ' Elements'

do ielem=1,nelem
   write(gid_msh,*) ielem, (element_coord(inode,ielem),inode=1,4),' 1'
enddo  

write (gid_msh,*) ' End elements'   
      
close(gid_msh)

write (gid_res,*) 'GiD Post Results File 1.0'
write (gid_res,*) 'GaussPoints "Group1" ElemType  Quadrilateral'
write (gid_res,*) 'Number Of Gauss Points: 1'
write (gid_res,*) 'Natural Coordinates: Internal'
write (gid_res,*) 'End GaussPoints'


End Subroutine OutputMesh


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

       SUBROUTINE OutputRes 
 
!outputs results files for both ParaView (as .csv files) and GiD (as post.res file)     
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

implicit none

integer :: i 
real (kind=DP) ::  ux, uy, umod, dispx, dispy, dispz, disp(ndimn)
character(len=10) :: dt_file, formt

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!ParaView output
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   
node_file = nnode + itimestep_sph 
stress_point_file = nstress + itimestep_sph 
free_surf_file = ntotal + itimestep_sph

formt = '(I6.6)'
write(dt_file,formt) itimestep_sph

!output specified information on nodes 
!-------------------------------------------------------------------------------

open(unit=node_file,file='nodes'//'.csv.'//trim(dt_file))

write(node_file,'(A)', advance = 'no') 'x-coord,'
write(node_file,'(A)', advance = 'no') 'y-coord,'!, &
if (vel_out(1) == 1) write(node_file,'(A)', advance = 'no') 'x-vel,' 
if (vel_out(2) == 1) write(node_file,'(A)', advance = 'no') 'y-vel,' 
if (stress_out(1) == 1) write(node_file,'(A)', advance = 'no') 'sxx,' 
if (stress_out(2) == 1) write(node_file,'(A)', advance = 'no') 'syy,' 
if (stress_out(3) == 1) write(node_file,'(A)', advance = 'no') 'sxy,' 
if (stress_out(4) == 1) write(node_file,'(A)', advance = 'no') 'szz,' 
if (strain_out == 1) write(node_file,'(A)', advance = 'no') 'strain,' 
if (disp_out == 1) write(node_file,'(A)', advance = 'no') 'disp,' 
if (rho_out == 1) write(node_file,'(A)', advance = 'no') 'density,' 
if (sml_out == 1) write(node_file,'(A)', advance = 'no') 'sml,' 
write(node_file,*) ' '
do i = 1,nnode
        write(node_file,'(F16.8)', advance = 'no') x(1,i)
        write(node_file,'(A)',advance = 'no') ','
        write(node_file,'(F16.8)', advance = 'no') x(2,i)
        write(node_file,'(A)',advance = 'no') ','
        if (vel_out(1) == 1) then
                write(node_file,'(F16.8)', advance = 'no') vel(1,i)
                write(node_file,'(A)',advance = 'no') ','
        end if
        if (vel_out(2) == 1) then
                write(node_file,'(F16.8)', advance = 'no') vel(2,i)
                write(node_file,'(A)',advance = 'no') ','
        end if
        if (stress_out(1) == 1) then
                write(node_file,'(F16.8)', advance = 'no') stress(1,i)
                write(node_file,'(A)',advance = 'no') ','
        end if
        if (stress_out(2) == 1) then
                write(node_file,'(F16.8)', advance = 'no') stress(2,i)
                write(node_file,'(A)',advance = 'no') ','
        end if
        if (stress_out(3) == 1) then
                write(node_file,'(F16.8)', advance = 'no') stress(3,i)
                write(node_file,'(A)',advance = 'no') ','
        end if
        if (stress_out(4) == 1) then
                write(node_file,'(F16.8)', advance = 'no') stress(4,i)
                write(node_file,'(A)',advance = 'no') ','
        end if
        if (strain_out == 1) then
                write(node_file,'(F16.8)', advance = 'no') internal_vars(1,i)
                write(node_file,'(A)',advance = 'no') ','
        end if
        if (disp_out == 1) then
                write(node_file,'(F16.8)', advance = 'no') sqrt(displ(1,i)**2 + displ(2,i)**2) 
                write(node_file,'(A)',advance = 'no') ','
        end if
        if (rho_out == 1) then
                write(node_file,'(F16.8)', advance = 'no') rho(i)
                write(node_file,'(A)',advance = 'no') ','
        end if
        if (sml_out == 1) then
                write(node_file,'(F16.8)', advance = 'no') hsml(i)
                write(node_file,'(A)',advance = 'no') ','
        end if
        write(node_file,*) ' '
end do

 close(unit=node_file)
 
!output specified information on stress-points
!-------------------------------------------------------------------------------

open(unit=stress_point_file,file='stress_points'//'.csv.'//trim(dt_file))

write(stress_point_file,'(A)', advance = 'no') 'x-coord,'
write(stress_point_file,'(A)', advance = 'no') 'y-coord,'!, &
if (vel_out(1) == 1) write(stress_point_file,'(A)', advance = 'no') 'x-vel,' 
if (vel_out(2) == 1) write(stress_point_file,'(A)', advance = 'no') 'y-vel,' 
if (stress_out(1) == 1) write(stress_point_file,'(A)', advance = 'no') 'sxx,' 
if (stress_out(2) == 1) write(stress_point_file,'(A)', advance = 'no') 'syy,' 
if (stress_out(3) == 1) write(stress_point_file,'(A)', advance = 'no') 'sxy,' 
if (stress_out(4) == 1) write(stress_point_file,'(A)', advance = 'no') 'szz,' 
if (strain_out == 1) write(stress_point_file,'(A)', advance = 'no') 'strain,' 
if (rho_out == 1) write(stress_point_file,'(A)', advance = 'no') 'density,' 
if (sml_out == 1) write(stress_point_file,'(A)', advance = 'no') 'sml,' 
write(stress_point_file,*) ' '
do i = nnode+1,ntotal
        write(stress_point_file,'(F16.8)', advance = 'no') x(1,i)
        write(stress_point_file,'(A)',advance = 'no') ','
        write(stress_point_file,'(F16.8)', advance = 'no') x(2,i)
        write(stress_point_file,'(A)',advance = 'no') ','
        if (vel_out(1) == 1) then
                write(stress_point_file,'(F16.8)', advance = 'no') vel(1,i)
                write(stress_point_file,'(A)',advance = 'no') ','
        end if
        if (vel_out(2) == 1) then
                write(stress_point_file,'(F16.8)', advance = 'no') vel(2,i)
                write(stress_point_file,'(A)',advance = 'no') ','
        end if
        if (stress_out(1) == 1) then
                write(stress_point_file,'(F16.8)', advance = 'no') stress(1,i)
                write(stress_point_file,'(A)',advance = 'no') ','
        end if
        if (stress_out(2) == 1) then
                write(stress_point_file,'(F16.8)', advance = 'no') stress(2,i)
                write(stress_point_file,'(A)',advance = 'no') ','
        end if
        if (stress_out(3) == 1) then
                write(stress_point_file,'(F16.8)', advance = 'no') stress(3,i)
                write(stress_point_file,'(A)',advance = 'no') ','
        end if
        if (stress_out(4) == 1) then
                write(stress_point_file,'(F16.8)', advance = 'no') stress(4,i)
                write(stress_point_file,'(A)',advance = 'no') ','
        end if
        if (strain_out == 1) then
                write(stress_point_file,'(F16.8)', advance = 'no') internal_vars(1,i)
                write(stress_point_file,'(A)',advance = 'no') ','
        end if       
        if (rho_out == 1) then
                write(stress_point_file,'(F16.8)', advance = 'no') rho(i)
                write(stress_point_file,'(A)',advance = 'no') ','
        end if
        if (sml_out == 1) then
                write(stress_point_file,'(F16.8)', advance = 'no') hsml(i)
                write(stress_point_file,'(A)',advance = 'no') ','
        end if
        write(stress_point_file,*) ' '
end do

 close(unit=stress_point_file)



!output nodes that lie on the free surface
!-------------------------------------------------------------------------------

open(unit=free_surf_file,file='surface_points'//'.csv.'//trim(dt_file))

write(free_surf_file,*) 'x-coord,','y-coord,'
do i = 1,nnode
        if (bc_or_not(i) == 2) then
                 write(free_surf_file,*) x(1,i), ',',x(2,i)
        end if
end do
 
 close(unit=free_surf_file)

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!GiD output
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!displacement always written to file to allow results deformation visualisation
!-------------------------------------------------------------------------------

write(gid_res,*) 'Result "disp" "disp" ',time_sph,' Vector OnNodes "where" '
write(gid_res,*) ' Values '
do i = 1,nnode
        if (update_x) then
                disp(:) = x(:,i) - x00(:,i)
                dispx   = disp(1)
                if (ndimn==1) then
                        dispy = 0.0
                else
                        dispy = disp(2)
                end if
	else
                dispx = displ(1,i)
	        dispy = displ(2,i)
	end if
        write(gid_res,*) i, dispx, dispy, '0.'
end do
write(gid_res,*) ' End Values '

!write velocity vector to GiD file (always included by default)
!-------------------------------------------------------------------------------

write(gid_res,*) 'Result "vel" "veloc" ',time_sph,' Vector OnNodes "where" '
write(gid_res,*) ' Values '
do i = 1, nnode
        ux = vel(1,i) 
        uy = vel(2,i) 
        umod = sqrt(vel(1,i)**2 + vel(2,i)**2)
        write(gid_res,*)  i,ux,uy,umod
end do
write(gid_res,*) ' End Values ' 

!write sigma_xx to GiD file
!-------------------------------------------------------------------------------

if (stress_out(1) == 1) then  
        write(gid_res,*) 'Result "sigmaxx" "sxx" ',time_sph,' Vector OnNodes "where" '
        write(gid_res,*) ' Values '
        do i = 1,nnode
                write(gid_res,*) i, stress(1,i), ' 0.   0.   '
        end do
        write(gid_res,*) ' End Values ' 
end if

!write sigma_yy to GiD file
!-------------------------------------------------------------------------------

if (stress_out(2) == 1) then
        write(gid_res,*) 'Result "sigmayy" "syy" ',time_sph,' Vector OnNodes "where" '
        write(gid_res,*) ' Values '
        do i = 1, nnode
                write(gid_res,*) i, stress(2,i),' 0.   0.   '
        end do
        write(gid_res,*) ' End Values ' 
end if

!write sigma_xy to GiD file
!-------------------------------------------------------------------------------

if (stress_out(3) == 1) then
         write(gid_res,*) 'Result "sigmaxy" "sxy" ',time_sph,' Vector OnNodes "where" '
         write(gid_res,*) ' Values '
         do i = 1,nnode
                 write(gid_res,*) i, stress(3,i),' 0.   0.   '
         end do
         write(gid_res,*) ' End Values ' 
end if

!write strain to GiD file
!-------------------------------------------------------------------------------

if (strain_out == 1) then
        write(gid_res,*) 'Result "plastic strain" "esp" ',time_sph,' Vector OnNodes "where" '
	write(gid_res,*) ' Values  ' 
	do i = 1,nnode
	         write(gid_res,*) i, internal_vars(1,i),' 0.  0.   '
	enddo 
	write(gid_res,*) ' End Values '
end if
 
        END SUBROUTINE OutputRes


END MODULE SPH_material_2018

  
