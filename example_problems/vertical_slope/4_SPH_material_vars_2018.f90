
!-------------------------------------------------------------

MODULE SPH_material_vars_2018

!------------------------------------------------------------

use variable_types

implicit none

real (kind=DP), allocatable:: x_1(:,:), rho_1(:), mass_1(:), hsml_1(:), vel_1(:,:), stress_1(:,:)
integer, allocatable:: itype_1(:), matno_1(:)

real (kind=DP), allocatable:: x_s  (:,:), rho_s(:), mass_s(:), hsml_s(:), stress_s(:,:), vel_s(:,:)
integer, allocatable:: itype_s(:), matno_s(:)

real (kind=DP), allocatable::  x_dummy(:,:), rho_dummy(:), mass_dummy(:), hsml_dummy(:), stress_dummy(:,:)
integer, allocatable::  itype_dummy(:), vel_dummy(:,:)

real (kind=DP), allocatable :: x_out(:,:)
	                                   
!   Variable for BCs

integer :: no_bc_vars                        !  number of prescribed vars
real (kind=DP), allocatable:: bc_list    (:,:)   !  info presc. values                   numpre_TG * nprerTG(=5)
integer :: no_segments_bc
real (kind=DP), allocatable:: Segment_BCs (:,:)   !  info on the segment where BCs are applied
integer :: no_nodes_bc
real (kind=DP), allocatable:: Nodal_BCs   (:,:)   !  info on the prescribed BC on nodes

integer :: ic_grav, tcurve_grav
real (kind=DP), allocatable :: cgrav(:)
real (kind=DP) :: ft_grav, dampingTG

integer :: ic_tcurve

END MODULE SPH_material_vars_2018
