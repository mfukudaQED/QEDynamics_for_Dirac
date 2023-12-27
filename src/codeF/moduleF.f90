!===================================================================
! Last Change: 12-Oct-2015.
!
! 2014.06.07
! All modules are summerized in this file.

! - module param_AM
! - module prop_param
! - module prop_tmp
! - module AA_nucatt_pg
! - module DIRRCI_coef
! - module PTcoef
! - module param_KRCI
!===================================================================
module param_AM
!===================================================================
  use Precision
  implicit none

!  integer :: direcB ! 1->x, 2->y, 3->z
!  complex(kind=dp) :: Bmag
  real(kind=dp),save :: Bmag(3)
!  integer :: direcB2 ! 1->x, 2->y, 3->z
!  real(kind=dp) :: Bmag2
  real(kind=dp),save :: Eelec(3)
  real(kind=dp),save :: magdip(3)
  real(kind=dp),save :: omega_A0M

end module param_AM

!===================================================================
module prop_param
!===================================================================
  use Precision
  implicit none

  integer,save :: meshx, meshy, meshz, mesh, meshyz
  real(kind=dp),save :: xori, xend, yori, yend, zori, zend
  real(kind=dp),save :: dx, dy, dz
  real(kind=dp),save :: x0, y0, z0
  integer,allocatable :: unitn(:)
  character(LEN=80),allocatable :: filename(:)

end module prop_param

!===================================================================
module prop_tmp
!===================================================================
  use Precision
  implicit none

  complex(kind=dp),allocatable :: sum0d(:,:,:), sum0c(:,:,:)
  complex(kind=dp),allocatable :: sum0s(:,:,:,:), sum0t(:,:,:,:), sum0z(:,:,:,:)
  complex(kind=dp),allocatable :: sum0j(:,:,:,:)
  complex(kind=dp),allocatable :: sum0tA(:,:,:,:), sum0A(:,:,:,:)
  complex(kind=dp),allocatable :: sum0tau(:,:,:,:,:)
  complex(kind=dp),allocatable :: sum0divj(:,:,:)
  complex(kind=dp),allocatable :: sum0tauA(:,:,:,:,:)
  complex(kind=dp),allocatable :: sum0tauAM(:,:,:,:,:)
  complex(kind=dp),allocatable :: sum0t_AM(:,:,:,:)
  complex(kind=dp),allocatable :: sum0pi(:,:,:,:),sum0orbang(:,:,:,:),sum0rots(:,:,:,:)
  real(kind=dp),allocatable :: sum0pol(:,:,:,:)
  real(kind=dp),allocatable :: sum0ef(:,:,:,:)

  complex(kind=dp),allocatable :: tmp_d(:,:,:,:,:), tmp_c(:,:,:,:,:)
  complex(kind=dp),allocatable :: tmp_s(:,:,:,:,:,:), tmp_t(:,:,:,:,:,:), tmp_z(:,:,:,:,:,:)
  complex(kind=dp),allocatable :: tmp_j(:,:,:,:,:,:)
  complex(kind=dp),allocatable :: tmp2_j_Qmat(:,:,:)
  complex(kind=dp),allocatable :: tmp_t_Arad_Qmat(:,:,:)
  complex(kind=dp),allocatable :: tmp_tau(:,:,:,:,:,:,:)
  complex(kind=dp),allocatable :: tmp_divj(:,:,:,:,:)
  complex(kind=dp),allocatable :: tmp_E_Qmat(:,:,:,:,:,:)
  complex(kind=dp),allocatable :: tmp_tauAM(:,:,:,:,:,:,:)
  complex(kind=dp),allocatable :: tmp_t_AM(:,:,:,:,:,:)
  complex(kind=dp),allocatable :: tmp_pi(:,:,:,:,:,:),tmp_orbang(:,:,:,:,:,:),tmp_rots(:,:,:,:,:,:)
end module prop_tmp

!===================================================================
module AA_nucatt_pg
!===================================================================
  use Precision
  real(kind=dp),allocatable :: nucatt_pg(:,:,:,:)
  complex(kind=dp),allocatable :: calEnpqm(:,:,:,:)

end module AA_nucatt_pg

!===================================================================
module DIRRCI_coef
!===================================================================
  use Precision

  !--- electron CI coefficients for natural orbital ---
  complex(kind=dp),allocatable,save :: c_natu(:,:,:) ! complex coefficients (4, NBS_L or NBS_S,2*NBS)
  real(kind=dp),allocatable,save :: occup_natu(:) ! occupation of natural orb (norb)
  real(kind=dp),save :: thresh_occdet 

end module DIRRCI_coef

!===================================================================
module PTcoef
!===================================================================
  use Precision
  implicit none
  
!!$  !--- electron solution coefficients ---
!!$  complex(kind=dp), allocatable :: c_Lapt(:,:), c_Lbpt(:,:), c_Sapt(:,:), c_Sbpt(:,:) ! complex coefficients (NBS_L or NBS_S,2*NBS)
!!$
!!$  !--- positron solution coefficients ---
!!$  complex(kind=dp), allocatable :: d_Lapt(:,:), d_Lbpt(:,:), d_Sapt(:,:), d_Sbpt(:,:) ! complex coefficients (NBS_L or NBS_S,2*NBS)

  !--- electron solution coefficients ---
  complex(kind=dp), allocatable :: c_psi(:,:,:) ! complex coefficients (4, NBS_L or NBS_S,2*NBS)
  !--- electron energy ---
  real(kind=dp), allocatable :: eneQ0(:),eneQ1(:),eneQ2(:) !(2*NBS)
  integer :: NEL
  double precision,allocatable :: occup(:)
!  integer :: direcB ! 1->x, 2->y, 3->z
!  complex(kind=dp) :: Bmag
!  complex(kind=dp) :: Bmag(3)
!  integer :: direcB2 ! 1->x, 2->y, 3->z
!  complex(kind=dp) :: Bmag2

end module PTcoef

!===================================================================
module param_KRCI
!===================================================================
  implicit none

  integer :: NBS_Eg, NBS_Pg, NBS_Eu, NBS_Pu ! number of molecular Electron and Positron orbitals

end module param_KRCI
