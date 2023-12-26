!===================================================================
module Precision
!===================================================================
  implicit none

!  integer, parameter :: qp = kind(1.q0)
  integer, parameter :: dp = kind(1.d0)
!  integer, parameter :: sp = kind(1.0)

end module Precision

!===================================================================
module DiracOutput
!===================================================================
  use Precision
  implicit none
  
  character(len=200),save :: FILESFOLDER
  character(len=120),save :: FIFI1,FIFI2,FIFI3,FIFI4,FIFI5
  integer,save :: SYMM
  integer,save :: KBAL !KBAL->1  ,URKBAL->0
  integer,save :: NROOT
  real(kind=dp),save :: U1,U2,U3,U4,U5,U6,UU

  integer,save :: NBS  ! total number of molecular orbitals
  integer,save :: NBS_E, NBS_P !number of molecular Electron and Positron orbitals
  integer,save :: NBS_TOTAL !NBS_E + NBS_P
  integer,save :: NBS00  ! total number of primitive gaussians (NBS_L+NBS_S)
  integer,save :: NBS0  !  number of primitive gaussians (max(NBS_L,NBS_S))
  integer,save :: NBS_LL,NBS_SS,NBS_LS,NBS_LT,NBS_ST,NBS_00

  integer,save :: NAT ! number of atoms
  real(kind=dp),allocatable,save :: xc(:),yc(:),zc(:) ! position of atom
  integer,allocatable,save :: cn(:) ! charge of atom

  !--- primitive gaussian ---
  integer,save :: NBS_L, NBS_S  ! number of primitive gaussians for large and small components
  real(kind=dp),allocatable,save :: aa_L(:),xx_L(:),yy_L(:),zz_L(:) ! position of primitive gaussian
  integer,allocatable,save :: nx_L(:),ny_L(:),nz_L(:) ! angular momentum of primitive gaussian
  real(kind=dp),allocatable,save :: aa_S(:),xx_S(:),yy_S(:),zz_S(:) ! position of primitive gaussian
  integer,allocatable,save :: nx_S(:),ny_S(:),nz_S(:) ! angular momentum of primitive gaussian
  real(kind=dp),allocatable,save :: c_Lnorm(:), c_Snorm(:) ! normalization constant of primitive gaussian
  integer,allocatable,save :: atnum_L(:),atnum_S(:) ! atom number to p.g. number
  integer,save :: nmax(3) ! max(nk_L+2,nk_S+2) k = x,y,z

  !--- electron solution coefficients ---
  complex(kind=dp),allocatable,save :: c_La(:,:), c_Lb(:,:), c_Sa(:,:), c_Sb(:,:) ! complex coefficients (NBS_L or NBS_S,NBS)

  !--- positron solution coefficients ---
  complex(kind=dp),allocatable,save :: d_La(:,:), d_Lb(:,:), d_Sa(:,:), d_Sb(:,:) ! complex coefficients (NBS_L or NBS_S,NBS)

  !---orbital energy---
  !double precision,allocatable,save :: e_eig(:),p_eig(:)
  real(kind=dp),allocatable,save :: e_eig(:),p_eig(:)

end module DiracOutput

!===================================================================
module DefineTypes
!===================================================================
  use Precision
  implicit none

  integer,parameter :: NMAX_PG = 3000  ! Assume NMAX_PG primitive gaussian functions at maximum for both large and small components.
  
  type primitive_gaussian
     ! Large (1,2) components
     real(kind=dp) :: aL(NMAX_PG)                  ! exponents 
     real(kind=dp) :: xL(NMAX_PG),yL(NMAX_PG),zL(NMAX_PG)      ! positions of p.g.
     integer      :: nxL(NMAX_PG),nyL(NMAX_PG),nzL(NMAX_PG)   ! angular momentum of p.g.
     ! Small (3,4) components
     real(kind=dp) :: aS(NMAX_PG)
     real(kind=dp) :: xS(NMAX_PG),yS(NMAX_PG),zS(NMAX_PG)
     integer      :: nxS(NMAX_PG),nyS(NMAX_PG),nzS(NMAX_PG)
  end type primitive_gaussian

  type nonzeromat4legs
     integer :: a,b,c,d ! specify non-zero components
     ! assume that there 4 spinor legs a,b,c,d
     complex(kind=dp) :: val
  end type nonzeromat4legs

end module DefineTypes

!===================================================================
module Constants
!===================================================================
  use Precision
  use DefineTypes
  implicit none

  real(kind=dp),save :: PI
  real(kind=dp),parameter :: CCC= 137.035999679_dp  !light speed
  real(kind=dp),parameter :: Ze = -1._dp  ! factor for electron charge
  complex(kind=dp),parameter :: IU = (0._dp,1._dp)    ! imaginary unit

  complex(kind=dp) :: Gam0(4,4),Gam(3,4,4) ! 3 denotes gamma1,gamma2,gamma3
  complex(kind=dp) :: Gam5(4,4)   ! gamma_5 = -gamma^5
  complex(kind=dp) :: Sigma(3,4,4) ! 4x4 Sigma matrix
  complex(kind=dp) :: Gam0k(3,4,4) ! gamma^0 gamma^k = alpha
!  complex(kind=dp) :: GamJJ(4,4,4,4) ! sum_{k=1}^3 [gamma^0 gamma^k]_{alpha beta} [gamma^0 gamma^k]_{gamma delta}
!                                     ! used in calc_intthetajj_mat
  integer,parameter :: Ngamjj = 24 ! number of non-zero components in GamJJ
  type(nonzeromat4legs) :: GamJJ(Ngamjj)
  
  real(kind=dp),parameter :: Rg_OVER_Msolar= 5.59e13_dp  ! to convert mass in solar mass to r_g in a.u.
  integer :: NOPT_GR
  integer :: NOPT_GRw  ! if 1, print out h_Qmat and stop.
  integer :: NOPT_GR_T, NOPT_GR_S, NOPT_GR_M, NOPT_GR_V

  !-----------------------------------------
  ! for discretization of photon momentum
  !-----------------------------------------
!  real(kind=dp),parameter :: P0MAX = 40._dp
  real(kind=dp) :: P0MAX
!!$  integer,parameter :: NP0   = 2  ! # of norm of photon momentum (photon energy)
!!$  integer,parameter :: NPTH  = 3  ! # in theta direction (0 <= theta <= pi)
!!$  integer,parameter :: NPPHI = 4  ! # in phi direction (0<= phi < 2pi)
!!$  integer,parameter :: Nph = NP0*NPTH*NPPHI*2 ! j=1~Nph
  integer :: NP0     ! # of norm of photon momentum (photon energy)
  integer :: NPTH    ! # in theta direction (0 <= theta <= pi)
  integer :: NPPHI   ! # in phi direction (0<= phi < 2pi)
  integer :: Nph  ! j=1~Nph

  integer :: NOCC  ! 1~nocc electron orbitals are occupied.

  !--------------------------------------------
  ! for numerical integration of calJ and calL
  !--------------------------------------------
  integer,parameter :: NB_Lret = 2  ! 3D integration (dr)
  integer,parameter :: NB_Jret = 2  ! 3D integration (dr)
  integer,parameter :: NB_L2   = 2  ! 2D integration (ds)
  integer,parameter :: NB_J2   = 2  ! 2D integration (ds)

  !-----------------------------------------
  ! for time evolution
  !-----------------------------------------
  real(kind=dp) :: DeltaT ! time step 

  !  complex(kind=dp) :: ALPHA_COH(Nph)   ! <coherent|hat(a)|coherent>
  complex(kind=dp), allocatable :: ALPHA_COH(:)   ! <coherent|hat(a)|coherent>
  
end module Constants

!===================================================================
module IntegralStorage
!===================================================================
  use Precision
  implicit none

  logical :: there_is_twoele
!  character(LEN=80),parameter :: file_twoele = 'twoele.dat'  ! 200
  character(LEN=80) :: file_twoele  ! 200
  real(kind=dp),parameter :: TH_twoele = 1.e-16_dp  ! do not store two electron integral smaller than this value.
  integer :: N_twoele ! number of two-electron integral larger than TH_twoele

  character(LEN=80),parameter :: file_twoele_bin = 'twoele.bin'  ! 300

  !------------------------------------------------------------------------------------------------
  ! these variables are now used only in subroutine setQmat_intIJJ_nz(intIJJ_Qmat_nz), which may not be used in future.
  logical :: there_is_intIJJ
  character(LEN=80) :: file_intIJJ  ! 204
  real(kind=dp),parameter :: TH_intIJJ = 1.e-16_dp  ! do not store IJJ integral smaller than this value.
  integer :: N_intIJJ ! number of IJJ integral larger than TH_intIJJ
  !------------------------------------------------------------------------------------------------

  logical :: there_is_intKJJ
  character(LEN=80) :: file_intKJJ  ! 204
  real(kind=dp),parameter :: TH_intKJJ = 1.e-16_dp  ! do not store KJJ integral smaller than this value.
  integer :: N_intKJJ ! number of KJJ integral larger than TH_intKJJ
  real(kind=dp) :: AlphaJJ_MAX ! Upper limit of alpha integration
  integer :: N_alphaJJ ! Number of grid spacing for alpha integration
  real(kind=dp) :: Diff_u_t ! u-t

!!$  logical :: there_is_intret
!!$  character(LEN=80),parameter :: file_intret = 'intret.dat'  ! 201
!!$  integer,parameter :: NRET = 7  ! number of time points to do integration for retardation

  logical :: there_is_nucele
  character(LEN=80),parameter :: file_nucele = 'nucele.dat'  ! 202  <--- fort.13

  logical :: there_is_twonuc
  character(LEN=80),parameter :: file_twonuc = 'twonuc.dat'  ! 203 <--- fort.14


end module IntegralStorage

!===================================================================
module NucBasis
!===================================================================
  use Precision
  implicit none

  integer,parameter :: Z_N = 1 ! nuclear charge for H
  integer,parameter :: N_PGN = 15 ! number of basis set function for expansion function for nuclear field.
  real(kind=dp) :: vecR_N(3,N_PGN)  ! center position of 1s gaussian
  integer,parameter :: NBS_PHI = 7 ! number of expansion functions for nuclear field
  integer :: NBS_N  ! =NBS_PHI for boson, =2xNBS_PHI for fermion
  real(kind=dp) :: L_well ! size of the single well potential
  real(kind=dp) :: L_pg ! range for p.g. (larger than L_well)

  ! parameters below are set in set_NucBasis in sub_int_nuc.f90
  integer :: NNUC ! number of nucleus
  real(kind=dp) :: m_N  ! nuclear mass  
  real(kind=dp) :: alpha_N  ! exponent for nucleus
!  real(kind=dp) :: R_N 
  integer :: NUCTYPE ! boson(0) of fermion(1)
  complex(kind=dp) :: c_nuc(NBS_PHI,N_PGN) ! coefficients to build plane waves from primitive gaussians

end module NucBasis

!===================================================================
module System_Medium
!===================================================================
  use Precision
  implicit none

  real(kind=dp) :: L_Med ! size of medium (M) (assumed to be cube)
  real(kind=dp) :: R_Sys ! size of system (A) (assumed to be ball)

end module System_Medium

!===================================================================
module Interface_Mod
!===================================================================
  Interface

     subroutine setQmat_twoele_nz(twoele_Qmat_nz)
       use DefineTypes
       type(nonzeromat4legs),allocatable,intent(out) :: twoele_Qmat_nz(:)
     end subroutine setQmat_twoele_nz

     subroutine setQmat_intIJJ_nz(intIJJ_Qmat_nz)
       use DefineTypes
       type(nonzeromat4legs),allocatable :: intIJJ_Qmat_nz(:)
     end subroutine setQmat_intIJJ_nz

     subroutine setQmat_intKJJ_nz(intKJJ_Qmat_nz)
       use DefineTypes
       type(nonzeromat4legs),allocatable :: intKJJ_Qmat_nz(:)
     end subroutine setQmat_intKJJ_nz

  end interface
end module Interface_mod
