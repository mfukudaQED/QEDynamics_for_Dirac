!=====================================================
MODULE EDM_CALCULATION
!=====================================================
  use Precision
  implicit none

  complex(kind=dp),allocatable,save :: intPsi_overlap(:,:,:)
  complex(kind=dp),allocatable,save :: intPsi_Lap(:,:,:)
  complex(kind=dp),allocatable,save :: intPsi_x_moment(:,:,:), intPsi_y_moment(:,:,:), intPsi_z_moment(:,:,:)
  complex(kind=dp),allocatable,save :: intPsi_xdy_ydx(:,:,:), intPsi_ydz_zdy(:,:,:), intPsi_zdx_xdz(:,:,:)
  complex(kind=dp),allocatable,save :: intPsi_ef_x(:,:,:),intPsi_ef_y(:,:,:),intPsi_ef_z(:,:,:)
  complex(kind=dp),allocatable,save :: psi_local(:,:,:,:), dpsi_local(:,:,:,:,:)
  complex(kind=dp),allocatable,save :: d2psi_local(:,:,:,:,:,:), ddpsi_local(:,:,:,:,:,:)
  complex(kind=dp),allocatable,save :: density(:)
  complex(kind=dp),allocatable,save :: spin(:,:)
  complex(kind=dp),allocatable,save :: spin_small(:,:)
  complex(kind=dp),allocatable,save :: grad_N(:,:)
  complex(kind=dp),allocatable,save :: je(:,:)
  complex(kind=dp),allocatable,save :: torq(:,:)
  complex(kind=dp),allocatable,save :: div_torq(:)
  complex(kind=dp),allocatable,save :: torq_AM(:,:)
  complex(kind=dp),allocatable,save :: zeta_force(:,:)
  complex(kind=dp),allocatable,save :: zeta_potential(:)
  complex(kind=dp),allocatable,save :: torqEDM_ele(:,:)
  complex(kind=dp),allocatable,save :: torqEDM_mag(:,:) 
  complex(kind=dp),allocatable,save :: EeffEDM_nuc(:)
  complex(kind=dp),allocatable,save :: chirality_density(:)
  complex(kind=dp),allocatable,save :: tau(:,:)

  real(kind=dp),allocatable,save :: vecA_M(:,:), vecE_nuc(:,:)
  real(kind=dp),allocatable,save :: x1(:), x2(:), x3(:)
  real(kind=dp),allocatable,save :: xA(:,:,:), yA(:,:,:), zA(:,:,:)

  character(LEN=8),save :: today

END MODULE EDM_CALCULATION

!=====================================================
MODULE BINOMIAL
!=====================================================
  implicit none

  integer,parameter :: nCr(0:16,0:16) =                                               &
    & reshape(                                                                        &
    &   (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,                                          &
    &     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,                                          &
    &     1,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,                                          &
    &     1,3,3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,                                          &
    &     1,4,6,4,1,0,0,0,0,0,0,0,0,0,0,0,0,                                          &
    &     1,5,10,10,5,1,0,0,0,0,0,0,0,0,0,0,0,                                        &
    &     1,6,15,20,15,6,1,0,0,0,0,0,0,0,0,0,0,                                       &
    &     1,7,21,35,35,21,7,1,0,0,0,0,0,0,0,0,0,                                      &
    &     1,8,28,56,70,56,28,8,1,0,0,0,0,0,0,0,0,                                     &
    &     1,9,36,84,126,126,84,36,9,1,0,0,0,0,0,0,0,                                  &
    &     1,10,45,120,210,252,210,120,45,10,1,0,0,0,0,0,0,                            &
    &     1,11,55,165,330,462,462,330,165,55,11,1,0,0,0,0,0,                          &
    &     1,12,66,220,495,792,924,792,495,220,66,12,1,0,0,0,0,                        &
    &     1,13,78,286,715,1287,1716,1716,1287,715,286,78,13,1,0,0,0,                  &
    &     1,14,91,364,1001,2002,3003,3432,3003,2002,1001,364,91,14,1,0,0,             &
    &     1,15,105,455,1365,3003,5005,6435,6435,5005,3003,1365,455,105,15,1,0,        &
    &     1,16,120,560,1820,4368,8008,11440,12870,11440,8008,4368,1820,560,120,16,1/) &
    & ,shape(nCr))
  integer,parameter :: d_factorial(-1:19) = (/1,1,1,2,3,8,15,48,105,384,945,3840,10395,46080,135135,&
                                            & 645120,2027025,10321920,34459425,185794560,654729075/)

END MODULE BINOMIAL

!=====================================================
MODULE MONTH_TO_DAY
!=====================================================
  implicit none

  integer,save :: days(12) = (/31,28,31,30,31,31,31,31,30,31,30,31/)

END MODULE MONTH_TO_DAY

