! Last Change:28-Sep-2015.
!====================================================================================
! 2013.11.07 
! - function intSpin_Qmat(k,nn,mm)
! - subroutine calc_intdensmom_Qmat(nn,mm,intdensmom_Qmat)
! - subroutine calc_intOrbang_Qmat(nn,mm,intOrbang_Qmat)
! - subroutine calc_intOrbang_AM_Qmat(nn,mm,intOrbang_AM_Qmat)
! - function inttau_Qmat(k,l,nn,mm)
! - function inttorq_Qmat(k,nn,mm)
! - function intjmoment_Qmat(k,l,nn,mm)
! - function intFj_Qmat(j,k,nn,mm)
! - function intHrel_Qmat(nn,mm)
! - subroutine setQmat_jAM(jAM_Qmat)
! - subroutine setQmat_jAMnuc(jAMnuc_Qmat)
! - subroutine setQmat_intrhoA0M(intrhoA0M_Qmat)
!====================================================================================

!============================================================
function intSpin_Qmat(k,nn,mm)
! t : time = timestep*DeltaT
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants  ! use Nph
  implicit none

  integer,intent(in) :: k
  complex(kind=dp) :: intSpin_mat
  complex(kind=dp) :: intSpin_Qmat
  integer,intent(in) :: nn,mm

  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b

  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  intSpin_Qmat = intSpin_mat(k,n,a,m,b)

  return
end function intSpin_Qmat

!============================================================
function intMagHyp_Qmat(posR,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants  ! use Nph
  implicit none

  complex(kind=dp) :: intMagHyp_mat
  complex(kind=dp) :: intMagHyp_Qmat
  integer,intent(in) :: nn,mm
  real(kind=dp),intent(in) :: posR(3)

  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b

  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  intMagHyp_Qmat = intMagHyp_mat(posR,n,a,m,b)

  return
end function intMagHyp_Qmat

!============================================================
subroutine calc_intdensmom_Qmat(nn,mm,intdensmom_Qmat)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  implicit none

  complex(kind=dp),intent(out) :: intdensmom_Qmat(3)
  integer,intent(in) :: nn,mm

  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  integer :: k
  character(LEN=1) :: a,b

  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  call calc_intdensmom_mat(n,a,m,b,intdensmom_Qmat)

!  write(*,*)'intdensmom_Qmat',nn,mm,intdensmom_Qmat(1),intdensmom_Qmat(2),intdensmom_Qmat(3)

  return
end subroutine calc_intdensmom_Qmat

!============================================================
subroutine calc_intOrbang_Qmat(nn,mm,intOrbang_Qmat)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  implicit none

  complex(kind=dp),intent(out) :: intOrbang_Qmat(3)
  integer,intent(in) :: nn,mm

  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  integer :: k
  character(LEN=1) :: a,b

  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  call calc_intOrbang_mat(n,a,m,b,intOrbang_Qmat)

!  write(*,*)'intOrbang_Qmat',nn,mm,intOrbang_Qmat(1),intOrbang_Qmat(2),intOrbang_Qmat(3)

  return
end subroutine calc_intOrbang_Qmat

!============================================================
subroutine calc_intOrbang_AM_Qmat(nn,mm,intOrbang_AM_Qmat)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  implicit none

  complex(kind=dp),intent(out) :: intOrbang_AM_Qmat(3)
  integer,intent(in) :: nn,mm

  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  integer :: k
  character(LEN=1) :: a,b

  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  call calc_intOrbang_AM_mat(n,a,m,b,intOrbang_AM_Qmat)

!  write(*,*)'intOrbang_Qmat',nn,mm,intOrbang_Qmat(1),intOrbang_Qmat(2),intOrbang_Qmat(3)

  return
end subroutine calc_intOrbang_AM_Qmat

!============================================================
function inttau_Qmat(k,l,nn,mm)
! t : time = timestep*DeltaT
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants  ! use Nph
  implicit none

  integer,intent(in) :: k,l
  complex(kind=dp) :: inttau_mat
  complex(kind=dp) :: inttau_Qmat
  integer,intent(in) :: nn,mm

  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b

  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  inttau_Qmat = inttau_mat(k,l,n,a,m,b)

  return
end function inttau_Qmat

!============================================================
function inttorq_Qmat(k,nn,mm)
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none

  integer,intent(in) :: nn,mm
  integer,intent(in) :: k

  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  complex(kind=dp) :: inttorq_Qmat
  complex(kind=dp) :: inttau_Qmat

  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  if(k==1) then
    inttorq_Qmat = -inttau_Qmat(2,3,nn,mm) + inttau_Qmat(3,2,nn,mm)
  else if(k==2) then
    inttorq_Qmat =  inttau_Qmat(1,3,nn,mm) - inttau_Qmat(3,1,nn,mm)
  else if(k==3) then
    inttorq_Qmat = -inttau_Qmat(1,2,nn,mm) + inttau_Qmat(2,1,nn,mm)
  end if

end function inttorq_Qmat 

!============================================================
function intjmoment_Qmat(k,l,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants  ! use Nph
  implicit none

  integer,intent(in) :: k,l
  complex(kind=dp) :: intjmoment_mat
  complex(kind=dp) :: intjmoment_Qmat
  integer,intent(in) :: nn,mm

  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b

  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  intjmoment_Qmat = intjmoment_mat(k,l,n,a,m,b)

  return
end function intjmoment_Qmat

!============================================================
function intFj_Qmat(j,k,nn,mm)
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none

  integer,intent(in) :: nn,mm
  integer,intent(in) :: j,k

  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  real(kind=dp) :: P0j
  real(kind=dp) :: vecPj(3) ! photon momentum vector in Cartesian coord.
  complex(kind=dp) :: vecpolj(3) ! (circular) polarization vector (depends on sig)
  complex(kind=dp) :: intF_mat
  complex(kind=dp) :: intFj_Qmat

    call calc_mode_pj(j,P0j,vecPj,vecpolj)
    call index_from_Qmat(nn,n,a)
    call index_from_Qmat(mm,m,b)

    intFj_Qmat = intF_mat(k,vecPj,n,a,m,b)

end function intFj_Qmat 

!============================================================
function intHrel_Qmat(nn,mm)
! t : time = timestep*DeltaT
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants  ! use Nph
  implicit none

  complex(kind=dp) :: intHrel_mat
  complex(kind=dp) :: intHrel_Qmat
  integer,intent(in) :: nn,mm

  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b

  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  intHrel_Qmat = intHrel_mat(n,a,m,b)

  return
end function intHrel_Qmat

!======================================================
subroutine setQmat_jAM(jAM_Qmat)
! 1 ~ 2*NBS : "+"
! 2*NBS +1 ~ 2*NBS : "-"
! in "+" and "-", 1,1bar,2,2bar, ... respectively
!======================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=8),intent(out) :: jAM_Qmat(4*NBS,4*NBS)

  complex(kind=8) :: intHrel_Qmat
  integer :: nn,mm

  do nn=1,4*NBS
     do mm=1,4*NBS

        jAM_Qmat(nn,mm) = intHrel_Qmat(nn,mm)
        
     end do
  end do

end subroutine setQmat_jAM

!======================================================
subroutine setQmat_jAMnuc(jAMnuc_Qmat)
! 1 ~ 2*NBS : "+"
! 2*NBS +1 ~ 2*NBS : "-"
! in "+" and "-", 1,1bar,2,2bar, ... respectively
!======================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=8),intent(out) :: jAMnuc_Qmat(4*NBS,4*NBS)

  complex(kind=8) :: intMagHyp_mat !function
  integer :: nn,mm

  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b

  integer :: l
  real(kind=dp) :: posR(3)

  do l=1,NAT
     posR(1)=xc(l); posR(2)=yc(l); posR(3)=zc(l)
     do nn=1,4*NBS
        do mm=1,4*NBS
           call index_from_Qmat(nn,n,a)
           call index_from_Qmat(mm,m,b)
           jAMnuc_Qmat(nn,mm) = - intMagHyp_mat(posR,n,a,m,b)
!           write(*,*)nn,mm,jAMnuc_Qmat(nn,mm)
        end do
     end do
  end do
!  stop

end subroutine setQmat_jAMnuc

!======================================================
subroutine setQmat_intrhoA0M(intrhoA0M_Qmat)
! 1 ~ 2*NBS : "+"
! 2*NBS +1 ~ 2*NBS : "-"
! in "+" and "-", 1,1bar,2,2bar, ... respectively
!======================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=8),intent(out) :: intrhoA0M_Qmat(4*NBS,4*NBS)

  complex(kind=8) :: intrhoA0M_mat !function
  integer :: nn,mm

  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b

!        write(*,*)'A0M'
  do nn=1,4*NBS
     do mm=1,4*NBS
        call index_from_Qmat(nn,n,a)
        call index_from_Qmat(mm,m,b)
        intrhoA0M_Qmat(nn,mm) = intrhoA0M_mat(n,a,m,b)
!        write(*,*)nn,mm,intrhoA0M_Qmat(nn,mm)
     end do
  end do

end subroutine setQmat_intrhoA0M
