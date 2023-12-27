! Last Change:08-Jun-2014.
!  function Aradj_vec(j,i,x,y,z)     <coherent|Arad^i(t)|coherent>
!  function Arad_vec(t,i,x,y,z)
!  function tau_Arad_mat(t,k,l,x,y,z,p,a,q,b)
!  function t_Arad_mat(t,i,x,y,z,p,a,q,b)
!  function tau_Arad_Qmat(t,k,l,x,y,z,nn,mm)
!  function t_Arad_Qmat(t,i,x,y,z,nn,mm)
!  subroutine calc_t_Arad_Qmat(tmp2_j_Qmat,tmp_Arad_vec,tmp_t_Arad_Qmat)
!  subroutine calc_tau_Arad_Qmat(tmp2_j_Qmat,tmp_Arad_vec,tmp_tau_Arad_Qmat)
!!$!======================================================================
function Aradj_vec(j,i,x,y,z)
! <coherent|Arad^i(t)|coherent>
! t : time = timestep*DeltaT
! p (photon momentum) specified in spherical coordinate (p0,th,phi)
! sig (+/-) specifies rotating direction of circular polarization (which depends on th and phi)
!============================================================
  use Constants  ! use P0MAX,NP0,NPTH,NPPHI,Nph,ALPHA_COH
  implicit none

  complex(kind=8) :: Aradj_vec
  integer,intent(in) :: i
  integer,intent(in) :: j ! collective index which expresees photon momentum and polarization
  integer :: k
  real(kind=8),intent(in) :: x,y,z

  real(kind=8) :: p0,th,phi ! photon momentum vector in spherical coord.
  character(LEN=1) :: sig ! polarization (+ or -)  
  integer :: ip0,ith,iphi
  real(kind=8) :: dp0,dth,dphi ! mesh width

  real(kind=8) :: vecP(3) ! photon momentum vector in Cartesian coord.
  real(kind=8) :: vecPj(3) ! photon momentum vector in Cartesian coord.
  complex(kind=8) :: vec_pol(3) ! (circular) polarization vector (depends on sig)
  complex(kind=8) :: vecpolj(3) ! (circular) polarization vector (depends on sig)
  complex(kind=8) :: tmpA
  real(kind=8) :: DeltaPj,P0j

  Aradj_vec = (0d0,0d0)

    call calc_mode_pj(j,P0j,vecPj,vecpolj)

    Aradj_vec = ALPHA_COH(j)*vecpolj(i)*cdexp(IU*(vecPj(1)*x +vecPj(2)*y +vecPj(3)*z))

!    write(*,*)'ALPHA_COH(',j,')=',ALPHA_COH(j)
!    write(*,*)'tmpA:',tmpA
!    write(*,*)'DeltaPj(',j,')=',DeltaPj
!    write(*,*) j, Aradj_vec

  return
end function Aradj_vec

!============================================================
function Arad_vec(t,i,x,y,z)
! 17-Oct-2011
! <coherent|Arad^i(t)|coherent>
! t : time = timestep*DeltaT
! p (photon momentum) specified in spherical coordinate (p0,th,phi)
! sig (+/-) specifies rotating direction of circular polarization (which depends on th and phi)
!============================================================
  use Constants  ! use P0MAX,NP0,NPTH,NPPHI,Nph,ALPHA_COH
  implicit none

  complex(kind=8) :: Arad_vec
  complex(kind=8) :: Aradj_vec
  complex(kind=8) :: Ajt_vec
  integer,intent(in) :: i
  integer :: j ! collective index which expresees photon momentum and polarization
  real(kind=8),intent(in) :: x,y,z
  real(kind=8),intent(in) :: t  ! time
  real(kind=8) :: DeltaPj,P0j
  Arad_vec = (0.d0,0.d0)

  do j=1,Nph
    call calc_DeltaPj_and_P0j(j,DeltaPj,P0j)
    Ajt_vec = Aradj_vec(j,i,x,y,z)*cdexp(-IU*CCC*P0j*t)*dsqrt(CCC/P0j)/(2.d0*PI)
    Arad_vec = Arad_vec + Ajt_vec + dconjg(Ajt_vec)
  end do

end function Arad_vec

!============================================================
function tau_Arad_mat(t,k,l,x,y,z,p,a,q,b)
!  17-Oct-2011
! t : time = timestep*DeltaT
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================

  use Constants  ! use CCC, Nph
  use DefineTypes
  use DiracOutput
  implicit none

  complex(kind=8) :: tau_Arad_mat
  real(kind=8),intent(in) :: t  ! time
  integer,intent(in) :: k,l
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b
  real(kind=8),intent(in) :: x,y,z

  complex(kind=8) :: j_mat,tau_mat,Arad_vec

!   tau_Arad_mat = tau_mat(k,l,x,y,z,p,a,q,b) + j_mat(l,x,y,z,p,a,q,b)*Arad_vec(t,k,x,y,z)/CCC
   tau_Arad_mat = j_mat(l,x,y,z,p,a,q,b)*Arad_vec(t,k,x,y,z)/CCC

  return
end function tau_Arad_mat

!============================================================
function t_Arad_mat(t,i,x,y,z,p,a,q,b)
! 17-Oct-2011
! t : time = timestep*DeltaT
! i-th component of t(ab)_pq, i=1,3
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use DefineTypes
  use DiracOutput
  use Constants  ! use Nph
  implicit none
 
  complex(kind=8) :: t_Arad_mat
  real(kind=8),intent(in) :: t  ! time
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  complex(kind=8) :: tau_Arad_mat

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in t_mat."
     stop
  end if

  if(i.eq.1) then
     t_Arad_mat = -tau_Arad_mat(t,2,3,x,y,z,p,a,q,b) +tau_Arad_mat(t,3,2,x,y,z,p,a,q,b)
  elseif(i.eq.2) then
     t_Arad_mat =  tau_Arad_mat(t,1,3,x,y,z,p,a,q,b) -tau_Arad_mat(t,3,1,x,y,z,p,a,q,b)
  elseif(i.eq.3) then
     t_Arad_mat = -tau_Arad_mat(t,1,2,x,y,z,p,a,q,b) +tau_Arad_mat(t,2,1,x,y,z,p,a,q,b)
  else
     write(*,*) "i should be 1-3 in t_mat."
     stop
  end if
  
  return
end function t_Arad_mat

!============================================================
function tau_Arad_Qmat(t,k,l,x,y,z,nn,mm)
! 4-Nov-2011
! t : time = timestep*DeltaT
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use DiracOutput
  use Constants  ! use Nph
  implicit none

  complex(kind=8) :: tau_Arad_Qmat
  real(kind=8),intent(in) :: t  ! time = timestep * DeltaT
  integer,intent(in) :: k,l
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=8) :: tau_Arad_mat
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b

  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  tau_Arad_Qmat = tau_Arad_mat(t,k,l,x,y,z,n,a,m,b)

  return
end function tau_Arad_Qmat

!============================================================
function t_Arad_Qmat(t,i,x,y,z,nn,mm)
! 4-Nov-2011
! t : time = timestep*DeltaT
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use DiracOutput
  use Constants  ! use Nph
  implicit none

  complex(kind=8) :: t_Arad_Qmat
  real(kind=8),intent(in) :: t  ! time = timestep * DeltaT
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=8) :: t_Arad_mat
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  
  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  t_Arad_Qmat = t_Arad_mat(t,i,x,y,z,n,a,m,b)

  return
end function t_Arad_Qmat

!============================================================================
subroutine calc_t_Arad_Qmat(tmp2_j_Qmat,tmp_Arad_vec,tmp_t_Arad_Qmat)
!============================================================================

  use DiracOutput
  use Constants
  implicit none

  integer k,l,nn,mm
  complex(kind=8),intent(in) :: tmp2_j_Qmat(3,4*NBS,4*NBS)
  complex(kind=8),intent(in) :: tmp_Arad_vec(3)
  complex(kind=8),intent(out) :: tmp_t_Arad_Qmat(3,4*NBS,4*NBS)
  complex(kind=8) :: tmp_tau_Arad_Qmat(3,3,4*NBS,4*NBS)

  do k=1,3
    do l=1,3
      do nn=1,4*NBS
        do mm=1,4*NBS
          tmp_tau_Arad_Qmat(k,l,nn,mm) = tmp2_j_Qmat(l,nn,mm)*tmp_Arad_vec(k)/CCC
        end do
      end do
    end do
  end do

  do nn=1,4*NBS
    do mm=1,4*NBS
      tmp_t_Arad_Qmat(1,nn,mm) = -tmp_tau_Arad_Qmat(2,3,nn,mm) +tmp_tau_Arad_Qmat(3,2,nn,mm)
      tmp_t_Arad_Qmat(2,nn,mm) =  tmp_tau_Arad_Qmat(1,3,nn,mm) -tmp_tau_Arad_Qmat(3,1,nn,mm)
      tmp_t_Arad_Qmat(3,nn,mm) = -tmp_tau_Arad_Qmat(1,2,nn,mm) +tmp_tau_Arad_Qmat(2,1,nn,mm)
    end do
  end do

end subroutine calc_t_Arad_Qmat
!============================================================================

!============================================================================
subroutine calc_tau_Arad_Qmat(tmp2_j_Qmat,tmp_Arad_vec,tmp_tau_Arad_Qmat)
!============================================================================

  use DiracOutput
  use Constants
  implicit none

  integer k,l,nn,mm
  complex(kind=8),intent(in) :: tmp2_j_Qmat(3,4*NBS,4*NBS)
  complex(kind=8),intent(in) :: tmp_Arad_vec(3)
  complex(kind=8),intent(out) :: tmp_tau_Arad_Qmat(3,3,4*NBS,4*NBS)

  do k=1,3
    do l=1,3
      do nn=1,4*NBS
        do mm=1,4*NBS
          tmp_tau_Arad_Qmat(k,l,nn,mm) = tmp2_j_Qmat(l,nn,mm)*tmp_Arad_vec(k)/CCC
        end do
      end do
    end do
  end do

end subroutine calc_tau_Arad_Qmat
!============================================================================
