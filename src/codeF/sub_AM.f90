! Last Change: 08-Jun-2014.
!====================================================================================
! 2013.11.07
! - subroutine read_param_AM
! - function func_AMvec(t,k,x,y,z)
! - function tau_AM_mat(t,k,l,x,y,z,p,a,q,b)
! - function t_AM_mat(t,i,x,y,z,p,a,q,b)
! - function tau_AM_Qmat(t,k,l,x,y,z,nn,mm)
! - function t_AM_Qmat(t,i,x,y,z,nn,mm)
!====================================================================================

!============================================================================
subroutine read_param_AM
!============================================================================
  use Precision
  use param_AM
  implicit none

  open(unit=110,file='param_AM.ini')
  read(110,*) Bmag(1)
  read(110,*) Bmag(2)
  read(110,*) Bmag(3)
  read(110,*) Eelec(1)
  read(110,*) Eelec(2)
  read(110,*) Eelec(3)
  read(110,*) magdip(1)
  read(110,*) magdip(2)
  read(110,*) magdip(3)
  read(110,*) omega_A0M
  close(unit=110)

  write(*,*)"# Bmag :",Bmag(1),Bmag(2),Bmag(3)
  write(*,*)"# Eelec :",Eelec(1),Eelec(2),Eelec(3)
  write(*,*)"# magdip :",magdip(1),magdip(2),magdip(3)
  write(*,*)"# omega_A0M :" , omega_A0M

end subroutine read_param_AM
!============================================================================

  
!!$!============================================================================
!!$subroutine set_AMvec(x,y,z,B,AMvec)
!!$  ! B = rot(AM) = (0,0,B)
!!$  ! AM = (-1/2*B*y,1/2*B*x,0)
!!$!============================================================================
!!$  use Precision
!!$  implicit none
!!$
!!$  real(kind=dp), intent(in) :: x,y,z
!!$  complex(kind=dp), intent(in) :: B
!!$  complex(kind=dp), intent(out) :: AMvec(3)
!!$
!!$  AMvec(1) = - 0.5d0 *B *y
!!$  AMvec(2) =   0.5d0 *B *x
!!$  AMvec(3) =   0.d0
!!$
!!$end subroutine set_AMvec
!!$!============================================================================

!============================================================================
function func_AMvec(t,k,x,y,z)
  ! BM = rot(AM)
!============================================================================
  use Precision
  use param_AM
  implicit none

  real(kind=8),intent(in) :: t  ! time
  integer, intent(in) :: k
  real(kind=dp), intent(in) :: x,y,z
!  complex(kind=dp) :: Bmag
  complex(kind=dp) :: func_AMvec

!!$  ! You have to change intHrel_mat and intHrel_pt_Qmat in sub_int_phys.f90 too.
!!$  if(direcB==1) then
!!$    !-------------------------------
!!$    ! BM = rot(AM) = (Bmag,0,0)
!!$    ! AM = (0,-1/2*Bmag*z,1/2*Bmag*y)
!!$!    Bmag = (1.d0,0.d0)
!!$    if(k==1) then
!!$      func_AMvec =   0.d0
!!$    else if(k==2) then
!!$      func_AMvec = - 0.5d0 *Bmag *z
!!$    else if(k==3) then
!!$      func_AMvec =   0.5d0 *Bmag *y
!!$    end if
!!$    !-------------------------------
!!$  else if(direcB==2) then
!!$    !-------------------------------
!!$    ! BM = rot(AM) = (0,Bmag,0)
!!$    ! AM = (1/2*Bmag*z,0,-1/2*Bmag*x)
!!$!    Bmag = (1.d0,0.d0)
!!$    if(k==1) then
!!$      func_AMvec =   0.5d0 *Bmag *z
!!$    else if(k==2) then
!!$      func_AMvec =   0.d0
!!$    else if(k==3) then
!!$      func_AMvec = - 0.5d0 *Bmag *x
!!$    end if
!!$    !-------------------------------
!!$  else if(direcB==3) then
!!$    !-------------------------------
!!$    ! BM = rot(AM) = (0,0,Bmag)
!!$    ! AM = (-1/2*B*y,1/2*B*x,0)
!!$!    Bmag = (1.d0,0.d0)
!!$    if(k==1) then
!!$      func_AMvec = - 0.5d0 *Bmag *y
!!$    else if(k==2) then
!!$      func_AMvec =   0.5d0 *Bmag *x
!!$    else if(k==3) then
!!$      func_AMvec =   0.d0
!!$    end if
!!$    !-------------------------------
!!$  else 
!!$    write(*,*)'# direcB should be 1-3.'
!!$  end if

    !-------------------------------
    ! BM = rot(AM) = (Bmag(1),Bmag(2),Bmag(3))
    ! AM = 0.5d0*\vec{BM} \times \vec{r}
    !    = (1/2*Bmagy*z -1/2*Bmagz*y, 1/2*Bmagz*x -1/2*Bmagx*z, 1/2*Bmagx*y -1/2*Bmagy*x)
    if(k==1) then
      func_AMvec = 0.5d0*Bmag(2)*z -0.5d0*Bmag(3)*y
    else if(k==2) then
      func_AMvec = 0.5d0*Bmag(3)*x -0.5d0*Bmag(1)*z
    else if(k==3) then
      func_AMvec = 0.5d0*Bmag(1)*y -0.5d0*Bmag(2)*x
    end if
    !-------------------------------

!!$    !-------------------------------
!!$    ! BM = rot(AM) = (Bmag,0,0)
!!$    ! AM = (0,0,Bmag*y)
!!$!    Bmag = (1.d0,0.d0)
!!$    if(k==1) then
!!$      func_AMvec = Bmag(1) *x
!!$    else if(k==2) then
!!$      func_AMvec = - Bmag(1) *y
!!$    else if(k==3) then
!!$      func_AMvec =   0.d0
!!$    end if
!!$    !-------------------------------
    return
end function func_AMvec
!============================================================================

!============================================================
function tau_AM_mat(t,k,l,x,y,z,p,a,q,b)
! t : time = timestep*DeltaT
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================

  use Constants  ! use CCC, Nph
  use DefineTypes
  use DiracOutput
  implicit none

  complex(kind=8) :: tau_AM_mat
  real(kind=8),intent(in) :: t  ! time
  integer,intent(in) :: k,l
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b
  real(kind=8),intent(in) :: x,y,z

  complex(kind=8) :: j_mat,tau_mat,func_AMvec

  tau_AM_mat = j_mat(l,x,y,z,p,a,q,b)*func_AMvec(t,k,x,y,z)/CCC

  return
end function tau_AM_mat

!============================================================
function t_AM_mat(t,i,x,y,z,p,a,q,b)
! t : time = timestep*DeltaT
! i-th component of t(ab)_pq, i=1,3
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use DefineTypes
  use DiracOutput
  use Constants  ! use Nph
  implicit none
 
  complex(kind=8) :: t_AM_mat
  real(kind=8),intent(in) :: t  ! time
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  complex(kind=8) :: tau_AM_mat

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in t_mat."
     stop
  end if

  if(i.eq.1) then
     t_AM_mat = -tau_AM_mat(t,2,3,x,y,z,p,a,q,b) +tau_AM_mat(t,3,2,x,y,z,p,a,q,b)
  elseif(i.eq.2) then
     t_AM_mat =  tau_AM_mat(t,1,3,x,y,z,p,a,q,b) -tau_AM_mat(t,3,1,x,y,z,p,a,q,b)
  elseif(i.eq.3) then
     t_AM_mat = -tau_AM_mat(t,1,2,x,y,z,p,a,q,b) +tau_AM_mat(t,2,1,x,y,z,p,a,q,b)
  else
     write(*,*) "i should be 1-3 in t_mat."
     stop
  end if
  
  return
end function t_AM_mat


!============================================================
function tau_AM_Qmat(t,k,l,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=8) :: tau_AM_Qmat
  real(kind=8),intent(in) :: t  ! time
  integer,intent(in) :: k,l
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=8) :: tau_AM_mat
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  
  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  tau_AM_Qmat = tau_AM_mat(t,k,l,x,y,z,n,a,m,b)

  return
end function tau_AM_Qmat

!============================================================
function t_AM_Qmat(t,i,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=8) :: t_AM_Qmat
  complex(kind=8) :: t_AM_mat
  real(kind=8),intent(in) :: t  ! time
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=8) :: t_mat
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  
  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  t_AM_Qmat = t_AM_mat(t,i,x,y,z,n,a,m,b)

  return
end function t_AM_Qmat
