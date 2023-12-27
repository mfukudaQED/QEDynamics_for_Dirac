! Last Change:09-Jul-2015.
!====================================================================================
! 2013.9.20 
! - function Hso_nuc_Qmat(x,y,z,nn,mm)
! - function Hso_nuc_mat(x,y,z,p,a,q,b)
! - function pi_Qmat(i,x,y,z,nn,mm)
! - function pi_mat(i,x,y,z,p,a,q,b)
! - subroutine calc_pi_pq(i,x,y,z,NL,NS,pg,c_p,c_q,pi_pq_i)
! - function orbang_Qmat(i,x,y,z,nn,mm)
! - function orbang_mat(i,x,y,z,p,a,q,b)
! - function rots_Qmat(i,x,y,z,nn,mm)
! - function rots_mat(i,x,y,z,p,a,q,b)
! - subroutine calc_rots_pq(i,x,y,z,NL,NS,pg,c_p,c_q,rots_pq_i)
! - function js_Qmat(i,j,x,y,z,nn,mm)
! - function js_mat(i,j,x,y,z,p,a,q,b)
! - subroutine calc_js_pq(i,j,x,y,z,NL,NS,pg,c_p,c_q,js_pq_ij)

! 2012.1.25
! -function chiral_Qmat(x,y,z,nn,mm)
! -function chiral_ele(x,y,z,NEL)
! -function chiral_mat(x,y,z,p,a,q,b)
! -subroutine calc_chiral_pq(x,y,z,NL,NS,pg,c_p,c_q,chiral_pq)
!
! -subroutine calc_AA_pq(ip,iq,posRA,NL,NS,pg,c_p,c_q,AA_pq)
! -subroutine calc_AA_pq2(ip,iq,NL,NS,pg,c_p,c_q,NBS0,nucatt_pg,AA_pq)
! -subroutine save_gauss_int(ip,iq,i,j,posRA,NL,NS,pg,overlap_normpg,nucatt_normpg,ef_normpg)
! -subroutine calc_spol_pq(i,x,y,z,NL,NS,pg,c_p,c_q,spol_pq_i)
!====================================================================================
!============================================================
function Hso_ele_Qmat(x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=8) :: Hso_ele_Qmat
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=8) :: Hso_ele_mat
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  
  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  Hso_ele_Qmat = Hso_ele_mat(x,y,z,n,a,m,b)

  return
end function Hso_ele_Qmat

!============================================================
function Hso_ele_mat(x,y,z,p,a,q,b,r,c,s,d)
! Hso_ele = \frac{i Z_ee \hbar^2}{4m^2c^2} \epsilon_{ijk} \psi^\dagger {E}^i \Sigma^k \partial_j \psi 
! \hat{E}^k_{\rm ele} = -\sum_{N_A} \frac{Z_{N_A}e (\vec{R}_{N_A}-\vec{r})^k}{{|\vec{R}_{N_A}-\vec{r}|^3}} 
! contribution from eleleus
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use DefineTypes
  use DiracOutput
  implicit none
 
  complex(kind=dp) :: Hso_ele_mat
  real(kind=dp),intent(in) :: x,y,z
  integer,intent(in) :: p,q,r,s
  character(LEN=1),intent(in) :: a,b,c,d

  integer :: i,k
  real(kind=dp) :: vecR(3)
  complex(kind=dp) :: Eele(3)
  complex(kind=dp) :: tmp_Eele(3)
  complex(kind=8) :: js_mat !function

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_r(4,NMAX_PG),c_s(4,NMAX_PG)

  vecR(1) = x; vecR(2) = y; vecR(3) = z

  call copy_DiracOutput(r,c,s,d,pg,c_r,c_s)  ! set pg,c_p,c_q
  Eele(:) = (0._dp,0._dp)

  do i=1,4
     call calc_intE_pq(i,i,vecR,NBS_L,NBS_S,pg,c_r,c_s,tmp_Eele) 
     do k=1,3
        Eele(k) = Eele(k) +tmp_Eele(k)
     end do
  end do

  Hso_ele_mat = Eele(1)*(js_mat(2,3,x,y,z,p,a,q,b)-js_mat(3,2,x,y,z,p,a,q,b)) &
              &+Eele(2)*(js_mat(3,1,x,y,z,p,a,q,b)-js_mat(1,3,x,y,z,p,a,q,b)) &
              &+Eele(3)*(js_mat(1,2,x,y,z,p,a,q,b)-js_mat(2,1,x,y,z,p,a,q,b))
  
  return
end function Hso_ele_mat

!============================================================
function Hso_nuc_Qmat(x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=8) :: Hso_nuc_Qmat
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=8) :: Hso_nuc_mat
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  
  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  Hso_nuc_Qmat = Hso_nuc_mat(x,y,z,n,a,m,b)

  return
end function Hso_nuc_Qmat


!============================================================
function Hso_nuc_mat(x,y,z,p,a,q,b)
! Hso_nuc = \frac{i Z_ee \hbar^2}{4m^2c^2} \epsilon_{ijk} \psi^\dagger {E}^i \Sigma^k \partial_j \psi 
! \hat{E}^k_{\rm nuc} = -\sum_{N_A} \frac{Z_{N_A}e (\vec{R}_{N_A}-\vec{r})^k}{{|\vec{R}_{N_A}-\vec{r}|^3}} 
! contribution from nucleus
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use DefineTypes
  use DiracOutput
  implicit none
 
  complex(kind=dp) :: Hso_nuc_mat
  real(kind=dp),intent(in) :: x,y,z
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  integer :: i,k,iN
  real(kind=dp) :: posRA(3)
  real(kind=dp) :: vecR(3)
  real(kind=dp) :: vecD(3) ! vecR-vecRA
  real(kind=dp) :: norm !  |vecR-vecRA|
  real(kind=dp) :: Enuc(3)
  complex(kind=8) :: js_mat !function

  vecR(1) = x; vecR(2) = y; vecR(3) = z

  Enuc(:) = 0._dp
  do k=1,3
     do iN=1,NAT
        posRA(1) = xc(iN); posRA(2) = yc(iN); posRA(3) = zc(iN)
        vecD(:) = vecR(:)-posRA(:)
        norm = sqrt(vecD(1)**2+vecD(2)**2+vecD(3)**2)
        Enuc(k) = Enuc(k) +cn(iN)*vecD(k)/norm**3
!        write(*,*) iN,posRA(:),vecD(:),norm
     end do
  end do

  Hso_nuc_mat = Enuc(1)*(js_mat(2,3,x,y,z,p,a,q,b)-js_mat(3,2,x,y,z,p,a,q,b)) &
              &+Enuc(2)*(js_mat(3,1,x,y,z,p,a,q,b)-js_mat(1,3,x,y,z,p,a,q,b)) &
              &+Enuc(3)*(js_mat(1,2,x,y,z,p,a,q,b)-js_mat(2,1,x,y,z,p,a,q,b))
  
  return
end function Hso_nuc_mat

!============================================================
function pi_Qmat(i,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=8) :: pi_Qmat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=8) :: pi_mat
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  
  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  pi_Qmat = pi_mat(i,x,y,z,n,a,m,b)

  return
end function pi_Qmat

!============================================================
function pi_mat(i,x,y,z,p,a,q,b)
! i-th component of j(ab)_pq , i=1-3
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use DefineTypes
  use DiracOutput
  implicit none
 
  complex(kind=8) :: pi_mat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=8) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in pi_mat."
     stop
  end if

  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  call calc_pi_pq(i,x,y,z,NBS_L,NBS_S,pg,c_p,c_q,pi_mat)
  
  return
end function pi_mat

!============================================================
subroutine calc_pi_pq(i,x,y,z,NL,NS,pg,c_p,c_q,pi_pq_i)
! calculate i-th component of j_pq
! output : pi^i_pq (bilinear for charge current density) (i=1-3)
!============================================================
  use DefineTypes
  use Constants
  implicit none
  
  integer,intent(in) :: i   !pi^i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=8),intent(out) :: pi_pq_i

  complex(kind=8) :: psi_p(4),psi_q(4)
  complex(kind=8) :: dpsi_p(3,4),dpsi_q(3,4)  ! 3:x,y,z
  integer :: j

  call calc_psi(x,y,z,NL,NS,pg,c_p,psi_p)
  call calc_psi(x,y,z,NL,NS,pg,c_q,psi_q)
  call calc_dpsi(x,y,z,NL,NS,pg,c_p,dpsi_p)
  call calc_dpsi(x,y,z,NL,NS,pg,c_q,dpsi_q)

  pi_pq_i = (0.d0,0.d0)
  do j=1,4
     pi_pq_i = pi_pq_i + conjg(psi_p(j))*dpsi_q(i,j) -conjg(dpsi_p(i,j))*psi_q(j) 
  end do
  pi_pq_i = - pi_pq_i *IU*0.5d0
  return
end subroutine calc_pi_pq

!============================================================
function orbang_Qmat(i,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=8) :: orbang_Qmat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=8) :: orbang_mat
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  
  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  orbang_Qmat = orbang_mat(i,x,y,z,n,a,m,b)

  return
end function orbang_Qmat

!============================================================
function orbang_mat(i,x,y,z,p,a,q,b)
! i-th component of le_pq, i=1,3
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use DefineTypes
  use DiracOutput
  implicit none
 
  complex(kind=8) :: orbang_mat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=8) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  complex(kind=8) :: pi_mat !function

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in orbang_mat."
     stop
  end if

  if(i.eq.1) then
     orbang_mat = y *pi_mat(3,x,y,z,p,a,q,b) - z *pi_mat(2,x,y,z,p,a,q,b)
  elseif(i.eq.2) then
     orbang_mat = z *pi_mat(1,x,y,z,p,a,q,b) - x *pi_mat(3,x,y,z,p,a,q,b)
  elseif(i.eq.3) then
     orbang_mat = x *pi_mat(2,x,y,z,p,a,q,b) - y *pi_mat(1,x,y,z,p,a,q,b)
  else
     write(*,*) "i should be 1-3 in orbang_mat."
     stop
  end if
  
  return
end function orbang_mat

!============================================================
function rots_Qmat(i,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=8) :: rots_Qmat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=8) :: rots_mat
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  
  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  rots_Qmat = rots_mat(i,x,y,z,n,a,m,b)

  return
end function rots_Qmat

!============================================================
function rots_mat(i,x,y,z,p,a,q,b)
! i-th component of j(ab)_pq , i=1-3
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use DefineTypes
  use DiracOutput
  implicit none
 
  complex(kind=8) :: rots_mat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=8) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in rots_mat."
     stop
  end if

  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  call calc_rots_pq(i,x,y,z,NBS_L,NBS_S,pg,c_p,c_q,rots_mat)
  
  return
end function rots_mat

!============================================================
subroutine calc_rots_pq(i,x,y,z,NL,NS,pg,c_p,c_q,rots_pq_i)
! calculate i-th component of rots_pq
!============================================================
  use DefineTypes
  use Constants
  implicit none
  
  integer,intent(in) :: i   !rots^i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=8),intent(out) :: rots_pq_i

  complex(kind=8) :: psi_p(4),psi_q(4)
  complex(kind=8) :: dpsi_p(3,4),dpsi_q(3,4)  ! 3:x,y,z
  integer :: j,k

  call calc_psi(x,y,z,NL,NS,pg,c_p,psi_p)
  call calc_psi(x,y,z,NL,NS,pg,c_q,psi_q)
  call calc_dpsi(x,y,z,NL,NS,pg,c_p,dpsi_p)
  call calc_dpsi(x,y,z,NL,NS,pg,c_q,dpsi_q)

  rots_pq_i = (0.d0,0.d0)

  if(i==1) then
     do j=1,4
        do k=1,4
           rots_pq_i = rots_pq_i +conjg(dpsi_p(2,j))*Sigma(3,j,k)*psi_q(k)& 
                               & +conjg(psi_p(j))*Sigma(3,j,k)*dpsi_q(2,k)&
                               & -conjg(dpsi_p(3,j))*Sigma(2,j,k)*psi_q(k)&
                               & -conjg(psi_p(j))*Sigma(2,j,k)*dpsi_q(3,k)
        end do
     end do
  else if(i==2) then
     do j=1,4
        do k=1,4
           rots_pq_i = rots_pq_i +conjg(dpsi_p(3,j))*Sigma(1,j,k)*psi_q(k)& 
                               & +conjg(psi_p(j))*Sigma(1,j,k)*dpsi_q(3,k)&
                               & -conjg(dpsi_p(1,j))*Sigma(3,j,k)*psi_q(k)&
                               & -conjg(psi_p(j))*Sigma(3,j,k)*dpsi_q(1,k)
        end do
     end do
  else if(i==3) then
     do j=1,4
        do k=1,4
           rots_pq_i = rots_pq_i +conjg(dpsi_p(1,j))*Sigma(2,j,k)*psi_q(k)& 
                               & +conjg(psi_p(j))*Sigma(2,j,k)*dpsi_q(1,k)&
                               & -conjg(dpsi_p(2,j))*Sigma(1,j,k)*psi_q(k)&
                               & -conjg(psi_p(j))*Sigma(1,j,k)*dpsi_q(2,k)
        end do
     end do
  else
     write(*,*) "i should be 1-3 in calc_rots_pq."
     stop
  end if

  rots_pq_i = rots_pq_i/2.d0
  return
end subroutine calc_rots_pq


!============================================================
function js_Qmat(i,j,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
! spin current = -(i/2)  \psi^\dagger \Sigma^i \partial_j \psi
! non Hermite
!============================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=8) :: js_Qmat
  integer,intent(in) :: i,j
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=8) :: js_mat
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  
  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  js_Qmat = js_mat(i,j,x,y,z,n,a,m,b)

  return
end function js_Qmat

!============================================================
function js_mat(i,j,x,y,z,p,a,q,b)
! i-th component of j(ab)_pq , i=1-3
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
! spin current = -(i/2) \psi^\dagger \Sigma^i \partial_j \psi
! non Hermite
!============================================================
  use DefineTypes
  use DiracOutput
  implicit none
 
  complex(kind=8) :: js_mat
  integer,intent(in) :: i,j
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=8) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in js_mat."
     stop
  end if
  if(j.lt.1 .or. j.gt.3) then
     write(*,*) "j should be 1-3 in js_mat."
     stop
  end if

  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  call calc_js_pq(i,j,x,y,z,NBS_L,NBS_S,pg,c_p,c_q,js_mat)
  
  return
end function js_mat

!============================================================
subroutine calc_js_pq(i,j,x,y,z,NL,NS,pg,c_p,c_q,js_pq_ij)
! calculate ij-th component of js_pq
! spin current = -(i/2)  \psi^\dagger \Sigma^i \partial_j \psi
!============================================================
  use DefineTypes
  use Constants
  implicit none
  
  integer,intent(in) :: i,j   !js^ij
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=8),intent(out) :: js_pq_ij

  complex(kind=8) :: psi_p(4),psi_q(4)
  complex(kind=8) :: dpsi_p(3,4),dpsi_q(3,4)  ! 3:x,y,z
  integer :: k,l

  call calc_psi(x,y,z,NL,NS,pg,c_p,psi_p)
  call calc_psi(x,y,z,NL,NS,pg,c_q,psi_q)
  call calc_dpsi(x,y,z,NL,NS,pg,c_p,dpsi_p)
  call calc_dpsi(x,y,z,NL,NS,pg,c_q,dpsi_q)

  js_pq_ij = (0.d0,0.d0)

     do k=1,4
        do l=1,4
           js_pq_ij = js_pq_ij +conjg(psi_p(k))*Sigma(i,k,l)*dpsi_q(j,l)
        end do
     end do


  js_pq_ij = - IU*js_pq_ij*0.5d0
  return
end subroutine calc_js_pq


!============================================================
function chiral_Qmat(x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use DiracOutput
  implicit none

  complex(kind=8) :: chiral_Qmat
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=8) :: chiral_mat
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  
  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  chiral_Qmat = chiral_mat(x,y,z,n,a,m,b)

  return
end function chiral_Qmat

!============================================================
function chiral_ele(x,y,z,NEL)
! component of electron chiral density
!============================================================
  implicit none

  complex(kind=8) :: chiral_ele
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NEL

  complex(kind=8) :: chiral_mat
  complex(kind=8) :: sum
  integer :: j
  
  sum = (0.d0,0.d0)
  if(mod(NEL,2).eq.0) then
     do j=1,NEL/2  ! assume Kramers-Restricted
        sum = sum +2.d0*chiral_mat(x,y,z,j,"+",j,"+")  ! Kramers pair
     end do
  else
     do j=1,(NEL-1)/2  
        sum = sum +2.d0*chiral_mat(x,y,z,j,"+",j,"+")  ! Kramers pair
     end do
     j=(NEL+1)/2
     sum = sum +chiral_mat(x,y,z,j,"+",j,"+")  
  end if

  chiral_ele = sum
  return
end function chiral_ele

!============================================================
function chiral_mat(x,y,z,p,a,q,b)
! component of chiral(ab)_pq 
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use DefineTypes
  use DiracOutput
  implicit none
 
  complex(kind=8) :: chiral_mat
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=8) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  call calc_chiral_pq(x,y,z,NBS_L,NBS_S,pg,c_p,c_q,chiral_mat)
  
  return
end function chiral_mat

!============================================================
subroutine calc_chiral_pq(x,y,z,NL,NS,pg,c_p,c_q,chiral_pq)
! calculate i-th component of chiral_pq
! attention!!! Gam(1,m2,m3)*Sigma(1,m3,m4) = Gam(2,m2,m3)*Sigma(2,m3,m4) = Gam(3,m2,m3)*Sigma(3,m3,m4) (sum over m3) 
!============================================================
  use DefineTypes
  use Constants
  implicit none
  
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=8),intent(out) :: chiral_pq

  complex(kind=8) :: psi_p(4),psi_q(4)
  integer :: m1,m2,m3,m4 ! spinor indice
  integer :: j,k,l

  call calc_psi(x,y,z,NL,NS,pg,c_p,psi_p)
  call calc_psi(x,y,z,NL,NS,pg,c_q,psi_q)

  chiral_pq = (0.d0,0.d0)
  do m1=1,4
     do m2=1,4
        do m3=1,4
           do m4=1,4
              chiral_pq = chiral_pq &
                   +conjg(psi_p(m1))*Gam0(m1,m2)*Gam(1,m2,m3)*Sigma(1,m3,m4)*psi_q(m4)
           end do
        end do
     end do
  end do
  chiral_pq = (CCC/2.d0)*chiral_pq
  return
end subroutine calc_chiral_pq

!=========================================================================
subroutine calc_AA_pq(ip,iq,posRA,NL,NS,pg,c_p,c_q,AA_pq)
! Int (psi^+)_ip (psi)_iq /|r-R|
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none

  integer,intent(in) :: ip,iq ! spinor indice
  real(kind=dp),intent(in) :: posRA(3)  ! position of nucleus
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=dp),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  
  complex(kind=dp),intent(out) :: AA_pq
  
  integer :: NLS
  integer :: i,j,k
  real(kind=dp) :: overlap,ef(3) ! --> dummy variables here 
  real(kind=dp) :: nucatt
  real(kind=dp) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=dp) :: alphaA, alphaB ! exponents
  integer :: nx,nbarx,ny,nbary,nz,nbarz
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg


  AA_pq = (0._dp,0._dp)

  do i=1,NLS(ip,NL,NS)
     do j=1,NLS(iq,NL,NS)
        call set_pg(ip,i,NL,NS,pg,posA,alphaA,nx,ny,nz)
        call set_pg(iq,j,NL,NS,pg,posB,alphaB,nbarx,nbary,nbarz)
        call gauss_int(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,posRA,overlap,nucatt,ef)
        f_norm_pg = norm_pg(alphaA,nx,ny,nz)*norm_pg(alphaB,nbarx,nbary,nbarz)
        AA_pq = AA_pq + conjg(c_p(ip,i))*c_q(iq,j)*nucatt*f_norm_pg
     end do
  end do
  return
end subroutine calc_AA_pq

!=========================================================================
subroutine calc_AA_pq2(ip,iq,NL,NS,pg,c_p,c_q,NBS0,nucatt_pg,AA_pq)
! Int (psi^+)_ip (psi)_iq /|r-R|
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none

  integer,intent(in) :: ip,iq ! spinor indice
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=dp),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  real(kind=dp),intent(in) :: nucatt_pg(4,4,NBS0,NBS0)
  integer,intent(in) :: NBS0
  complex(kind=dp),intent(out) :: AA_pq
  
  integer :: NLS
  integer :: i,j,k


  AA_pq = (0._dp,0._dp)

  do i=1,NLS(ip,NL,NS)
     do j=1,NLS(iq,NL,NS)
        AA_pq = AA_pq + conjg(c_p(ip,i))*c_q(iq,j)*nucatt_pg(ip,iq,i,j)
     end do
  end do
  return
end subroutine calc_AA_pq2

!============================================================
subroutine save_gauss_int(ip,iq,i,j,posRA,NL,NS,pg,overlap_normpg,nucatt_normpg,ef_normpg)
!============================================================
  use Precision
  use DefineTypes
  implicit none

  integer,intent(in) :: ip,iq ! spinor indice
  integer,intent(in) :: i,j
  real(kind=dp),intent(in) :: posRA(3)  ! position of nucleus
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  
  real(kind=dp),intent(out) :: overlap_normpg,ef_normpg(3) ! --> dummy variables here 
  real(kind=dp),intent(out) :: nucatt_normpg

  integer :: NLS
  integer :: k
  real(kind=dp) :: overlap,ef(3) ! --> dummy variables here 
  real(kind=dp) :: nucatt
  real(kind=dp) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=dp) :: alphaA, alphaB ! exponents
  integer :: nx,nbarx,ny,nbary,nz,nbarz
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg

        call set_pg(ip,i,NL,NS,pg,posA,alphaA,nx,ny,nz)
        call set_pg(iq,j,NL,NS,pg,posB,alphaB,nbarx,nbary,nbarz)
        call gauss_int(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,posRA,overlap,nucatt,ef)
        f_norm_pg = norm_pg(alphaA,nx,ny,nz)*norm_pg(alphaB,nbarx,nbary,nbarz)
        overlap_normpg = overlap*f_norm_pg
        nucatt_normpg = nucatt*f_norm_pg
        do k=1,3
          ef_normpg(k) = ef(k)*f_norm_pg
        end do
  return
end subroutine save_gauss_int

!============================================================
subroutine calc_spol_pq(i,x,y,z,NL,NS,pg,c_p,c_q,spol_pq_i)
! calculate i-th component of spol_pq
! (1/2) psi gam0 sigma^k psi
!============================================================
  use DefineTypes
  use Constants
  implicit none
  
  integer,intent(in) :: i   !s^i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=8),intent(out) :: spol_pq_i

  complex(kind=8) :: psi_p(4),psi_q(4)
  integer :: j,k,l

  call calc_psi(x,y,z,NL,NS,pg,c_p,psi_p)
  call calc_psi(x,y,z,NL,NS,pg,c_q,psi_q)

  spol_pq_i = (0.d0,0.d0)
  do j=1,4
     do k=1,4
        do l=1,4
           spol_pq_i = spol_pq_i +conjg(psi_p(j))*Gam0(j,k)*Sigma(i,k,l)*psi_q(l)
        end do
     end do
  end do
  spol_pq_i = spol_pq_i/2.d0
  return
end subroutine calc_spol_pq

!============================================================
subroutine calc_gamk_pq(i,x,y,z,NL,NS,pg,c_p,c_q,gamk_pq_i)
! calculate i-th component of gamk_pq
! (1/2) psi gam0 sigma^k psi
!============================================================
  use DefineTypes
  use Constants
  implicit none
  
  integer,intent(in) :: i   !s^i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=8),intent(out) :: gamk_pq_i

  complex(kind=8) :: psi_p(4),psi_q(4)
  integer :: j,k,l

  call calc_psi(x,y,z,NL,NS,pg,c_p,psi_p)
  call calc_psi(x,y,z,NL,NS,pg,c_q,psi_q)

  gamk_pq_i = (0.d0,0.d0)
  do j=1,4
     do k=1,4
        gamk_pq_i = gamk_pq_i +conjg(psi_p(j))*Gam(i,j,k)*psi_q(k)
     end do
  end do
  gamk_pq_i = gamk_pq_i
  return
end subroutine calc_gamk_pq


!============================================================
subroutine calc_d2psi(x,y,z,NL,NS,pg,c,d2psi)
! output : d2psi/dxdx,d2psidxdy,d2psidxdz,...
!============================================================
  use DefineTypes
  implicit none
  
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(in) :: c(4,NMAX_PG)  ! NMAX_PG is defined in DefineTypes
  
  complex(kind=8),intent(out) :: d2psi(3,3,4) ! d_k d_l psi

  integer :: nx,ny,nz
  real(kind=8) :: posA(3),alpha
  real(kind=8) :: func_d2pg  ! derivative of primitive gaussian function (normalized)
  integer :: j,k,l

  d2psi(1:3,1:3,1:4) = (0.d0,0.d0)
  do l=1,3
     do k=1,3
        do j=1,NL
           posA(1)=pg%xL(j); posA(2)=pg%yL(j); posA(3)=pg%zL(j); alpha=pg%aL(j); nx=pg%nxL(j); ny=pg%nyL(j); nz=pg%nzL(j)
           d2psi(k,l,1) = d2psi(k,l,1) +c(1,j)*func_d2pg(k,l,x,y,z,posA,alpha,nx,ny,nz)
           d2psi(k,l,2) = d2psi(k,l,2) +c(2,j)*func_d2pg(k,l,x,y,z,posA,alpha,nx,ny,nz)
        end do
        do j=1,NS
           posA(1)=pg%xS(j); posA(2)=pg%yS(j); posA(3)=pg%zS(j); alpha=pg%aS(j); nx=pg%nxS(j); ny=pg%nyS(j); nz=pg%nzS(j)
           d2psi(k,l,3) = d2psi(k,l,3) +c(3,j)*func_d2pg(k,l,x,y,z,posA,alpha,nx,ny,nz)
           d2psi(k,l,4) = d2psi(k,l,4) +c(4,j)*func_d2pg(k,l,x,y,z,posA,alpha,nx,ny,nz)
        end do
     end do
  end do
  return
end subroutine calc_d2psi
