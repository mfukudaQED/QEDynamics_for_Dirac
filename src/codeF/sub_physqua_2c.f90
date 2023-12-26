!============================================================
function zetapot2c_pt_Qmat(x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef
  implicit none

  complex(kind=8) :: zetapot2c_pt_Qmat
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg
  
  call copy_DiracOutput_pg(pg)
  call calc_zetapot2c_pq(x,y,z,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),zetapot2c_pt_Qmat)

  return

end function zetapot2c_pt_Qmat

!============================================================
subroutine calc_zetapot2c_pq(x,y,z,NL,NS,pg,c_p,c_q,zetapot2c_pq)
! calculate zeta potential 2c ver.
!============================================================
  use Precision
  use DefineTypes
  use Constants
  implicit none
  
  real(kind=dp),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=dp),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=dp),intent(out) :: zetapot2c_pq

  complex(kind=dp) :: psi_p(4),psi_q(4)
  complex(kind=dp) :: dpsi_p(3,4),dpsi_q(3,4)  ! 3:x,y,z

  call calc_psi(x,y,z,NL,NS,pg,c_p,psi_p)
  call calc_psi(x,y,z,NL,NS,pg,c_q,psi_q)
  call calc_dpsi(x,y,z,NL,NS,pg,c_p,dpsi_p)
  call calc_dpsi(x,y,z,NL,NS,pg,c_q,dpsi_q)

  zetapot2c_pq = (0.d0,0.d0)

  zetapot2c_pq = zetapot2c_pq +    conjg(dpsi_p(1,2))*psi_q(1) +    conjg(dpsi_p(1,1))*psi_q(2)& 
                            & + IU*conjg(dpsi_p(2,2))*psi_q(1) - IU*conjg(dpsi_p(2,1))*psi_q(2)&
                            & +    conjg(dpsi_p(3,1))*psi_q(1) -    conjg(dpsi_p(3,2))*psi_q(2)&
                            & -    conjg(psi_q(2))*dpsi_p(1,1) -    conjg(psi_q(1))*dpsi_p(1,2)&
                            & - IU*conjg(psi_q(2))*dpsi_p(2,1) + IU*conjg(psi_q(1))*dpsi_p(2,2)&
                            & -    conjg(psi_q(1))*dpsi_p(3,1) +    conjg(psi_q(2))*dpsi_p(3,2)
  zetapot2c_pq = zetapot2c_pq*IU*0.25_dp

  return
end subroutine calc_zetapot2c_pq

!============================================================
function paudd2c_pt_Qmat(k,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef
  implicit none

  integer,intent(in) :: k !1~3
  complex(kind=8) :: paudd2c_pt_Qmat
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg
  
  call copy_DiracOutput_pg(pg)
  call calc_paudd2c_pq(k,x,y,z,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),paudd2c_pt_Qmat)

  return

end function paudd2c_pt_Qmat

!============================================================
subroutine calc_paudd2c_pq(k,x,y,z,NL,NS,pg,c_p,c_q,paudd2c_pq)
! calculate (ihbar^2/4/m) [ psi_L^+ pauli^k d_m d_m psi_L - (d_m d_m psi_L)^+ pauli^k psi_L]
!============================================================
  use Precision
  use DefineTypes
  use Constants
  implicit none
  
  integer,intent(in) :: k !1~3
  real(kind=dp),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=dp),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=dp),intent(out) :: paudd2c_pq

  complex(kind=dp) :: psi_p(4),psi_q(4)
  complex(kind=dp) :: d2psi_p(3,3,4),d2psi_q(3,3,4)  ! 3:x,y,z
  complex(kind=dp) :: sumd2psi_p(4), sumd2psi_q(4)
  integer :: i,j

  call calc_psi(x,y,z,NL,NS,pg,c_p,psi_p)
  call calc_psi(x,y,z,NL,NS,pg,c_q,psi_q)
  call calc_d2psi(x,y,z,NL,NS,pg,c_p,d2psi_p)
  call calc_d2psi(x,y,z,NL,NS,pg,c_q,d2psi_q)

  sumd2psi_p(1:4) = (0.d0,0.d0)
  sumd2psi_q(1:4) = (0.d0,0.d0)
  do i=1,4
     do j=1,3
        sumd2psi_p(i) = sumd2psi_p(i) + d2psi_p(j,j,i)
        sumd2psi_q(i) = sumd2psi_q(i) + d2psi_q(j,j,i)
     end do
  end do

  paudd2c_pq = (0.d0,0.d0)

  if (k.eq.1) then
     paudd2c_pq =      conjg(psi_q(2))*sumd2psi_p(1) +    conjg(psi_q(1))*sumd2psi_p(2)&
                & -    conjg(sumd2psi_p(2))*psi_q(1) -    conjg(sumd2psi_p(1))*psi_q(2) 
  elseif(k.eq.2) then
     paudd2c_pq =   IU*conjg(psi_q(2))*sumd2psi_p(1) - IU*conjg(psi_q(1))*sumd2psi_p(2)&
                & - IU*conjg(sumd2psi_p(2))*psi_q(1) + IU*conjg(sumd2psi_p(1))*psi_q(2)
  elseif(k.eq.3) then
     paudd2c_pq =      conjg(psi_q(1))*sumd2psi_p(1) -    conjg(psi_q(2))*sumd2psi_p(2)&
                & -    conjg(sumd2psi_p(1))*psi_q(1) +    conjg(sumd2psi_p(2))*psi_q(2)
  end if
                            
  paudd2c_pq = paudd2c_pq*IU*0.25_dp

  return
end subroutine calc_paudd2c_pq

!============================================================
function intN2c_pt_Qmat(nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef
  implicit none
 
  complex(kind=dp) :: intN2c_pt_Qmat
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: intN2c_pq
  integer :: i
  
  call copy_DiracOutput_pg(pg)
  
  intN2c_pt_Qmat = (0._dp,0._dp)
  do i=1,2
     call calc_intN2c_pq(i,i,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intN2c_pq) 
!     write(*,*) i,intN2c_pq
     intN2c_pt_Qmat = intN2c_pt_Qmat +intN2c_pq
  end do
  return
end function intN2c_pt_Qmat
!=========================================================================
subroutine calc_intN2c_pq(ip,iq,NL,NS,pg,c_p,c_q,intN2c_pq)
! int (psi^+)_ip (psi)_iq 
!=========================================================================  
  use Precision
  use DefineTypes
  use Constants
  implicit none

  integer,intent(in) :: ip,iq ! spinor indice
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=dp),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  
  complex(kind=dp),intent(out) :: intN2c_pq
  
  integer :: NLS
  integer :: i,j
  real(kind=dp) :: overlap,ene_kin
  real(kind=dp) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=dp) :: alphaA, alphaB ! exponents
  integer :: nx,nbarx,ny,nbary,nz,nbarz
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg

  intN2c_pq = (0._dp,0._dp)

  do i=1,NL
     do j=1,NL
        call set_pg(ip,i,NL,NS,pg,posA,alphaA,nx,ny,nz)
        call set_pg(iq,j,NL,NS,pg,posB,alphaB,nbarx,nbary,nbarz)
        call gauss_int_overlap(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,overlap)
        call gauss_int_KE(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,ene_kin)
        f_norm_pg = norm_pg(alphaA,nx,ny,nz)*norm_pg(alphaB,nbarx,nbary,nbarz)
        intN2c_pq = intN2c_pq + conjg(c_p(ip,i))*c_q(iq,j)*overlap*f_norm_pg
        intN2c_pq = intN2c_pq + conjg(c_p(ip,i))*c_q(iq,j)*ene_kin*f_norm_pg*0.5_dp/CCC/CCC
     end do
  end do
  return
end subroutine calc_intN2c_pq
