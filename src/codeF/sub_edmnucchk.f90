! Last Change:08-Sep-2015.

!============================================================
function intEeff_EDM_nuc2_pt_Qmat(nn,mm)
! nn,mm : labels for molecular orbitals (KP) . 1~2*NBS
! - (\gamma^0 -I) \Sigma^k E^k
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use PTcoef, only : c_psi
  implicit none

  integer,intent(in) :: nn,mm
  integer :: i,j,k,l ! p.g. indice
  integer :: NLS !function
  integer :: NL,NS
  integer :: a,b ! spinor indices 1,2->L, 3,4->S

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: tmp, tmptot
  complex(kind=dp) :: intEeff_EDM_nuc2_pt_Qmat
  real(kind=dp) :: posR(3)  
  complex(kind=dp) :: intE_pg !function

  call set_GammaMatrix
  call copy_DiracOutput_pg(pg)
  NL=NBS_L; NS=NBS_S

    intEeff_EDM_nuc2_pt_Qmat=(0._dp,0._dp)
    tmp=(0._dp,0._dp)
  do k=1,NAT !number of atoms
    posR(1)=xc(k) ! position of nucleus
    posR(2)=yc(k)
    posR(3)=zc(k)
!    write(*,*)'posR',k, posR(1), posR(2), posR(3)
    tmptot=(0._dp,0._dp)
    !m=1
!!$    !LL
!!$    a=1
!!$    do i=1,NLS(a,NL,NS)
!!$       do j=1,NLS(a,NL,NS)
!!$         tmp =  conjg(c_psi(2,i,nn))*c_psi(1,j,mm)& 
!!$              &+conjg(c_psi(1,i,nn))*c_psi(2,j,mm)
!!$         tmptot = tmptot + tmp*intE_pg(1,posR,i,j,a,a,NL,NS,pg)
!!$       end do
!!$    end do

    !SS
    a=3
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
         tmp =  conjg(c_psi(4,i,nn))*c_psi(3,j,mm)& 
              &+conjg(c_psi(3,i,nn))*c_psi(4,j,mm)
         tmptot = tmptot + tmp*intE_pg(1,posR,i,j,a,a,NL,NS,pg)
       end do
    end do

    !m=2
!!$    !LL
!!$    a=1
!!$    do i=1,NLS(a,NL,NS)
!!$       do j=1,NLS(a,NL,NS)
!!$         tmp =  IU*conjg(c_psi(2,i,nn))*c_psi(1,j,mm)& 
!!$              &-IU*conjg(c_psi(1,i,nn))*c_psi(2,j,mm)
!!$         tmptot = tmptot + tmp*intE_pg(2,posR,i,j,a,a,NL,NS,pg)
!!$       end do
!!$    end do

    !SS
    a=3
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
         tmp =  IU*conjg(c_psi(4,i,nn))*c_psi(3,j,mm)& 
              &-IU*conjg(c_psi(3,i,nn))*c_psi(4,j,mm)
         tmptot = tmptot + tmp*intE_pg(2,posR,i,j,a,a,NL,NS,pg)
       end do
    end do

!!$    !m=3
!!$    !LL
!!$    a=1
!!$    do i=1,NLS(a,NL,NS)
!!$       do j=1,NLS(a,NL,NS)
!!$         tmp =  conjg(c_psi(1,i,nn))*c_psi(1,j,mm)& 
!!$              &-conjg(c_psi(2,i,nn))*c_psi(2,j,mm)
!!$         tmptot = tmptot + tmp*intE_pg(3,posR,i,j,a,a,NL,NS,pg)
!!$       end do
!!$    end do

    !SS
    a=3
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
         tmp = +conjg(c_psi(3,i,nn))*c_psi(3,j,mm)& 
              &-conjg(c_psi(4,i,nn))*c_psi(4,j,mm)
         tmptot = tmptot + tmp*intE_pg(3,posR,i,j,a,a,NL,NS,pg)
       end do
    end do

    intEeff_EDM_nuc2_pt_Qmat = intEeff_EDM_nuc2_pt_Qmat + cn(k)*tmptot ! cn -> nuclear charge
!    write(*,*)'NAT=',k,cn(k), nn,mm, tmptot
  end do !k=1,NAT

  intEeff_EDM_nuc2_pt_Qmat = - intEeff_EDM_nuc2_pt_Qmat ! minus for intE_pg , added in 20140204
  intEeff_EDM_nuc2_pt_Qmat = - 2._dp*intEeff_EDM_nuc2_pt_Qmat !minus is added in 20140207
    
    return
end function intEeff_EDM_nuc2_pt_Qmat


!============================================================
function intEeff_EDM_nuc3_pt_Qmat(nn,mm)
! nn,mm : labels for molecular orbitals (KP) . 1~2*NBS
! \int \psi \Sigma^k E^k \psi d^3r
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use PTcoef, only : c_psi
  implicit none

  integer,intent(in) :: nn,mm
  integer :: i,j,k,l ! p.g. indice
  integer :: NLS !function
  integer :: NL,NS
  integer :: a,b ! spinor indices 1,2->L, 3,4->S

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: tmp, tmptot
  complex(kind=dp) :: intEeff_EDM_nuc3_pt_Qmat
  real(kind=dp) :: posR(3)  
  complex(kind=dp) :: intE_pg !function

  call set_GammaMatrix
  call copy_DiracOutput_pg(pg)
  NL=NBS_L; NS=NBS_S

    intEeff_EDM_nuc3_pt_Qmat=(0._dp,0._dp)
    tmp=(0._dp,0._dp)
  do k=1,NAT !number of atoms
    posR(1)=xc(k) ! position of nucleus
    posR(2)=yc(k)
    posR(3)=zc(k)
!    write(*,*)'posR',k, posR(1), posR(2), posR(3)
    tmptot=(0._dp,0._dp)
    !m=1
    !LL
    a=1
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
         tmp =  conjg(c_psi(2,i,nn))*c_psi(1,j,mm)& 
              &+conjg(c_psi(1,i,nn))*c_psi(2,j,mm)
         tmptot = tmptot + tmp*intE_pg(1,posR,i,j,a,a,NL,NS,pg)
       end do
    end do

    !SS
    a=3
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
         tmp = -conjg(c_psi(4,i,nn))*c_psi(3,j,mm)& 
              &-conjg(c_psi(3,i,nn))*c_psi(4,j,mm)
         tmptot = tmptot - tmp*intE_pg(1,posR,i,j,a,a,NL,NS,pg)
       end do
    end do

    !m=2
    !LL
    a=1
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
         tmp =  IU*conjg(c_psi(2,i,nn))*c_psi(1,j,mm)& 
              &-IU*conjg(c_psi(1,i,nn))*c_psi(2,j,mm)
         tmptot = tmptot + tmp*intE_pg(2,posR,i,j,a,a,NL,NS,pg)
       end do
    end do

    !SS
    a=3
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
         tmp = -IU*conjg(c_psi(4,i,nn))*c_psi(3,j,mm)& 
              &+IU*conjg(c_psi(3,i,nn))*c_psi(4,j,mm)
         tmptot = tmptot - tmp*intE_pg(2,posR,i,j,a,a,NL,NS,pg)
       end do
    end do

    !m=3
    !LL
    a=1
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
         tmp =  conjg(c_psi(1,i,nn))*c_psi(1,j,mm)& 
              &-conjg(c_psi(2,i,nn))*c_psi(2,j,mm)
         tmptot = tmptot + tmp*intE_pg(3,posR,i,j,a,a,NL,NS,pg)
       end do
    end do

    !SS
    a=3
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
         tmp = -conjg(c_psi(3,i,nn))*c_psi(3,j,mm)& 
              &+conjg(c_psi(4,i,nn))*c_psi(4,j,mm)
         tmptot = tmptot - tmp*intE_pg(3,posR,i,j,a,a,NL,NS,pg)
       end do
    end do

    intEeff_EDM_nuc3_pt_Qmat = intEeff_EDM_nuc3_pt_Qmat + cn(k)*tmptot ! cn -> nuclear charge
!    write(*,*)'NAT=',k,cn(k), nn,mm, tmptot
  end do !k=1,NAT

  intEeff_EDM_nuc3_pt_Qmat = - intEeff_EDM_nuc3_pt_Qmat
    
    return
end function intEeff_EDM_nuc3_pt_Qmat

