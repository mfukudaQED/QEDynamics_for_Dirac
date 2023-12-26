!==========================================================================
! SUBROUTINE
!==========================================================================
SUBROUTINE SET_PSI_LOCAL(actorb,occ,coef)
!==========================================================================
  use Precision
  use DiracOutput
  use prop_param
  use EDM_calculation
  implicit none

  integer,intent(in) :: actorb
  real(kind=dp),intent(in) :: occ(actorb)
  complex(kind=dp),intent(in) :: coef(NBS0,2*NBS,4)

  real(kind=dp), allocatable,dimension(:,:) :: gL,gS
  real(kind=dp), allocatable,dimension(:,:,:) :: dgL,dgS
  real(kind=dp), allocatable,dimension(:,:,:) :: d2gL,d2gS
  complex(kind=dp), allocatable,dimension(:,:,:) :: psi_nat
  complex(kind=dp), allocatable,dimension(:,:,:,:) :: dpsi_nat
  integer :: i,j,k,l1,l2,flag
  character(len=12) :: get_time

  k = 1
  flag = 0
  !$omp parallel do private(flag,i,j,l1,l2,gL,gS,dgL,dgS,d2gL,d2gS,psi_nat,dpsi_nat)
  do i = 1,meshx
    if (flag==0) then
      allocate(gL(meshz,NBS_L))
      allocate(gS(meshz,NBS_S))
      allocate(dgL(meshz,NBS_L,3))
      allocate(dgS(meshz,NBS_S,3))
      allocate(d2gL(meshz,NBS_L,6))
      allocate(d2gS(meshz,NBS_S,6))
      allocate(psi_nat(meshz,actorb,4))
      allocate(dpsi_nat(meshz,actorb,4,3))
      flag = 1
    end if
    do j = 1,meshy

      l1 = (i -1) *meshyz +(j -1)*meshz +1
      l2 = l1 +meshz -1

      call calc_gaussian_local(i,j,gL,gS,dgL,dgS,d2gL,d2gS)
      call calc_psi_local(l1,l2,actorb,coef,occ,gL,gS,psi_nat)
      call calc_dpsi_local(l1,l2,actorb,coef,occ,dgL,dgS,psi_nat,dpsi_nat)
      call calc_d2psi_local(l1,l2,actorb,coef,occ,d2gL,d2gS,psi_nat,dpsi_nat)

    end do
    write(*,'(i10,a,i10,2a)') k, ' of ', meshx, ' completed at : ', get_time()
    k = k + 1
  end do
  !$omp end parallel do

  return
END SUBROUTINE SET_PSI_LOCAL

!==========================================================================
SUBROUTINE SET_PG_ORBITAL_ANGULAR_MOMENTUM
!==========================================================================
  use Precision
  use DiracOutput
  use prop_param
  use EDM_calculation
  implicit none

  integer :: i,n

  xA(1:meshx,-2:-1,1:NAT) = 0._dp
  yA(1:meshy,-2:-1,1:NAT) = 0._dp
  zA(1:meshz,-2:-1,1:NAT) = 0._dp
  xA(1:meshx,0,1:NAT) = 1._dp
  yA(1:meshy,0,1:NAT) = 1._dp
  zA(1:meshz,0,1:NAT) = 1._dp

  do i = 1,NAT
    xA(1:meshx,1,i) = x1(1:meshx)-xc(i)
    yA(1:meshy,1,i) = x2(1:meshy)-yc(i)
    zA(1:meshz,1,i) = x3(1:meshz)-zc(i)
    xA(1:meshx,2,i) = xA(1:meshx,1,i)**2
    yA(1:meshy,2,i) = yA(1:meshy,1,i)**2
    zA(1:meshz,2,i) = zA(1:meshz,1,i)**2
    do n = 3,nmax(1)
      xA(1:meshx,n,i) = xA(1:meshx,n-1,i)*xA(1:meshx,1,i)
    end do
    do n = 3,nmax(2)
      yA(1:meshy,n,i) = yA(1:meshy,n-1,i)*yA(1:meshy,1,i)
    end do
    do n = 3,nmax(3)
      zA(1:meshz,n,i) = zA(1:meshz,n-1,i)*zA(1:meshz,1,i)
    end do
  end do

END SUBROUTINE SET_PG_ORBITAL_ANGULAR_MOMENTUM

!==========================================================================
SUBROUTINE CALC_PG_COMPONENT(i,j,k,inat,alpha,nx,ny,nz,exp_rA2,d_x,d_y,d_z)
!==========================================================================
  use Precision
  use DiracOutput
  use EDM_calculation
  implicit none

  integer,intent(in) :: i,j,k,inat
  real(kind=dp),intent(in) :: alpha
  integer,intent(in) :: nx,ny,nz

  real(kind=dp),intent(out) :: exp_rA2,d_x,d_y,d_z

  exp_rA2 = exp(-alpha*(xA(i,2,inat)+yA(j,2,inat)+zA(k,2,inat)))
  d_x = nx*xA(i,nx-1,inat) -2*alpha*xA(i,nx+1,inat)
  d_y = ny*yA(j,ny-1,inat) -2*alpha*yA(j,ny+1,inat)
  d_z = nz*zA(k,nz-1,inat) -2*alpha*zA(k,nz+1,inat)

END SUBROUTINE CALC_PG_COMPONENT

!==========================================================================
SUBROUTINE CALC_D2PG_COMPONENT(i,j,k,inat,alpha,nx,ny,nz,d_xx,d_yy,d_zz)
!==========================================================================
  use Precision
  use DiracOutput
  use EDM_calculation
  implicit none

  integer,intent(in) :: i,j,k,inat
  real(kind=dp),intent(in) :: alpha
  integer,intent(in) :: nx,ny,nz

  real(kind=dp),intent(out) :: d_xx,d_yy,d_zz

  d_xx = nx*(nx-1)*xA(i,nx-2,inat) -2*alpha*(2*nx+1)*xA(i,nx,inat) +4*alpha**2*xA(i,nx+2,inat)
  d_yy = ny*(ny-1)*yA(j,ny-2,inat) -2*alpha*(2*ny+1)*yA(j,ny,inat) +4*alpha**2*yA(j,ny+2,inat)
  d_zz = nz*(nz-1)*zA(k,nz-2,inat) -2*alpha*(2*nz+1)*zA(k,nz,inat) +4*alpha**2*zA(k,nz+2,inat)

END SUBROUTINE CALC_D2PG_COMPONENT

!==========================================================================
SUBROUTINE CALC_GAUSSIAN_LOCAL(i,j,gL,gS,dgL,dgS,d2gL,d2gS)
!==========================================================================
  use Precision
  use DiracOutput
  use prop_param
  use EDM_calculation
  implicit none

  integer,intent(in) :: i,j
  real(kind=dp),intent(out) :: gL(meshz,NBS_L),gS(meshz,NBS_S)
  real(kind=dp),intent(out) :: dgL(meshz,NBS_L,3),dgS(meshz,NBS_S,3)
  real(kind=dp),intent(out) :: d2gL(meshz,NBS_L,6),d2gS(meshz,NBS_S,6)

  real(kind=dp) :: exp_rA2,d_x,d_y,d_z,d_xx,d_yy,d_zz
  integer :: k,m

  do m = 1,NBS_L
    do k = 1,meshz
      call calc_pg_component(i,j,k,atnum_L(m),aa_L(m),nx_L(m),ny_L(m),nz_L(m),exp_rA2,d_x,d_y,d_z)
      gL (k,m) = c_Lnorm(m) *xA(i,nx_L(m),atnum_L(m))*yA(j,ny_L(m),atnum_L(m))*zA(k,nz_L(m),atnum_L(m))*exp_rA2
      dgL(k,m,1) = c_Lnorm(m) *d_x*yA(j,ny_L(m),atnum_L(m))*zA(k,nz_L(m),atnum_L(m))*exp_rA2
      dgL(k,m,2) = c_Lnorm(m) *xA(i,nx_L(m),atnum_L(m))*d_y*zA(k,nz_L(m),atnum_L(m))*exp_rA2
      dgL(k,m,3) = c_Lnorm(m) *xA(i,nx_L(m),atnum_L(m))*yA(j,ny_L(m),atnum_L(m))*d_z*exp_rA2
     call calc_d2pg_component(i,j,k,atnum_L(m),aa_L(m),nx_L(m),ny_L(m),nz_L(m),d_xx,d_yy,d_zz)
      d2gL(k,m,1) = c_Lnorm(m) *d_xx*yA(j,ny_L(m),atnum_L(m))*zA(k,nz_L(m),atnum_L(m))*exp_rA2 !d_xx
      d2gL(k,m,2) = c_Lnorm(m) *xA(i,nx_L(m),atnum_L(m))*d_yy*zA(k,nz_L(m),atnum_L(m))*exp_rA2 !d_yy
      d2gL(k,m,3) = c_Lnorm(m) *xA(i,nx_L(m),atnum_L(m))*yA(j,ny_L(m),atnum_L(m))*d_zz*exp_rA2 !d_zz
      d2gL(k,m,4) = c_Lnorm(m) *d_x*d_y*zA(k,nz_L(m),atnum_L(m))*exp_rA2 !d_x *d_y
      d2gL(k,m,5) = c_Lnorm(m) *xA(i,nx_L(m),atnum_L(m))*d_y*d_z*exp_rA2 !d_y *d_z
      d2gL(k,m,6) = c_Lnorm(m) *d_x*yA(j,ny_L(m),atnum_L(m))*d_z*exp_rA2 !d_z *d_x
    end do
  end do
  do m = 1,NBS_S
    do k = 1,meshz
      call calc_pg_component(i,j,k,atnum_S(m),aa_S(m),nx_S(m),ny_S(m),nz_S(m),exp_rA2,d_x,d_y,d_z)
      gS (k,m) = c_Snorm(m) *xA(i,nx_S(m),atnum_S(m))*yA(j,ny_S(m),atnum_S(m))*zA(k,nz_S(m),atnum_S(m))*exp_rA2
      dgS(k,m,1) = c_Snorm(m) *d_x*yA(j,ny_S(m),atnum_S(m))*zA(k,nz_S(m),atnum_S(m))*exp_rA2
      dgS(k,m,2) = c_Snorm(m) *xA(i,nx_S(m),atnum_S(m))*d_y*zA(k,nz_S(m),atnum_S(m))*exp_rA2
      dgS(k,m,3) = c_Snorm(m) *xA(i,nx_S(m),atnum_S(m))*yA(j,ny_S(m),atnum_S(m))*d_z*exp_rA2
      call calc_d2pg_component(i,j,k,atnum_S(m),aa_S(m),nx_S(m),ny_S(m),nz_S(m),d_xx,d_yy,d_zz)
      d2gS(k,m,1) = c_Snorm(m) *d_xx*yA(j,ny_S(m),atnum_S(m))*zA(k,nz_S(m),atnum_S(m))*exp_rA2 !d_xx
      d2gS(k,m,2) = c_Snorm(m) *xA(i,nx_S(m),atnum_S(m))*d_yy*zA(k,nz_S(m),atnum_S(m))*exp_rA2 !d_yy
      d2gS(k,m,3) = c_Snorm(m) *xA(i,nx_S(m),atnum_S(m))*yA(j,ny_S(m),atnum_S(m))*d_zz*exp_rA2 !d_zz
      d2gS(k,m,4) = c_Snorm(m) *d_x*d_y*zA(k,nz_S(m),atnum_S(m))*exp_rA2 !d_x *d_y
      d2gS(k,m,5) = c_Snorm(m) *xA(i,nx_S(m),atnum_S(m))*d_y*d_z*exp_rA2 !d_y *d_z
      d2gS(k,m,6) = c_Snorm(m) *d_x*yA(j,ny_S(m),atnum_S(m))*d_z*exp_rA2 !d_z *d_x
    end do
  end do

  return
END SUBROUTINE CALC_GAUSSIAN_LOCAL

!==========================================================================
SUBROUTINE CALC_PSI_LOCAL(l1,l2,actorb,coef,occ,gL,gS,psi_nat)
!==========================================================================
  use Precision
  use DiracOutput
  use prop_param
  use EDM_calculation
  implicit none

  integer,intent(in) :: l1,l2,actorb
  complex(kind=dp),intent(in) :: coef(NBS0,2*NBS,4)
  real(kind=dp),intent(in) :: occ(actorb)
  real(kind=dp),intent(in) :: gL(meshz,NBS_L),gS(meshz,NBS_S)
  complex(kind=dp),intent(out) :: psi_nat(meshz,actorb,4)

  integer :: a,b

!{\psi}_a
  do a = 1,2		!{\psi}_{L \alpha}, {\psi}_{L \beta}
    psi_nat(1:meshz,1:actorb,a) = matmul(gL(1:meshz,1:NBS_L), coef(1:NBS_L,1:actorb,a))
  end do
  do a = 3,4		!{\psi}_{S \alpha}, {\psi}_{S \beta}
    psi_nat(1:meshz,1:actorb,a) = matmul(gS(1:meshz,1:NBS_S), coef(1:NBS_S,1:actorb,a))
  end do

  do b = 1,4
    do a = 1,b
      psi_local(l1:l2,1:actorb,a,b) = conjg(psi_nat(1:meshz,1:actorb,a))*psi_nat(1:meshz,1:actorb,b)		! {\psi}_a^{\ast} {\psi}_b
      psi_local(l1:l2,0,a,b) = matmul(psi_local(l1:l2,1:actorb,a,b), occ(1:actorb))
    end do
  end do
  do b = 1,3
    do a = b+1,4
      psi_local(l1:l2,0:actorb,a,b) = conjg(psi_local(l1:l2,0:actorb,b,a))		! {\psi}_a^{\ast} {\psi}_b = ( {\psi}_b^{\ast} {\psi}_a )^{\ast}
    end do
  end do

  return
END SUBROUTINE CALC_PSI_LOCAL

!==========================================================================
SUBROUTINE CALC_DPSI_LOCAL(l1,l2,actorb,coef,occ,dgL,dgS,psi_nat,dpsi_nat)
!==========================================================================
  use Precision
  use DiracOutput
  use prop_param
  use EDM_calculation
  implicit none

  integer,intent(in) :: l1,l2,actorb
  complex(kind=dp),intent(in) :: coef(NBS0,2*NBS,4)
  real(kind=dp),intent(in) :: occ(actorb)
  real(kind=dp),intent(in) :: dgL(meshz,NBS_L,3),dgS(meshz,NBS_S,3)
  complex(kind=dp),intent(in) :: psi_nat(meshz,actorb,4)
  complex(kind=dp),intent(out) :: dpsi_nat(meshz,actorb,4,3)

  integer :: a,b,l

!{\partial}_l {\psi}_a
  do l = 1,3
    do a = 1,2		!{\partial}_l {\psi}_{L \alpha}, {\partial}_l {\psi}_{L \beta}
      dpsi_nat(1:meshz,1:actorb,a,l) = matmul(dgL(1:meshz,1:NBS_L,l), coef(1:NBS_L,1:actorb,a))
    end do
    do a = 3,4		!{\partial}_l {\psi}_{S \alpha}, {\partial}_l {\psi}_{S \beta}
      dpsi_nat(1:meshz,1:actorb,a,l) = matmul(dgS(1:meshz,1:NBS_S,l), coef(1:NBS_S,1:actorb,a))
    end do
  end do

  do l = 1,3
    do b = 1,4
      do a = 1,4
        dpsi_local(l1:l2,1:actorb,a,b,l) = conjg(psi_nat(1:meshz,1:actorb,a))*dpsi_nat(1:meshz,1:actorb,b,l)		!{\psi}_a^{\ast} {\partial}_l {\psi}_b
        dpsi_local(l1:l2,0,a,b,l) = matmul(dpsi_local(l1:l2,1:actorb,a,b,l), occ(1:actorb))
      end do
    end do
  end do

  return
END SUBROUTINE CALC_DPSI_LOCAL

!==========================================================================
SUBROUTINE CALC_D2PSI_LOCAL(l1,l2,actorb,coef,occ,d2gL,d2gS,psi_nat,dpsi_nat)
!==========================================================================
  use Precision
  use DiracOutput
  use prop_param
  use EDM_calculation
  implicit none

  integer,intent(in) :: l1,l2,actorb
  complex(kind=dp),intent(in) :: coef(NBS0,2*NBS,4)
  real(kind=dp),intent(in) :: occ(actorb)
  real(kind=dp),intent(in) :: d2gL(meshz,NBS_L,6),d2gS(meshz,NBS_S,6)

  complex(kind=dp),intent(in) :: psi_nat(meshz,actorb,4)
  complex(kind=dp),intent(in) :: dpsi_nat(meshz,actorb,4,3)
  complex(kind=dp) :: d2psi_nat(meshz,actorb,4,6)
  integer :: a,b,l,m

  do l = 1,6
    do a = 1,2
      d2psi_nat(1:meshz,1:actorb,a,l) = matmul(d2gL(1:meshz,1:NBS_L,l), coef(1:NBS_L,1:actorb,a))
    end do
    do a = 3,4
      d2psi_nat(1:meshz,1:actorb,a,l) = matmul(d2gS(1:meshz,1:NBS_S,l), coef(1:NBS_S,1:actorb,a))
    end do
  end do

  do b = 1,4
    do a = 1,4
      d2psi_local(l1:l2,1:actorb,a,b,1,1) = conjg(psi_nat(1:meshz,1:actorb,a))*d2psi_nat(1:meshz,1:actorb,b,1)		!{\psi}_a^{\ast} {\partial}_1 {\partial}_1 {\psi}_b
      d2psi_local(l1:l2,1:actorb,a,b,2,2) = conjg(psi_nat(1:meshz,1:actorb,a))*d2psi_nat(1:meshz,1:actorb,b,2)		!{\psi}_a^{\ast} {\partial}_2 {\partial}_2 {\psi}_b
      d2psi_local(l1:l2,1:actorb,a,b,3,3) = conjg(psi_nat(1:meshz,1:actorb,a))*d2psi_nat(1:meshz,1:actorb,b,3)		!{\psi}_a^{\ast} {\partial}_3 {\partial}_3 {\psi}_b
      d2psi_local(l1:l2,1:actorb,a,b,1,2) = conjg(psi_nat(1:meshz,1:actorb,a))*d2psi_nat(1:meshz,1:actorb,b,4)		!{\psi}_a^{\ast} {\partial}_1 {\partial}_2 {\psi}_b
      d2psi_local(l1:l2,1:actorb,a,b,2,3) = conjg(psi_nat(1:meshz,1:actorb,a))*d2psi_nat(1:meshz,1:actorb,b,5)		!{\psi}_a^{\ast} {\partial}_2 {\partial}_3 {\psi}_b
      d2psi_local(l1:l2,1:actorb,a,b,3,1) = conjg(psi_nat(1:meshz,1:actorb,a))*d2psi_nat(1:meshz,1:actorb,b,6)		!{\psi}_a^{\ast} {\partial}_3 {\partial}_1 {\psi}_b
      d2psi_local(l1:l2,1:actorb,a,b,2,1) = d2psi_local(l1:l2,1:actorb,a,b,1,2)		!{\psi}_a^{\ast} {\partial}_2 {\partial}_1 {\psi}_b = {\psi}_a^{\ast} {\partial}_1 {\partial}_2 {\psi}_b
      d2psi_local(l1:l2,1:actorb,a,b,3,2) = d2psi_local(l1:l2,1:actorb,a,b,2,3)		!{\psi}_a^{\ast} {\partial}_3 {\partial}_2 {\psi}_b = {\psi}_a^{\ast} {\partial}_2 {\partial}_3 {\psi}_b
      d2psi_local(l1:l2,1:actorb,a,b,1,3) = d2psi_local(l1:l2,1:actorb,a,b,3,1)		!{\psi}_a^{\ast} {\partial}_1 {\partial}_3 {\psi}_b = {\psi}_a^{\ast} {\partial}_3 {\partial}_1 {\psi}_b
    end do
  end do

  do m = 1,3
    do l = 1,3
      do b = 1,4
        do a = 1,4
          ddpsi_local(l1:l2,1:actorb,a,b,l,m) = conjg(dpsi_nat(1:meshz,1:actorb,a,l))*dpsi_nat(1:meshz,1:actorb,b,m)		!( {\partial}_l {\psi}_a^{\ast} ) {\partial}_m {\psi}_b 
        end do
      end do
    end do
  end do

  do m = 1,3
    do l = 1,3
      do b = 1,4
        do a = 1,4
          ddpsi_local(l1:l2,0,a,b,l,m) = matmul(ddpsi_local(l1:l2,1:actorb,a,b,l,m), occ(1:actorb))
          d2psi_local(l1:l2,0,a,b,l,m) = matmul(d2psi_local(l1:l2,1:actorb,a,b,l,m), occ(1:actorb))
        end do
      end do
    end do
  end do

  return
END SUBROUTINE CALC_D2PSI_LOCAL

!==========================================================================
SUBROUTINE CALC_DENSITY(orb)
!{\rho} = {\psi}^\dagger {\psi}
!       = {\psi}_{L \alpha}^\ast {\psi}_{L \alpha} + {\psi}_{L \beta}^\ast {\psi}_{L \beta} 
!       + {\psi}_{S \alpha}^\ast {\psi}_{S \alpha} + {\psi}_{S \beta}^\ast {\psi}_{S \beta}
!==========================================================================
  use Precision
  use prop_param ! mesh
  use EDM_calculation
  implicit none

  integer,intent(in) :: orb

  density(1:mesh) = psi_local(1:mesh,orb,1,1) +psi_local(1:mesh,orb,2,2)&
                  &+psi_local(1:mesh,orb,3,3) +psi_local(1:mesh,orb,4,4)

  return
END SUBROUTINE CALC_DENSITY

!==========================================================================
SUBROUTINE CALC_SPIN(orb)
!s_e^i = {\psi}^\dag {\hbar}/{2} \Sigma^i {\psi}
!      = {\hbar}/{2} ( {{\psi}_L}^\dag \sigma^i {\psi}_L} + {{\psi}_S}^\dag \sigma^i {{\psi}_S} )
!==========================================================================
  use Precision
  use constants
  use prop_param ! mesh
  use EDM_calculation
  implicit none

  integer,intent(in) :: orb

  spin_small(1:mesh,1) =  psi_local(1:mesh,orb,3,4) +psi_local(1:mesh,orb,4,3)
  spin_small(1:mesh,2) =  IU*(-psi_local(1:mesh,orb,3,4) +psi_local(1:mesh,orb,4,3))
  spin_small(1:mesh,3) =  psi_local(1:mesh,orb,3,3) -psi_local(1:mesh,orb,4,4)
  spin(1:mesh,1) =  psi_local(1:mesh,orb,1,2) +psi_local(1:mesh,orb,2,1) +spin_small(1:mesh,1)
  spin(1:mesh,2) =  IU*(-psi_local(1:mesh,orb,1,2) +psi_local(1:mesh,orb,2,1)) +spin_small(1:mesh,2)
  spin(1:mesh,3) =  psi_local(1:mesh,orb,1,1) -psi_local(1:mesh,orb,2,2) +spin_small(1:mesh,3)

  spin_small = 0.5_dp *spin_small
  spin = 0.5_dp *spin

  return
END SUBROUTINE CALC_SPIN

!==========================================================================
SUBROUTINE CALC_GRAD_N(orb)
!\partial_i {\rho}
!= \partial_i ( {\psi}^\dagger {\psi} )
!= {\psi}_{L \alpha}^\ast \partial_i {\psi}_{L \alpha} + {\psi}_{L \beta}^\ast \partial_i {\psi}_{L \beta} 
!+ {\psi}_{S \alpha}^\ast \partial_i {\psi}_{S \alpha} + {\psi}_{S \beta}^\ast \partial_i {\psi}_{S \beta} + c.c.
!==========================================================================
  use Precision
  use prop_param ! mesh
  use EDM_calculation
  implicit none

  integer,intent(in) :: orb

  grad_N(1:mesh,1) = dpsi_local(1:mesh,orb,1,1,1) +dpsi_local(1:mesh,orb,2,2,1)&
                   &+dpsi_local(1:mesh,orb,3,3,1) +dpsi_local(1:mesh,orb,4,4,1)
  grad_N(1:mesh,2) = dpsi_local(1:mesh,orb,1,1,2) +dpsi_local(1:mesh,orb,2,2,2)&
                   &+dpsi_local(1:mesh,orb,3,3,2) +dpsi_local(1:mesh,orb,4,4,2)
  grad_N(1:mesh,3) = dpsi_local(1:mesh,orb,1,1,3) +dpsi_local(1:mesh,orb,2,2,3)&
                   &+dpsi_local(1:mesh,orb,3,3,3) +dpsi_local(1:mesh,orb,4,4,3)

  grad_N = grad_N + conjg(grad_N)

  return
END SUBROUTINE CALC_GRAD_N

!==========================================================================
SUBROUTINE CALC_JE(orb)
!j_e^i = Z_eec {\psi}^\dagger \gamma^0 \gamma^i {\psi}
!	   = Z_eec {{\psi}_L}^\dag \sigma^i {{\psi}_S} + h.c.
!==========================================================================
  use Precision
  use Constants
  use prop_param ! mesh
  use EDM_calculation
  implicit none

  integer,intent(in) :: orb

  je(1:mesh,1) = psi_local(1:mesh,orb,1,4) +psi_local(1:mesh,orb,2,3)
  je(1:mesh,2) = IU*(-psi_local(1:mesh,orb,1,4) +psi_local(1:mesh,orb,2,3))
  je(1:mesh,3) = psi_local(1:mesh,orb,1,3) -psi_local(1:mesh,orb,2,4)

  je = Ze *CCC *je
  je = je + conjg(je)

  return
END SUBROUTINE CALC_JE

!==========================================================================
SUBROUTINE CALC_SPIN_TORQ(orb)
!t^i = - \epsilon_{ijk} {\tau}^{A jk} 
!    = - \epsilon_{ijk} {\tau}^{\Pi jk}
!{\tau}^{\Pi kl} = {i\hbar c}/{2} {\psi}^\dag \gamma^0 \gamma^l \partial_k {\psi} + h.c. 
!				 = {i\hbar c}/{2} \left[ {{\psi}_S}^\dag \sigma^l \partial_k {\psi}_L 
!				 + {{\psi}_L}^\dag \sigma^l \partial_k {\psi}_S \right] + h.c. 
!==========================================================================
  use Precision
  use Constants
  use prop_param ! mesh
  use EDM_calculation
  implicit none

  integer,intent(in) :: orb

  !{\tau}^{\Pi 11}
  tau(1:mesh,1) = dpsi_local(1:mesh,orb,1,4,1) +dpsi_local(1:mesh,orb,2,3,1)&
                &+dpsi_local(1:mesh,orb,3,2,1) +dpsi_local(1:mesh,orb,4,1,1)

  !{\tau}^{\Pi 12}
  tau(1:mesh,2) = IU*(-dpsi_local(1:mesh,orb,1,4,1) +dpsi_local(1:mesh,orb,2,3,1)&
                     &-dpsi_local(1:mesh,orb,3,2,1) +dpsi_local(1:mesh,orb,4,1,1))

  !{\tau}^{\Pi 13}
  tau(1:mesh,3) = dpsi_local(1:mesh,orb,1,3,1) -dpsi_local(1:mesh,orb,2,4,1)&
                &+dpsi_local(1:mesh,orb,3,1,1) -dpsi_local(1:mesh,orb,4,2,1)

  !{\tau}^{\Pi 21}
  tau(1:mesh,4) = dpsi_local(1:mesh,orb,1,4,2) +dpsi_local(1:mesh,orb,2,3,2)&
                &+dpsi_local(1:mesh,orb,3,2,2) +dpsi_local(1:mesh,orb,4,1,2)

  !{\tau}^{\Pi 22}
  tau(1:mesh,5) = IU*(-dpsi_local(1:mesh,orb,1,4,2) +dpsi_local(1:mesh,orb,2,3,2)&
                     &-dpsi_local(1:mesh,orb,3,2,2) +dpsi_local(1:mesh,orb,4,1,2))

  !{\tau}^{\Pi 23}
  tau(1:mesh,6) = dpsi_local(1:mesh,orb,1,3,2) -dpsi_local(1:mesh,orb,2,4,2)&
                &+dpsi_local(1:mesh,orb,3,1,2) -dpsi_local(1:mesh,orb,4,2,2)

  !{\tau}^{\Pi 31}
  tau(1:mesh,7) = dpsi_local(1:mesh,orb,1,4,3) +dpsi_local(1:mesh,orb,2,3,3)&
                &+dpsi_local(1:mesh,orb,3,2,3) +dpsi_local(1:mesh,orb,4,1,3)

  !{\tau}^{\Pi 32}
  tau(1:mesh,8) = IU*(-dpsi_local(1:mesh,orb,1,4,3) +dpsi_local(1:mesh,orb,2,3,3)&
                     &-dpsi_local(1:mesh,orb,3,2,3) +dpsi_local(1:mesh,orb,4,1,3))

  !{\tau}^{\Pi 33}
  tau(1:mesh,9) = dpsi_local(1:mesh,orb,1,3,3) -dpsi_local(1:mesh,orb,2,4,3)&
                &+dpsi_local(1:mesh,orb,3,1,3) -dpsi_local(1:mesh,orb,4,2,3)

  tau = 0.5_dp *IU *CCC *tau
  tau = tau + conjg(tau)

  torq(1:mesh,1) = -tau(1:mesh,6) +tau(1:mesh,8)
  torq(1:mesh,2) = -tau(1:mesh,7) +tau(1:mesh,3)
  torq(1:mesh,3) = -tau(1:mesh,2) +tau(1:mesh,4)

  return
END SUBROUTINE CALC_SPIN_TORQ

!==========================================================================
SUBROUTINE CALC_SPIN_TORQ_AM
!t_{e A}^i = -Z_e e {\epsilon}_{ijk} {\psi}^{\dagger} {\gamma}^0 {\gamma}^k A^j_M {\psi} 
!		   = - {1}/{c} {\epsilon}_{ijk} ( Z_e e c {\psi}^{\dagger} {\gamma}^0 {\gamma}^k {\psi} ) A^j_M 
!		   = - {1}/{c} {\epsilon}_{ijk} j^k_e A^j_M
! You have to set je & vecA_M in advance.
!==========================================================================
  use Precision
  use Constants
  use EDM_calculation
  use prop_param ! mesh
  use param_AM
  implicit none

  torq_AM(1:mesh,1) = -vecA_M(1:mesh,2)*je(1:mesh,3) +vecA_M(1:mesh,3)*je(1:mesh,2)
  torq_AM(1:mesh,2) = -vecA_M(1:mesh,3)*je(1:mesh,1) +vecA_M(1:mesh,1)*je(1:mesh,3)
  torq_AM(1:mesh,3) = -vecA_M(1:mesh,1)*je(1:mesh,2) +vecA_M(1:mesh,2)*je(1:mesh,1)

  torq_AM = torq_AM/CCC

  return
END SUBROUTINE CALC_SPIN_TORQ_AM

!==========================================================================
SUBROUTINE CALC_ZETA_POTENTIAL(orb)
!{\phi}_5 = {\hbar}/{2 Z_e e} {j}_5^0 
!= {\hbar c}/{2} {\psi}^\dag \gamma_5 {\psi}
!= {\hbar c}/{2} {{\psi}_L}^\dag {{\psi}_S} + h.c.
!==========================================================================
  use Precision
  use Constants
  use prop_param ! mesh
  use EDM_calculation
  implicit none

  integer,intent(in) :: orb

  chirality_density(1:mesh) = psi_local(1:mesh,orb,1,3) +psi_local(1:mesh,orb,2,4)
  chirality_density = chirality_density + conjg(chirality_density)

  zeta_potential = 0.5_dp *CCC *chirality_density

  return
END SUBROUTINE CALC_ZETA_POTENTIAL

!==========================================================================
SUBROUTINE CALC_ZETA_FORCE(orb)
!{\zeta}_e^i = - \partial_i {\phi}_5
!==========================================================================
  use Precision
  use Constants
  use prop_param ! mesh
  use EDM_calculation
  implicit none

  integer,intent(in) :: orb

  zeta_force(1:mesh,1) = dpsi_local(1:mesh,orb,1,3,1) +dpsi_local(1:mesh,orb,2,4,1)&
                       &+dpsi_local(1:mesh,orb,3,1,1) +dpsi_local(1:mesh,orb,4,2,1)
  zeta_force(1:mesh,2) = dpsi_local(1:mesh,orb,1,3,2) +dpsi_local(1:mesh,orb,2,4,2)&
                       &+dpsi_local(1:mesh,orb,3,1,2) +dpsi_local(1:mesh,orb,4,2,2)
  zeta_force(1:mesh,3) = dpsi_local(1:mesh,orb,1,3,3) +dpsi_local(1:mesh,orb,2,4,3)&
                      & +dpsi_local(1:mesh,orb,3,1,3) +dpsi_local(1:mesh,orb,4,2,3)

  zeta_force = -0.5_dp *CCC *zeta_force
  zeta_force = zeta_force + conjg(zeta_force)

  return
END SUBROUTINE CALC_ZETA_FORCE

!==========================================================================
SUBROUTINE CALC_TORQEDM_ELE(orb)
!t_{\rm EDM}^{E k} 
!= d_e \epsilon_{ijk} {\psi}^\dag \gamma^0 \Sigma^i E^j {\psi}
!= d_e \epsilon_{ijk} ( {{\psi}_L}^\dag \sigma^i {{\psi}_L} - {{\psi}_S}^\dag \sigma^i {\psi}_S ) E^j
!
!E^j = Eelec(j): external electric field for initial state
!==========================================================================
  use Precision
  use Constants
  use prop_param ! mesh
  use param_AM
  use EDM_calculation
  implicit none

  integer,intent(in) :: orb

  complex(kind=dp) :: work_vec(mesh,3)

  work_vec(1:mesh,1) = psi_local(1:mesh,orb,1,2) +psi_local(1:mesh,orb,2,1)&
                     &-psi_local(1:mesh,orb,3,4) -psi_local(1:mesh,orb,4,3)
  work_vec(1:mesh,2) = IU*(-psi_local(1:mesh,orb,1,2) +psi_local(1:mesh,orb,2,1)&
                          &+psi_local(1:mesh,orb,3,4) -psi_local(1:mesh,orb,4,3))
  work_vec(1:mesh,3) = psi_local(1:mesh,orb,1,1) -psi_local(1:mesh,orb,2,2)&
                     &-psi_local(1:mesh,orb,3,3) +psi_local(1:mesh,orb,4,4)

  torqEDM_ele(1:mesh,1) = work_vec(1:mesh,2)*Eelec(3) -work_vec(1:mesh,3)*Eelec(2)
  torqEDM_ele(1:mesh,2) = work_vec(1:mesh,3)*Eelec(1) -work_vec(1:mesh,1)*Eelec(3)
  torqEDM_ele(1:mesh,3) = work_vec(1:mesh,1)*Eelec(2) -work_vec(1:mesh,2)*Eelec(1)

  return
END SUBROUTINE CALC_TORQEDM_ELE

!==========================================================================
SUBROUTINE CALC_TORQEDM_MAG(orb)
!t_{\rm EDM}^{B k} 
!= d_e i \epsilon_{ijk} {\psi}^\dag \gamma^i B^j {\psi}
!= d_e i \epsilon_{ijk} ( -{{\psi}_S}^\dag \sigma^i {{\psi}_L} + {{\psi}_L}^\dag \sigma^i {\psi}_S ) B^j
!
!B^j = Bmag(j): external magnetic field for initial state
!==========================================================================
  use Precision
  use Constants
  use prop_param ! mesh
  use param_AM
  use EDM_calculation
  implicit none

  integer,intent(in) :: orb

  complex(kind=dp) :: work_vec(mesh,3)

  work_vec(1:mesh,1) = psi_local(1:mesh,orb,1,4) +psi_local(1:mesh,orb,2,3)&
                     &-psi_local(1:mesh,orb,3,2) -psi_local(1:mesh,orb,4,1)
  work_vec(1:mesh,2) = IU*(-psi_local(1:mesh,orb,1,4) +psi_local(1:mesh,orb,2,3)&
                          &+psi_local(1:mesh,orb,3,2) -psi_local(1:mesh,orb,4,1))
  work_vec(1:mesh,3) = psi_local(1:mesh,orb,1,3) -psi_local(1:mesh,orb,2,4)&
                     &-psi_local(1:mesh,orb,3,1) +psi_local(1:mesh,orb,4,2)

  torqEDM_mag(1:mesh,1) = work_vec(1:mesh,2)*Bmag(3) -work_vec(1:mesh,3)*Bmag(2)
  torqEDM_mag(1:mesh,2) = work_vec(1:mesh,3)*Bmag(1) -work_vec(1:mesh,1)*Bmag(3)
  torqEDM_mag(1:mesh,3) = work_vec(1:mesh,1)*Bmag(2) -work_vec(1:mesh,2)*Bmag(1)

  torqEDM_mag = IU *torqEDM_mag

  return
END SUBROUTINE CALC_TORQEDM_MAG

!==========================================================================
SUBROUTINE SET_VECA_M
!\vec{A}_{ext} = {1}/{2} \vec{B}_{ext} \times \vec{r}
!==========================================================================
  use Precision
  use prop_param
  use param_AM
  use EDM_calculation
  implicit none

  integer :: i,j,k,l

  l = 1
  do i = 1,meshx
    do j = 1,meshy
      do k = 1,meshz

        vecA_M(l,1) = Bmag(2)*x3(k) -Bmag(3)*x2(j)
        vecA_M(l,2) = Bmag(3)*x1(i) -Bmag(1)*x3(k)
        vecA_M(l,3) = Bmag(1)*x2(j) -Bmag(2)*x1(i)
        l = l + 1

      end do
    end do
  end do

  vecA_M = 0.5_dp *vecA_M

  return
END SUBROUTINE SET_VECA_M

!==========================================================================
SUBROUTINE CALC_EEFFEDM_NUC
!{\cal E}_{\rm eff}^{\rm nuc} 
!\approx \int \left\langle \Phi \left|: 2 \hat{\psi}_S^\dagger \vec{\sigma} \cdot \hat{\vec{E}}^{\rm nuc} \hat{\psi}_S :\right| \Phi \right\rangle d^3\vec{r}
!
! EeffEDM_nuc is integrand: 2 \hat{\psi}_S^\dagger \vec{\sigma} \cdot \hat{\vec{E}}^{\rm nuc} \hat{\psi}_S
!==========================================================================
  use Precision
  use prop_param ! mesh
  use EDM_calculation
  implicit none

  EeffEDM_nuc(1:mesh) = spin_small(1:mesh,1)*vecE_nuc(1:mesh,1)&
                      &+spin_small(1:mesh,2)*vecE_nuc(1:mesh,2)&
                      &+spin_small(1:mesh,3)*vecE_nuc(1:mesh,3)

  EeffEDM_nuc = -4._dp*EeffEDM_nuc

  return
END SUBROUTINE CALC_EEFFEDM_NUC

!==========================================================================
SUBROUTINE PRINT_QUANTITIES(filepath,frmt,orb)
!==========================================================================
  use Precision
  use prop_param
  use EDM_calculation
  implicit none

  character(*),intent(in) :: filepath,frmt
  integer,intent(in) :: orb

  real(kind=dp) :: x,y,z
  integer :: i,j,k,l,m
  character(len=3) :: orb_num

  write(orb_num,'(i3.3)') orb

  open(unit=31,file=filepath//'/'//orb_num//'_dens.dat')
  open(unit=32,file=filepath//'/'//orb_num//'_spin.dat')
  open(unit=33,file=filepath//'/'//orb_num//'_spin_small.dat')
  open(unit=34,file=filepath//'/'//orb_num//'_grad_N.dat')
  open(unit=35,file=filepath//'/'//orb_num//'_je.dat')
  open(unit=36,file=filepath//'/'//orb_num//'_torq.dat')
  open(unit=37,file=filepath//'/'//orb_num//'_torq_AM.dat')
  open(unit=38,file=filepath//'/'//orb_num//'_zeta_force.dat')
  open(unit=40,file=filepath//'/'//orb_num//'_tz.dat')
  open(unit=42,file=filepath//'/'//orb_num//'_zeta_potential.dat')
  open(unit=43,file=filepath//'/'//orb_num//'_torqEDM_ele.dat')
  open(unit=44,file=filepath//'/'//orb_num//'_torqEDM_mag.dat')
  open(unit=45,file=filepath//'/'//orb_num//'_EeffEDM_nuc.dat')
  open(unit=107,file=filepath//'/'//orb_num//'_chirality_density.dat')

  do i = 1,meshx
    do j = 1,meshy
      do k = 1,meshz

        l = (i -1) *meshyz +(j -1) *meshz +k

        write(31,frmt) x1(i), x2(j), x3(k), density(l)
        write(32,frmt) x1(i), x2(j), x3(k), (spin(l,m), m = 1,3)
        write(33,frmt) x1(i), x2(j), x3(k), (spin_small(l,m), m = 1,3)
        write(34,frmt) x1(i), x2(j), x3(k), (grad_N(l,m), m = 1,3)
        write(35,frmt) x1(i), x2(j), x3(k), (je(l,m), m = 1,3)
        write(36,frmt) x1(i), x2(j), x3(k), (torq(l,m), m = 1,3)
        write(37,frmt) x1(i), x2(j), x3(k), (torq_AM(l,m), m = 1,3)
        write(38,frmt) x1(i), x2(j), x3(k), (zeta_force(l,m), m = 1,3)
        write(39,frmt) x1(i), x2(j), x3(k), (torq(l,m) +torq_AM(l,m), m = 1,3)
        write(40,frmt) x1(i), x2(j), x3(k), (torq(l,m) +zeta_force(l,m), m = 1,3)
        write(42,frmt) x1(i), x2(j), x3(k), zeta_potential(l)
        write(43,frmt) x1(i), x2(j), x3(k), (torqEDM_ele(l,m), m = 1,3)
        write(44,frmt) x1(i), x2(j), x3(k), (torqEDM_mag(l,m), m = 1,3)
        write(45,frmt) x1(i), x2(j), x3(k), EeffEDM_nuc(l)
        write(107,frmt) x1(i), x2(j), x3(k), (chirality_density(l))

      end do ! k = 1,meshz
      write(31,*)
      write(32,*)
      write(33,*)
      write(34,*)
      write(35,*)
      write(36,*)
      write(37,*)
      write(38,*)
      write(39,*)
      write(40,*)
      write(42,*)
      write(43,*)
      write(44,*)
      write(45,*)
      write(107,*)

    end do ! j = 1,meshy
    write(31,*)
    write(32,*)
    write(33,*)
    write(34,*)
    write(35,*)
    write(36,*)
    write(37,*)
    write(38,*)
    write(39,*)
    write(40,*)
    write(42,*)
    write(43,*)
    write(44,*)
    write(45,*)
    write(107,*)

  end do ! i = 1,meshx
  close(unit=31)
  close(unit=32)
  close(unit=33)
  close(unit=34)
  close(unit=35)
  close(unit=36)
  close(unit=37)
  close(unit=38)
  close(unit=39)
  close(unit=40)
  close(unit=42)
  close(unit=43)
  close(unit=44)
  close(unit=45)
  close(unit=107)

  return
END SUBROUTINE PRINT_QUANTITIES
