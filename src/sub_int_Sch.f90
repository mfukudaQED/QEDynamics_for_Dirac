!============================================================
function intSsch_jk_mat(j,k,p,a,q,b)
! S^{jk}_{n^a m^b} = -hbar c int s^j psi^dagger_{n^a}(s) gamma^0 gamma^k psi_{m^b}(s) d^3s
!
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!
! 140220 
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none
 
  complex(kind=dp) :: intSsch_jk_mat
  integer,intent(in) :: p,q
  integer,intent(in) :: j,k
  character(LEN=1),intent(in) :: a,b
  
  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  
  complex(kind=dp) :: intM_pq_14(3),intM_pq_23(3),intM_pq_32(3),intM_pq_41(3)
  complex(kind=dp) :: intM_pq_13(3),intM_pq_24(3),intM_pq_31(3),intM_pq_42(3)

  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q

  ! beta = gamma^0 = diag(1,1,-1,-1)
  
  call calc_intMschi_pq(1,4,NBS_L,NBS_S,pg,c_p,c_q,intM_pq_14) 
  call calc_intMschi_pq(2,3,NBS_L,NBS_S,pg,c_p,c_q,intM_pq_23) 
  call calc_intMschi_pq(3,2,NBS_L,NBS_S,pg,c_p,c_q,intM_pq_32) 
  call calc_intMschi_pq(4,1,NBS_L,NBS_S,pg,c_p,c_q,intM_pq_41) 
  call calc_intMschi_pq(1,3,NBS_L,NBS_S,pg,c_p,c_q,intM_pq_13) 
  call calc_intMschi_pq(2,4,NBS_L,NBS_S,pg,c_p,c_q,intM_pq_24) 
  call calc_intMschi_pq(3,1,NBS_L,NBS_S,pg,c_p,c_q,intM_pq_31) 
  call calc_intMschi_pq(4,2,NBS_L,NBS_S,pg,c_p,c_q,intM_pq_42) 

  if(k==1) then
    intSsch_jk_mat = intM_pq_14(j)+intM_pq_23(j)+intM_pq_32(j)+intM_pq_41(j)
  else if(k==2) then
    intSsch_jk_mat = -IU*(intM_pq_14(j)+intM_pq_32(j)) +IU*(intM_pq_23(j)+intM_pq_41(j))
  else if(k==3) then
    intSsch_jk_mat = (intM_pq_13(j)+intM_pq_31(j)) -(intM_pq_24(j)+intM_pq_42(j))
  end if

  intSsch_jk_mat = - IU*CCC*intSsch_jk_mat

  return
end function intSsch_jk_mat

!============================================================
function intMsch_i_mat(i,p,a,q,b)
! M^j_{n^a m^b} = mc^2 int s^j psi^dagger_{n^a}(s) gamma^0 psi_{m^b}(s) d^3s
!
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!
! 140220 
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none
 
  complex(kind=dp) :: intMsch_i_mat
  integer,intent(in) :: p,q
  integer,intent(in) :: i
  character(LEN=1),intent(in) :: a,b
  
  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  complex(kind=dp) :: intM_pq_11(3),intM_pq_22(3),intM_pq_33(3),intM_pq_44(3)
  
  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q

  ! beta = gamma^0 = diag(1,1,-1,-1)

  call calc_intMschi_pq(1,1,NBS_L,NBS_S,pg,c_p,c_q,intM_pq_11) 
  call calc_intMschi_pq(2,2,NBS_L,NBS_S,pg,c_p,c_q,intM_pq_22) 
  call calc_intMschi_pq(3,3,NBS_L,NBS_S,pg,c_p,c_q,intM_pq_33) 
  call calc_intMschi_pq(4,4,NBS_L,NBS_S,pg,c_p,c_q,intM_pq_44) 

  intMsch_i_mat = (intM_pq_11(i) +intM_pq_22(i) -intM_pq_33(i) -intM_pq_44(i))*CCC**2

  return
end function intMsch_i_mat


!=========================================================================
subroutine calc_intMschi_pq(ip,iq,NL,NS,pg,c_p,c_q,intM_pq)
! 
! 140220
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none

  integer,intent(in) :: ip,iq ! spinor indice
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=dp),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  
  complex(kind=dp),intent(out) :: intM_pq(3)
  
  integer :: NLS
  integer :: i,j,k
  real(kind=dp) :: overlap
  real(kind=dp) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=dp) :: posC(3)
  real(kind=dp) :: alphaA, alphaB ! exponents
  integer :: nx,nbarx,ny,nbary,nz,nbarz
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg

  real(kind=8) :: moment(3)

  intM_pq(:) = (0._dp,0._dp)

  posC(:)=0._dp

  do i=1,NLS(ip,NL,NS)
     do j=1,NLS(iq,NL,NS)
        call set_pg(ip,i,NL,NS,pg,posA,alphaA,nx,ny,nz)
        call set_pg(iq,j,NL,NS,pg,posB,alphaB,nbarx,nbary,nbarz)
        call gauss_int_moment(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,posC,moment)
        f_norm_pg = norm_pg(alphaA,nx,ny,nz)*norm_pg(alphaB,nbarx,nbary,nbarz)
        do k=1,3
          intM_pq(k) = intM_pq(k) + conjg(c_p(ip,i))*c_q(iq,j)*moment(k)*f_norm_pg
        end do
     end do
  end do
  return
end subroutine calc_intMschi_pq


!============================================================
function inthsch2_mat(rad,mass,p,a,q,b)
! Integration of one-electron terms for Schwarzshild spacetime
! Tsch+Ssch+Msch+V
! for given radius (R) and mass
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!
! rad is given in units of Schwarzshild radius (rad = R/r_g)
! mass is given in units of Solar mass
!
! 140220
!============================================================
  use Precision
  use Constants ! 140115 modified to use Rg_OVER_Msolar= 5.59e13_dp
  implicit none

  complex(kind=dp) :: inthsch2_mat
  real(kind=dp),intent(in) :: rad,mass
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b
  
  complex(kind=dp) :: intTsch_mat(3),intSsch_mat,intM_mat,intV_mat

  complex(kind=dp) :: intTsch_jkl_mat(7) ! (T^111,T^122,T^133,T^212,T^221,T^313,T^331)
  complex(kind=dp) :: intTsch_kl_mat(3)  ! (T^11, T^22, T^33)
  real(kind=dp) :: c1,c2,c3,c4,c5,c6,c7,c8
  real(kind=dp) :: TR,SR,r_g
  complex(kind=dp) :: T11,T22,T33
  complex(kind=dp) :: T111,T122,T133,T212,T221,T313,T331

  integer :: i,j
  complex(kind=dp) :: intMsch_i_mat !function
  complex(kind=dp) :: intSsch_jk_mat

  r_g = Rg_OVER_Msolar*mass ! convert mass in solar mass to r_g in a.u.
!  TR = sqrt(1._dp -1._dp/rad)
!  SR = (1._dp/rad -(3._dp/4._dp)/rad**2)/r_g
!  write(*,*) rad,mass
!  write(*,*) r_g,TR,SR
!  stop
 
  c1 = 1._dp -1._dp/rad
  c2 = c1/rad/r_g
  c3 = 1._dp/rad**2/r_g
  c4 = 1._dp -0.5_dp/rad
  c5 = c4/rad/r_g
  c6 = (1._dp/rad -(3._dp/4._dp)/rad**2)/r_g !S1(R)
  c7 = -(1._dp/rad -(3._dp/2._dp)/rad**2)/r_g**2/rad ! -S2(R)/R
  c8 = -c5/2._dp/rad/r_g
!  call calc_Tsch_kl_mat(p,a,q,b,intTsch_kl_mat)
  call calc_Tsch_mat(p,a,q,b,intTsch_kl_mat)
  T11 = intTsch_kl_mat(1)
  T22 = intTsch_kl_mat(2)
  T33 = intTsch_kl_mat(3)
  call calc_Tsch_jkl_mat(p,a,q,b,intTsch_jkl_mat)
  T111 = intTsch_jkl_mat(1)
  T122 = intTsch_jkl_mat(2)
  T133 = intTsch_jkl_mat(3)
  T212 = intTsch_jkl_mat(4)
  T221 = intTsch_jkl_mat(5)
  T313 = intTsch_jkl_mat(6)
  T331 = intTsch_jkl_mat(7)

  ! NOPT_GR: 1 for including inhomogeneity to 1st order, 2 for homogeneous only.
  ! NOPT_GR_T: 1 for including T, 2 for not including T
  ! NOPT_GR_S: 1 for including S, 2 for not including S
  ! NOPT_GR_M: 1 for including M, 2 for not including M
  ! NOPT_GR_V: 1 for including V, 2 for not including V

  if(NOPT_GR.eq.1) then
     intTsch_mat(1) =  c1*T11 +c2*(T212+T313) +c3*T111
     intTsch_mat(2) =  c4*T22 -c5*T221 +0.5_dp*c3*T122
     intTsch_mat(3) =  c4*T33 -c5*T331 +0.5_dp*c3*T133
  elseif(NOPT_GR.eq.2) then
     intTsch_mat(1) =  c1*T11 !+c2*(T212+T313) +c3*T111
     intTsch_mat(2) =  c4*T22 !-c5*T221 +0.5_dp*c3*T122
     intTsch_mat(3) =  c4*T33 !-c5*T331 +0.5_dp*c3*T133
  else
     write(*,*) "Check NOPT_GR."
     stop
  end if
  
  inthsch2_mat = (0._dp,0._dp)

  if(NOPT_GR_T.eq.1) then
    inthsch2_mat = inthsch2_mat + intTsch_mat(1)+intTsch_mat(2)+intTsch_mat(3)
  end if

  if(NOPT_GR_S.eq.1) then
    if(NOPT_GR.eq.1) then
      inthsch2_mat = inthsch2_mat + c6*intSsch_mat(1,p,a,q,b) + c7*intSsch_jk_mat(1,1,p,a,q,b) + c8*intSsch_jk_mat(2,2,p,a,q,b)
    else if(NOPT_GR.eq.2) then
      inthsch2_mat = inthsch2_mat + c6*intSsch_mat(1,p,a,q,b) !+ c7*intSsch_jk_mat(1,1,p,a,q,b) + c8*intSsch_jk_mat(2,2,p,a,q,b)
    else
       write(*,*) "Check NOPT_GR."
       stop
    end if
  end if

  if(NOPT_GR_M.eq.1) then
    if(NOPT_GR.eq.1) then
      inthsch2_mat = inthsch2_mat + c4*intM_mat(p,a,q,b) +0.5_dp*c3*intMsch_i_mat(1,p,a,q,b) 
    else if(NOPT_GR.eq.2) then
      inthsch2_mat = inthsch2_mat + c4*intM_mat(p,a,q,b) !+0.5_dp*c3*intMsch_i_mat(1,p,a,q,b) 
    else
       write(*,*) "Check NOPT_GR."
       stop
    end if
  end if

  if(NOPT_GR_V.eq.1) then
    inthsch2_mat = inthsch2_mat + intV_mat(p,a,q,b)
  end if

  
  return
end function inthsch2_mat

!!$!============================================================
!!$function inthsch2_mat(rad,mass,p,a,q,b)
!!$! Integration of one-electron terms for Schwarzshild spacetime
!!$! Tsch+Ssch+Msch+V
!!$! for given radius (R) and mass
!!$! p,q : labels for molecular orbitals.
!!$! a,b : labels for electron/positron ("+"/"-")
!!$!
!!$! rad is given in units of Schwarzshild radius (rad = R/r_g)
!!$! mass is given in units of Solar mass
!!$!
!!$! 140115
!!$! 
!!$!============================================================
!!$  use Precision
!!$  use Constants ! 140115 modified to use Rg_OVER_Msolar= 5.59e13_dp
!!$  implicit none
!!$
!!$  complex(kind=dp) :: inthsch2_mat
!!$  real(kind=dp),intent(in) :: rad,mass
!!$  integer,intent(in) :: p,q
!!$  character(LEN=1),intent(in) :: a,b
!!$  
!!$  complex(kind=dp) :: intTsch_mat(3),intSsch_mat,intM_mat,intV_mat
!!$
!!$  complex(kind=dp) :: intTsch_jkl_mat(7) ! (T^111,T^123,T^132,T^212,T^231,T^313,T^321)
!!$  complex(kind=dp) :: intTsch_kl_mat(3)  ! (T^11, T^23, T^32)
!!$  real(kind=dp) :: c1,c2,c3,c4,c5
!!$  real(kind=dp) :: TR,SR,r_g
!!$  complex(kind=dp) :: T11,T23,T32
!!$  complex(kind=dp) :: T111,T123,T132,T212,T231,T313,T321
!!$
!!$  r_g = Rg_OVER_Msolar*mass ! convert mass in solar mass to r_g in a.u.
!!$!  TR = sqrt(1._dp -1._dp/rad)
!!$!  SR = (1._dp/rad -(3._dp/4._dp)/rad**2)/r_g
!!$!  write(*,*) rad,mass
!!$!  write(*,*) r_g,TR,SR
!!$!  stop
!!$ 
!!$  c1 = 1._dp -1._dp/rad
!!$  c2 = c1/rad/r_g
!!$  c3 = 1._dp/rad**2/r_g
!!$  c4 = 1._dp -0.5_dp/rad
!!$  c5 = c4/rad/r_g
!!$  call calc_Tsch_kl_mat(p,a,q,b,intTsch_kl_mat)
!!$  T11 = intTsch_kl_mat(1)
!!$  T23 = intTsch_kl_mat(2)
!!$  T32 = intTsch_kl_mat(3)
!!$  call calc_Tsch_jkl_mat(p,a,q,b,intTsch_jkl_mat)
!!$  T111 = intTsch_jkl_mat(1)
!!$  T123 = intTsch_jkl_mat(2)
!!$  T132 = intTsch_jkl_mat(3)
!!$  T212 = intTsch_jkl_mat(4)
!!$  T231 = intTsch_jkl_mat(5)
!!$  T313 = intTsch_jkl_mat(6)
!!$  T321 = intTsch_jkl_mat(7)
!!$
!!$  intTsch_mat(1) =  c1*T11 +c2*(T212+T313) +c3*T111
!!$  intTsch_mat(2) = -c4*T23 +c5*T321 -0.5_dp*c3*T123
!!$  intTsch_mat(3) =  c4*T32 -c5*T231 +0.5_dp*c3*T132
!!$
!!$  intTsch_mat(1) =  c1*T11 !+c2*(T212+T313) +c3*T111
!!$  intTsch_mat(2) = -c4*T23 !+c5*T321 -0.5_dp*c3*T123
!!$  intTsch_mat(3) =  c4*T32 !-c5*T231 +0.5_dp*c3*T132
!!$
!!$  intTsch_mat(1) =  c1*T11 !+c2*(T212+T313) +c3*T111
!!$  intTsch_mat(2) =  c4*T23 !+c5*T321 -0.5_dp*c3*T123
!!$  intTsch_mat(3) =  c4*T32 !-c5*T231 +0.5_dp*c3*T132
!!$
!!$
!!$  
!!$  inthsch2_mat = intTsch_mat(1)+intTsch_mat(2)+intTsch_mat(3) &
!!$       & +c4*intM_mat(p,a,q,b) +intV_mat(p,a,q,b)
!!$
!!$
!!$!  call calc_Tsch_mat(p,a,q,b,intTsch_mat)
!!$!  inthsch_mat = TR**2*intTsch_mat(1) +TR*(intTsch_mat(2)+intTsch_mat(3)) &
!!$!       & +SR*intSsch_mat(1,p,a,q,b) +TR*intM_mat(p,a,q,b) +intV_mat(p,a,q,b)
!!$  
!!$  return
!!$end function inthsch2_mat





!============================================================
subroutine calc_Tsch_jkl_mat(p,a,q,b,intTsch_jkl_mat)
! Integration for kinetic energy for Schwarzshild spacetime
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!
! intTsch_jkl_mat(7) = (T^111,T^123,T^132,T^212,T^231,T^313,T^321)
!
! 140121
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none
 
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b
  complex(kind=dp),intent(out) :: intTsch_jkl_mat(7)  ! (T^111,T^122,T^133,T^212,T^221,T^313,T^331)
  ! T^111 : intTsch_jkl_mat(1)
  ! T^122 : intTsch_jkl_mat(2)
  ! T^133 : intTsch_jkl_mat(3)
  ! T^212 : intTsch_jkl_mat(4)
  ! T^221 : intTsch_jkl_mat(5)
  ! T^313 : intTsch_jkl_mat(6)
  ! T^331 : intTsch_jkl_mat(7)
  
  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  complex(kind=dp) :: intTjl_pq(3,3), intTjl_pq_ab(3,3,4,4)
  integer :: j,k,l
  integer :: alpha,beta
  complex(kind=dp) :: sum
  
  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  
  ! compute int x^j psi^+_alpha (d/dx^l psi)_beta  for j,l=1,2,3
  do alpha=1,4
     do beta=1,4
        call calc_intTjl_pq(alpha,beta,NBS_L,NBS_S,pg,c_p,c_q,intTjl_pq) 
        do j=1,3
           do l=1,3
              intTjl_pq_ab(j,l,alpha,beta) = intTjl_pq(j,l)
           end do
        end do
     end do
  end do
  
  ! compute sum_{alpha,beta} (-i hbar c) int x^j psi^+_alpha [gamma^0 gamma^k]_{alpha,beta} (d/dx^l psi)_beta  
  ! for (j,k,l) = (111),(122),(133),(212),(221),(313),(331)

  ! T^111
  j=1; k=1; l=1
  sum = (0._dp,0._dp)
  do alpha=1,4
     do beta=1,4
        sum = sum + Gam0k(k,alpha,beta)*intTjl_pq_ab(j,l,alpha,beta)
     end do
  end do
  intTsch_jkl_mat(1) = sum*(-IU*CCC)
  
  ! T^122
  j=1; k=2; l=2
  sum = (0._dp,0._dp)
  do alpha=1,4
     do beta=1,4
        sum = sum + Gam0k(k,alpha,beta)*intTjl_pq_ab(j,l,alpha,beta)
     end do
  end do
  intTsch_jkl_mat(2) = sum*(-IU*CCC)
  
  ! T^133
  j=1; k=3; l=3
  sum = (0._dp,0._dp)
  do alpha=1,4
     do beta=1,4
        sum = sum + Gam0k(k,alpha,beta)*intTjl_pq_ab(j,l,alpha,beta)
     end do
  end do
  intTsch_jkl_mat(3) = sum*(-IU*CCC)

  ! T^212
  j=2; k=1; l=2
  sum = (0._dp,0._dp)
  do alpha=1,4
     do beta=1,4
        sum = sum + Gam0k(k,alpha,beta)*intTjl_pq_ab(j,l,alpha,beta)
     end do
  end do
  intTsch_jkl_mat(4) = sum*(-IU*CCC)

  ! T^221
  j=2; k=2; l=1
  sum = (0._dp,0._dp)
  do alpha=1,4
     do beta=1,4
        sum = sum + Gam0k(k,alpha,beta)*intTjl_pq_ab(j,l,alpha,beta)
     end do
  end do
  intTsch_jkl_mat(5) = sum*(-IU*CCC)

  ! T^313
  j=3; k=1; l=3
  sum = (0._dp,0._dp)
  do alpha=1,4
     do beta=1,4
        sum = sum + Gam0k(k,alpha,beta)*intTjl_pq_ab(j,l,alpha,beta)
     end do
  end do
  intTsch_jkl_mat(6) = sum*(-IU*CCC)

  ! T^331
  j=3; k=3; l=1
  sum = (0._dp,0._dp)
  do alpha=1,4
     do beta=1,4
        sum = sum + Gam0k(k,alpha,beta)*intTjl_pq_ab(j,l,alpha,beta)
     end do
  end do
  intTsch_jkl_mat(7) = sum*(-IU*CCC)

  return
end subroutine calc_Tsch_jkl_mat

!=========================================================================
subroutine calc_intTjl_pq(ip,iq,NL,NS,pg,c_p,c_q,intTjl_pq)
! Int x^j (psi^+)_ip (d/dx^l psi)_iq for j,l=1-3
!
! 140115
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none
  
  integer,intent(in) :: ip,iq ! spinor indice
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=dp),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  
  complex(kind=dp),intent(out) :: intTjl_pq(3,3)
  
  integer :: NLS
  integer :: i,j
  integer :: jj,ll !j,l=1-3
  real(kind=dp) :: momgrad(3,3)
  real(kind=dp) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=dp) :: alphaA, alphaB ! exponents
  integer :: nx,nbarx,ny,nbary,nz,nbarz
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg

  intTjl_pq(1:3,1:3) = (0._dp,0._dp)

  do i=1,NLS(ip,NL,NS)
     do j=1,NLS(iq,NL,NS)
        call set_pg(ip,i,NL,NS,pg,posA,alphaA,nx,ny,nz)
        call set_pg(iq,j,NL,NS,pg,posB,alphaB,nbarx,nbary,nbarz)
!        call gauss_int_relKE(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,grad) <--- calc_intT_pq
        call gauss_int_momgrad(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,momgrad)
        f_norm_pg = norm_pg(alphaA,nx,ny,nz)*norm_pg(alphaB,nbarx,nbary,nbarz)
        do jj=1,3
           do ll=1,3
              intTjl_pq(jj,ll) = intTjl_pq(jj,ll) + conjg(c_p(ip,i))*c_q(iq,j)*momgrad(jj,ll)*f_norm_pg
           end do
        end do
     end do
  end do
  return
end subroutine calc_intTjl_pq



!============================================================
subroutine calc_Tsch_kl_mat(p,a,q,b,intTsch_kl_mat)
! Integration for kinetic energy for Schwarzshild spacetime
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!
! intTsch_kl_mat(3) = (T^11, T^23, T^32)
!
! 140121
! 
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none
 
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b
  complex(kind=dp),intent(out) :: intTsch_kl_mat(3)  ! (T^11, T^23, T^32)
  
  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  complex(kind=dp) :: intT_pq(3),intT_pq_ab(3,4,4)
  integer :: l
  integer :: alpha,beta
  complex(kind=dp) :: sum

  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q

  ! compute int psi^+_alpha (d/dx^l psi)_beta  for l=1,2,3
  do alpha=1,4
     do beta=1,4
        call calc_intT_pq(alpha,beta,NBS_L,NBS_S,pg,c_p,c_q,intT_pq) 
        do l=1,3
           intT_pq_ab(l,alpha,beta) = intT_pq(l)
        end do
     end do
  end do
  
  ! compute sum_{alpha,beta} int (-i hbar c) psi^+_alpha [gamma^0 gamma^k]_{alpha,beta} (d/dx^l psi)_beta  
  ! for (k,l) = (1,1), (2,3), (3,2)

  ! T^11 
  sum = (0._dp,0._dp)
  do alpha=1,4
     do beta=1,4
        sum = sum + Gam0k(1,alpha,beta)*intT_pq_ab(1,alpha,beta)
     End do
  end do
  intTsch_kl_mat(1) = sum*(-IU*CCC)

  ! T^23
  sum = (0._dp,0._dp)
  do alpha=1,4
     do beta=1,4
        sum = sum + Gam0k(2,alpha,beta)*intT_pq_ab(3,alpha,beta)
!        sum = sum + Gam0k(2,alpha,beta)*intT_pq_ab(2,alpha,beta)  !T^22
     end do
  end do
  intTsch_kl_mat(2) = sum*(-IU*CCC)

  ! T^32
  sum = (0._dp,0._dp)
  do alpha=1,4
     do beta=1,4
        sum = sum + Gam0k(3,alpha,beta)*intT_pq_ab(2,alpha,beta)
!        sum = sum + Gam0k(3,alpha,beta)*intT_pq_ab(3,alpha,beta) !T^33
     end do
  end do
  intTsch_kl_mat(3) = sum*(-IU*CCC)
  
  return
end subroutine calc_Tsch_kl_mat



!====================================================================================
subroutine gauss_int_momgrad(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,momgrad)
!
! [Input]
! posA : position of primitive gaussian A (unnormalized) 
! alphaA : exponent of A
! nx,ny,nz : type of p.g. A 
! posB : position of primitive gaussian B (unnormalized) 
! alphaB : exponenet of B
! nbarx,nbary,nbarz : type of p.g. B
! 
! [output]
! momgrad(j,l)  : int x^j g grad_l g (j,l=1-3)
!
! 140115
!====================================================================================
  use Precision
  implicit none

  real(kind=dp),intent(in) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=dp),intent(in) :: alphaA, alphaB ! exponents
  integer,intent(in) :: nx,nbarx,ny,nbary,nz,nbarz

  real(kind=dp),intent(out) :: momgrad(3,3)

  real(kind=dp) :: grad(3),grad_x(3),grad_y(3),grad_z(3)
  integer :: j,l

  call gauss_int_relKE(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,grad) ! grad^l(vec{A},alphaA,vec{n}; vec{B},alphaB,vec{nbar})
  call gauss_int_relKE(posA,alphaA,nx+1,ny,nz,posB,alphaB,nbarx,nbary,nbarz,grad_x) ! grad^l(vec{A},alphaA,vec{n}+(1,0,0); vec{B},alphaB,vec{nbar})
  call gauss_int_relKE(posA,alphaA,nx,ny+1,nz,posB,alphaB,nbarx,nbary,nbarz,grad_y) ! grad^l(vec{A},alphaA,vec{n}+(0,1,0); vec{B},alphaB,vec{nbar})
  call gauss_int_relKE(posA,alphaA,nx,ny,nz+1,posB,alphaB,nbarx,nbary,nbarz,grad_z) ! grad^l(vec{A},alphaA,vec{n}+(0,0,1); vec{B},alphaB,vec{nbar})
  
  j=1
  do l=1,3
     momgrad(j,l) = grad_x(l) + posA(j)*grad(l)
  end do
  
  j=2
  do l=1,3
     momgrad(j,l) = grad_y(l) + posA(j)*grad(l)
  end do
  
  j=3
  do l=1,3
     momgrad(j,l) = grad_z(l) + posA(j)*grad(l)
  end do

  return
end subroutine gauss_int_momgrad



!==========================================================================



!======================================================
subroutine setQmat_hsch(rad,mass,hsch_Qmat)
! 1 ~ 2*NBS : "+"
! 2*NBS +1 ~ 2*NBS : "-"
! in "+" and "-", 1,1bar,2,2bar, ... respectively
!
! rad is given in units of Schwarzshild radius
! mass is given in units of Solar mass
!
! 140116
!======================================================
  use Precision
  use DiracOutput
  implicit none

  real(kind=dp),intent(in) :: rad,mass
  complex(kind=dp),intent(out) :: hsch_Qmat(4*NBS,4*NBS)

  complex(kind=dp) :: inthsch_mat,inthsch2_mat
  integer :: nn,mm  ! index for 1~4*NBS
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b

  do nn=1,4*NBS
     do mm=1,4*NBS
        
        call index_from_Qmat(nn,n,a)
        call index_from_Qmat(mm,m,b)

!        hsch_Qmat(nn,mm) = inthsch_mat(rad,mass,n,a,m,b)
        hsch_Qmat(nn,mm) = inthsch2_mat(rad,mass,n,a,m,b) !140116
        
     end do
  end do

end subroutine setQmat_hsch

!============================================================
function inthsch_mat(rad,mass,p,a,q,b)
! Integration of one-electron terms for Schwarzshild spacetime
! Tsch+Ssch+Msch+V
! for given radius and mass
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!
! rad is given in units of Schwarzshild radius
! mass is given in units of Solar mass
!
! 130119
!============================================================
  use Precision
  use Constants ! 140115 modified to use Rg_OVER_Msolar= 5.59e13_dp
  implicit none

  complex(kind=dp) :: inthsch_mat
  real(kind=dp),intent(in) :: rad,mass
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  complex(kind=dp) :: intTsch_mat(3),intSsch_mat,intM_mat,intV_mat
  real(kind=dp) :: TR,SR,r_g

!  r_g = 5.59e13_dp*mass ! convert mass in solar mass to r_g in a.u.
  r_g = Rg_OVER_Msolar*mass ! convert mass in solar mass to r_g in a.u.
  TR = sqrt(1._dp -1._dp/rad)
  SR = (1._dp/rad -(3._dp/4._dp)/rad**2)/r_g
!  write(*,*) rad,mass
!  write(*,*) r_g,TR,SR
!  stop
  call calc_Tsch_mat(p,a,q,b,intTsch_mat)
  inthsch_mat = TR**2*intTsch_mat(1) +TR*(intTsch_mat(2)+intTsch_mat(3)) &
       & +SR*intSsch_mat(1,p,a,q,b) +TR*intM_mat(p,a,q,b) +intV_mat(p,a,q,b)
  
  return
end function inthsch_mat


!============================================================
function intSsch_mat(k,p,a,q,b)
! Integration of (-i hbar c) psi^+_{pa} gamma^0 gamma^k psi_{qb} 
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
! k : vector index
!
! used for S^(Sch)_{pa,qb} 
! 130119
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none
 
  complex(kind=dp) :: intSsch_mat
  integer,intent(in) :: k
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b
  
  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  complex(kind=dp) :: intM_pq(4,4)
  integer :: alpha,beta
  complex(kind=dp) :: sum
  

  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q

  ! compute int psi^+_alpha psi_beta
  do alpha=1,4
     do beta=1,4
        call calc_intN_pq(alpha,beta,NBS_L,NBS_S,pg,c_p,c_q,intM_pq(alpha,beta)) 
     end do
  end do

  ! compute int psi^+_alpha gamma^0 gamma^k psi_beta
  sum = (0._dp,0._dp)
  do alpha=1,4
     do beta=1,4
        sum = sum + intM_pq(alpha,beta)*Gam0k(k,alpha,beta)
     end do
  end do
  
  intSsch_mat = sum*(-IU*CCC)

  return
end function intSsch_mat

!============================================================
subroutine calc_Tsch_mat(p,a,q,b,intTsch_mat)
! Integration of kinetic energy for Schwarzshild spacetime
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!
! specially for computation on x-axis
! 130119
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none
 
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b
  complex(kind=dp),intent(out) :: intTsch_mat(3)
  
  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  complex(kind=dp) :: intT_pq(3),intT_pq_ab(3,4,4)
  integer :: k
  integer :: alpha,beta
  complex(kind=dp) :: sum

  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q

  ! compute int psi^+_alpha (d/dx_k psi)_beta  for k=1,2,3
  do alpha=1,4
     do beta=1,4
        call calc_intT_pq(alpha,beta,NBS_L,NBS_S,pg,c_p,c_q,intT_pq) 
        do k=1,3
           intT_pq_ab(k,alpha,beta) = intT_pq(k)
        end do
     end do
  end do
  
  ! T^(Sch)_1
  sum = (0._dp,0._dp)
  do alpha=1,4
     do beta=1,4
        sum = sum + intT_pq_ab(1,alpha,beta)*Gam0k(1,alpha,beta)
     end do
  end do
  intTsch_mat(1) = sum*(-IU*CCC)

  ! T^(Sch)_2
  sum = (0._dp,0._dp)
  do alpha=1,4
     do beta=1,4
!        sum = sum + (-intT_pq_ab(3,alpha,beta))*Gam0k(2,alpha,beta)
        sum = sum + intT_pq_ab(2,alpha,beta)*Gam0k(2,alpha,beta)
     end do
  end do
  intTsch_mat(2) = sum*(-IU*CCC)

  ! T^(Sch)_3
  sum = (0._dp,0._dp)
  do alpha=1,4
     do beta=1,4
!        sum = sum + intT_pq_ab(2,alpha,beta)*Gam0k(3,alpha,beta)
        sum = sum + intT_pq_ab(3,alpha,beta)*Gam0k(3,alpha,beta)
     end do
  end do
  intTsch_mat(3) = sum*(-IU*CCC)

  return
end subroutine calc_Tsch_mat


