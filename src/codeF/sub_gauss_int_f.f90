! Last Change:27-Nov-2014.
!====================================================================================
! 2013.11.07
! - subroutine gauss_int_xigradj(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,xigradj)
! - function inttwoelegrad1_pg(m,i,j,k,l,s_i,s_j,s_k,s_l,NL,NS,pg)
! - function inttwoelegrad_pg(m,i,j,k,l,s_i,s_j,s_k,s_l,NL,NS,pg)
! - function intEeff_EDM_ob_pg(i,j,s_i,s_j,NL,NS,pg)
! - subroutine gauss_int_2nd_moment(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,posC,2moment)
! - subroutine calc_2nd_moment_integral(i,j,k,alphaP,vecPC,int_2mom)
!====================================================================================

!====================================================================================
subroutine gauss_int_xigradj(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,xigradj)
! [input]
! posA : position of primitive gaussian A (unnormalized) 
! alphaA : exponent of A
! nx,ny,nz : type of p.g. A 
! posB : position of primitive gaussian B (unnormalized) 
! alphaB : exponenet of B
! nbarx,nbary,nbarz : type of p.g. B
! 
! [output]
! xigradj(3,3)  : g*x^i \partial_j g
!
!====================================================================================
  implicit none

  real(kind=8) :: PI
  !  integer :: nx,nbarx,ny,nbary,nz,nbarz ! angular momentum
  
  real(kind=8),intent(in) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=8),intent(in) :: alphaA, alphaB ! exponents
  integer,intent(in) :: nx,nbarx,ny,nbary,nz,nbarz

  real(kind=8),intent(out) :: xigradj(3,3)

  real(kind=8) :: c_PAB ! prefactor when two gaussians are combined (E_IJ)
  real(kind=8) :: alphaP, posP(3) ! exponents and position of combined gaussian
  real(kind=8) :: vecPA(3), vecPB(3) ! P-A, P-B
  integer :: i,j,k,l
  
  real(kind=8) :: grad(3), grad1(3,3)

  PI = atan(1.d0)*4.d0

  call gauss_int_relKE(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,grad)
  call gauss_int_relKE(posA,alphaA,nx+1,ny,nz,posB,alphaB,nbarx,nbary,nbarz,grad1(1,:))
  call gauss_int_relKE(posA,alphaA,nx,ny+1,nz,posB,alphaB,nbarx,nbary,nbarz,grad1(2,:))
  call gauss_int_relKE(posA,alphaA,nx,ny,nz+1,posB,alphaB,nbarx,nbary,nbarz,grad1(3,:))

  do j=1,3
    do i=1,3
      xigradj(i,j) = grad1(i,j) + posA(i)*grad(j)
    end do
  end do
  ! YOU HAVE TO CONSIDER norm_pg

end subroutine gauss_int_xigradj

!=========================================================================
function inttwoelegrad1_pg(m,i,j,k,l,s_i,s_j,s_k,s_l,NL,NS,pg)
! (\partial_mNM|PQ) = Int dr ds \partial_m[(g^si(r)^+)_i] (g^sj(r))_j (1/|r-s|)(g^sk(s)^+)_k (g^sl(s))_l 
! 4c integration for primitive gaussians.
! Note that g depends on large or small components
! g^si = g^L if si=1,2
!      = g^S if si=3,4
! 2013.11.07 written by fukuda
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none

  real(kind=dp) :: inttwoelegrad1_pg

  integer,intent(in) :: i,j,k,l ! p.g. indice
  integer,intent(in) :: m ! vec component 1-3
  integer,intent(in) :: s_i,s_j,s_k,s_l ! spinor indice
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  
  integer :: NLS
  real(kind=dp) :: twoelea,twoeleb
  real(kind=dp) :: posA(3), posB(3), posC(3), posD(3) ! position of center of gaussian 
  real(kind=dp) :: alphaA, alphaB, alphaC, alphaD ! exponents
  integer :: n1(3),nbar1(3),n2(3),nbar2(3)  
  integer :: n1a(3), n1b(3)
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg_a, f_norm_pg_b

  call set_pg(s_i,i,NL,NS,pg,posA,alphaA,n1(1),n1(2),n1(3))
  call set_pg(s_j,j,NL,NS,pg,posB,alphaB,nbar1(1),nbar1(2),nbar1(3))
  call set_pg(s_k,k,NL,NS,pg,posC,alphaC,n2(1),n2(2),n2(3))
  call set_pg(s_l,l,NL,NS,pg,posD,alphaD,nbar2(1),nbar2(2),nbar2(3))

  n1a(:)=n1(:)
  n1b(:)=n1(:)
  if(m.lt.1 .or. m.gt.3) then
     write(*,*) "m should be 1-3 in inttwoelegrad1_pg."
     stop
  end if
  n1a(m)=n1(m)-1
  n1b(m)=n1(m)+1
  
  if(n1(m).eq.0) then
    twoelea=0.d0; f_norm_pg_a=0.d0

  else

  call gauss_int_twoele(posA,alphaA,n1a,posB,alphaB,nbar1,posC,alphaC,n2,posD,alphaD,nbar2,twoelea)  ! (AB|CD)
  f_norm_pg_a = norm_pg(alphaA,n1a(1),n1a(2),n1a(3))*norm_pg(alphaB,nbar1(1),nbar1(2),nbar1(3)) &
       *norm_pg(alphaC,n2(1),n2(2),n2(3))*norm_pg(alphaD,nbar2(1),nbar2(2),nbar2(3))
  end if

  call gauss_int_twoele(posA,alphaA,n1b,posB,alphaB,nbar1,posC,alphaC,n2,posD,alphaD,nbar2,twoeleb)  ! (AB|CD)
  f_norm_pg_b = norm_pg(alphaA,n1b(1),n1b(2),n1b(3))*norm_pg(alphaB,nbar1(1),nbar1(2),nbar1(3)) &
       *norm_pg(alphaC,n2(1),n2(2),n2(3))*norm_pg(alphaD,nbar2(1),nbar2(2),nbar2(3))

  inttwoelegrad1_pg = n1(m)*twoelea*f_norm_pg_a - 2.d0*alphaA*twoeleb*f_norm_pg_b

!  write(*,*) 'inttwoelegrad1_pg',inttwoelegrad1_pg

  return
end function inttwoelegrad1_pg

!=========================================================================
function inttwoelegrad_pg(m,i,j,k,l,s_i,s_j,s_k,s_l,NL,NS,pg)
! (NM|^|PQ)^m = Int dr ds (g^si(r)^+)_i (g^sj(r))_j ((r-s)^k/|r-s|^3)(g^sk(s)^+)_k (g^sl(s))_l 
!             = - Int dr ds \partial_m[(g^si(r)^+)_i (g^sj(r))_j] (1/|r-s|)(g^sk(s)^+)_k (g^sl(s))_l 
! 4c integration for primitive gaussians.
! Note that g depends on large or small components
! g^si = g^L if si=1,2
!      = g^S if si=3,4
! 2013.11.07 
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none

  real(kind=dp) :: inttwoelegrad_pg

  integer,intent(in) :: i,j,k,l ! p.g. indice
  integer,intent(in) :: m ! vec component 1-3
  integer,intent(in) :: s_i,s_j,s_k,s_l ! spinor indice
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg

  real(kind=dp) :: inttwoelegrad1_pg !function

  if(i==j) then
    inttwoelegrad_pg = - 2.d0*inttwoelegrad1_pg(m,i,j,k,l,s_i,s_j,s_k,s_l,NL,NS,pg)
  else
    inttwoelegrad_pg = - inttwoelegrad1_pg(m,i,j,k,l,s_i,s_j,s_k,s_l,NL,NS,pg) &
                       - inttwoelegrad1_pg(m,j,i,k,l,s_j,s_i,s_k,s_l,NL,NS,pg) 
  end if

!  write(*,*) 'inttwoelegrad_pg',inttwoelegrad_pg

  return
end function inttwoelegrad_pg

!=========================================================================
function intEeff_EDM_ob_pg(i,j,s_i,s_j,NL,NS,pg)
! Int ds (g^si(r)^+)_i (nabla^2)  (g^sj(r))_j 
! Note that g depends on large or small components
! g^si = g^L if si=1,2
!      = g^S if si=3,4
!=========================================================================
  use Precision
  use DefineTypes
  implicit none

  real(kind=dp) :: intEeff_EDM_ob_pg

  integer,intent(in) :: i,j ! p.g. indice
  integer,intent(in) :: s_i,s_j ! spinor indice
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg

  real(kind=8) :: ene_kin
  real(kind=dp) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=dp) :: alphaA, alphaB ! exponents
  integer :: nx,nbarx,ny,nbary,nz,nbarz
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg

  call set_pg(s_i,i,NL,NS,pg,posA,alphaA,nx,ny,nz)
  call set_pg(s_j,j,NL,NS,pg,posB,alphaB,nbarx,nbary,nbarz)
  call gauss_int_KE(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,ene_kin) ! nonrelativistic KE (for electron)
  ! normalizatoin does not depende on p.g.
  f_norm_pg = norm_pg(alphaA,nx,ny,nz)*norm_pg(alphaB,nbarx,nbary,nbarz)

  intEeff_EDM_ob_pg = - 2._dp* f_norm_pg*ene_kin

    return
end function intEeff_EDM_ob_pg

!===============================================================================================
subroutine gauss_int_2nd_moment(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,posC,moment2)
! [input]
! posA : position of primitive gaussian A (unnormalized) 
! alphaA : exponent of A
! nx,ny,nz : type of p.g. A 
! posB : position of primitive gaussian B (unnormalized) 
! alphaB : exponenet of B
! nbarx,nbary,nbarz : type of p.g. B
! posC : position related to 2nd_moment
! 
! [output]
! moment2 : 2nd moment integral  (A|xC^2|B)  (xC = x-Cx) , (A|yC^2|B), (A|zC^2|B), (A|xCyC|B),...  p.223 of M&D
!===============================================================================================
  implicit none

  real(kind=8) :: PI
  
  real(kind=8),intent(in) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=8),intent(in) :: alphaA, alphaB ! exponents
  real(kind=8),intent(in) :: posC(3) 
  integer,intent(in) :: nx,nbarx,ny,nbary,nz,nbarz

  real(kind=8),intent(out) :: moment2(6)

  real(kind=8) :: c_PAB ! prefactor when two gaussians are combined
  real(kind=8) :: alphaP, posP(3) ! exponents and position of combined gaussian
  real(kind=8) :: vecPA(3), vecPB(3) ! P-A, P-B
  real(kind=8) :: vecPC(3) ! P-C
  real(kind=8) :: T

  integer :: i,j,k,l
  
  real(kind=8), allocatable, dimension(:,:,:) :: d,e,f
  
  integer :: n_sum
  integer :: nx_sum,ny_sum,nz_sum  ! upper bounds of N,L,M
  
  real(kind=8) :: int_2mom(6), dijk
  real(kind=8) :: sum(6)
  real(kind=8) :: norm

  PI = atan(1.d0)*4.d0
  
  nx_sum = nx+nbarx
  ny_sum = ny+nbary
  nz_sum = nz+nbarz
  n_sum = nx_sum+ny_sum+nz_sum

  call calc_gaussP(posA,alphaA,posB,alphaB,c_PAB,posP,alphaP,vecPA,vecPB)
  
  norm = 1.d0 ! unnormalized

  allocate(d(0:nx,0:nbarx,0:nx_sum))
  allocate(e(0:ny,0:nbary,0:ny_sum))
  allocate(f(0:nz,0:nbarz,0:nz_sum))

  call calc_d(alphaP,vecPA(1),vecPB(1),nx,nbarx,d)
  call calc_d(alphaP,vecPA(2),vecPB(2),ny,nbary,e)
  call calc_d(alphaP,vecPA(3),vecPB(3),nz,nbarz,f)

  call calc_PC_and_T(posC,posP,alphaP,vecPC,T)  ! compute vecPC. we do not use T
!  write(*,*) vecPC
!  stop

  !-------------------------------------
  ! 2nd moment integral
  !-------------------------------------
  do l=1,6
     sum(l)=0.d0
  end do
  do i=0,nx_sum
     do j=0,ny_sum
        do k=0,nz_sum
           call calc_2nd_moment_integral(i,j,k,alphaP,vecPC,int_2mom)
           dijk = d(nx,nbarx,i)*e(ny,nbary,j)*f(nz,nbarz,k)
           do l=1,6
              sum(l) = sum(l) +dijk*int_2mom(l)
           end do
        end do
     end do
  end do

  do l=1,6
     moment2(l) = c_PAB*norm**2*sum(l)
!     write(*,*) l,moment2(l)
  end do
!  stop

  return
end subroutine gauss_int_2nd_moment

!========================================================================
subroutine calc_2nd_moment_integral(i,j,k,alphaP,vecPC,int_2mom)
! calculation of 
! int_2mom(1-3): [NLM|xC^2], [NLM|yC^2], [NLM|zC^2] in McMurchie & Davidson Eq. (3.7) for 2nd moment integral
! int_2mom(4-5): [NLM|xCyC], [NLM|yCzC], [NLM|zCxC] in McMurchie & Davidson Eq. (3.8) for 2nd moment integral
!========================================================================
  implicit none
  integer,intent(in) :: i,j,k 
  real(kind=8),intent(in) :: alphaP,vecPC(3)
  real(kind=8),intent(out) :: int_2mom(6)
  real(kind=8) :: PI

  PI = atan(1.d0)*4.d0

  ! [NLM|xC^2]
  if((j.eq.0).and.(k.eq.0)) then
     if(i.eq.0) then
        int_2mom(1) = (vecPC(1)**2 + 1.0d0/(2.d0*alphaP))*(PI/alphaP)**1.5d0
     elseif(i.eq.1) then
        int_2mom(1) = 2.d0*vecPC(1)*(PI/alphaP)**1.5d0
     elseif(i.eq.1) then
        int_2mom(1) = 2.d0*(PI/alphaP)**1.5d0
     else
        int_2mom(1) = 0.d0
     end if
  else
     int_2mom(1) = 0.d0
  end if

  ! [NLM|yC^2]
  if((i.eq.0).and.(k.eq.0)) then
     if(j.eq.0) then
        int_2mom(2) = (vecPC(2)**2 + 1.0d0/(2.d0*alphaP))*(PI/alphaP)**1.5d0
     elseif(j.eq.1) then
        int_2mom(2) = 2.d0*vecPC(2)*(PI/alphaP)**1.5d0
     elseif(j.eq.2) then
        int_2mom(2) = 2.d0*(PI/alphaP)**1.5d0
     else
        int_2mom(2) = 0.d0
     end if
  else
     int_2mom(2) = 0.d0
  end if

  ! [NLM|zC^2]
  if((i.eq.0).and.(j.eq.0)) then
     if(k.eq.0) then
        int_2mom(3) = (vecPC(3)**2 + 1.0d0/(2.d0*alphaP))*(PI/alphaP)**1.5d0
     elseif(k.eq.1) then
        int_2mom(3) = 2.d0*vecPC(3)*(PI/alphaP)**1.5d0
     elseif(k.eq.2) then
        int_2mom(3) = 2.d0*(PI/alphaP)**1.5d0
     else
        int_2mom(3) = 0.d0
     end if
  else
     int_2mom(3) = 0.d0
  end if

  ! [NLM|xCyC]
  if(k.eq.0) then
     if((i.eq.0).and.(j.eq.0)) then
        int_2mom(4) = (vecPC(1)*vecPC(2))*(PI/alphaP)**1.5d0
     elseif((i.eq.1).and.(j.eq.0)) then
        int_2mom(4) = vecPC(2)*(PI/alphaP)**1.5d0
     elseif((i.eq.0).and.(j.eq.1)) then
        int_2mom(4) = vecPC(1)*(PI/alphaP)**1.5d0
     elseif((i.eq.1).and.(j.eq.1)) then
        int_2mom(4) = (PI/alphaP)**1.5d0
     else
        int_2mom(4) = 0.d0
     end if
  else
     int_2mom(4) = 0.d0
  end if

  ! [NLM|yCzC]
  if(i.eq.0) then
     if((j.eq.0).and.(k.eq.0)) then
        int_2mom(5) = (vecPC(2)*vecPC(3))*(PI/alphaP)**1.5d0
     elseif((j.eq.1).and.(k.eq.0)) then
        int_2mom(5) = vecPC(3)*(PI/alphaP)**1.5d0
     elseif((j.eq.0).and.(k.eq.1)) then
        int_2mom(5) = vecPC(2)*(PI/alphaP)**1.5d0
     elseif((j.eq.1).and.(k.eq.1)) then
        int_2mom(5) = (PI/alphaP)**1.5d0
     else
        int_2mom(5) = 0.d0
     end if
  else
     int_2mom(5) = 0.d0
  end if

  ! [NLM|zCxC]
  if(j.eq.0) then
     if((k.eq.0).and.(i.eq.0)) then
        int_2mom(6) = (vecPC(3)*vecPC(1))*(PI/alphaP)**1.5d0
     elseif((k.eq.1).and.(i.eq.0)) then
        int_2mom(6) = vecPC(1)*(PI/alphaP)**1.5d0
     elseif((k.eq.0).and.(i.eq.1)) then
        int_2mom(6) = vecPC(3)*(PI/alphaP)**1.5d0
     elseif((k.eq.1).and.(i.eq.1)) then
        int_2mom(6) = (PI/alphaP)**1.5d0
     else
        int_2mom(6) = 0.d0
     end if
  else
     int_2mom(6) = 0.d0
  end if

  return
end subroutine calc_2nd_moment_integral


!==================================================================
function func_d2pg(i,j,x,y,z,posA,alpha,nx,ny,nz)
! i(=1,2,3=x,y,z) derivative of normalized gaussian function
! center at (Ax,Ay,Az)
!==================================================================
  implicit none

  real(kind=8) :: func_d2pg
  integer,intent(in) :: i,j
  real(kind=8),intent(in) :: x,y,z,posA(3),alpha
  integer,intent(in) :: nx,ny,nz
  real(kind=8) :: norm_pg,func_d2unpg

  func_d2pg = norm_pg(alpha,nx,ny,nz) *func_d2unpg(i,j,x,y,z,posA,alpha,nx,ny,nz)

  return
end function func_d2pg

!==================================================================
function func_d2unpg(i,j,x,y,z,posA,alpha,nx,ny,nz)
! i,j(=1,2,3=x,y,z) derivative of UNnormalized gaussian function
! center at (Ax,Ay,Az)
!==================================================================
  implicit none

  real(kind=8) :: func_d2unpg
  integer,intent(in) :: i,j
  real(kind=8),intent(in) :: x,y,z,posA(3),alpha
  integer,intent(in) :: nx,ny,nz
  real(kind=8) :: func_unpg
  real(kind=8) :: func_dunpg

  if(i.eq.1) then !d^2/dxdx^j
     func_d2unpg = nx*func_dunpg(j,x,y,z,posA,alpha,nx-1,ny,nz) -2.d0*alpha*func_dunpg(j,x,y,z,posA,alpha,nx+1,ny,nz)
  elseif(i.eq.2) then ! d/dy
     func_d2unpg = ny*func_dunpg(j,x,y,z,posA,alpha,nx,ny-1,nz) -2.d0*alpha*func_dunpg(j,x,y,z,posA,alpha,nx,ny+1,nz)
  elseif(i.eq.3) then ! d/dz
     func_d2unpg = nz*func_dunpg(j,x,y,z,posA,alpha,nx,ny,nz-1) -2.d0*alpha*func_dunpg(j,x,y,z,posA,alpha,nx,ny,nz+1)
  else
     stop "i should be 1 or 2 or 3 in fund_dunpg."
  end if

  return
end function func_d2unpg
