!================================================================
! subroutines to compute expectation values of operators of 
! physical quantities from density matrices.
!
! 2013.6.22
!
!================================================================

!==============================================================================================
!subroutine calc_energy_offset_BO(calE,h_Qmat,twoele_Qmat,calFj0_Qmat,alpha,energy_offset)
subroutine calc_energy_offset_BO(calE,calE_NO0,h_Qmat,twoele_Qmat,energy_offset)
!  calculate offset of (integrated) energy with BO approximation
!  <Phi|H|Phi> - <Phi|:H:|Phi>
!
!  130711
!==============================================================================================
  use Precision
  use DiracOutput
  use Constants
  implicit none

  complex(kind=dp),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ = <Phi|e^dagger_P e_Q|Phi> (t=0)
  complex(kind=dp),intent(in) :: calE_NO0(4*NBS,4*NBS)  ! <Phi|:e^dagger_P e_Q :|Phi> (t=0)
  complex(kind=dp),intent(in) :: h_Qmat(4*NBS,4*NBS)  ! h_PQ
  complex(kind=dp),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)
!  complex(kind=dp),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
!  complex(kind=dp),intent(in) :: alpha(Nph)   ! <coherent|hat(a)|coherent>

  complex(kind=dp),intent(out) :: energy_offset

!  real(kind=dp) :: dpj ! tilde(Delta p)_j
!  real(kind=dp) :: p0j ! p_j

  integer :: pp,qq,rr,ss,nn,mm,j
  complex(kind=dp) :: sum
  complex(kind=dp) :: sumee
  

  sum = (0._dp,0._dp)
  ! h_NM *(<Phi|calE_NM|Phi> -<Phi|:calE_NM:|Phi> )
  sum = (0._dp,0._dp)
  do nn=1,4*NBS
     do mm=1,4*NBS
        sum = sum +h_Qmat(nn,mm)*(calE(nn,mm)-calE_NO0(nn,mm))
     end do
  end do

  sumee = (0._dp,0._dp)
  ! (1/2)(NQ|MP)*(<Phi|e^dagger_N e^dagger_M e_P e_Q |Phi> -<Phi|:e^dagger_N e^dagger_M e_P e_Q:|Phi> )
  do nn=1,4*NBS
     do mm=1,4*NBS
        do pp=1,4*NBS
           do qq=1,4*NBS
              sumee = sumee +twoele_Qmat(nn,qq,mm,pp)* (calE(nn,qq)*calE(mm,pp)-calE(nn,pp)*calE(mm,qq)) &
                   &  -twoele_Qmat(nn,qq,mm,pp)* (calE_NO0(nn,qq)*calE_NO0(mm,pp)-calE_NO0(nn,pp)*calE_NO0(mm,qq))
           end do
        end do
     end do
  end do

  energy_offset = sum +sumee/2._dp

  return
end subroutine calc_energy_offset_BO



!==============================================================================================
!subroutine calc_energy_eae_BO(t,calE,calEa,h_Qmat,calFj0_Qmat,energy)
subroutine calc_energy_eae_BO(t,calE,calEa,h_Qmat,twoele_Qmat,calFj0_Qmat,energy)
!  calculate (integrated) energy with BO approximation
!  using contribution of eae
!  A_A is neglected. Two-electron term is neglected.
!  made from calc_derivatives_BO_N
!
!  130711 
!==============================================================================================
  use Precision
  use DiracOutput
  use Constants
  implicit none

  real(kind=dp),intent(in) :: t  ! time
  complex(kind=dp),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  complex(kind=dp),intent(in) :: calEa(4*NBS,4*NBS,Nph) ! calEa_PQj
  complex(kind=dp),intent(in) :: h_Qmat(4*NBS,4*NBS)  ! h_PQ
  complex(kind=dp),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
  complex(kind=dp),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)

  complex(kind=dp),intent(out) :: energy

  real(kind=dp) :: dp0_th_phi ! dp0*dtheta*dphi
  real(kind=dp) :: p0j,thj,phij ! p_j,theta_j,phi_j

  integer :: i,j,k,l
  integer :: pp,qq,rr,ss,nn,mm
  complex(kind=dp) :: sum

  sum = (0._dp,0._dp)

  ! h_NM * calE_NM
  sum = (0._dp,0._dp)
  do nn=1,4*NBS
     do mm=1,4*NBS
        sum = sum +h_Qmat(nn,mm)*calE(nn,mm)
     end do
  end do

  ! -1/sqrt(2 pi^2 c)* int d^3p/sqrt(2*p0) 
  ! ( F^k_NM(p)e^k(p) exp(-i c p0 t) calEa_NMj + F^*k_MN(p)e^*k(p) exp(i c p0 t) calEa^*_MNj
  ! be careful about the order of subscripts.
  do nn=1,4*NBS
     do mm=1,4*NBS
        do j=1,Nph
           call get_pj(j,p0j,thj,phij,dp0_th_phi)
           sum = sum -dp0_th_phi*p0j**2*sin(thj)/( 2._dp*PI*sqrt(CCC)*sqrt(p0j) ) &
                & *( calFj0_Qmat(nn,mm,j)*exp(-IU*CCC*p0j*t)*calEa(nn,mm,j) &
                &      + conjg(calFj0_Qmat(mm,nn,j))*exp(IU*CCC*p0j*t)*conjg(calEa(mm,nn,j)) )
        end do
     end do
  end do

  ! (1/2)(NQ|MP) calE_NMPQ  --> calE_NMPQ is approximated to be calE_NQ*calE_MP -calE_NP*calE_MQ 
  do nn=1,4*NBS
     do mm=1,4*NBS
        do pp=1,4*NBS
           do qq=1,4*NBS
              sum = sum +twoele_Qmat(nn,qq,mm,pp)*(calE(nn,qq)*calE(mm,pp)-calE(nn,pp)*calE(mm,qq))/2._dp
!              write(*,"(4i6,2es16.6)") nn,mm,pp,qq,twoele_Qmat(nn,qq,mm,pp)
           end do
        end do
     end do
  end do

  energy = sum

  return
end subroutine calc_energy_eae_BO


!==============================================================================================
!subroutine calc_energy_BO(t,calE,h_Qmat,calFj0_Qmat,alpha,energy)
subroutine calc_energy_BO(t,calE,h_Qmat,twoele_Qmat,calFj0_Qmat,alpha,energy)
!  calculate (integrated) energy with BO approximation
!  do not use contribution of eae (do not consider coherent state)
!  A_A is neglected. Two-electron term is neglected.
!  made from calc_derivatives_BO_N
!
!  130630 
!==============================================================================================
  use Precision
  use DiracOutput
  use Constants
  implicit none

  real(kind=dp),intent(in) :: t  ! time
  complex(kind=dp),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  complex(kind=dp),intent(in) :: h_Qmat(4*NBS,4*NBS)  ! h_PQ
  complex(kind=dp),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)
  complex(kind=dp),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
  complex(kind=dp),intent(in) :: alpha(Nph)   ! <coherent|hat(a)|coherent>

  complex(kind=dp),intent(out) :: energy

  real(kind=dp) :: dpj ! tilde(Delta p)_j
  real(kind=dp) :: p0j ! p_j

  integer :: pp,qq,rr,ss,nn,mm,j
  complex(kind=dp) :: sum

  sum = (0._dp,0._dp)

  ! h_NM * calE_NM
  sum = (0._dp,0._dp)
  do nn=1,4*NBS
     do mm=1,4*NBS
        sum = sum +h_Qmat(nn,mm)*calE(nn,mm)
     end do
  end do
!  write(*,*) sum

  ! contribution from coherent state
  do nn=1,4*NBS
     do mm=1,4*NBS
        do j=1,Nph
           call calc_DeltaPj_and_P0j(j,dpj,p0j) ! sub_int.f90
           sum = sum +( calFj0_Qmat(nn,mm,j)       *exp(-IU*CCC*p0j*t)*alpha(j) &
                & +conjg(calFj0_Qmat(mm,nn,j))*exp(+IU*CCC*p0j*t)*conjg(alpha(j)) )*calE(nn,mm)/(2._dp*PI*sqrt(p0j*CCC))
        end do
     end do
  end do

  ! (1/2)(NQ|MP) calE_NMPQ  --> calE_NMPQ is approximated to be calE_NQ*calE_MP -calE_NP*calE_MQ 
  do nn=1,4*NBS
     do mm=1,4*NBS
        do pp=1,4*NBS
           do qq=1,4*NBS
              sum = sum +twoele_Qmat(nn,qq,mm,pp)*(calE(nn,qq)*calE(mm,pp)-calE(nn,pp)*calE(mm,qq))/2._dp
!              write(*,"(4i6,2es16.6)") nn,mm,pp,qq,twoele_Qmat(nn,qq,mm,pp)
           end do
        end do
     end do
  end do
!  write(*,*) sum

  energy = sum

  return
end subroutine calc_energy_BO

!================================================================
function dAraddt_4E(k,t,vecR,alpha)
! k : vector index 1~3
! t : time
! vecR : position
! alpha : <coherent|hat(a)|coherent>
! (-1/c)dArad^k/dt 
! Note that (-1/c) is included
! 121106
!================================================================
  use Precision
  use DiracOutput
  use Constants
  implicit none
  
  real(kind=dp) :: dAraddt_4E
  integer,intent(in) :: k  ! k : vector index 1~3
  real(kind=dp),intent(in) :: t   ! t : time
  real(kind=dp),intent(in) :: vecR(3)   ! vecR : position
  complex(kind=dp),intent(in) :: alpha(Nph)   ! <coherent|hat(a)|coherent>

  integer :: j
  real(kind=dp) :: p0j,vecPj(3) ! p^0_j, vec{p}_j
  complex(kind=dp) :: vecpolj(3) ! (circular) polarization vector (depends on sig)  e^k_j
  complex(kind=dp) :: sum
  real(kind=dp) :: PdotR ! vec{p}.vec{r}

  sum=(0._dp,0._dp)
  do j=1,Nph
     call calc_mode_pj(j,p0j,vecPj,vecpolj) ! in sub_int.f90
     PdotR = dot_product(vecPj,vecR)
!     write(*,*) k,j,vecPj,vecR,PdotR
     sum = sum + sqrt(p0j)*alpha(j)*vecpolj(k)*exp(-IU*CCC*p0j*t)*exp(IU*PdotR)
  end do
  
!  write(*,*) k,IU*sqrt(CCC)*(sum-conjg(sum))/2._dp/PI ,dp
     
  dAraddt_4E = real(IU*sqrt(CCC)*(sum-conjg(sum))/2._dp/PI ,dp)

  return
end function dAraddt_4E


!================================================================
subroutine calc_pol_BO(vecR,calE,E_Qmat,pol)
! Polarization density (for BO approximation)
! 121102
!================================================================
  use Precision
  use DiracOutput
  use Constants
  implicit none
  
!  real(kind=dp),intent(in) :: x,y,z
  real(kind=dp),intent(in) :: vecR(3)
  complex(kind=dp),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  complex(kind=dp),intent(in) :: E_Qmat(3,4*NBS,4*NBS)  ! E^k_NM
  real(kind=dp), intent(out) :: pol(3)
  integer :: k,nn,mm
  complex(kind=dp) :: sum
  real(kind=8) :: posRA(3),pol_nuc
  integer :: i,iN
  real(kind=dp) :: vecD(3) ! vecR-vecRA
  real(kind=dp) :: norm !  |vecR-vecRA|

  do k=1,3
     ! contribution from electron 
     sum = (0._dp,0._dp)
     do nn=1,4*NBS
        do mm=1,4*NBS
           sum = sum -E_Qmat(k,nn,mm)*calE(nn,mm)
        end do
     end do
!     write(*,*) k,sum
     pol(k) = real(sum,dp)

     ! contribution from nucleus
     pol_nuc = 0._dp
     do iN=1,NAT
        posRA(1) = xc(iN); posRA(2) = yc(iN); posRA(3) = zc(iN)
        vecD(:) = vecR(:)-posRA(:)
        norm = sqrt(vecD(1)**2+vecD(2)**2+vecD(3)**2)
        pol_nuc = pol_nuc +cn(iN)*vecD(k)/norm**3
!        write(*,*) iN,posRA(:),vecD(:),norm
     end do
     pol_nuc = -pol_nuc/(4._dp*PI)

!     write(*,*) k,pol(k),pol_nuc
     pol(k) = pol(k) +pol_nuc
  end do
  
  return
end subroutine calc_pol_BO


!================================================================
subroutine calc_rhonuc_tot(calC,rhonuc_tot)
! total nuclear charge density 
! (integrated over whole space -> trace of d.m.)
!================================================================
  use Precision
  use NucBasis  ! for NBS_N,Z_N
  implicit none
  
  complex(kind=dp),intent(in) :: calC(NBS_N,NBS_N)  ! calC_ij
  real(kind=dp), intent(out) :: rhonuc_tot
  integer :: nn
  complex(kind=dp) :: sum

  sum = (0._dp,0._dp)
  do nn=1,NBS_N
     sum = sum +calC(nn,nn)
  end do

  rhonuc_tot = real(sum,dp)
  
  return
end subroutine calc_rhonuc_tot

!================================================================
subroutine calc_rhonuc(x,y,z,calC,rhonuc)
! nuclear charge density (should be real)
!================================================================
  use Precision
  use NucBasis  ! for NBS_N
  implicit none
  
  real(kind=dp),intent(in) :: x,y,z
  complex(kind=dp),intent(in) :: calC(NBS_N,NBS_N)  ! calC_ij
  real(kind=dp), intent(out) :: rhonuc
  integer :: nn,mm
  complex(kind=dp) :: sum
  complex(kind=dp) :: rho_nuc_mat

  sum = (0._dp,0._dp)
  do nn=1,NBS_N
     do mm=1,NBS_N
        sum = sum +rho_nuc_mat(x,y,z,nn,mm)*calC(nn,mm)
     end do
  end do

  rhonuc = real(sum,dp)
  
  return
end subroutine calc_rhonuc

!================================================================
subroutine calc_rhoe_tot(calE,rhoe_tot)
! total electron charge density (integrated over whole space)
! (should be conserved)
! trace of density matrix
!================================================================
  use Precision
  use DiracOutput
  use Constants ! for Ze=-1
  implicit none
  
  complex(kind=dp),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  real(kind=dp), intent(out) :: rhoe_tot
  integer :: nn
  complex(kind=dp) :: sum

  sum = (0._dp,0._dp)
  do nn=1,4*NBS
     sum = sum +calE(nn,nn)
  end do
  
  rhoe_tot = Ze*real(sum,dp)
  
  return
end subroutine calc_rhoe_tot

!================================================================
subroutine calc_rhoe(x,y,z,calE,rhoe)
! electron charge density
! (should be real)
!================================================================
  use Precision
  use DiracOutput
  implicit none
  
  real(kind=dp),intent(in) :: x,y,z
  complex(kind=dp),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  real(kind=dp), intent(out) :: rhoe
  integer :: nn,mm
  complex(kind=dp) :: sum
  complex(kind=dp) :: rho_Qmat 

  sum = (0._dp,0._dp)
  do nn=1,4*NBS
     do mm=1,4*NBS
        sum = sum +rho_Qmat(x,y,z,nn,mm)*calE(nn,mm)
     end do
  end do
  
  rhoe = real(sum,dp)
  
  return
end subroutine calc_rhoe

!================================================================
subroutine calc_tauS_eigvec(x,y,z,calE,eig,vec)
! symmetrize the stress tensor 
! eigenvalues and eigenvectors of real symmetric electronic stress tensor 
!
!================================================================
  use Precision
  use DiracOutput
  implicit none
  
  real(kind=dp),intent(in) :: x,y,z
  complex(kind=dp),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  real(kind=dp), intent(out) :: eig(3),vec(3,3)
  complex(kind=dp) :: ctau_kl, ctau(3,3)
  real(kind=dp) :: tau_kl, tauT_kl, tauS(3,3)
  integer :: k,l

  do k=1,3
     do l=1,3
        call calc_tau_exval(k,l,x,y,z,calE,ctau_kl)
        ctau(k,l) = ctau_kl
     end do
  end do

  ! calculate Hernitian part
  ! tauPI = tauS + tauA, S^T = S, A^T = -A
  ! PI = S+A, PI^T = S-A, S=(PI+PI^T)/2
  do k=1,3
     do l=1,3
        tau_kl  = real(ctau(k,l),dp)
        tauT_kl = real(ctau(l,k),dp)
        tauS(k,l) = (tau_kl+tauT_kl)/2._dp
!        write(*,*) k,l,tauS(k,l)
     end do
  end do
 
  call calc_eigen(tauS,eig,vec)  ! defined in sub_eigsym.f90
!  call calc_eig_ch(3,tauS,eig,vec)

!  do k=1,3
!     write(*,'(4es16.6)') eig(k),vec(1,k),vec(2,k),vec(3,k) 
!  end do

  return
end subroutine calc_tauS_eigvec

!================================================================
subroutine  calc_tau_exval(k,l,x,y,z,calE,tau_kl)
! expectation value of electronic stress tensor
! this should be real matrix (but not necessarily symmetric)
! but we compute this as complex at this stage.
!
! <tau^PI>^kl = (tau^kl_PQ + A^k*j^l_PQ/c) * calE_PQ
!================================================================
  use Precision
  use DiracOutput
  use Constants ! for CCC
  implicit none
  
  integer,intent(in) :: k,l
  real(kind=dp),intent(in) :: x,y,z
  complex(kind=dp),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  complex(kind=dp),intent(out) :: tau_kl
  integer :: nn,mm
  complex(kind=dp) :: sum
  complex(kind=dp) :: tau_Qmat 

  sum = (0._dp,0._dp)
  do nn=1,4*NBS
     do mm=1,4*NBS
        sum = sum +tau_Qmat(k,l,x,y,z,nn,mm)*calE(nn,mm)
     end do
  end do
  
  tau_kl = sum

  return
end subroutine calc_tau_exval



!=======================================================================================
! the subroutine below is for Hermite matrix.
! may not be used for stress tensor.
!=======================================================================================

!!$!==================================================
!!$subroutine calc_eig_ch(n,mat_ch,eig,vec)
!!$! -calculate eigen values (eig) and eigenvectors (vec) of n by n 
!!$! complex Hermite matrix (mat_ch).
!!$! -sorted like eig(1)<eig(2)<eig(3)< ...
!!$! -it seems that vec is normalized (norm=1)
!!$!
!!$! input: 
!!$!   n: dimension
!!$!   mat_ch(n,n): complex Hermite matrix
!!$! output:
!!$!   eig(i): eigenvalues (real) (i=1~n)
!!$!   vec(n,i): eigenvectors of eig(i) (complex) 
!!$!================================================== 
!!$  use Precision
!!$  implicit none
!!$  
!!$  integer, intent(in) :: n
!!$  complex(kind=dp), intent(in) :: mat_ch(n,n)
!!$
!!$  real(kind=dp), intent(out) :: eig(n)
!!$  complex(kind=dp), intent(out) :: vec(n,n)
!!$
!!$!  integer, parameter :: nm = 4
!!$  integer :: nm  ! we only treat nm=n case
!!$!  integer, parameter :: n = 4
!!$  integer, parameter :: matz = 1 ! 0 for only eigenvalue
!!$  integer :: ierr
!!$!  real(kind=dp) :: ar(nm,n),ai(nm,n)   ! input matrix
!!$  real(kind=dp) :: ar(n,n),ai(n,n)   ! input matrix
!!$  real(kind=dp) :: w(n)  ! eignevalues
!!$!  real(kind=dp) :: zr(nm,n),zi(nm,n)   ! eigenvectors
!!$  real(kind=dp) :: zr(n,n),zi(n,n)   ! eigenvectors
!!$  real(kind=dp) :: fv1(n),fv2(n),fm1(2,n)   ! temporary alleys
!!$
!!$  complex(kind=dp),parameter :: IU=(0._dp,1._dp)
!!$  integer :: i,j
!!$
!!$  nm = n
!!$  do i=1,n
!!$     do j=1,n
!!$        ar(i,j) = real(mat_ch(i,j))
!!$        ai(i,j) = aimag(mat_ch(i,j))
!!$     end do
!!$  end do
!!$   
!!$  call ch(nm,n,ar,ai,w,matz,zr,zi,fv1,fv2,fm1,ierr)
!!$  ! w(i) : i th eigenvector w(1)<w(2)<w(3)< ...
!!$  ! zr(j,i), zi(j,i) : j th component of eigenvector corresponding to w(i)
!!$  
!!$  do i=1,n
!!$     eig(i) = w(i)
!!$  end do
!!$
!!$  do i=1,n
!!$     do j=1,n
!!$        vec(i,j) = zr(i,j) +IU*zi(i,j)
!!$     end do
!!$  end do
!!$
!!$  return
!!$end subroutine calc_eig_ch
