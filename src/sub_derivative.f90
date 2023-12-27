!================================================================================
! subroutines for derivatives of Rigged QED
!
!================================================================================


!========================================================================================================
subroutine calc_derivatives_NonBO_OPQ(use_exchange,t,calE,calEa,calC,calCa, & 
     & TM_Qmat,Tnuc_mat,twoele_Qmat,calFj0_Qmat,nucele_Qmat,twonuc_mat, &
     & dcalE,dcalEa,dcalC,dcalCa)
!  derivatives for Non-BO calculation with Arad interaction
!  140205
!========================================================================================================
  use Precision
  use DiracOutput
  use Constants
  use NucBasis
  implicit none

  logical,intent(in) :: use_exchange  !.true. --> include exchange terms in diff. eq. , .false. --> without exchange terms
  real(kind=dp),intent(in) :: t  ! time
  complex(kind=dp),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  complex(kind=dp),intent(in) :: calEa(4*NBS,4*NBS,Nph) ! calEa_PQj
  complex(kind=dp),intent(in) :: calC(NBS_N,NBS_N)  ! calC_klj
  complex(kind=dp),intent(in) :: calCa(NBS_N,NBS_N,Nph)  ! calC_klj
  complex(kind=dp),intent(in) :: TM_Qmat(4*NBS,4*NBS)  ! T_PQ + M_PQ
  complex(kind=dp),intent(in) :: Tnuc_mat(NBS_N,NBS_N)  ! T_akl
  complex(kind=dp),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)
  complex(kind=dp),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
  complex(kind=dp),intent(in) :: nucele_Qmat(NBS_N,NBS_N,4*NBS,4*NBS)  ! (ij|PQ)
  complex(kind=dp),intent(in) :: twonuc_mat(NBS_N,NBS_N,NBS_N,NBS_N)  ! (ij|kl)

  complex(kind=dp),intent(out) :: dcalE(4*NBS,4*NBS)  ! d(calE_PQ)/dt *DeltaT
  complex(kind=dp),intent(out) :: dcalEa(4*NBS,4*NBS,Nph) ! d(calEa_PQj)/dt *DeltaT
  complex(kind=dp),intent(out) :: dcalC(NBS_N,NBS_N)  ! d(calC_ij)/dt *DeltaT
  complex(kind=dp),intent(out) :: dcalCa(NBS_N,NBS_N,Nph)  ! d(calC_ij)/dt *DeltaT

  real(kind=dp) :: dp0_th_phi ! dp0*dtheta*dphi
  real(kind=dp) :: p0j,thj,phij ! p_j,theta_j,phi_j

  integer :: i,j,k,l,p,q,r,s,n,m
  integer :: ii,jj,kk,ll,pp,qq,rr,ss,nn,mm
  complex(kind=dp) :: sum

  complex(kind=dp) :: O_Qmat(4*NBS,4*NBS), P_Qmat(4*NBS,4*NBS,Nph), Q_Qmat(4*NBS,4*NBS,Nph)
  complex(kind=dp) :: I_Qmat(4*NBS,4*NBS)
  complex(kind=dp) :: O_Nmat(NBS_N,NBS_N), P_Nmat(NBS_N,NBS_N,Nph), Q_Nmat(NBS_N,NBS_N,Nph)
  complex(kind=dp) :: I_Nmat(NBS_N,NBS_N)

  sum = (0._dp,0._dp)
  O_Qmat(:,:) = (0._dp,0._dp)
  P_Qmat(:,:,:) = (0._dp,0._dp)
  Q_Qmat(:,:,:) = (0._dp,0._dp)
  O_Nmat(:,:) = (0._dp,0._dp)
  P_Nmat(:,:,:) = (0._dp,0._dp)
  Q_Nmat(:,:,:) = (0._dp,0._dp)
  
  !--------------------------------------------
  !
  ! Electron part: O_Qmat, P_Qmat, Q_Qmat
  !
  !--------------------------------------------

  !---------------------------------------
  !   calculation of O_NM
  !---------------------------------------
  !----- (T+M)_MR calE_NR -----

  do nn=1,4*NBS
     do mm=1,4*NBS
        sum = (0._dp,0._dp)
        do rr=1,4*NBS
           sum = sum +TM_Qmat(mm,rr)*calE(nn,rr)
        end do
        O_Qmat(nn,mm) = sum/IU
     end do
  end do

  !----- (MR|PQ) calE_NR calE_PQ -----
  I_Qmat(:,:) = (0._dp,0._dp)
  do mm=1,4*NBS
     do rr=1,4*NBS
        sum = (0._dp,0._dp)
        do pp=1,4*NBS
           do qq=1,4*NBS
              sum = sum + twoele_Qmat(mm,rr,pp,qq)*calE(pp,qq)
           end do
        end do
        I_Qmat(mm,rr) = sum
     end do
  end do
  
  ! add to O^Phi_NM
  do nn=1,4*NBS
     do mm=1,4*NBS
        sum = (0._dp,0._dp)
        do rr=1,4*NBS
           sum = sum + I_Qmat(mm,rr)*calE(nn,rr)
        end do
        O_Qmat(nn,mm) = O_Qmat(nn,mm) + sum/IU
     end do
  end do
  
  if(use_exchange) then
     !-----  -(MR|PQ) calE_NQ calE_PR -----
     I_Qmat(:,:) = (0._dp,0._dp)
     do mm=1,4*NBS
        do qq=1,4*NBS
           sum = (0._dp,0._dp)
           do pp=1,4*NBS
              do rr=1,4*NBS
                 sum = sum +twoele_Qmat(mm,rr,pp,qq)*calE(pp,rr)
              end do
           end do
           I_Qmat(mm,qq) = sum
        end do
     end do
     
     ! add to O^Phi_NM
     do nn=1,4*NBS
        do mm=1,4*NBS
           sum = (0._dp,0._dp)
           do qq=1,4*NBS
              sum = sum + I_Qmat(mm,qq)*calE(nn,qq)
           end do
           O_Qmat(nn,mm) = O_Qmat(nn,mm) - sum/IU  ! note the minus sign here.
        end do
     end do
     
  end if

  !----- (MR|ij) calE_NR calC_ij -----
  I_Qmat(:,:) = (0._dp,0._dp)
  do mm=1,4*NBS
     do rr=1,4*NBS
        sum = (0._dp,0._dp)
        do ii=1,NBS_N
           do jj=1,NBS_N
              sum = sum + nucele_Qmat(ii,jj,mm,rr)*calC(ii,jj)
           end do
        end do
        I_Qmat(mm,rr) = sum
     end do
  end do

  ! add to O^Phi_NM
  do nn=1,4*NBS
     do mm=1,4*NBS
        sum = (0._dp,0._dp)
        do rr=1,4*NBS
           sum = sum + I_Qmat(mm,rr)*calE(nn,rr)
        end do
        O_Qmat(nn,mm) = O_Qmat(nn,mm) + sum/IU
     end do
  end do
 
  !-----  -1/sqrt(2 pi^2 c)* int d^3p/sqrt(2*p0) 
  !        ( F^k_MR(p)e^k(p) exp(-i c p0 t) calEa_NRj + F^*k_RM(p)e^*k(p) exp(i c p0 t) calEa^*_RNj -----
  do nn=1,4*NBS
     do mm=1,4*NBS

        sum = (0._dp,0._dp)
        do j=1,Nph
           call get_pj(j,p0j,thj,phij,dp0_th_phi)
           do rr=1,4*NBS
              sum = sum -dp0_th_phi*p0j**2*sin(thj)/(2._dp*PI*sqrt(CCC)*sqrt(p0j)) *( calFj0_Qmat(mm,rr,j)*exp(-IU*CCC*p0j*t)*calEa(nn,rr,j) &
                   & + conjg(calFj0_Qmat(rr,mm,j))*exp(IU*CCC*p0j*t)*conjg(calEa(rr,nn,j)) )
              sum = sum -calFj0_Qmat(mm,rr,j)
           end do
        end do

        O_Qmat(nn,mm) = O_Qmat(nn,mm) + sum/IU

     end do
  end do

  !---------------------------------------
  !   calculation of P^N_NM (n=0)
  !---------------------------------------
  !----- (T+M)_MR calEa_NRj -----
  do nn=1,4*NBS
     do mm=1,4*NBS
        do j=1,Nph
           sum = (0._dp,0._dp)
           do rr=1,4*NBS
              sum = sum +TM_Qmat(mm,rr)*calEa(nn,rr,j)
           end do
           P_Qmat(nn,mm,j) = sum/IU
        end do
     end do
  end do

  do nn=1,4*NBS
     do mm=1,4*NBS
        do j=1,Nph

           call get_pj(j,p0j,thj,phij,dp0_th_phi)
           sum = (0._dp,0._dp)
           do rr=1,4*NBS
              sum = sum -1._dp/(2._dp*PI*sqrt(CCC)*sqrt(p0j))*conjg(calFj0_Qmat(rr,mm,j))*exp(IU*CCC*p0j*t)*calE(nn,rr)
           end do
           P_Qmat(nn,mm,j) = P_Qmat(nn,mm,j) + sum/IU

        end do
     end do
  end do

  !---------------------------------------
  !   calculation of Q^N_NM (n=0)
  !---------------------------------------
  !----- (T+M)_MR calEa^*_RNj -----
  do nn=1,4*NBS
     do mm=1,4*NBS
        do j=1,Nph
           sum = (0._dp,0._dp)
           do rr=1,4*NBS
              sum = sum +TM_Qmat(mm,rr)*conjg(calEa(rr,nn,j))
           end do
           Q_Qmat(nn,mm,j) = sum/IU
        end do
     end do
  end do

  !---------------------------------------
  ! calculation of derivatives(*DeltaT)
  !---------------------------------------
  !----- d(calE_NM)/dt = O_NM + O_NM^dagger -----
  do nn=1,4*NBS
     do mm=1,4*NBS
        dcalE(nn,mm) = (O_Qmat(nn,mm) +conjg(O_Qmat(mm,nn)))*DeltaT
     end do
  end do

  !----- d(calEa_NMj)/dt = (Q_NMj)^dagger + P_NMj -----
  do nn=1,4*NBS
     do mm=1,4*NBS
        do j=1,Nph
           dcalEa(nn,mm,j) = ( conjg(Q_Qmat(mm,nn,j)) +P_Qmat(nn,mm,j) )*DeltaT  
        end do
     end do
  end do

  !----------------------------------------------
  !
  ! Atomic nucleus part: O_Nmat, P_Nmat, Q_Nmat
  !
  !----------------------------------------------
  
  !-----------------------------------------
  ! calculation of O_aij 
  !-----------------------------------------
  !----- T_ajm calC_aim -----
  do ii=1,NBS_N
     do jj=1,NBS_N
        sum = (0._dp,0._dp)
        do mm = 1,NBS_N
           sum = sum + Tnuc_mat(jj,mm)*calC(ii,mm)
        end do
        O_Nmat(ii,jj) = sum/IU
     end do
  end do

  !----- (jm|PQ) calC_im calE_PQ -----
  I_Nmat(:,:) = (0._dp,0._dp)
  do jj=1,NBS_N
     do mm=1,NBS_N
        sum = (0._dp,0._dp)
        do pp=1,4*NBS
           do qq=1,4*NBS
              sum = sum + nucele_Qmat(jj,mm,pp,qq)*calE(pp,qq)
           end do
        end do
        I_Nmat(jj,mm) = sum
     end do
  end do

  do ii=1,NBS_N
     do jj=1,NBS_N
        sum = (0._dp,0._dp)
        do mm=1,NBS_N
           sum = sum +I_Nmat(jj,mm)*calC(ii,mm)
        end do
        O_Nmat(ii,jj) = O_Nmat(ii,jj) +sum/IU
     end do
  end do

  !----- (jm|kl) calC_im calC_kl ------
  I_Nmat(:,:) = (0._dp,0._dp)
  do jj=1,NBS_N
     do mm=1,NBS_N
        sum = (0._dp,0._dp)
        do kk=1,NBS_N
           do ll=1,NBS_N
              sum = sum + twonuc_mat(jj,mm,kk,ll)*calC(kk,ll)
           end do
        end do
        I_Nmat(jj,mm) = sum
     end do
  end do

  do ii=1,NBS_N
     do jj=1,NBS_N
        sum = (0._dp,0._dp)
        do mm=1,NBS_N
           sum = sum + I_Nmat(jj,mm)*calC(ii,mm)
        end do
        O_Nmat(ii,jj) = O_Nmat(ii,jj) +sum/IU
     end do
  end do
  
  if(use_exchange) then
     !----- (jm|kl) calC_il calC_km ------
     I_Nmat(:,:) = (0._dp,0._dp)
     do jj=1,NBS_N
        do ll=1,NBS_N
           sum = (0._dp,0._dp)
           do kk=1,NBS_N
              do mm=1,NBS_N
                 sum = sum + twonuc_mat(jj,mm,kk,ll)*calC(kk,mm)
              end do
           end do
           I_Nmat(jj,ll) = sum
        end do
     end do
     
     do ii=1,NBS_N
        do jj=1,NBS_N
           sum = (0._dp,0._dp)
           do ll=1,NBS_N
              sum = sum + I_Nmat(jj,ll)*calC(ii,ll)
           end do
           if(NUCTYPE.eq.0) then !boson
              O_Nmat(ii,jj) = O_Nmat(ii,jj) +sum/IU
           elseif(NUCTYPE.eq.1) then !fermion
              O_Nmat(ii,jj) = O_Nmat(ii,jj) -sum/IU
           else
              write(*,*) "NUCTYPE should be 0 or 1."
              stop
           end if
        end do
     end do
  end if

  
  !-----------------------------------------
  ! calculation of P_aiPj 
  !-----------------------------------------
  !----- T_ajm calC_aiPm -----
  do ii=1,NBS_N
     do jj=1,NBS_N
        do j=1,Nph
           sum = (0._dp,0._dp)
           do mm = 1,NBS_N
              sum = sum + Tnuc_mat(jj,mm)*calCa(ii,mm,j)
           end do
        end do
        P_Nmat(ii,jj,j) = sum/IU
     end do
  end do
  

  !-----------------------------------------
  ! calculation of Q_aiPj 
  !-----------------------------------------
  !----- T_ajm (calC_amPi)^* -----
  do ii=1,NBS_N
     do jj=1,NBS_N
        do j=1,Nph
           sum = (0._dp,0._dp)
           do mm = 1,NBS_N
              sum = sum + Tnuc_mat(jj,mm)*conjg(calCa(mm,ii,j))
           end do
        end do
        Q_Nmat(ii,jj,j) = sum/IU
     end do
  end do
  
  !---------------------------------------
  ! calculation of derivatives(*DeltaT)
  !---------------------------------------
  !----- d(calC_ij)/dt = (O_ji)^* + O_ij -----
  do ii=1,NBS_N
     do jj=1,NBS_N
        dcalC(ii,jj) = (conjg(O_Nmat(jj,ii)) + O_Nmat(ii,jj))*DeltaT
     end do
  end do

  !----- d(calCa_iPj)/dt = (Q_jPi)^* + P_iPj -----
  do ii=1,NBS_N
     do jj=1,NBS_N
        do j=1,Nph
           dcalCa(ii,jj,j) = (conjg(Q_Nmat(jj,ii,j)) + P_Nmat(ii,jj,j))*DeltaT
        end do
     end do
  end do
  
  return
end subroutine calc_derivatives_NonBO_OPQ



!========================================================================================================
subroutine calc_derivatives_NonBOwithEx(use_exchange,t,calE,calC, & 
     & TM_Qmat,Tnuc_mat,twoele_Qmat,nucele_Qmat,twonuc_mat,calFj0_Qmat,alpha, &
     & dcalE,dcalC)
!  derivatives for Non-BO calculation
!  added use_exchange (120712)
!========================================================================================================
  use Precision
  use DiracOutput
  use Constants
  use NucBasis
  implicit none

  logical,intent(in) :: use_exchange  !.true. --> include exchange terms in diff. eq. , .false. --> without exchange terms
  real(kind=dp),intent(in) :: t  ! time
  complex(kind=dp),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  complex(kind=dp),intent(in) :: calC(NBS_N,NBS_N)  ! calC_ij
  complex(kind=dp),intent(in) :: TM_Qmat(4*NBS,4*NBS)  ! T_PQ + M_PQ
  complex(kind=dp),intent(in) :: Tnuc_mat(NBS_N,NBS_N)  ! T_Nij
  complex(kind=dp),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)
  complex(kind=dp),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
  complex(kind=dp),intent(in) :: alpha(Nph)   ! <coherent|hat(a)|coherent>
  complex(kind=dp),intent(in) :: nucele_Qmat(NBS_N,NBS_N,4*NBS,4*NBS)  ! (ij|PQ)
  complex(kind=dp),intent(in) :: twonuc_mat(NBS_N,NBS_N,NBS_N,NBS_N)  ! (ij|kl)

  complex(kind=dp),intent(out) :: dcalE(4*NBS,4*NBS)  ! d(calE_PQ)/dt *DeltaT
  complex(kind=dp),intent(out) :: dcalC(NBS_N,NBS_N)  ! d(calC_ij)/dt *DeltaT

  real(kind=dp) :: dpj ! tilde(Delta p)_j
  real(kind=dp) :: p0j ! p_j

  integer :: i,j,k,l,p,q,r,s,n,m
  complex(kind=dp) :: sum_te,sum_ne,sum_tn
  complex(kind=dp) :: sumE,sumC
  complex(kind=dp) :: dcalF

  complex(kind=dp) :: I2_Qmat(4*NBS,4*NBS), I4_Qmat(4*NBS,4*NBS) 
  complex(kind=dp) :: I_Qmat(4*NBS,4*NBS)  ! I_PQ = (I_1+I_2+I_3+I_4)_PQ = (TM +I2 +I4)_PQ
  complex(kind=dp) :: I4_mat(NBS_N,NBS_N) 
  complex(kind=dp) :: I_mat(NBS_N,NBS_N)  ! I_ij = (I_1+I_2+I_3+I_4)_ij = (T +I4)_ij

  I2_Qmat  = (0._dp,0._dp)
  I4_Qmat  = (0._dp,0._dp)

  ! calculation of I4
  do p=1,4*NBS
     do q=1,4*NBS

        sum_te = (0._dp,0._dp)
        do r=1,4*NBS
           do s=1,4*NBS
              sum_te = sum_te + twoele_Qmat(p,q,r,s)*calE(r,s)
           end do
        end do
        sum_ne = (0._dp,0._dp)
        do i=1,NBS_N
           do j=1,NBS_N
              sum_ne = sum_ne + nucele_Qmat(i,j,p,q)*calC(i,j)  ! (PQ|ij)=(ij|PQ)
           end do
        end do
        I4_Qmat(p,q) = sum_te + sum_ne

     end do
  end do

!!$  ! calculation of I2 (Arad only)
!!$  do p=1,4*NBS
!!$     do q=1,4*NBS
!!$
!!$        sum_calF=(0._dp,0._dp)
!!$        do j=1,Nph
!!$           call calc_DeltaPj_and_P0j(j,dpj,p0j)
!!$           
!!$           dcalF = calFj0_Qmat(nn,mm,j)       *exp(-IU*CCC*p0j*t)*alpha(j) &
!!$                & +conjg(calFj0_Qmat(mm,nn,j))*exp(+IU*CCC*p0j*t)*conjg(alpha(j))
!!$           
!!$           sum_calF = sum_calF -dpj*dcalF
!!$
!!$        end do
!!$        
!!$        I2_Qmat(p,q) = sum_calF
!!$
!!$     end do
!!$  end do

  I_Qmat(:,:) = TM_Qmat(:,:) + I2_Qmat(:,:) + I4_Qmat(:,:)

  ! calculation of derivatives(*DeltaT)
  do p=1,4*NBS
     do q=p,4*NBS
        sumE = (0._dp,0._dp)
        ! non-exchange terms
        do r=1,4*NBS
           sumE  = sumE  -I_Qmat(r,p)*calE(r,q)  +I_Qmat(q,r)*calE(p,r)
        end do

        if(use_exchange) then
           ! exchange terms
           do r=1,4*NBS
              do n=1,4*NBS
                 do m=1,4*NBS
                    sumE = sumE +twoele_Qmat(r,p,n,m)*calE(r,m)*calE(n,q) -twoele_Qmat(q,r,n,m)*calE(p,m)*calE(n,r) 
                 end do
              end do
           end do
        end if
        
        dcalE(p,q)  = -IU*sumE*DeltaT
     end do
  end do
  
  do p=2,4*NBS
     do q=1,p-1
        dcalE(p,q)  = conjg(dcalE(q,p))
     end do
  end do
  
  !--------------------------------------------------------------------------

  I4_mat  = (0._dp,0._dp)

  ! calculation of I4
  do i=1,NBS_N
     do j=1,NBS_N

        sum_ne = (0._dp,0._dp)
        do p=1,4*NBS
           do q=1,4*NBS
              sum_ne = sum_ne + nucele_Qmat(i,j,p,q)*calE(p,q)
           end do
        end do
        sum_tn = (0._dp,0._dp)
        do k=1,NBS_N
           do l=1,NBS_N
              sum_tn = sum_tn + twonuc_mat(i,j,k,l)*calC(k,l)  
           end do
        end do
        I4_mat(i,j) = sum_tn + sum_ne

     end do
  end do

  I_mat(:,:) = Tnuc_mat(:,:) + I4_mat(:,:)

  do i=1,NBS_N
     do j=i,NBS_N

        sumC = (0._dp,0._dp)
        ! non-exchange terms
        do k=1,NBS_N
           sumC  = sumC  -I_mat(k,i)*calC(k,j)  +I_mat(j,k)*calC(i,k)
        end do

        if(use_exchange) then
           ! exchange terms
           if(NUCTYPE.eq.0) then ! boson
              do k=1,NBS_N
                 do n=1,NBS_N
                    do m=1,NBS_N
                       sumC = sumC -twonuc_mat(k,i,n,m)*calC(k,m)*calC(n,j) +twonuc_mat(j,k,n,m)*calC(i,m)*calC(n,k)
                    end do
                 end do
              end do
           elseif(NUCTYPE.eq.1) then ! fermion
              do k=1,NBS_N
                 do n=1,NBS_N
                    do m=1,NBS_N
                       sumC = sumC +twonuc_mat(k,i,n,m)*calC(k,m)*calC(n,j) -twonuc_mat(j,k,n,m)*calC(i,m)*calC(n,k)
                    end do
                 end do
              end do
           else
              write(*,*) "NUCTYPE should be 0 or 1."
              stop
           end if
        end if

        dcalC(i,j)  = -IU*sumC*DeltaT
     end do
  end do
  
  do i=2,NBS_N
     do j=1,i-1
        dcalC(i,j)  = conjg(dcalC(j,i))
     end do
  end do

  return
end subroutine calc_derivatives_NonBOwithEx

!==============================================================================================
subroutine calc_derivatives_BOwithEx(Use_exchange,t,calE,calEV,calDV, &
     & h_Qmat,twoele_Qmat,calFj0_Qmat,alpha,dcalE,dcalEV,dcalDV)
!  derivatives for BO approximation
!  added use_exchange (120712)
!==============================================================================================
  use Precision
  use DiracOutput
  use Constants
  implicit none

  logical,intent(in) :: use_exchange  !.true. --> include exchange terms in diff. eq. , .false. --> without exchange terms
  real(kind=dp),intent(in) :: t  ! time
  complex(kind=dp),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  complex(kind=dp),intent(in) :: calEV(4*NBS,4*NBS)  ! calEV_PQ
  complex(kind=dp),intent(in) :: calDV(4*NBS,4*NBS)  ! calDV_PQ
  complex(kind=dp),intent(in) :: h_Qmat(4*NBS,4*NBS)  ! h_PQ
  complex(kind=dp),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)
  complex(kind=dp),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
  complex(kind=dp),intent(in) :: alpha(Nph)   ! <coherent|hat(a)|coherent>

  complex(kind=dp),intent(out) :: dcalE(4*NBS,4*NBS)  ! d(calE_PQ)/dt *DeltaT
  complex(kind=dp),intent(out) :: dcalEV(4*NBS,4*NBS)  ! d(calEV_PQ)/dt *DeltaT
  complex(kind=dp),intent(out) :: dcalDV(4*NBS,4*NBS)  ! d(calDV_PQ)/dt *DeltaT

  real(kind=dp) :: dpj ! tilde(Delta p)_j
  real(kind=dp) :: p0j ! p_j

  integer :: i,j,k,l
  integer :: pp,qq,rr,ss,nn,mm
  complex(kind=dp) :: sum_te,sum_calF,sum,sumEV,sumDV
  complex(kind=dp) :: dcalF

  complex(kind=dp) :: I2_Qmat(4*NBS,4*NBS), I4_Qmat(4*NBS,4*NBS) 
  complex(kind=dp) :: I_Qmat(4*NBS,4*NBS)  ! I_PQ = (I_1+I_2+I_3+I_4)_PQ = (TM +I2 +I4)_PQ

  complex(kind=dp) :: I4V_Qmat(4*NBS,4*NBS),IV_Qmat(4*NBS,4*NBS)    ! for VEVs
  complex(kind=dp) :: sum_teV

  I2_Qmat  = (0._dp,0._dp)
  I4_Qmat  = (0._dp,0._dp)
  I4V_Qmat  = (0._dp,0._dp)

  sum_calF=(0._dp,0._dp); sum_te=(0._dp,0._dp)

  ! calculation of I4 (h_PQ + (PQ|RS)*calE_RS ) and I4V (h_PQ + (PQ|RS)*calEV_RS )
  !$omp parallel do default(none), shared(I4_Qmat,I4V_Qmat,twoele_Qmat,h_Qmat,calE,calEV,NBS), private(nn,mm,rr,ss,sum_te,sum_teV)
  do nn=1,4*NBS
     do mm=1,4*NBS

        sum_te = (0._dp,0._dp)
        sum_teV = (0._dp,0._dp)
        do rr=1,4*NBS
           do ss=1,4*NBS
              sum_te = sum_te + twoele_Qmat(nn,mm,rr,ss)*calE(rr,ss)
              sum_teV = sum_teV + twoele_Qmat(nn,mm,rr,ss)*calEV(rr,ss)
           end do
        end do

        I4_Qmat(nn,mm) = sum_te + h_Qmat(nn,mm)
        I4V_Qmat(nn,mm) = sum_teV + h_Qmat(nn,mm)

     end do
  end do
  !$omp end parallel do

  ! calculation of I2 (Arad only)
  ! 120828 modified to use continuous-mode coherent state
  !$omp parallel do default(none), shared(I2_Qmat,calFj0_Qmat,calE,alpha,t,NBS,Nph,PI),&
                            !$omp& private(nn,mm,j,dpj,p0j,sum_calF,dcalF)
  do nn=1,4*NBS
     do mm=1,4*NBS

        sum_calF=(0._dp,0._dp)
        do j=1,Nph
           call calc_DeltaPj_and_P0j(j,dpj,p0j) ! sub_int.f90
           ! we do not use dpj in our modified version 120828
!           write(*,*) j,dpj,p0j,alpha(j)
           
           dcalF = calFj0_Qmat(nn,mm,j)       *exp(-IU*CCC*p0j*t)*alpha(j) &
                & +conjg(calFj0_Qmat(mm,nn,j))*exp(+IU*CCC*p0j*t)*conjg(alpha(j))
           
!           sum_calF = sum_calF -dpj*dcalF
           sum_calF = sum_calF -dcalF/(2._dp*PI*sqrt(p0j*CCC))

        end do
!        stop

        I2_Qmat(nn,mm) = sum_calF

     end do
  end do
  !$omp end parallel do


  I_Qmat(:,:) = I2_Qmat(:,:) + I4_Qmat(:,:)
  IV_Qmat(:,:) = I2_Qmat(:,:) + I4V_Qmat(:,:)


  ! calculation of derivatives(*DeltaT)
  do pp=1,4*NBS
     do qq=pp,4*NBS
        sum = (0._dp,0._dp);sumEV = (0._dp,0._dp);sumDV = (0._dp,0._dp)
        
        ! non-exchange terms
        do rr=1,4*NBS
           sum   = sum   -I_Qmat(rr,pp) *calE(rr,qq)  +I_Qmat(qq,rr) *calE(pp,rr)
           sumEV = sumEV -IV_Qmat(rr,pp)*calEV(rr,qq) +IV_Qmat(qq,rr)*calEV(pp,rr)
           sumDV = sumDV +IV_Qmat(pp,rr)*calDV(rr,qq) -IV_Qmat(rr,qq)*calDV(pp,rr)
        end do
        
        if(use_exchange) then
           ! exchange terms
           do rr=1,4*NBS
              do nn=1,4*NBS
                 do mm=1,4*NBS
                    sum   = sum   +twoele_Qmat(rr,pp,nn,mm)*calE(rr,mm) *calE(nn,qq)  -twoele_Qmat(qq,rr,nn,mm)*calE(pp,mm) *calE(nn,rr) 
                    sumEV = sumEV +twoele_Qmat(rr,pp,nn,mm)*calEV(rr,mm)*calEV(nn,qq) -twoele_Qmat(qq,rr,nn,mm)*calEV(pp,mm)*calEV(nn,rr) 
                    sumDV = sumDV -twoele_Qmat(pp,rr,nn,mm)*calEV(nn,rr)*calDV(mm,qq) +twoele_Qmat(rr,qq,nn,mm)*calDV(pp,nn)*calEV(rr,mm) 
                 end do
              end do
           end do
        end if

        dcalE(pp,qq)  = -IU*sum*DeltaT
        dcalEV(pp,qq) = -IU*sumEV*DeltaT
        dcalDV(pp,qq) = -IU*sumDV*DeltaT
     end do
  end do
     
  do pp=2,4*NBS
     do qq=1,pp-1
        dcalE(pp,qq)  = conjg(dcalE(qq,pp))
        dcalEV(pp,qq) = conjg(dcalEV(qq,pp))
        dcalDV(pp,qq) = conjg(dcalDV(qq,pp))
     end do
  end do
  
  return
end subroutine calc_derivatives_BOwithEx


!==============================================================================================
subroutine calc_L_Qmat(t,calE,calFj0_Qmat,L_Qmat)
!  To check logarithmic divergence
!
!  L_NM = IU*dp0_th_phi*p0_j*sinth_j* calF_MRj* calF^*_SRj* calE_NS
!
!  130124
!==============================================================================================
  use Precision
  use DiracOutput
  use Constants
  implicit none

  real(kind=dp),intent(in) :: t  ! time
  complex(kind=dp),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  complex(kind=dp),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
  
  complex(kind=dp),intent(out) :: L_Qmat(4*NBS,4*NBS)  

  real(kind=dp) :: dp0_th_phi ! dp0*dtheta*dphi
  real(kind=dp) :: p0j,thj,phij ! p_j,theta_j,phi_j

  integer :: i,j,k,l
  integer :: pp,qq,rr,ss,nn,mm
  complex(kind=dp) :: sum,sum1,sum2

!  complex(kind=dp) :: FF_Qmat(4*NBS,4*NBS)  ! sum_j sum_R (calF_MRj* calF^*_SRj)
  complex(kind=dp) :: FcalE1_Qmat(4*NBS,4*NBS,Nph)  ! sum_ss   
  complex(kind=dp) :: FcalE2_Qmat(4*NBS,4*NBS,Nph)  ! sum_ss   


!!$  do j=1,Nph
!!$     call get_pj(j,p0j,thj,phij,dp0_th_phi)
!!$     write(10,"(1i6,3es16.6)",advance="no") j,p0j,thj,phij
!!$!     write(10,"(1i6,5es16.6)") j,p0j,thj,phij,calFj0_Qmat(3,4,j)    
!!$     do mm=1,4*NBS
!!$        do ss=1,4*NBS
!!$           write(10,"(2es16.6)",advance="no") calFj0_Qmat(mm,ss,j)
!!$        end do
!!$     end do
!!$     write(10,*)
!!$  end do
!!$  stop

  do nn=1,4*NBS
     do rr=1,4*NBS
        do j=1,Nph
           
           sum1 =(0._dp,0._dp)
           sum2 =(0._dp,0._dp)
           do ss=1,4*NBS  ! sum over r'^e'
              sum1 = sum1 + conjg(calFj0_Qmat(ss,rr,j))*calE(nn,ss)  ! propto calE[p]_NR
              sum2 = sum2 + calFj0_Qmat(ss,nn,j)*calE(ss,rr)
           end do
           FcalE1_Qmat(nn,rr,j) = sum1
           FcalE2_Qmat(nn,rr,j) = sum2

!           call get_pj(j,p0j,thj,phij,dp0_th_phi)
!           write(11,"(3i6,2es16.6)") nn,rr,j,sum1/(2._dp*PI*sqrt(p0j))
        end do
     end do
  end do
  
  do nn=1,4*NBS
     do mm=1,4*NBS
        
        sum = (0._dp,0._dp)
        do j=1,Nph
           call get_pj(j,p0j,thj,phij,dp0_th_phi)
           do rr=1,4*NBS
              sum = sum + 0.5_dp*dp0_th_phi*p0j*sin(thj)* &
                   & ( calFj0_Qmat(mm,rr,j)*exp(-IU*CCC*p0j*t)*FcalE1_Qmat(nn,rr,j) &
                   & -conjg(calFj0_Qmat(rr,mm,j))*exp(IU*CCC*p0j*t)*FcalE2_Qmat(nn,rr,j) )
           end do
        end do
        L_Qmat(nn,mm) = sum

     end do
  end do

!!$  do nn=1,4*NBS
!!$     do rr=1,4*NBS
!!$        do j=1,Nph
!!$           
!!$           sum1 =(0._dp,0._dp)
!!$           sum2 =(0._dp,0._dp)
!!$           do ss=1,4*NBS  ! sum over r'^e'
!!$              sum1 = sum1 + conjg(calFj0_Qmat(ss,rr,j))*calE(nn,ss)
!!$              sum2 = sum2 + conjg(calFj0_Qmat(nn,ss,j))*calE(ss,rr)
!!$           end do
!!$           FcalE1_Qmat(nn,rr,j) = sum1
!!$           FcalE2_Qmat(nn,rr,j) = sum2
!!$           
!!$        end do
!!$     end do
!!$  end do
!!$  
!!$  do nn=1,4*NBS
!!$     do mm=1,4*NBS
!!$        
!!$        sum = (0._dp,0._dp)
!!$        do j=1,Nph
!!$           do rr=1,4*NBS
!!$              call get_pj(j,p0j,thj,phij,dp0_th_phi)
!!$              sum = sum + 0.5_dp*dp0_th_phi*p0j*sin(thj)*calFj0_Qmat(mm,rr,j)* &
!!$                   & (exp(-IU*CCC*p0j*t)*FcalE1_Qmat(nn,rr,j) -exp(IU*CCC*p0j*t)*FcalE2_Qmat(nn,rr,j))
!!$           end do
!!$        end do
!!$        L_Qmat(nn,mm) = sum
!!$
!!$     end do
!!$  end do

! we used incorrect expressions below
!!$  do mm=1,4*NBS
!!$     do ss=1,4*NBS
!!$        
!!$        sum = (0._dp,0._dp)
!!$        do rr=1,4*NBS
!!$           do j=1,Nph
!!$              call get_pj(j,p0j,thj,phij,dp0_th_phi)
!!$              sum = sum + IU*dp0_th_phi*p0j*sin(thj)*sin(CCC*p0j*t)*calFj0_Qmat(mm,rr,j)*conjg(calFj0_Qmat(ss,rr,j))
!!$!              write(10,"(1i6,4es16.6)") j,p0j,thj,phij,sin(CCC*p0j*t)
!!$!              write(10,"(2i6,6es16.6)") rr,j,calFj0_Qmat(mm,rr,j),calFj0_Qmat(ss,rr,j),calFj0_Qmat(mm,rr,j)*conjg(calFj0_Qmat(ss,rr,j))
!!$           end do
!!$!           write(*,*) sum
!!$        end do
!!$        FF_Qmat(mm,ss) = sum
!!$!        write(*,"(2i6,3es16.6)") mm,ss,sum,dp0_th_phi
!!$!        stop
!!$
!!$     end do
!!$  end do
!!$!  stop
!!$
!!$  do nn=1,4*NBS
!!$     do mm=1,4*NBS
!!$
!!$        sum =(0._dp,0._dp)
!!$        do ss=1,4*NBS
!!$           sum = sum + FF_Qmat(mm,ss)*calE(nn,ss)
!!$        end do
!!$        L_Qmat(nn,mm) = sum
!!$        
!!$     end do
!!$  end do

  return
end subroutine calc_L_Qmat



!==============================================================================================
subroutine calc_derivatives_BO_N(t,calE,calEa,h_Qmat,calFj0_Qmat,dcalE,dcalEa)
!  Derivatives for BO approximation including those for <e^dagger a e>
!  neglecting the term which is nonlinear in calE (which is indicated by "N")
!
!  130513 
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

  complex(kind=dp),intent(out) :: dcalE(4*NBS,4*NBS)  ! d(calE_PQ)/dt *DeltaT
  complex(kind=dp),intent(out) :: dcalEa(4*NBS,4*NBS,Nph) ! d(calEa_PQj)/dt *DeltaT

  real(kind=dp) :: dp0_th_phi ! dp0*dtheta*dphi
  real(kind=dp) :: p0j,thj,phij ! p_j,theta_j,phi_j

  integer :: i,j,k,l
  integer :: pp,qq,rr,ss,nn,mm
  complex(kind=dp) :: sum

  complex(kind=dp) :: ON_Qmat(4*NBS,4*NBS), PN_Qmat(4*NBS,4*NBS,Nph), QN_Qmat(4*NBS,4*NBS,Nph)

  sum = (0._dp,0._dp)
  ON_Qmat(:,:) = (0._dp,0._dp)
  PN_Qmat(:,:,:) = (0._dp,0._dp)
  QN_Qmat(:,:,:) = (0._dp,0._dp)

  !---------------------------------------
  !   calculation of O^N_NM
  !---------------------------------------
  ! h_MR calE_NR
  do nn=1,4*NBS
     do mm=1,4*NBS
        sum = (0._dp,0._dp)
        do rr=1,4*NBS
           sum = sum +h_Qmat(mm,rr)*calE(nn,rr)
        end do
        ON_Qmat(nn,mm) = sum/IU
     end do
  end do

  ! -1/sqrt(2 pi^2 c)* int d^3p/sqrt(2*p0) 
  ! ( F^k_MR(p)e^k(p) exp(-i c p0 t) calEa_NRj + F^*k_RM(p)e^*k(p) exp(i c p0 t) calEa^*_RNj
!  write(*,*) "F term in deriv"
!  write(10,*) "F term in deriv"
  do nn=1,4*NBS
     do mm=1,4*NBS

        sum = (0._dp,0._dp)
        do j=1,Nph
           call get_pj(j,p0j,thj,phij,dp0_th_phi)
           do rr=1,4*NBS
              sum = sum -dp0_th_phi*p0j**2*sin(thj)/(2._dp*PI*sqrt(CCC)*sqrt(p0j)) *( calFj0_Qmat(mm,rr,j)*exp(-IU*CCC*p0j*t)*calEa(nn,rr,j) &
                   & + conjg(calFj0_Qmat(rr,mm,j))*exp(IU*CCC*p0j*t)*conjg(calEa(rr,nn,j)) )
           end do
        end do

!!$        write(*,"(2i6,4es16.6)") nn,mm,sum/IU
!!$        write(10,"(2i6,4es16.6)") nn,mm,sum/IU

        ON_Qmat(nn,mm) = ON_Qmat(nn,mm) + sum/IU
!        ON_Qmat(nn,mm) = sum/IU

     end do
  end do

  !---------------------------------------
  !   calculation of P^N_NM (n=0)
  !---------------------------------------
  ! h_MR calEa_NRj
  do nn=1,4*NBS
     do mm=1,4*NBS
        do j=1,Nph
           sum = (0._dp,0._dp)
           do rr=1,4*NBS
              sum = sum +h_Qmat(mm,rr)*calEa(nn,rr,j)
           end do
           PN_Qmat(nn,mm,j) = sum/IU
        end do
     end do
  end do

!  write(*,*) "P in deriv"
!  write(10,*) "P in deriv"
  do nn=1,4*NBS
     do mm=1,4*NBS
        do j=1,Nph

           call get_pj(j,p0j,thj,phij,dp0_th_phi)
           sum = (0._dp,0._dp)
           do rr=1,4*NBS
              sum = sum -1._dp/(2._dp*PI*sqrt(CCC)*sqrt(p0j))*conjg(calFj0_Qmat(rr,mm,j))*exp(IU*CCC*p0j*t)*calE(nn,rr)
           end do
!           write(9,"(3i6,2es16.6)") nn,mm,j,sum
           PN_Qmat(nn,mm,j) = PN_Qmat(nn,mm,j) + sum/IU
!           write(*,"(2i6,4es16.6)") nn,mm,sum/IU
!           write(10,"(2i6,4es16.6)") nn,mm,sum/IU
!           PN_Qmat(nn,mm,j) = sum/IU

        end do
     end do
  end do

  
  !---------------------------------------
  !   calculation of Q^N_NM (n=0)
  !---------------------------------------
  ! h_MR calEa^*_RNj
  do nn=1,4*NBS
     do mm=1,4*NBS
        do j=1,Nph
           sum = (0._dp,0._dp)
           do rr=1,4*NBS
              sum = sum +h_Qmat(mm,rr)*conjg(calEa(rr,nn,j))
           end do
           QN_Qmat(nn,mm,j) = sum/IU
        end do
     end do
  end do

  !---------------------------------------
  ! calculation of derivatives(*DeltaT)
  !---------------------------------------
  ! d(calE_NM)/dt = O_NM + O_NM^dagger
!  write(*,*) "dcalE in deriv"
!  write(10,*) "dcalE in deriv"
  do nn=1,4*NBS
     do mm=1,4*NBS
        dcalE(nn,mm) = (ON_Qmat(nn,mm) +conjg(ON_Qmat(mm,nn)))*DeltaT
!        write(*,"(2i6,4es16.6)") nn,mm,dcalE(nn,mm)
!        write(10,"(2i6,4es16.6)") nn,mm,dcalE(nn,mm)
     end do
  end do

  ! d(calEa_NMj)/dt = (Q_NMj)^dagger + P_NMj
  do nn=1,4*NBS
     do mm=1,4*NBS
        do j=1,Nph
           dcalEa(nn,mm,j) = ( conjg(QN_Qmat(mm,nn,j)) +PN_Qmat(nn,mm,j) )*DeltaT  
        end do
     end do
  end do
  
  return
end subroutine calc_derivatives_BO_N


!==============================================================================================
!subroutine calc_derivatives_BO_Phi_test(t,calE,h_Qmat,calFj0_Qmat,twoele_Qmat_nz,dcalE)
subroutine calc_derivatives_BO_Phi_test(t,calE,h_Qmat,calFj0_Qmat,twoele_Qmat,dcalE)
!  Derivatives for BO approximation 
!  to test new implimentation of two-electron integral.
!  130314
!==============================================================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use IntegralStorage
  implicit none

  real(kind=dp),intent(in) :: t  ! time
  complex(kind=dp),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  complex(kind=dp),intent(in) :: h_Qmat(4*NBS,4*NBS)  ! h_PQ
  complex(kind=dp),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
!  type(nonzeromat4legs),intent(in) :: twoele_Qmat_nz(N_twoele)  
  complex(kind=dp),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)

  complex(kind=dp),intent(out) :: dcalE(4*NBS,4*NBS)  ! d(calE_PQ)/dt *DeltaT

  real(kind=dp) :: dp0_th_phi ! dp0*dtheta*dphi
  real(kind=dp) :: p0j,thj,phij ! p_j,theta_j,phi_j

  integer :: i,j,k,l
  integer :: pp,qq,rr,ss,nn,mm
  complex(kind=dp) :: sum

  complex(kind=dp) :: ON_Qmat(4*NBS,4*NBS)
!  complex(kind=dp) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS),I4_Qmat(4*NBS,4*NBS)
  complex(kind=dp) :: I4_Qmat(4*NBS,4*NBS)

  sum = (0._dp,0._dp)
  ON_Qmat(:,:) = (0._dp,0._dp)
  I4_Qmat(:,:) = (0._dp,0._dp)
!  twoele_Qmat(:,:,:,:) = (0._dp,0._dp)

  !---------------------------------------
  !   calculation of O^N_NM
  !---------------------------------------
  ! h_MR calE_NR
  do nn=1,4*NBS
     do mm=1,4*NBS
        sum = (0._dp,0._dp)
        do rr=1,4*NBS
           sum = sum +h_Qmat(mm,rr)*calE(nn,rr)
        end do
        ON_Qmat(nn,mm) = sum/IU
     end do
  end do

  ! -1/sqrt(2 pi^2 c)* int d^3p/sqrt(2*p0) 
  ! ( F^k_MR(p)e^k(p) exp(-i c p0 t) calEa_NRj + F^*k_RM(p)e^*k(p) exp(i c p0 t) calEa^*_RNj
!  write(*,*) "F term in deriv"
!  write(10,*) "F term in deriv"
!!$  do nn=1,4*NBS
!!$     do mm=1,4*NBS
!!$
!!$        sum = (0._dp,0._dp)
!!$        do j=1,Nph
!!$           call get_pj(j,p0j,thj,phij,dp0_th_phi)
!!$           do rr=1,4*NBS
!!$              sum = sum -dp0_th_phi*p0j**2*sin(thj)/(2._dp*PI*sqrt(CCC)*sqrt(p0j)) *( calFj0_Qmat(mm,rr,j)*exp(-IU*CCC*p0j*t)*calEa(nn,rr,j) &
!!$                   & + conjg(calFj0_Qmat(rr,mm,j))*exp(IU*CCC*p0j*t)*conjg(calEa(rr,nn,j)) )
!!$           end do
!!$        end do
!!$
!!$        ON_Qmat(nn,mm) = ON_Qmat(nn,mm) + sum/IU
!!$     end do
!!$  end do


!!$  do i=1,N_twoele
!!$     nn = twoele_Qmat_nz(i)%a
!!$     mm = twoele_Qmat_nz(i)%b
!!$     pp = twoele_Qmat_nz(i)%c
!!$     qq = twoele_Qmat_nz(i)%d
!!$     twoele_Qmat(nn,mm,pp,qq) = twoele_Qmat_nz(i)%val
!!$  end do

  do mm=1,4*NBS
     do rr=1,4*NBS

        sum = (0._dp,0._dp)
        do pp=1,4*NBS
           do qq=1,4*NBS
              sum = sum + twoele_Qmat(mm,rr,pp,qq)*calE(pp,qq)
           end do
        end do
        I4_Qmat(mm,rr) = sum

     end do
  end do

  do nn=1,4*NBS
     do mm=1,4*NBS
        sum = (0._dp,0._dp)
        do rr=1,4*NBS
           sum = sum + I4_Qmat(mm,rr)*calE(nn,rr)
        end do

        ON_Qmat(nn,mm) = ON_Qmat(nn,mm) + sum/IU
     end do
  end do



!!$  do nn=1,4*NBS
!!$     do mm=1,4*NBS
!!$        sum = (0._dp,0._dp)
!!$        do i=1,N_twoele
!!$           rr = twoele_Qmat_nz(i)%b
!!$           pp = twoele_Qmat_nz(i)%c
!!$           qq = twoele_Qmat_nz(i)%d
!!$           if(twoele_Qmat_nz(i)%a.eq.mm) then
!!$              sum = sum +twoele_Qmat_nz(i)%val*calE(nn,rr)*calE(pp,qq)
!!$           end if
!!$        end do
!!$
!!$        ON_Qmat(nn,mm) = ON_Qmat(nn,mm) + sum/IU
!!$     end do
!!$  end do

     


  !---------------------------------------
  ! calculation of derivatives(*DeltaT)
  !---------------------------------------
  ! d(calE_NM)/dt = O_NM + O_NM^dagger
!  write(*,*) "dcalE in deriv"
!  write(10,*) "dcalE in deriv"
  do nn=1,4*NBS
     do mm=1,4*NBS
        dcalE(nn,mm) = (ON_Qmat(nn,mm) +conjg(ON_Qmat(mm,nn)))*DeltaT
!        write(*,"(2i6,4es16.6)") nn,mm,dcalE(nn,mm)
!        write(10,"(2i6,4es16.6)") nn,mm,dcalE(nn,mm)
     end do
  end do

  return
end subroutine calc_derivatives_BO_Phi_test


!==============================================================================================
!subroutine calc_derivatives_BO_Phi(use_retarded,use_exchange,t,calE,calE_past,h_Qmat,&
!     & calFj0_Qmat,alpha,twoele_Qmat,intIJJ_Qmat,dcalE)
subroutine calc_derivatives_BO_O(use_retarded,use_exchange,t,calE,calE_past,h_Qmat,&
     & calFj0_Qmat,alpha,twoele_Qmat,intIJJ_Qmat,dcalE)  ! 131227
!  Derivatives for BO approximation 
!  to include IJJ
!  131227
!==============================================================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use IntegralStorage
  implicit none

  logical,intent(in) :: use_retarded  ! .true. --> include retarded potential, .false. --> ignore retarded potential
  logical,intent(in) :: use_exchange  !.true. --> include exchange terms in diff. eq. , .false. --> without exchange terms
  real(kind=dp),intent(in) :: t  ! time
  complex(kind=dp),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  complex(kind=dp),intent(in) :: calE_past(4*NBS,4*NBS)  ! calE_PQ at one time step before
  complex(kind=dp),intent(in) :: h_Qmat(4*NBS,4*NBS)  ! h_PQ
  complex(kind=dp),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
  complex(kind=dp),intent(in) :: alpha(Nph)   ! <coherent|hat(a)|coherent>
  complex(kind=dp),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)
  complex(kind=dp),intent(in) :: intIJJ_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! I_JJ

  complex(kind=dp),intent(out) :: dcalE(4*NBS,4*NBS)  ! d(calE_PQ)/dt *DeltaT

  real(kind=dp) :: dp0_th_phi ! dp0*dtheta*dphi
  real(kind=dp) :: p0j,thj,phij ! p_j,theta_j,phi_j
  real(kind=dp) :: dpj ! tilde(Delta p)_j
  complex(kind=dp) :: dcalF

  integer :: i,j,k,l
  integer :: pp,qq,rr,ss,nn,mm
  complex(kind=dp) :: sum

  complex(kind=dp) :: O_Qmat(4*NBS,4*NBS)
  complex(kind=dp) :: I_Qmat(4*NBS,4*NBS)

  O_Qmat(:,:) = (0._dp,0._dp)

  !---------------------------------------
  !   calculation of O^Phi_NM
  !---------------------------------------

  !----------h_MR calE_NR---------------
  do nn=1,4*NBS
     do mm=1,4*NBS
        sum = (0._dp,0._dp)
        do rr=1,4*NBS
           sum = sum +h_Qmat(mm,rr)*calE(nn,rr)
        end do
        O_Qmat(nn,mm) = sum/IU
     end do
  end do
  !---------------------------------------

  !--------(MR|PQ) calE_NR calE_PQ-----------
  I_Qmat(:,:) = (0._dp,0._dp)
  do mm=1,4*NBS
     do rr=1,4*NBS

        sum = (0._dp,0._dp)
        do pp=1,4*NBS
           do qq=1,4*NBS
              sum = sum + twoele_Qmat(mm,rr,pp,qq)*calE(pp,qq)
           end do
        end do
        I_Qmat(mm,rr) = sum

     end do
  end do

  ! add to O^Phi_NM
  do nn=1,4*NBS
     do mm=1,4*NBS
        sum = (0._dp,0._dp)
        do rr=1,4*NBS
           sum = sum + I_Qmat(mm,rr)*calE(nn,rr)
        end do
!        write(*,*) nn,mm,sum
        O_Qmat(nn,mm) = O_Qmat(nn,mm) + sum/IU
     end do
  end do
  !---------------------------------------

  if(use_exchange) then
     !-----  -(MR|PQ) calE_NQ calE_PR -----
     I_Qmat(:,:) = (0._dp,0._dp)
     do mm=1,4*NBS
        do qq=1,4*NBS
           sum = (0._dp,0._dp)
           do pp=1,4*NBS
              do rr=1,4*NBS
                 sum = sum +twoele_Qmat(mm,rr,pp,qq)*calE(pp,rr)
              end do
           end do
           I_Qmat(mm,qq) = sum
        end do
     end do
     
     ! add to O^Phi_NM
     do nn=1,4*NBS
        do mm=1,4*NBS
           sum = (0._dp,0._dp)
           do qq=1,4*NBS
              sum = sum + I_Qmat(mm,qq)*calE(nn,qq)
           end do
           O_Qmat(nn,mm) = O_Qmat(nn,mm) - sum/IU  ! note the minus sign here.
        end do
     end do
     
  end if



  !--- Arad term using coherent state ---
  I_Qmat(:,:) = (0._dp,0._dp)
  ! compute sum over j (coherent mode)
  do mm=1,4*NBS
     do rr=1,4*NBS
        sum = (0._dp,0._dp)
        do j=1,Nph
           call calc_DeltaPj_and_P0j(j,dpj,p0j) ! sub_int.f90 (we use only p0j below)
           dcalF = calFj0_Qmat(mm,rr,j)       *exp(-IU*CCC*p0j*t)*alpha(j) &
                & +conjg(calFj0_Qmat(rr,mm,j))*exp(+IU*CCC*p0j*t)*conjg(alpha(j))
           sum = sum -dcalF/(2._dp*PI*sqrt(p0j*CCC))
        end do
        I_Qmat(mm,rr) = sum
     end do
  end do

  ! add to O^Phi_NM
  do nn=1,4*NBS
     do mm=1,4*NBS
        sum = (0._dp,0._dp)
        do rr=1,4*NBS
           sum = sum + I_Qmat(mm,rr)*calE(nn,rr)
        end do
        O_Qmat(nn,mm) = O_Qmat(nn,mm) + sum/IU
     end do
  end do
  

  if(use_retarded) then
     !---------retarded potential term : IJJ_MR calE_NR-------------------
     ! construct I_jT(t)_MR
     I_Qmat(:,:) = (0._dp,0._dp)
     do mm=1,4*NBS
        do rr=1,4*NBS
           
           sum = (0._dp,0._dp)
           do pp=1,4*NBS
              do qq=1,4*NBS
                 !if(intIJJ_Qmat(mm,rr,pp,qq).ne.(0._dp,0._dp)) then
                 !   write(10,"(4i6,4es16.6)") mm,rr,pp,qq,intIJJ_Qmat(mm,rr,pp,qq),calE_past(pp,qq)
                 !end if
                 !if(calE_past(pp,qq).ne.(0._dp,0._dp)) then
                 !   write(10,"(4i6,4es16.6)") mm,rr,pp,qq,intIJJ_Qmat(mm,rr,pp,qq),calE_past(pp,qq)
                 !end if
                 sum = sum + intIJJ_Qmat(mm,rr,pp,qq)*calE_past(pp,qq)
              end do
           end do
           I_Qmat(mm,rr) = sum/(-CCC**3*PI) *t    ! approximated int du
           
!!$           if( abs(I_Qmat(mm,rr)).gt.1.e-5_dp) then
!!$              !           write(*,"(2i6,2es16.6)") mm,rr,I_Qmat(mm,rr)
!!$              !           write(9,"(2i6,2es16.6)") mm,rr,I_Qmat(mm,rr)
!!$           end if
           
        end do
     end do
     
     do nn=1,4*NBS
        do mm=1,4*NBS
           sum = (0._dp,0._dp)
           do rr=1,4*NBS
              sum = sum + I_Qmat(mm,rr)*calE(nn,rr)
           end do
           O_Qmat(nn,mm) = O_Qmat(nn,mm) + sum/IU
        end do
     end do
     !---------------------------------------
  end if

  !---------------------------------------
  ! calculation of derivatives(*DeltaT)
  !---------------------------------------
  ! d(calE_NM)/dt = O_NM + O_NM^dagger
!  write(*,*) "dcalE in deriv"
!  write(10,*) "dcalE in deriv"
  do nn=1,4*NBS
     do mm=1,4*NBS
        dcalE(nn,mm) = (O_Qmat(nn,mm) +conjg(O_Qmat(mm,nn)))*DeltaT
!        write(*,"(2i6,4es16.6)") nn,mm,dcalE(nn,mm)
!        write(10,"(2i6,4es16.6)") nn,mm,dcalE(nn,mm)
     end do
  end do

  return
end subroutine calc_derivatives_BO_O



!==============================================================================================
!subroutine calc_derivatives_BO_OPQ(use_exchange,t,calE,calEa,h_Qmat,calFj0_Qmat,twoele_Qmat,dcalE,dcalEa)
!subroutine calc_derivatives_BO_OPQ(use_exchange,t,calE,calEa,h_Qmat,VT_Qmat,M_Qmat,calFj0_Qmat,twoele_Qmat,dcalE,dcalEa,Zme) !131123
!subroutine calc_derivatives_BO_OPQ(use_exchange,t,Zme_in,calE,calEa,VT_Qmat,M_Qmat,calFj0_Qmat,twoele_Qmat,dcalE,dcalEa,Zme) !131123
subroutine calc_derivatives_BO_OPQ(use_exchange,t,Zme_in,calE,calEa,VT_Qmat,M_Qmat,calFj0_Qmat,twoele_Qmat,dcalE,dcalEa) !131123
!  Derivatives for BO approximation including those for <e^dagger a e>
!
!  131125 
!==============================================================================================
  use Precision
  use DiracOutput
  use Constants
  implicit none

  logical,intent(in) :: use_exchange  !.true. --> include exchange terms in diff. eq. , .false. --> without exchange terms
  real(kind=dp),intent(in) :: t  ! time
  complex(kind=dp),intent(in) :: Zme_in ! renormalization factor to be used for time evolution
  
  complex(kind=dp),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  complex(kind=dp),intent(in) :: calEa(4*NBS,4*NBS,Nph) ! calEa_PQj
!  complex(kind=dp),intent(in) :: h_Qmat(4*NBS,4*NBS)  ! h_PQ
  complex(kind=dp),intent(in) :: VT_Qmat(4*NBS,4*NBS)  ! VT_PQ
  complex(kind=dp),intent(in) :: M_Qmat(4*NBS,4*NBS)  ! M_PQ
  complex(kind=dp),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
  complex(kind=dp),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)

  complex(kind=dp),intent(out) :: dcalE(4*NBS,4*NBS)  ! d(calE_PQ)/dt *DeltaT
  complex(kind=dp),intent(out) :: dcalEa(4*NBS,4*NBS,Nph) ! d(calEa_PQj)/dt *DeltaT
!  complex(kind=dp),intent(out) :: Zme ! new renormalization factor calculated after time evolution

  complex(kind=dp) :: h_Qmat(4*NBS,4*NBS)  ! h_PQ

  real(kind=dp) :: dp0_th_phi ! dp0*dtheta*dphi
  real(kind=dp) :: p0j,thj,phij ! p_j,theta_j,phi_j

  integer :: i,j,k,l
  integer :: pp,qq,rr,ss,nn,mm
  complex(kind=dp) :: sum

  complex(kind=dp) :: O_Qmat(4*NBS,4*NBS), P_Qmat(4*NBS,4*NBS,Nph), Q_Qmat(4*NBS,4*NBS,Nph)
  complex(kind=dp) :: I_Qmat(4*NBS,4*NBS)
  
!  complex(kind=dp) :: zme_te,zme_VT,zme_O,zme_M
  

  sum = (0._dp,0._dp)
  O_Qmat(:,:) = (0._dp,0._dp)
  P_Qmat(:,:,:) = (0._dp,0._dp)
  Q_Qmat(:,:,:) = (0._dp,0._dp)

  !---------------------------------------
  !   calculation of O^N_NM
  !---------------------------------------
  !----- h_MR calE_NR -----
  h_Qmat = VT_Qmat + M_Qmat/Zme_in ! 131123 added

  do nn=1,4*NBS
     do mm=1,4*NBS
        sum = (0._dp,0._dp)
        do rr=1,4*NBS
           sum = sum +h_Qmat(mm,rr)*calE(nn,rr)
        end do
        O_Qmat(nn,mm) = sum/IU
     end do
  end do

  !----- (MR|PQ) calE_NR calE_PQ -----
  I_Qmat(:,:) = (0._dp,0._dp)
  do mm=1,4*NBS
     do rr=1,4*NBS
        sum = (0._dp,0._dp)
        do pp=1,4*NBS
           do qq=1,4*NBS
              sum = sum + twoele_Qmat(mm,rr,pp,qq)*calE(pp,qq)
           end do
        end do
        I_Qmat(mm,rr) = sum
     end do
  end do
  
  ! add to O^Phi_NM
  do nn=1,4*NBS
     do mm=1,4*NBS
        sum = (0._dp,0._dp)
        do rr=1,4*NBS
           sum = sum + I_Qmat(mm,rr)*calE(nn,rr)
        end do
        !        write(*,*) nn,mm,sum
        O_Qmat(nn,mm) = O_Qmat(nn,mm) + sum/IU
     end do
  end do
  
  if(use_exchange) then
     !-----  -(MR|PQ) calE_NQ calE_PR -----
     I_Qmat(:,:) = (0._dp,0._dp)
     do mm=1,4*NBS
        do qq=1,4*NBS
           sum = (0._dp,0._dp)
           do pp=1,4*NBS
              do rr=1,4*NBS
                 sum = sum +twoele_Qmat(mm,rr,pp,qq)*calE(pp,rr)
              end do
           end do
           I_Qmat(mm,qq) = sum
        end do
     end do
     
     ! add to O^Phi_NM
     do nn=1,4*NBS
        do mm=1,4*NBS
           sum = (0._dp,0._dp)
           do qq=1,4*NBS
              sum = sum + I_Qmat(mm,qq)*calE(nn,qq)
           end do
           O_Qmat(nn,mm) = O_Qmat(nn,mm) - sum/IU  ! note the minus sign here.
        end do
     end do
     
  end if



  
  !-----  -1/sqrt(2 pi^2 c)* int d^3p/sqrt(2*p0) 
  !        ( F^k_MR(p)e^k(p) exp(-i c p0 t) calEa_NRj + F^*k_RM(p)e^*k(p) exp(i c p0 t) calEa^*_RNj -----
  do nn=1,4*NBS
     do mm=1,4*NBS

        sum = (0._dp,0._dp)
        do j=1,Nph
           call get_pj(j,p0j,thj,phij,dp0_th_phi)
           do rr=1,4*NBS
              sum = sum -dp0_th_phi*p0j**2*sin(thj)/(2._dp*PI*sqrt(CCC)*sqrt(p0j)) *( calFj0_Qmat(mm,rr,j)*exp(-IU*CCC*p0j*t)*calEa(nn,rr,j) &
                   & + conjg(calFj0_Qmat(rr,mm,j))*exp(IU*CCC*p0j*t)*conjg(calEa(rr,nn,j)) )
           end do
        end do

        O_Qmat(nn,mm) = O_Qmat(nn,mm) + sum/IU

     end do
  end do


  !---------------------------------------
  !   calculation of P^N_NM (n=0)
  !---------------------------------------
  !----- h_MR calEa_NRj -----
  do nn=1,4*NBS
     do mm=1,4*NBS
        do j=1,Nph
           sum = (0._dp,0._dp)
           do rr=1,4*NBS
              sum = sum +h_Qmat(mm,rr)*calEa(nn,rr,j)
           end do
           P_Qmat(nn,mm,j) = sum/IU
        end do
     end do
  end do

  do nn=1,4*NBS
     do mm=1,4*NBS
        do j=1,Nph

           call get_pj(j,p0j,thj,phij,dp0_th_phi)
           sum = (0._dp,0._dp)
           do rr=1,4*NBS
              sum = sum -1._dp/(2._dp*PI*sqrt(CCC)*sqrt(p0j))*conjg(calFj0_Qmat(rr,mm,j))*exp(IU*CCC*p0j*t)*calE(nn,rr)
           end do
           P_Qmat(nn,mm,j) = P_Qmat(nn,mm,j) + sum/IU

        end do
     end do
  end do

  
  !---------------------------------------
  !   calculation of Q^N_NM (n=0)
  !---------------------------------------
  !----- h_MR calEa^*_RNj -----
  do nn=1,4*NBS
     do mm=1,4*NBS
        do j=1,Nph
           sum = (0._dp,0._dp)
           do rr=1,4*NBS
              sum = sum +h_Qmat(mm,rr)*conjg(calEa(rr,nn,j))
           end do
           Q_Qmat(nn,mm,j) = sum/IU
        end do
     end do
  end do

  !---------------------------------------
  ! calculation of derivatives(*DeltaT)
  !---------------------------------------
  !----- d(calE_NM)/dt = O_NM + O_NM^dagger -----
  do nn=1,4*NBS
     do mm=1,4*NBS
        dcalE(nn,mm) = (O_Qmat(nn,mm) +conjg(O_Qmat(mm,nn)))*DeltaT
     end do
  end do

  !----- d(calEa_NMj)/dt = (Q_NMj)^dagger + P_NMj -----
  do nn=1,4*NBS
     do mm=1,4*NBS
        do j=1,Nph
           dcalEa(nn,mm,j) = ( conjg(Q_Qmat(mm,nn,j)) +P_Qmat(nn,mm,j) )*DeltaT  
        end do
     end do
  end do


!!$  !---------------------------------------
!!$  !   calculation of Zme
!!$  !---------------------------------------
!!$
!!$  !----- I*calO_NN -----
!!$  sum = (0._dp,0._dp)
!!$  do nn=1,4*NBS
!!$     sum = sum +O_Qmat(nn,nn)
!!$  end do
!!$  zme_O = IU*sum
!!$
!!$  !----- VT_NM calE_NM -----
!!$  sum = (0._dp,0._dp)
!!$  do nn=1,4*NBS
!!$     do mm=1,4*NBS
!!$        sum = sum +VT_Qmat(nn,mm)*calE(nn,mm)
!!$     end do
!!$  end do
!!$  zme_VT = sum
!!$
!!$  !----- (NM|PQ) (calE_NM calE_PQ -calE_NQ calE_PM ) -----
!!$  !----- (NM|PQ) (calE_NM calE_PQ ) -----
!!$  sum = (0._dp,0._dp)
!!$  do nn=1,4*NBS
!!$     do mm=1,4*NBS
!!$        do pp=1,4*NBS
!!$           do qq=1,4*NBS
!!$              sum = sum + twoele_Qmat(nn,mm,pp,qq)*calE(nn,mm)*calE(pp,qq)
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$  if(use_exchange) then
!!$     !----- (NM|PQ) (-calE_NQ calE_PM ) -----
!!$     do nn=1,4*NBS
!!$        do mm=1,4*NBS
!!$           do pp=1,4*NBS
!!$              do qq=1,4*NBS
!!$                 sum = sum - twoele_Qmat(nn,mm,pp,qq)*calE(nn,qq)*calE(pp,mm)
!!$              end do
!!$           end do
!!$        end do
!!$     end do
!!$  end if
!!$  zme_te = sum
!!$  
!!$  !----- M_NM calE_NM -----
!!$  sum = (0._dp,0._dp)
!!$  do nn=1,4*NBS
!!$     do mm=1,4*NBS
!!$        sum = sum +M_Qmat(nn,mm)*calE(nn,mm)
!!$     end do
!!$  end do
!!$  zme_M = sum
!!$
!!$  Zme = (zme_O -zme_VT -zme_te)/zme_M
!!$  
  
  return
end subroutine calc_derivatives_BO_OPQ





!==============================================================================================
subroutine mass_renormalize(use_exchange,t,calE,calEa,VT_Qmat,M_Qmat,calFj0_Qmat,twoele_Qmat,Zme) 
!  Renormalize mass
!
!  131125 
!==============================================================================================
  use Precision
  use DiracOutput
  use Constants
  implicit none

  logical,intent(in) :: use_exchange  !.true. --> include exchange terms in diff. eq. , .false. --> without exchange terms
  real(kind=dp),intent(in) :: t  ! time
  
  complex(kind=dp),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  complex(kind=dp),intent(in) :: calEa(4*NBS,4*NBS,Nph) ! calEa_PQj
!  complex(kind=dp),intent(in) :: h_Qmat(4*NBS,4*NBS)  ! h_PQ
  complex(kind=dp),intent(in) :: VT_Qmat(4*NBS,4*NBS)  ! VT_PQ
  complex(kind=dp),intent(in) :: M_Qmat(4*NBS,4*NBS)  ! M_PQ
  complex(kind=dp),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
  complex(kind=dp),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)

  complex(kind=dp),intent(out) :: Zme ! new renormalization factor calculated after time evolution

  complex(kind=dp) :: h_Qmat(4*NBS,4*NBS)  ! h_PQ

  real(kind=dp) :: dp0_th_phi ! dp0*dtheta*dphi
  real(kind=dp) :: p0j,thj,phij ! p_j,theta_j,phi_j

  integer :: i,j,k,l
  integer :: pp,qq,rr,ss,nn,mm
  complex(kind=dp) :: sum

  complex(kind=dp) :: O_Qmat(4*NBS,4*NBS)
  complex(kind=dp) :: I_Qmat(4*NBS,4*NBS)
  
  complex(kind=dp) :: zme_te,zme_VT,zme_O,zme_M
  

  sum = (0._dp,0._dp)
  O_Qmat(:,:) = (0._dp,0._dp)

  !-----  -1/sqrt(2 pi^2 c)* int d^3p/sqrt(2*p0) 
  !        ( F^k_MR(p)e^k(p) exp(-i c p0 t) calEa_NRj + F^*k_RM(p)e^*k(p) exp(i c p0 t) calEa^*_RNj -----
  do nn=1,4*NBS
     do mm=1,4*NBS

        sum = (0._dp,0._dp)
        do j=1,Nph
           call get_pj(j,p0j,thj,phij,dp0_th_phi)
           do rr=1,4*NBS
              sum = sum -dp0_th_phi*p0j**2*sin(thj)/(2._dp*PI*sqrt(CCC)*sqrt(p0j)) *( calFj0_Qmat(mm,rr,j)*exp(-IU*CCC*p0j*t)*calEa(nn,rr,j) &
                   & + conjg(calFj0_Qmat(rr,mm,j))*exp(IU*CCC*p0j*t)*conjg(calEa(rr,nn,j)) )
           end do
        end do

!        O_Qmat(nn,mm) = O_Qmat(nn,mm) + sum/IU
        O_Qmat(nn,mm) = O_Qmat(nn,mm) + sum

     end do
  end do


  !---------------------------------------
  !   calculation of Zme
  !---------------------------------------

  !----- I*calO_NN -----
  sum = (0._dp,0._dp)
  do nn=1,4*NBS
     sum = sum +O_Qmat(nn,nn)
  end do
!  zme_O = IU*sum
  zme_O = sum

  !----- M_NM calE_NM -----
  sum = (0._dp,0._dp)
  do nn=1,4*NBS
     do mm=1,4*NBS
        sum = sum +M_Qmat(nn,mm)*calE(nn,mm)
     end do
  end do
  zme_M = sum

!  Zme = (zme_O -zme_VT -zme_te)/zme_M
  Zme = zme_O/zme_M +1._dp
  
  return
end subroutine mass_renormalize

!==============================================================================================
subroutine calc_derivatives_BO_OPQ_notRen(Use_exchange,t,calE,calEa,h_Qmat,calFj0_Qmat,twoele_Qmat,dcalE,dcalEa)
!  Derivatives for BO approximation including those for <e^dagger a e>
!
!  131125 
!==============================================================================================
  use Precision
  use DiracOutput
  use Constants
  implicit none

  logical,intent(in) :: use_exchange  !.true. --> include exchange terms in diff. eq. , .false. --> without exchange terms
  real(kind=dp),intent(in) :: t  ! time
  complex(kind=dp),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  complex(kind=dp),intent(in) :: calEa(4*NBS,4*NBS,Nph) ! calEa_PQj
  complex(kind=dp),intent(in) :: h_Qmat(4*NBS,4*NBS)  ! h_PQ
  complex(kind=dp),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
  complex(kind=dp),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)

  complex(kind=dp),intent(out) :: dcalE(4*NBS,4*NBS)  ! d(calE_PQ)/dt *DeltaT
  complex(kind=dp),intent(out) :: dcalEa(4*NBS,4*NBS,Nph) ! d(calEa_PQj)/dt *DeltaT

  real(kind=dp) :: dp0_th_phi ! dp0*dtheta*dphi
  real(kind=dp) :: p0j,thj,phij ! p_j,theta_j,phi_j

  integer :: i,j,k,l
  integer :: pp,qq,rr,ss,nn,mm
  complex(kind=dp) :: sum

  complex(kind=dp) :: O_Qmat(4*NBS,4*NBS), P_Qmat(4*NBS,4*NBS,Nph), Q_Qmat(4*NBS,4*NBS,Nph)
  complex(kind=dp) :: I_Qmat(4*NBS,4*NBS)
  
  sum = (0._dp,0._dp)
  O_Qmat(:,:) = (0._dp,0._dp)
  P_Qmat(:,:,:) = (0._dp,0._dp)
  Q_Qmat(:,:,:) = (0._dp,0._dp)

  !---------------------------------------
  !   calculation of O^N_NM
  !---------------------------------------
  !----- h_MR calE_NR -----
  do nn=1,4*NBS
     do mm=1,4*NBS
        sum = (0._dp,0._dp)
        do rr=1,4*NBS
           sum = sum +h_Qmat(mm,rr)*calE(nn,rr)
        end do
        O_Qmat(nn,mm) = sum/IU
     end do
  end do

  !----- (MR|PQ) calE_NR calE_PQ -----
  I_Qmat(:,:) = (0._dp,0._dp)
  do mm=1,4*NBS
     do rr=1,4*NBS
        sum = (0._dp,0._dp)
        do pp=1,4*NBS
           do qq=1,4*NBS
              sum = sum + twoele_Qmat(mm,rr,pp,qq)*calE(pp,qq)
           end do
        end do
        I_Qmat(mm,rr) = sum
     end do
  end do
  
  ! add to O^Phi_NM
  do nn=1,4*NBS
     do mm=1,4*NBS
        sum = (0._dp,0._dp)
        do rr=1,4*NBS
           sum = sum + I_Qmat(mm,rr)*calE(nn,rr)
        end do
        !        write(*,*) nn,mm,sum
        O_Qmat(nn,mm) = O_Qmat(nn,mm) + sum/IU
     end do
  end do
  
  if(use_exchange) then
     !-----  -(MR|PQ) calE_NQ calE_PR -----
     I_Qmat(:,:) = (0._dp,0._dp)
     do mm=1,4*NBS
        do qq=1,4*NBS
           sum = (0._dp,0._dp)
           do pp=1,4*NBS
              do rr=1,4*NBS
                 sum = sum +twoele_Qmat(mm,rr,pp,qq)*calE(pp,rr)
              end do
           end do
           I_Qmat(mm,qq) = sum
        end do
     end do
     
     ! add to O^Phi_NM
     do nn=1,4*NBS
        do mm=1,4*NBS
           sum = (0._dp,0._dp)
           do qq=1,4*NBS
              sum = sum + I_Qmat(mm,qq)*calE(nn,qq)
           end do
           O_Qmat(nn,mm) = O_Qmat(nn,mm) - sum/IU  ! note the minus sign here.
        end do
     end do
     
  end if
  
  !-----  -1/sqrt(2 pi^2 c)* int d^3p/sqrt(2*p0) 
  !        ( F^k_MR(p)e^k(p) exp(-i c p0 t) calEa_NRj + F^*k_RM(p)e^*k(p) exp(i c p0 t) calEa^*_RNj -----
  do nn=1,4*NBS
     do mm=1,4*NBS

        sum = (0._dp,0._dp)
        do j=1,Nph
           call get_pj(j,p0j,thj,phij,dp0_th_phi)
           do rr=1,4*NBS
              sum = sum -dp0_th_phi*p0j**2*sin(thj)/(2._dp*PI*sqrt(CCC)*sqrt(p0j)) *( calFj0_Qmat(mm,rr,j)*exp(-IU*CCC*p0j*t)*calEa(nn,rr,j) &
                   & + conjg(calFj0_Qmat(rr,mm,j))*exp(IU*CCC*p0j*t)*conjg(calEa(rr,nn,j)) )
           end do
        end do

        O_Qmat(nn,mm) = O_Qmat(nn,mm) + sum/IU

     end do
  end do


  !---------------------------------------
  !   calculation of P^N_NM (n=0)
  !---------------------------------------
  !----- h_MR calEa_NRj -----
  do nn=1,4*NBS
     do mm=1,4*NBS
        do j=1,Nph
           sum = (0._dp,0._dp)
           do rr=1,4*NBS
              sum = sum +h_Qmat(mm,rr)*calEa(nn,rr,j)
           end do
           P_Qmat(nn,mm,j) = sum/IU
        end do
     end do
  end do

  do nn=1,4*NBS
     do mm=1,4*NBS
        do j=1,Nph

           call get_pj(j,p0j,thj,phij,dp0_th_phi)
           sum = (0._dp,0._dp)
           do rr=1,4*NBS
              sum = sum -1._dp/(2._dp*PI*sqrt(CCC)*sqrt(p0j))*conjg(calFj0_Qmat(rr,mm,j))*exp(IU*CCC*p0j*t)*calE(nn,rr)
           end do
           P_Qmat(nn,mm,j) = P_Qmat(nn,mm,j) + sum/IU

        end do
     end do
  end do

  
  !---------------------------------------
  !   calculation of Q^N_NM (n=0)
  !---------------------------------------
  !----- h_MR calEa^*_RNj -----
  do nn=1,4*NBS
     do mm=1,4*NBS
        do j=1,Nph
           sum = (0._dp,0._dp)
           do rr=1,4*NBS
              sum = sum +h_Qmat(mm,rr)*conjg(calEa(rr,nn,j))
           end do
           Q_Qmat(nn,mm,j) = sum/IU
        end do
     end do
  end do

  !---------------------------------------
  ! calculation of derivatives(*DeltaT)
  !---------------------------------------
  !----- d(calE_NM)/dt = O_NM + O_NM^dagger -----
  do nn=1,4*NBS
     do mm=1,4*NBS
        dcalE(nn,mm) = (O_Qmat(nn,mm) +conjg(O_Qmat(mm,nn)))*DeltaT
     end do
  end do

  !----- d(calEa_NMj)/dt = (Q_NMj)^dagger + P_NMj -----
  do nn=1,4*NBS
     do mm=1,4*NBS
        do j=1,Nph
           dcalEa(nn,mm,j) = ( conjg(Q_Qmat(mm,nn,j)) +P_Qmat(nn,mm,j) )*DeltaT  
        end do
     end do
  end do
  
  return
end subroutine calc_derivatives_BO_OPQ_notRen

!------------------------------------------------------------------------------------------
!
!  subroutines below are not use now.
!
!------------------------------------------------------------------------------------------

!========================================================================================================
subroutine calc_derivatives_NonBO(t,calE,calEV,calDV,calC,calCV,calBV, & 
     & TM_Qmat,Tnuc_mat,twoele_Qmat,nucele_Qmat,twonuc_mat,calFj0_Qmat,alpha, &
     & dcalE,dcalEV,dcalDV,dcalC,dcalCV,dcalBV)
!========================================================================================================
  Use DiracOutput
  use Constants
  use NucBasis
  implicit none

  real(kind=8),intent(in) :: t  ! time
  complex(kind=8),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  complex(kind=8),intent(in) :: calEV(4*NBS,4*NBS)  ! calEV_PQ
  complex(kind=8),intent(in) :: calDV(4*NBS,4*NBS)  ! calDV_PQ
  complex(kind=8),intent(in) :: calC(NBS_N,NBS_N)  ! calC_ij
  complex(kind=8),intent(in) :: calCV(NBS_N,NBS_N)  ! calC_ij
  complex(kind=8),intent(in) :: calBV(NBS_N,NBS_N)  ! calC_ij
  complex(kind=8),intent(in) :: TM_Qmat(4*NBS,4*NBS)  ! T_PQ + M_PQ
  complex(kind=8),intent(in) :: Tnuc_mat(NBS_N,NBS_N)  ! T_Nij
  complex(kind=8),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)
  complex(kind=8),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
  complex(kind=8),intent(in) :: alpha(Nph)   ! <coherent|hat(a)|coherent>
  complex(kind=8),intent(in) :: nucele_Qmat(NBS_N,NBS_N,4*NBS,4*NBS)  ! (ij|PQ)
  complex(kind=8),intent(in) :: twonuc_mat(NBS_N,NBS_N,NBS_N,NBS_N)  ! (ij|kl)

  complex(kind=8),intent(out) :: dcalE(4*NBS,4*NBS)  ! d(calE_PQ)/dt *DeltaT
  complex(kind=8),intent(out) :: dcalEV(4*NBS,4*NBS)  ! d(calEV_PQ)/dt *DeltaT
  complex(kind=8),intent(out) :: dcalDV(4*NBS,4*NBS)  ! d(calDV_PQ)/dt *DeltaT
  complex(kind=8),intent(out) :: dcalC(NBS_N,NBS_N)  ! d(calC_ij)/dt *DeltaT
  complex(kind=8),intent(out) :: dcalCV(NBS_N,NBS_N)  ! d(calCV_ij)/dt *DeltaT
  complex(kind=8),intent(out) :: dcalBV(NBS_N,NBS_N)  ! d(calBV_ij)/dt *DeltaT


  real(kind=8) :: dpj ! tilde(Delta p)_j
  real(kind=8) :: p0j ! p_j

  integer :: i,j,k,l,p,q,r,s
  complex(kind=8) :: sum_te,sum_ne,sum_tn
  complex(kind=8) :: sumE,sumEV,sumDV,sumC,sumCV,sumBV
  complex(kind=8) :: dcalF

  complex(kind=8) :: I2_Qmat(4*NBS,4*NBS), I4_Qmat(4*NBS,4*NBS) 
  complex(kind=8) :: I_Qmat(4*NBS,4*NBS)  ! I_PQ = (I_1+I_2+I_3+I_4)_PQ = (TM +I2 +I4)_PQ
  complex(kind=8) :: I4_mat(NBS_N,NBS_N) 
  complex(kind=8) :: I_mat(NBS_N,NBS_N)  ! I_ij = (I_1+I_2+I_3+I_4)_ij = (T +I4)_ij

!  complex(kind=8) :: I4V_Qmat(4*NBS,4*NBS),IV_Qmat(4*NBS,4*NBS)    ! for VEVs
!  complex(kind=8) :: I4V_mat(NBS_N,NBS_N),IV_mat(NBS_N,NBS_N)
!  complex(kind=8) :: sum_teV,sum_ne,sum_tnV


  I2_Qmat  = (0.d0,0.d0)
  I4_Qmat  = (0.d0,0.d0)

  ! calculation of I4
  do p=1,4*NBS
     do q=1,4*NBS

        sum_te = (0.d0,0.d0)
        do r=1,4*NBS
           do s=1,4*NBS
              sum_te = sum_te + twoele_Qmat(p,q,r,s)*calE(r,s)
           end do
        end do
        sum_ne = (0.d0,0.d0)
        do i=1,NBS_N
           do j=1,NBS_N
              sum_ne = sum_ne + nucele_Qmat(i,j,p,q)*calC(i,j)  ! (PQ|ij)=(ij|PQ)
           end do
        end do
        I4_Qmat(p,q) = sum_te + sum_ne

     end do
  end do

!!$  ! calculation of I2 (Arad only)
!!$  do p=1,4*NBS
!!$     do q=1,4*NBS
!!$
!!$        sum_calF=(0.d0,0.d0)
!!$        do j=1,Nph
!!$           call calc_DeltaPj_and_P0j(j,dpj,p0j)
!!$           
!!$           dcalF = calFj0_Qmat(nn,mm,j)       *exp(-IU*CCC*p0j*t)*alpha(j) &
!!$                & +conjg(calFj0_Qmat(mm,nn,j))*exp(+IU*CCC*p0j*t)*conjg(alpha(j))
!!$           
!!$           sum_calF = sum_calF -dpj*dcalF
!!$
!!$        end do
!!$        
!!$        I2_Qmat(p,q) = sum_calF
!!$
!!$     end do
!!$  end do

  I_Qmat(:,:) = TM_Qmat(:,:) + I2_Qmat(:,:) + I4_Qmat(:,:)

  ! calculation of derivatives(*DeltaT)
  do p=1,4*NBS
     do q=p,4*NBS
        sumE = (0.d0,0.d0);sumEV = (0.d0,0.d0);sumDV = (0.d0,0.d0)
        do r=1,4*NBS
           sumE  = sumE  -I_Qmat(r,p)*calE(r,q)  +I_Qmat(q,r)*calE(p,r)
           sumEV = sumEV -I_Qmat(r,p)*calEV(r,q) +I_Qmat(q,r)*calEV(p,r)
           sumDV = sumDV +I_Qmat(p,r)*calDV(r,q) -I_Qmat(r,q)*calDV(p,r)
        end do
        dcalE(p,q)  = -IU*sumE *DeltaT
        dcalEV(p,q) = -IU*sumEV*DeltaT
        dcalDV(p,q) = -IU*sumDV*DeltaT
     end do
  end do
  
  do p=2,4*NBS
     do q=1,p-1
        dcalE(p,q)  = conjg(dcalE(q,p))
        dcalEV(p,q) = conjg(dcalEV(q,p))
        dcalDV(p,q) = conjg(dcalDV(q,p))
     end do
  end do
  
  !--------------------------------------------------------------------------

  I4_mat  = (0.d0,0.d0)

  ! calculation of I4
  do i=1,NBS_N
     do j=1,NBS_N

        sum_ne = (0.d0,0.d0)
        do p=1,4*NBS
           do q=1,4*NBS
              sum_ne = sum_ne + nucele_Qmat(i,j,p,q)*calE(p,q)
           end do
        end do
        sum_tn = (0.d0,0.d0)
        do k=1,NBS_N
           do l=1,NBS_N
              sum_tn = sum_tn + twonuc_mat(i,j,k,l)*calC(k,l)  
           end do
        end do
        I4_mat(i,j) = sum_tn + sum_ne

     end do
  end do

  I_mat(:,:) = Tnuc_mat(:,:) + I4_mat(:,:)

  do i=1,NBS_N
     do j=i,NBS_N
        sumC = (0.d0,0.d0); sumCV = (0.d0,0.d0); sumBV = (0.d0,0.d0);
        do k=1,NBS_N
           sumC  = sumC  -I_mat(k,i)*calC(k,j)  +I_mat(j,k)*calC(i,k)
           sumCV = sumCV -I_mat(k,i)*calCV(k,j) +I_mat(j,k)*calCV(i,k)
           sumBV = sumBV +I_mat(i,k)*calBV(k,j) -I_mat(k,j)*calBV(i,k)
        end do
        dcalC(i,j)  = -IU*sumC *DeltaT
        dcalCV(i,j) = -IU*sumCV*DeltaT
        dcalBV(i,j) = -IU*sumBV*DeltaT
     end do
  end do
  
  do i=2,NBS_N
     do j=1,i-1
        dcalC(i,j)  = conjg(dcalC(j,i))
        dcalCV(i,j) = conjg(dcalCV(j,i))
        dcalBV(i,j) = conjg(dcalBV(j,i))
     end do
  end do



  
  return
end subroutine calc_derivatives_NonBO







!==============================================================================================
subroutine calc_derivatives_BO(t,calE,calEV,calDV,h_Qmat,twoele_Qmat,calFj0_Qmat,alpha, &
     & dcalE,dcalEV,dcalDV)
! 
!==============================================================================================
  use DiracOutput
  use Constants
  implicit none

  real(kind=8),intent(in) :: t  ! time
  complex(kind=8),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  complex(kind=8),intent(in) :: calEV(4*NBS,4*NBS)  ! calEV_PQ
  complex(kind=8),intent(in) :: calDV(4*NBS,4*NBS)  ! calDV_PQ
  complex(kind=8),intent(in) :: h_Qmat(4*NBS,4*NBS)  ! h_PQ
  complex(kind=8),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)
  complex(kind=8),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
  complex(kind=8),intent(in) :: alpha(Nph)   ! <coherent|hat(a)|coherent>

  complex(kind=8),intent(out) :: dcalE(4*NBS,4*NBS)  ! d(calE_PQ)/dt *DeltaT
  complex(kind=8),intent(out) :: dcalEV(4*NBS,4*NBS)  ! d(calEV_PQ)/dt *DeltaT
  complex(kind=8),intent(out) :: dcalDV(4*NBS,4*NBS)  ! d(calDV_PQ)/dt *DeltaT

  real(kind=8) :: dpj ! tilde(Delta p)_j
  real(kind=8) :: p0j ! p_j

  integer :: i,j,k,l
  integer :: pp,qq,rr,ss,nn,mm
  complex(kind=8) :: sum_te,sum_calF,sum,sumEV,sumDV
  complex(kind=8) :: dcalF

  complex(kind=8) :: I2_Qmat(4*NBS,4*NBS), I4_Qmat(4*NBS,4*NBS) 
  complex(kind=8) :: I_Qmat(4*NBS,4*NBS)  ! I_PQ = (I_1+I_2+I_3+I_4)_PQ = (TM +I2 +I4)_PQ


  I2_Qmat  = (0.d0,0.d0)
  I4_Qmat  = (0.d0,0.d0)

  sum_calF=(0.d0,0.d0); sum_te=(0.d0,0.d0)

  ! calculation of I4 (h_PQ + (PQ|RS)*calE_RS )
  do nn=1,4*NBS
     do mm=1,4*NBS

        sum_te = (0.d0,0.d0)
        do rr=1,4*NBS
           do ss=1,4*NBS
              sum_te = sum_te + twoele_Qmat(nn,mm,rr,ss)*calE(rr,ss)
           end do
        end do

        I4_Qmat(nn,mm) = sum_te + h_Qmat(nn,mm)

     end do
  end do

  ! calculation of I2 (Arad only)
  do nn=1,4*NBS
     do mm=1,4*NBS

        sum_calF=(0.d0,0.d0)
        do j=1,Nph
           call calc_DeltaPj_and_P0j(j,dpj,p0j)
           
           dcalF = calFj0_Qmat(nn,mm,j)       *exp(-IU*CCC*p0j*t)*alpha(j) &
                & +conjg(calFj0_Qmat(mm,nn,j))*exp(+IU*CCC*p0j*t)*conjg(alpha(j))
           
           sum_calF = sum_calF -dpj*dcalF

        end do
        
        I2_Qmat(nn,mm) = sum_calF

     end do
  end do


  I_Qmat(:,:) = I2_Qmat(:,:) + I4_Qmat(:,:)


  ! calculation of derivatives(*DeltaT)
  do pp=1,4*NBS
     do qq=pp,4*NBS
        sum = (0.d0,0.d0);sumEV = (0.d0,0.d0);sumDV = (0.d0,0.d0)
        do rr=1,4*NBS
           sum   = sum   -I_Qmat(rr,pp)*calE(rr,qq)  +I_Qmat(qq,rr)*calE(pp,rr)
           sumEV = sumEV -I_Qmat(rr,pp)*calEV(rr,qq) +I_Qmat(qq,rr)*calEV(pp,rr)
           sumDV = sumDV +I_Qmat(pp,rr)*calDV(rr,qq) -I_Qmat(rr,qq)*calDV(pp,rr)
        end do
        dcalE(pp,qq)  = -IU*sum*DeltaT
        dcalEV(pp,qq) = -IU*sumEV*DeltaT
        dcalDV(pp,qq) = -IU*sumDV*DeltaT
     end do
  end do
  
  do pp=2,4*NBS
     do qq=1,pp-1
        dcalE(pp,qq)  = conjg(dcalE(qq,pp))
        dcalEV(pp,qq) = conjg(dcalEV(qq,pp))
        dcalDV(pp,qq) = conjg(dcalDV(qq,pp))
     end do
  end do

  return
end subroutine calc_derivatives_BO


















!====================================================================================
!subroutine calc_derivative(nn,mm,it,calE,h_Qmat,calFj0_Qmat,alpha,dcalEdt)
subroutine calc_derivative(nn,mm,it,calE,h_Qmat,twoele_Qmat,calFj0_Qmat,alpha,dcalEdt)
!====================================================================================
  use DiracOutput
  use Constants
  implicit none

  integer,intent(in) :: nn,mm ! compute (nn,mm) component of d(calE_NM)/dt
  integer,intent(in) :: it  ! it-th time step
  complex(kind=8),intent(in) :: calE(4*NBS,4*NBS)  ! calE_NM
  complex(kind=8),intent(in) :: h_Qmat(4*NBS,4*NBS)  ! h_NM
  complex(kind=8),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (NM|PQ)
  complex(kind=8),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_NMj(t=0)
  complex(kind=8),intent(in) :: alpha(Nph)   ! <coherent|hat(a)|coherent>

  complex(kind=8),intent(out) :: dcalEdt  ! d(calE_NM)/dt

  real(kind=8) :: dpj ! tilde(Delta p)_j
  real(kind=8) :: p0j ! p_j
  real(kind=8) :: t ! time = it*DeltaT

  integer :: i,j,k,l
  integer :: rr,pp,qq
  complex(kind=8) :: sum_h,sum_te,sum_calF,sum_dcalF
  complex(kind=8) :: dcalF

  sum_h =(0.d0,0.d0); sum_calF=(0.d0,0.d0); sum_te=(0.d0,0.d0)

  ! it starts from 0 -> t starts from 0
  t = DeltaT*it  ! DeltaT is defined in Constants and set in the main routine.
  

  ! h^*_NR *calE_RM - h_MR *calE_NR
  do rr=1,4*NBS
     sum_h = sum_h +conjg(h_Qmat(nn,rr))*calE(rr,mm) -h_Qmat(mm,rr)*calE(nn,rr)
!     write(*,"(3i6,6es16.6)") nn,mm,rr,h_Qmat(mm,rr),calE(nn,rr),sum_h
  end do

  ! two electron integral 
  do rr=1,4*NBS
     do pp=1,4*NBS
        do qq=1,4*NBS
           
           sum_te = sum_te &
                &  + conjg(twoele_Qmat(nn,rr,pp,qq))*calE(rr,mm)*calE(qq,pp) &
                &  - twoele_Qmat(mm,rr,pp,qq)       *calE(pp,qq)*calE(nn,rr)
           
        end do
     end do
  end do

!  write(*,"(2i6,4es16.6)") nn,mm,sum_h,sum_te

  
  ! Arad part
  do j=1,Nph
     call calc_DeltaPj_and_P0j(j,dpj,p0j)
     sum_dcalF=(0.d0,0.d0)
     do rr=1,4*NBS
        dcalF = conjg(calFj0_Qmat(nn,rr,j))*exp(+IU*CCC*p0j*t)*calE(rr,mm)*conjg(alpha(j)) &
             & +calFj0_Qmat(rr,nn,j)       *exp(-IU*CCC*p0j*t)*calE(rr,mm)*alpha(j) &
             & -calFj0_Qmat(mm,rr,j)       *exp(-IU*CCC*p0j*t)*calE(nn,rr)*alpha(j) &
             & -conjg(calFj0_Qmat(rr,mm,j))*exp(+IU*CCC*p0j*t)*calE(nn,rr)*conjg(alpha(j)) 
        sum_dcalF = sum_dcalF +dcalF
     end do
!     write(*,"(3i6,2es16.6,4es16.6)") nn,mm,j,dpj,p0j,calFj0_Qmat(nn,rr,j),sum_calF
!     write(*,"(3i6,2es16.6)") nn,mm,j,sum_calF
     sum_calF = sum_calF -dpj*sum_dcalF
!     write(*,"(3i6,4es16.6)") nn,mm,j,sum_dcalF,sum_calF
  end do
!  write(*,"(2i6,2es16.6)") nn,mm,sum_calF
  
!  stop
  dcalEdt = IU*(sum_h+sum_te+sum_calF)
!  dcalEdt = IU*(sum_h+sum_calF)
!  dcalEdt = IU*(sum_h)
!  dcalEdt = IU*(sum_calF)
!  write(*,"(2i6,4es16.6)") nn,mm,sum_calF,dcalEdt
!  dcalEdt = IU*(sum_h)
  
  return
end subroutine calc_derivative



