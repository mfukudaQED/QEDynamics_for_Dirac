! Last Change:08-Jun-2014.
!================================================================================
! subroutines for derivatives of Rigged QED
!
!================================================================================

!!========================================================================================================
!subroutine calc_derivatives_NonBOwithEx(use_exchange,t,calE,calC, & 
!     & TM_Qmat,Tnuc_mat,twoele_Qmat,nucele_Qmat,twonuc_mat,calFj0_Qmat,alpha, &
!     & dcalE,dcalC)
!!  derivatives for Non-BO calculation
!!  (120712)
!!========================================================================================================
!  use Precision
!  use DiracOutput
!  use Constants
!  use NucBasis
!  implicit none
!
!  logical,intent(in) :: use_exchange  !.true. --> include exchange terms in diff. eq. , .false. --> without exchange terms
!  real(kind=dp),intent(in) :: t  ! time
!  complex(kind=dp),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
!  complex(kind=dp),intent(in) :: calC(NBS_N,NBS_N)  ! calC_ij
!  complex(kind=dp),intent(in) :: TM_Qmat(4*NBS,4*NBS)  ! T_PQ + M_PQ
!  complex(kind=dp),intent(in) :: Tnuc_mat(NBS_N,NBS_N)  ! T_Nij
!  complex(kind=dp),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)
!  complex(kind=dp),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
!  complex(kind=dp),intent(in) :: alpha(Nph)   ! <coherent|hat(a)|coherent>
!  complex(kind=dp),intent(in) :: nucele_Qmat(NBS_N,NBS_N,4*NBS,4*NBS)  ! (ij|PQ)
!  complex(kind=dp),intent(in) :: twonuc_mat(NBS_N,NBS_N,NBS_N,NBS_N)  ! (ij|kl)
!
!  complex(kind=dp),intent(out) :: dcalE(4*NBS,4*NBS)  ! d(calE_PQ)/dt *DeltaT
!  complex(kind=dp),intent(out) :: dcalC(NBS_N,NBS_N)  ! d(calC_ij)/dt *DeltaT
!
!  real(kind=dp) :: dpj ! tilde(Delta p)_j
!  real(kind=dp) :: p0j ! p_j
!
!  integer :: i,j,k,l,p,q,r,s,n,m
!  complex(kind=dp) :: sum_te,sum_ne,sum_tn
!  complex(kind=dp) :: sumE,sumC
!  complex(kind=dp) :: dcalF
!
!  complex(kind=dp) :: I2_Qmat(4*NBS,4*NBS), I4_Qmat(4*NBS,4*NBS) 
!  complex(kind=dp) :: I_Qmat(4*NBS,4*NBS)  ! I_PQ = (I_1+I_2+I_3+I_4)_PQ = (TM +I2 +I4)_PQ
!  complex(kind=dp) :: I4_mat(NBS_N,NBS_N) 
!  complex(kind=dp) :: I_mat(NBS_N,NBS_N)  ! I_ij = (I_1+I_2+I_3+I_4)_ij = (T +I4)_ij
!
!  I2_Qmat  = (0._dp,0._dp)
!  I4_Qmat  = (0._dp,0._dp)
!
!  ! calculation of I4
!  do p=1,4*NBS
!     do q=1,4*NBS
!
!        sum_te = (0._dp,0._dp)
!        do r=1,4*NBS
!           do s=1,4*NBS
!              sum_te = sum_te + twoele_Qmat(p,q,r,s)*calE(r,s)
!           end do
!        end do
!        sum_ne = (0._dp,0._dp)
!        do i=1,NBS_N
!           do j=1,NBS_N
!              sum_ne = sum_ne + nucele_Qmat(i,j,p,q)*calC(i,j)  ! (PQ|ij)=(ij|PQ)
!           end do
!        end do
!        I4_Qmat(p,q) = sum_te + sum_ne
!
!     end do
!  end do
!
!!!$  ! calculation of I2 (Arad only)
!!!$  do p=1,4*NBS
!!!$     do q=1,4*NBS
!!!$
!!!$        sum_calF=(0._dp,0._dp)
!!!$        do j=1,Nph
!!!$           call calc_DeltaPj_and_P0j(j,dpj,p0j)
!!!$           
!!!$           dcalF = calFj0_Qmat(nn,mm,j)       *exp(-IU*CCC*p0j*t)*alpha(j) &
!!!$                & +conjg(calFj0_Qmat(mm,nn,j))*exp(+IU*CCC*p0j*t)*conjg(alpha(j))
!!!$           
!!!$           sum_calF = sum_calF -dpj*dcalF
!!!$
!!!$        end do
!!!$        
!!!$        I2_Qmat(p,q) = sum_calF
!!!$
!!!$     end do
!!!$  end do
!
!  I_Qmat(:,:) = TM_Qmat(:,:) + I2_Qmat(:,:) + I4_Qmat(:,:)
!
!  ! calculation of derivatives(*DeltaT)
!  do p=1,4*NBS
!     do q=p,4*NBS
!        sumE = (0._dp,0._dp)
!        ! non-exchange terms
!        do r=1,4*NBS
!           sumE  = sumE  -I_Qmat(r,p)*calE(r,q)  +I_Qmat(q,r)*calE(p,r)
!        end do
!
!        if(use_exchange) then
!           ! exchange terms
!           do r=1,4*NBS
!              do n=1,4*NBS
!                 do m=1,4*NBS
!                    sumE = sumE +twoele_Qmat(r,p,n,m)*calE(r,m)*calE(n,q) -twoele_Qmat(q,r,n,m)*calE(p,m)*calE(n,r) 
!                 end do
!              end do
!           end do
!        end if
!        
!        dcalE(p,q)  = -IU*sumE*DeltaT
!     end do
!  end do
!  
!  do p=2,4*NBS
!     do q=1,p-1
!        dcalE(p,q)  = conjg(dcalE(q,p))
!     end do
!  end do
!  
!  !--------------------------------------------------------------------------
!
!  I4_mat  = (0._dp,0._dp)
!
!  ! calculation of I4
!  do i=1,NBS_N
!     do j=1,NBS_N
!
!        sum_ne = (0._dp,0._dp)
!        do p=1,4*NBS
!           do q=1,4*NBS
!              sum_ne = sum_ne + nucele_Qmat(i,j,p,q)*calE(p,q)
!           end do
!        end do
!        sum_tn = (0._dp,0._dp)
!        do k=1,NBS_N
!           do l=1,NBS_N
!              sum_tn = sum_tn + twonuc_mat(i,j,k,l)*calC(k,l)  
!           end do
!        end do
!        I4_mat(i,j) = sum_tn + sum_ne
!
!     end do
!  end do
!
!  I_mat(:,:) = Tnuc_mat(:,:) + I4_mat(:,:)
!
!  do i=1,NBS_N
!     do j=i,NBS_N
!
!        sumC = (0._dp,0._dp)
!        ! non-exchange terms
!        do k=1,NBS_N
!           sumC  = sumC  -I_mat(k,i)*calC(k,j)  +I_mat(j,k)*calC(i,k)
!        end do
!
!        if(use_exchange) then
!           ! exchange terms
!           if(NUCTYPE.eq.0) then ! boson
!              do k=1,NBS_N
!                 do n=1,NBS_N
!                    do m=1,NBS_N
!                       sumC = sumC -twonuc_mat(k,i,n,m)*calC(k,m)*calC(n,j) +twonuc_mat(j,k,n,m)*calC(i,m)*calC(n,k)
!                    end do
!                 end do
!              end do
!           elseif(NUCTYPE.eq.1) then ! fermion
!              do k=1,NBS_N
!                 do n=1,NBS_N
!                    do m=1,NBS_N
!                       sumC = sumC +twonuc_mat(k,i,n,m)*calC(k,m)*calC(n,j) -twonuc_mat(j,k,n,m)*calC(i,m)*calC(n,k)
!                    end do
!                 end do
!              end do
!           else
!              write(*,*) "NUCTYPE should be 0 or 1."
!              stop
!           end if
!        end if
!
!        dcalC(i,j)  = -IU*sumC*DeltaT
!     end do
!  end do
!  
!  do i=2,NBS_N
!     do j=1,i-1
!        dcalC(i,j)  = conjg(dcalC(j,i))
!     end do
!  end do
!
!  return
!end subroutine calc_derivatives_NonBOwithEx

!!==============================================================================================
!subroutine calc_derivatives_BOwithEx(Use_exchange,t,calE,calEV,calDV, &
!     & h_Qmat,twoele_Qmat,calFj0_Qmat,alpha,dcalE,dcalEV,dcalDV,jAM_Qmat)
!!  derivatives for BO approximation
!!  (120712)
!!==============================================================================================
!  use Precision
!  use DiracOutput
!  use Constants
!  implicit none
!
!  logical,intent(in) :: use_exchange  !.true. --> include exchange terms in diff. eq. , .false. --> without exchange terms
!  real(kind=dp),intent(in) :: t  ! time
!  complex(kind=dp),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
!  complex(kind=dp),intent(in) :: calEV(4*NBS,4*NBS)  ! calEV_PQ
!  complex(kind=dp),intent(in) :: calDV(4*NBS,4*NBS)  ! calDV_PQ
!  complex(kind=dp),intent(in) :: h_Qmat(4*NBS,4*NBS)  ! h_PQ
!  complex(kind=dp),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)
!  complex(kind=dp),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
!  complex(kind=dp),intent(in) :: alpha(Nph)   ! <coherent|hat(a)|coherent>
!  complex(kind=dp),intent(in) :: jAM_Qmat(4*NBS,4*NBS)  ! jAM_PQ
!
!  complex(kind=dp),intent(out) :: dcalE(4*NBS,4*NBS)  ! d(calE_PQ)/dt *DeltaT
!  complex(kind=dp),intent(out) :: dcalEV(4*NBS,4*NBS)  ! d(calEV_PQ)/dt *DeltaT
!  complex(kind=dp),intent(out) :: dcalDV(4*NBS,4*NBS)  ! d(calDV_PQ)/dt *DeltaT
!
!  real(kind=dp) :: dpj ! tilde(Delta p)_j
!  real(kind=dp) :: p0j ! p_j
!
!  integer :: i,j,k,l
!  integer :: pp,qq,rr,ss,nn,mm
!  complex(kind=dp) :: sum_te,sum_calF,sum,sumEV,sumDV
!  complex(kind=dp) :: dcalF
!
!  complex(kind=dp) :: I2_Qmat(4*NBS,4*NBS), I4_Qmat(4*NBS,4*NBS) 
!  complex(kind=dp) :: I_Qmat(4*NBS,4*NBS)  ! I_PQ = (I_1+I_2+I_3+I_4)_PQ = (TM +I2 +I4)_PQ
!
!  complex(kind=dp) :: I4V_Qmat(4*NBS,4*NBS),IV_Qmat(4*NBS,4*NBS)    ! for VEVs
!  complex(kind=dp) :: sum_teV
!
!  I2_Qmat  = (0._dp,0._dp)
!  I4_Qmat  = (0._dp,0._dp)
!  I4V_Qmat  = (0._dp,0._dp)
!
!  sum_calF=(0._dp,0._dp); sum_te=(0._dp,0._dp)
!
!  ! calculation of I4 (h_PQ + (PQ|RS)*calE_RS ) and I4V (h_PQ + (PQ|RS)*calEV_RS )
!  !$omp parallel do default(none), shared(I4_Qmat,I4V_Qmat,twoele_Qmat,h_Qmat,calE,calEV,NBS), private(nn,mm,rr,ss,sum_te,sum_teV)
!  do nn=1,4*NBS
!     do mm=1,4*NBS
!
!        sum_te = (0._dp,0._dp)
!        sum_teV = (0._dp,0._dp)
!!!$ !notwoelecalc
!!!$        do rr=1,4*NBS
!!!$           do ss=1,4*NBS
!!!$              sum_te = sum_te + twoele_Qmat(nn,mm,rr,ss)*calE(rr,ss)
!!!$              sum_teV = sum_teV + twoele_Qmat(nn,mm,rr,ss)*calEV(rr,ss)
!!!$           end do
!!!$        end do
!
!        I4_Qmat(nn,mm) = sum_te + h_Qmat(nn,mm)
!        I4V_Qmat(nn,mm) = sum_teV + h_Qmat(nn,mm)
!
!     end do
!  end do
!  !$omp end parallel do
!
!  ! calculation of I2 (Arad only)
!  !$omp parallel do default(none), shared(I2_Qmat,calFj0_Qmat,calE,alpha,t,NBS,Nph,PI,jAM_Qmat),&
!                            !$omp& private(nn,mm,j,dpj,p0j,sum_calF,dcalF)
!  do nn=1,4*NBS
!     do mm=1,4*NBS
!
!        sum_calF=(0._dp,0._dp)
!        do j=1,Nph
!           call calc_DeltaPj_and_P0j(j,dpj,p0j) ! sub_int.f90
!           ! we do not use dpj in our modified version 120828
!!           write(*,*) j,dpj,p0j,alpha(j)
!           
!           dcalF = calFj0_Qmat(nn,mm,j)       *exp(-IU*CCC*p0j*t)*alpha(j) &
!                & +conjg(calFj0_Qmat(mm,nn,j))*exp(+IU*CCC*p0j*t)*conjg(alpha(j))
!           
!!           sum_calF = sum_calF -dpj*dcalF
!           sum_calF = sum_calF -dcalF/(2._dp*PI*sqrt(p0j*CCC))
!
!        end do
!!        stop
!
!!!!
!        I2_Qmat(nn,mm) = sum_calF
!!        I2_Qmat(nn,mm) = sum_calF + jAM_Qmat(nn,mm)
!!!!
!!!$        I2_Qmat(nn,mm) = sum_calF
!!!$
!     end do
!  end do
!  !$omp end parallel do
!
!
!  I_Qmat(:,:) = I2_Qmat(:,:) + I4_Qmat(:,:)
!  IV_Qmat(:,:) = I2_Qmat(:,:) + I4V_Qmat(:,:)
!
!
!  ! calculation of derivatives(*DeltaT)
!  do pp=1,4*NBS
!     do qq=pp,4*NBS
!        sum = (0._dp,0._dp);sumEV = (0._dp,0._dp);sumDV = (0._dp,0._dp)
!        
!        ! non-exchange terms
!        do rr=1,4*NBS
!           sum   = sum   -I_Qmat(rr,pp) *calE(rr,qq)  +I_Qmat(qq,rr) *calE(pp,rr)
!           sumEV = sumEV -IV_Qmat(rr,pp)*calEV(rr,qq) +IV_Qmat(qq,rr)*calEV(pp,rr)
!           sumDV = sumDV +IV_Qmat(pp,rr)*calDV(rr,qq) -IV_Qmat(rr,qq)*calDV(pp,rr)
!        end do
!        
!        if(use_exchange) then
!           ! exchange terms
!           do rr=1,4*NBS
!              do nn=1,4*NBS
!                 do mm=1,4*NBS
!                    sum   = sum   +twoele_Qmat(rr,pp,nn,mm)*calE(rr,mm) *calE(nn,qq)  -twoele_Qmat(qq,rr,nn,mm)*calE(pp,mm) *calE(nn,rr) 
!                    sumEV = sumEV +twoele_Qmat(rr,pp,nn,mm)*calEV(rr,mm)*calEV(nn,qq) -twoele_Qmat(qq,rr,nn,mm)*calEV(pp,mm)*calEV(nn,rr) 
!                    sumDV = sumDV -twoele_Qmat(pp,rr,nn,mm)*calEV(nn,rr)*calDV(mm,qq) +twoele_Qmat(rr,qq,nn,mm)*calDV(pp,nn)*calEV(rr,mm) 
!                 end do
!              end do
!           end do
!        end if
!
!        dcalE(pp,qq)  = -IU*sum*DeltaT
!        dcalEV(pp,qq) = -IU*sumEV*DeltaT
!        dcalDV(pp,qq) = -IU*sumDV*DeltaT
!     end do
!  end do
!     
!  do pp=2,4*NBS
!     do qq=1,pp-1
!        dcalE(pp,qq)  = conjg(dcalE(qq,pp))
!        dcalEV(pp,qq) = conjg(dcalEV(qq,pp))
!        dcalDV(pp,qq) = conjg(dcalDV(qq,pp))
!     end do
!  end do
!  
!  return
!end subroutine calc_derivatives_BOwithEx



!------------------------------------------------------------------------------------------
!
!  subroutines below are not use now.
!
!------------------------------------------------------------------------------------------

!!========================================================================================================
!subroutine calc_derivatives_NonBO(t,calE,calEV,calDV,calC,calCV,calBV, & 
!     & TM_Qmat,Tnuc_mat,twoele_Qmat,nucele_Qmat,twonuc_mat,calFj0_Qmat,alpha, &
!     & dcalE,dcalEV,dcalDV,dcalC,dcalCV,dcalBV)
!!========================================================================================================
!  Use DiracOutput
!  use Constants
!  use NucBasis
!  implicit none
!
!  real(kind=8),intent(in) :: t  ! time
!  complex(kind=8),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
!  complex(kind=8),intent(in) :: calEV(4*NBS,4*NBS)  ! calEV_PQ
!  complex(kind=8),intent(in) :: calDV(4*NBS,4*NBS)  ! calDV_PQ
!  complex(kind=8),intent(in) :: calC(NBS_N,NBS_N)  ! calC_ij
!  complex(kind=8),intent(in) :: calCV(NBS_N,NBS_N)  ! calC_ij
!  complex(kind=8),intent(in) :: calBV(NBS_N,NBS_N)  ! calC_ij
!  complex(kind=8),intent(in) :: TM_Qmat(4*NBS,4*NBS)  ! T_PQ + M_PQ
!  complex(kind=8),intent(in) :: Tnuc_mat(NBS_N,NBS_N)  ! T_Nij
!  complex(kind=8),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)
!  complex(kind=8),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
!  complex(kind=8),intent(in) :: alpha(Nph)   ! <coherent|hat(a)|coherent>
!  complex(kind=8),intent(in) :: nucele_Qmat(NBS_N,NBS_N,4*NBS,4*NBS)  ! (ij|PQ)
!  complex(kind=8),intent(in) :: twonuc_mat(NBS_N,NBS_N,NBS_N,NBS_N)  ! (ij|kl)
!
!  complex(kind=8),intent(out) :: dcalE(4*NBS,4*NBS)  ! d(calE_PQ)/dt *DeltaT
!  complex(kind=8),intent(out) :: dcalEV(4*NBS,4*NBS)  ! d(calEV_PQ)/dt *DeltaT
!  complex(kind=8),intent(out) :: dcalDV(4*NBS,4*NBS)  ! d(calDV_PQ)/dt *DeltaT
!  complex(kind=8),intent(out) :: dcalC(NBS_N,NBS_N)  ! d(calC_ij)/dt *DeltaT
!  complex(kind=8),intent(out) :: dcalCV(NBS_N,NBS_N)  ! d(calCV_ij)/dt *DeltaT
!  complex(kind=8),intent(out) :: dcalBV(NBS_N,NBS_N)  ! d(calBV_ij)/dt *DeltaT
!
!
!  real(kind=8) :: dpj ! tilde(Delta p)_j
!  real(kind=8) :: p0j ! p_j
!
!  integer :: i,j,k,l,p,q,r,s
!  complex(kind=8) :: sum_te,sum_ne,sum_tn
!  complex(kind=8) :: sumE,sumEV,sumDV,sumC,sumCV,sumBV
!  complex(kind=8) :: dcalF
!
!  complex(kind=8) :: I2_Qmat(4*NBS,4*NBS), I4_Qmat(4*NBS,4*NBS) 
!  complex(kind=8) :: I_Qmat(4*NBS,4*NBS)  ! I_PQ = (I_1+I_2+I_3+I_4)_PQ = (TM +I2 +I4)_PQ
!  complex(kind=8) :: I4_mat(NBS_N,NBS_N) 
!  complex(kind=8) :: I_mat(NBS_N,NBS_N)  ! I_ij = (I_1+I_2+I_3+I_4)_ij = (T +I4)_ij
!
!!  complex(kind=8) :: I4V_Qmat(4*NBS,4*NBS),IV_Qmat(4*NBS,4*NBS)    ! for VEVs
!!  complex(kind=8) :: I4V_mat(NBS_N,NBS_N),IV_mat(NBS_N,NBS_N)
!!  complex(kind=8) :: sum_teV,sum_ne,sum_tnV
!
!
!  I2_Qmat  = (0.d0,0.d0)
!  I4_Qmat  = (0.d0,0.d0)
!
!  ! calculation of I4
!  do p=1,4*NBS
!     do q=1,4*NBS
!
!        sum_te = (0.d0,0.d0)
!        do r=1,4*NBS
!           do s=1,4*NBS
!              sum_te = sum_te + twoele_Qmat(p,q,r,s)*calE(r,s)
!           end do
!        end do
!        sum_ne = (0.d0,0.d0)
!        do i=1,NBS_N
!           do j=1,NBS_N
!              sum_ne = sum_ne + nucele_Qmat(i,j,p,q)*calC(i,j)  ! (PQ|ij)=(ij|PQ)
!           end do
!        end do
!        I4_Qmat(p,q) = sum_te + sum_ne
!
!     end do
!  end do
!
!!!$  ! calculation of I2 (Arad only)
!!!$  do p=1,4*NBS
!!!$     do q=1,4*NBS
!!!$
!!!$        sum_calF=(0.d0,0.d0)
!!!$        do j=1,Nph
!!!$           call calc_DeltaPj_and_P0j(j,dpj,p0j)
!!!$           
!!!$           dcalF = calFj0_Qmat(nn,mm,j)       *exp(-IU*CCC*p0j*t)*alpha(j) &
!!!$                & +conjg(calFj0_Qmat(mm,nn,j))*exp(+IU*CCC*p0j*t)*conjg(alpha(j))
!!!$           
!!!$           sum_calF = sum_calF -dpj*dcalF
!!!$
!!!$        end do
!!!$        
!!!$        I2_Qmat(p,q) = sum_calF
!!!$
!!!$     end do
!!!$  end do
!
!  I_Qmat(:,:) = TM_Qmat(:,:) + I2_Qmat(:,:) + I4_Qmat(:,:)
!
!  ! calculation of derivatives(*DeltaT)
!  do p=1,4*NBS
!     do q=p,4*NBS
!        sumE = (0.d0,0.d0);sumEV = (0.d0,0.d0);sumDV = (0.d0,0.d0)
!        do r=1,4*NBS
!           sumE  = sumE  -I_Qmat(r,p)*calE(r,q)  +I_Qmat(q,r)*calE(p,r)
!           sumEV = sumEV -I_Qmat(r,p)*calEV(r,q) +I_Qmat(q,r)*calEV(p,r)
!           sumDV = sumDV +I_Qmat(p,r)*calDV(r,q) -I_Qmat(r,q)*calDV(p,r)
!        end do
!        dcalE(p,q)  = -IU*sumE *DeltaT
!        dcalEV(p,q) = -IU*sumEV*DeltaT
!        dcalDV(p,q) = -IU*sumDV*DeltaT
!     end do
!  end do
!  
!  do p=2,4*NBS
!     do q=1,p-1
!        dcalE(p,q)  = conjg(dcalE(q,p))
!        dcalEV(p,q) = conjg(dcalEV(q,p))
!        dcalDV(p,q) = conjg(dcalDV(q,p))
!     end do
!  end do
!  
!  !--------------------------------------------------------------------------
!
!  I4_mat  = (0.d0,0.d0)
!
!  ! calculation of I4
!  do i=1,NBS_N
!     do j=1,NBS_N
!
!        sum_ne = (0.d0,0.d0)
!        do p=1,4*NBS
!           do q=1,4*NBS
!              sum_ne = sum_ne + nucele_Qmat(i,j,p,q)*calE(p,q)
!           end do
!        end do
!        sum_tn = (0.d0,0.d0)
!        do k=1,NBS_N
!           do l=1,NBS_N
!              sum_tn = sum_tn + twonuc_mat(i,j,k,l)*calC(k,l)  
!           end do
!        end do
!        I4_mat(i,j) = sum_tn + sum_ne
!
!     end do
!  end do
!
!  I_mat(:,:) = Tnuc_mat(:,:) + I4_mat(:,:)
!
!  do i=1,NBS_N
!     do j=i,NBS_N
!        sumC = (0.d0,0.d0); sumCV = (0.d0,0.d0); sumBV = (0.d0,0.d0);
!        do k=1,NBS_N
!           sumC  = sumC  -I_mat(k,i)*calC(k,j)  +I_mat(j,k)*calC(i,k)
!           sumCV = sumCV -I_mat(k,i)*calCV(k,j) +I_mat(j,k)*calCV(i,k)
!           sumBV = sumBV +I_mat(i,k)*calBV(k,j) -I_mat(k,j)*calBV(i,k)
!        end do
!        dcalC(i,j)  = -IU*sumC *DeltaT
!        dcalCV(i,j) = -IU*sumCV*DeltaT
!        dcalBV(i,j) = -IU*sumBV*DeltaT
!     end do
!  end do
!  
!  do i=2,NBS_N
!     do j=1,i-1
!        dcalC(i,j)  = conjg(dcalC(j,i))
!        dcalCV(i,j) = conjg(dcalCV(j,i))
!        dcalBV(i,j) = conjg(dcalBV(j,i))
!     end do
!  end do
!
!
!
!  
!  return
!end subroutine calc_derivatives_NonBO






!!==============================================================================================
!subroutine calc_derivatives_BO(t,calE,calEV,calDV,h_Qmat,twoele_Qmat,calFj0_Qmat,alpha, &
!     & dcalE,dcalEV,dcalDV,jAM_Qmat)
!! 
!!  It seems to be incorrect calculation for VEV. (KI, 120711)
!!==============================================================================================
!  use DiracOutput
!  use Constants
!  implicit none
!
!  real(kind=8),intent(in) :: t  ! time
!  complex(kind=8),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
!  complex(kind=8),intent(in) :: calEV(4*NBS,4*NBS)  ! calEV_PQ
!  complex(kind=8),intent(in) :: calDV(4*NBS,4*NBS)  ! calDV_PQ
!  complex(kind=8),intent(in) :: h_Qmat(4*NBS,4*NBS)  ! h_PQ
!  complex(kind=8),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)
!  complex(kind=8),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
!  complex(kind=8),intent(in) :: alpha(Nph)   ! <coherent|hat(a)|coherent>
!  complex(kind=dp),intent(in) :: jAM_Qmat(4*NBS,4*NBS)  ! jAM_PQ
!
!  complex(kind=8),intent(out) :: dcalE(4*NBS,4*NBS)  ! d(calE_PQ)/dt *DeltaT
!  complex(kind=8),intent(out) :: dcalEV(4*NBS,4*NBS)  ! d(calEV_PQ)/dt *DeltaT
!  complex(kind=8),intent(out) :: dcalDV(4*NBS,4*NBS)  ! d(calDV_PQ)/dt *DeltaT
!
!  real(kind=8) :: dpj ! tilde(Delta p)_j
!  real(kind=8) :: p0j ! p_j
!
!  integer :: i,j,k,l
!  integer :: pp,qq,rr,ss,nn,mm
!  complex(kind=8) :: sum_te,sum_calF,sum,sumEV,sumDV
!  complex(kind=8) :: dcalF
!
!  complex(kind=8) :: I2_Qmat(4*NBS,4*NBS), I4_Qmat(4*NBS,4*NBS) 
!  complex(kind=8) :: I_Qmat(4*NBS,4*NBS)  ! I_PQ = (I_1+I_2+I_3+I_4)_PQ = (TM +I2 +I4)_PQ
!
!
!  I2_Qmat  = (0.d0,0.d0)
!  I4_Qmat  = (0.d0,0.d0)
!
!  sum_calF=(0.d0,0.d0); sum_te=(0.d0,0.d0)
!
!  ! calculation of I4 (h_PQ + (PQ|RS)*calE_RS )
!  do nn=1,4*NBS
!     do mm=1,4*NBS
!
!        sum_te = (0.d0,0.d0)
! !notwoelecalc
!        do rr=1,4*NBS
!           do ss=1,4*NBS
!              sum_te = sum_te + twoele_Qmat(nn,mm,rr,ss)*calE(rr,ss)
!           end do
!        end do
!
!        I4_Qmat(nn,mm) = sum_te + h_Qmat(nn,mm)
!
!     end do
!  end do
!
!  ! calculation of I2 (Arad only)
!  do nn=1,4*NBS
!     do mm=1,4*NBS
!
!        sum_calF=(0.d0,0.d0)
!!!!        do j=1,Nph
!!!!           call calc_DeltaPj_and_P0j(j,dpj,p0j)
!!!!           
!!!!           dcalF = calFj0_Qmat(nn,mm,j)       *exp(-IU*CCC*p0j*t)*alpha(j) &
!!!!                & +conjg(calFj0_Qmat(mm,nn,j))*exp(+IU*CCC*p0j*t)*conjg(alpha(j))
!!!!           
!!!!           sum_calF = sum_calF -dpj*dcalF
!!!!
!!!!        end do
!        
!!!! 
!!        I2_Qmat(nn,mm) = sum_calF
!        I2_Qmat(nn,mm) = sum_calF + jAM_Qmat(nn,mm)
!!!!
!
!     end do
!  end do
!
!
!  I_Qmat(:,:) = I2_Qmat(:,:) + I4_Qmat(:,:)
!
!
!  ! calculation of derivatives(*DeltaT)
!  do pp=1,4*NBS
!     do qq=pp,4*NBS
!        sum = (0.d0,0.d0);sumEV = (0.d0,0.d0);sumDV = (0.d0,0.d0)
!        do rr=1,4*NBS
!           sum   = sum   -I_Qmat(rr,pp)*calE(rr,qq)  +I_Qmat(qq,rr)*calE(pp,rr)
!           sumEV = sumEV -I_Qmat(rr,pp)*calEV(rr,qq) +I_Qmat(qq,rr)*calEV(pp,rr)
!           sumDV = sumDV +I_Qmat(pp,rr)*calDV(rr,qq) -I_Qmat(rr,qq)*calDV(pp,rr)
!        end do
!        dcalE(pp,qq)  = -IU*sum*DeltaT
!        dcalEV(pp,qq) = -IU*sumEV*DeltaT
!        dcalDV(pp,qq) = -IU*sumDV*DeltaT
!     end do
!  end do
!  
!  do pp=2,4*NBS
!     do qq=1,pp-1
!        dcalE(pp,qq)  = conjg(dcalE(qq,pp))
!        dcalEV(pp,qq) = conjg(dcalEV(qq,pp))
!        dcalDV(pp,qq) = conjg(dcalDV(qq,pp))
!     end do
!  end do
!
!  return
!end subroutine calc_derivatives_BO


!==============================================================================================
subroutine calc_derivatives_BO_withAMnuc(t,calE,calEV,calDV,h_Qmat,twoele_Qmat,calFj0_Qmat,alpha, &
     & dcalE,dcalEV,dcalDV,jAM_Qmat,jAMnuc_Qmat)
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
  complex(kind=dp),intent(in) :: jAM_Qmat(4*NBS,4*NBS)  ! jAM_PQ
  complex(kind=dp),intent(in) :: jAMnuc_Qmat(4*NBS,4*NBS)  ! jAM_PQ

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

!  write(*,*)'call subroutine calc_derivatives_BOwithAMnuc'

  I2_Qmat  = (0.d0,0.d0)
  I4_Qmat  = (0.d0,0.d0)

  sum_calF=(0.d0,0.d0); sum_te=(0.d0,0.d0)

  ! calculation of I4 (h_PQ + (PQ|RS)*calE_RS )
  do nn=1,4*NBS
     do mm=1,4*NBS

        sum_te = (0.d0,0.d0)
 !notwoelecalc
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
!!!        do j=1,Nph
!!!           call calc_DeltaPj_and_P0j(j,dpj,p0j)
!!!           
!!!           dcalF = calFj0_Qmat(nn,mm,j)       *exp(-IU*CCC*p0j*t)*alpha(j) &
!!!                & +conjg(calFj0_Qmat(mm,nn,j))*exp(+IU*CCC*p0j*t)*conjg(alpha(j))
!!!           
!!!           sum_calF = sum_calF -dpj*dcalF
!!!
!!!        end do
        
!!! 
!        I2_Qmat(nn,mm) = sum_calF
        I2_Qmat(nn,mm) = sum_calF + jAM_Qmat(nn,mm) + jAMnuc_Qmat(nn,mm)
!!!

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

!  write(*,*)'stop subroutine calc_derivatives_BOwithAMnuc'
!  stop
  return
end subroutine calc_derivatives_BO_withAMnuc

!==============================================================================================
subroutine calc_derivatives_BO_withA0M(t,calE,calEV,calDV,h_Qmat,twoele_Qmat,calFj0_Qmat,alpha, &
     & dcalE,dcalEV,dcalDV,intrhoA0M_Qmat)
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
  complex(kind=dp),intent(in) :: intrhoA0M_Qmat(4*NBS,4*NBS)  !

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
!!$ !notwoelecalc
!!$        do rr=1,4*NBS
!!$           do ss=1,4*NBS
!!$              sum_te = sum_te + twoele_Qmat(nn,mm,rr,ss)*calE(rr,ss)
!!$           end do
!!$        end do

!        I4_Qmat(nn,mm) = sum_te + h_Qmat(nn,mm)
        I4_Qmat(nn,mm) = sum_te + h_Qmat(nn,mm) + intrhoA0M_Qmat(nn,mm)

     end do
  end do

!!$  ! calculation of I2 (Arad only)
!!$  do nn=1,4*NBS
!!$     do mm=1,4*NBS
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
!!$!!! vecAM
!!$        I2_Qmat(nn,mm) = sum_calF
!!$!        I2_Qmat(nn,mm) = sum_calF + jAM_Qmat(nn,mm)
!!$!!!
!!$
!!$     end do
!!$  end do


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
end subroutine calc_derivatives_BO_withA0M

!==============================================================================================
subroutine calc_derivatives_BO_withA0Msin(t,calE,calEV,calDV,h_Qmat,twoele_Qmat,calFj0_Qmat,alpha, &
     & dcalE,dcalEV,dcalDV,intrhoA0M_Qmat)
!==============================================================================================
  use DiracOutput
  use Constants
  use param_AM
  implicit none

  real(kind=8),intent(in) :: t  ! time
  complex(kind=8),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  complex(kind=8),intent(in) :: calEV(4*NBS,4*NBS)  ! calEV_PQ
  complex(kind=8),intent(in) :: calDV(4*NBS,4*NBS)  ! calDV_PQ
  complex(kind=8),intent(in) :: h_Qmat(4*NBS,4*NBS)  ! h_PQ
  complex(kind=8),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)
  complex(kind=8),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
  complex(kind=8),intent(in) :: alpha(Nph)   ! <coherent|hat(a)|coherent>
  complex(kind=dp),intent(in) :: intrhoA0M_Qmat(4*NBS,4*NBS)  !

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
!!$ !notwoelecalc
!!$        do rr=1,4*NBS
!!$           do ss=1,4*NBS
!!$              sum_te = sum_te + twoele_Qmat(nn,mm,rr,ss)*calE(rr,ss)
!!$           end do
!!$        end do

!        I4_Qmat(nn,mm) = sum_te + h_Qmat(nn,mm)
        I4_Qmat(nn,mm) = sum_te + h_Qmat(nn,mm) + intrhoA0M_Qmat(nn,mm)*sin(omega_A0M*t)

     end do
  end do

!!$  ! calculation of I2 (Arad only)
!!$  do nn=1,4*NBS
!!$     do mm=1,4*NBS
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
!!$!!! vecAM
!!$        I2_Qmat(nn,mm) = sum_calF
!!$!        I2_Qmat(nn,mm) = sum_calF + jAM_Qmat(nn,mm)
!!$!!!
!!$
!!$     end do
!!$  end do


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
end subroutine calc_derivatives_BO_withA0Msin


!==============================================================================================
subroutine calc_derivatives_BOwithEx_withA0M(t,calE,calEV,calDV,h_Qmat,twoele_Qmat,calFj0_Qmat,alpha, &
     & dcalE,dcalEV,dcalDV,intrhoA0M_Qmat)
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
  complex(kind=dp),intent(in) :: intrhoA0M_Qmat(4*NBS,4*NBS)  !

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
 !notwoelecalc
!!$        do rr=1,4*NBS
!!$           do ss=1,4*NBS
!!$              sum_te = sum_te + twoele_Qmat(nn,mm,rr,ss)*calE(rr,ss)
!!$           end do
!!$        end do

!        I4_Qmat(nn,mm) = sum_te + h_Qmat(nn,mm)
        I4_Qmat(nn,mm) = sum_te + h_Qmat(nn,mm) + intrhoA0M_Qmat(nn,mm)

     end do
  end do

!!$  ! calculation of I2 (Arad only)
!!$  do nn=1,4*NBS
!!$     do mm=1,4*NBS
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
!!$!!! vecAM
!!$        I2_Qmat(nn,mm) = sum_calF
!!$!        I2_Qmat(nn,mm) = sum_calF + jAM_Qmat(nn,mm)
!!$!!!
!!$
!!$     end do
!!$  end do


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
end subroutine calc_derivatives_BOwithEx_withA0M

















!!====================================================================================
!!subroutine calc_derivative(nn,mm,it,calE,h_Qmat,calFj0_Qmat,alpha,dcalEdt)
!subroutine calc_derivative(nn,mm,it,calE,h_Qmat,twoele_Qmat,calFj0_Qmat,alpha,dcalEdt)
!!====================================================================================
!  use DiracOutput
!  use Constants
!  implicit none
!
!  integer,intent(in) :: nn,mm ! compute (nn,mm) component of d(calE_NM)/dt
!  integer,intent(in) :: it  ! it-th time step
!  complex(kind=8),intent(in) :: calE(4*NBS,4*NBS)  ! calE_NM
!  complex(kind=8),intent(in) :: h_Qmat(4*NBS,4*NBS)  ! h_NM
!  complex(kind=8),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (NM|PQ)
!  complex(kind=8),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_NMj(t=0)
!  complex(kind=8),intent(in) :: alpha(Nph)   ! <coherent|hat(a)|coherent>
!
!  complex(kind=8),intent(out) :: dcalEdt  ! d(calE_NM)/dt
!
!  real(kind=8) :: dpj ! tilde(Delta p)_j
!  real(kind=8) :: p0j ! p_j
!  real(kind=8) :: t ! time = it*DeltaT
!
!  integer :: i,j,k,l
!  integer :: rr,pp,qq
!  complex(kind=8) :: sum_h,sum_te,sum_calF,sum_dcalF
!  complex(kind=8) :: dcalF
!
!  sum_h =(0.d0,0.d0); sum_calF=(0.d0,0.d0); sum_te=(0.d0,0.d0)
!
!  ! it starts from 0 -> t starts from 0
!  t = DeltaT*it  ! DeltaT is defined in GammaMatrix and set in the main routine.
!  
!
!  ! h^*_NR *calE_RM - h_MR *calE_NR
!  do rr=1,4*NBS
!     sum_h = sum_h +conjg(h_Qmat(nn,rr))*calE(rr,mm) -h_Qmat(mm,rr)*calE(nn,rr)
!!     write(*,"(3i6,6es16.6)") nn,mm,rr,h_Qmat(mm,rr),calE(nn,rr),sum_h
!  end do
!
!  ! two electron integral 
!  do rr=1,4*NBS
!     do pp=1,4*NBS
!        do qq=1,4*NBS
!           
!           sum_te = sum_te &
!                &  + conjg(twoele_Qmat(nn,rr,pp,qq))*calE(rr,mm)*calE(qq,pp) &
!                &  - twoele_Qmat(mm,rr,pp,qq)       *calE(pp,qq)*calE(nn,rr)
!           
!        end do
!     end do
!  end do
!
!!  write(*,"(2i6,4es16.6)") nn,mm,sum_h,sum_te
!
!  
!  ! Arad part
!  do j=1,Nph
!     call calc_DeltaPj_and_P0j(j,dpj,p0j)
!     sum_dcalF=(0.d0,0.d0)
!     do rr=1,4*NBS
!        dcalF = conjg(calFj0_Qmat(nn,rr,j))*exp(+IU*CCC*p0j*t)*calE(rr,mm)*conjg(alpha(j)) &
!             & +calFj0_Qmat(rr,nn,j)       *exp(-IU*CCC*p0j*t)*calE(rr,mm)*alpha(j) &
!             & -calFj0_Qmat(mm,rr,j)       *exp(-IU*CCC*p0j*t)*calE(nn,rr)*alpha(j) &
!             & -conjg(calFj0_Qmat(rr,mm,j))*exp(+IU*CCC*p0j*t)*calE(nn,rr)*conjg(alpha(j)) 
!        sum_dcalF = sum_dcalF +dcalF
!     end do
!!     write(*,"(3i6,2es16.6,4es16.6)") nn,mm,j,dpj,p0j,calFj0_Qmat(nn,rr,j),sum_calF
!!     write(*,"(3i6,2es16.6)") nn,mm,j,sum_calF
!     sum_calF = sum_calF -dpj*sum_dcalF
!!     write(*,"(3i6,4es16.6)") nn,mm,j,sum_dcalF,sum_calF
!  end do
!!  write(*,"(2i6,2es16.6)") nn,mm,sum_calF
!  
!!  stop
!  dcalEdt = IU*(sum_h+sum_te+sum_calF)
!!  dcalEdt = IU*(sum_h+sum_calF)
!!  dcalEdt = IU*(sum_h)
!!  dcalEdt = IU*(sum_calF)
!!  write(*,"(2i6,4es16.6)") nn,mm,sum_calF,dcalEdt
!!  dcalEdt = IU*(sum_h)
!  
!  return
!end subroutine calc_derivative



