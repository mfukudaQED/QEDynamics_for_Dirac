!=====================================================================
!
! subroutines and functions for calculating I_JE 
!
!=====================================================================

!!$!=====================================================================================================
!!$subroutine gauss_int_thetaJE(alpha,posA,alphaA,n1,posB,alphaB,nbar1,posC,alphaC,n2,posD,alphaD,nbar2,thetaJE)
!!$! calculate integral for primitive gaussians (AB|theta_jE|CD)
!!$! 
!!$! alpha=0  -> 0
!!$! alpha/=0 -> use formula for [NLM|theta_jE|N'L'M']
!!$!
!!$! 140109
!!$!
!!$!=====================================================================================================
!!$  use Precision
!!$  use Constants  ! for IU and CCC
!!$  implicit none
!!$  
!!$  real(kind=dp),intent(in) :: alpha
!!$  real(kind=8),intent(in) :: posA(3),posB(3),posC(3),posD(3)
!!$  real(kind=8),intent(in) :: alphaA,alphaB,alphaC,alphaD
!!$  integer,intent(in) :: n1(3),nbar1(3),n2(3),nbar2(3)  
!!$  real(kind=8),intent(out) :: thetaJE(3) ! k=1,2,3
!!$
!!$  real(kind=8) :: PI
!!$  integer :: i,j,k
!!$  integer :: i1,j1,k1,i2,j2,k2
!!$  
!!$  real(kind=8), allocatable, dimension(:,:,:) :: d1,e1,f1,d2,e2,f2
!!$  real(kind=8) :: d1ijk,d2ijk
!!$
!!$  real(kind=8) :: posP(3),alphaP,c_PAB,vecPA(3),vecPB(3) ! mixing of A and B
!!$  real(kind=8) :: posQ(3),alphaQ,c_QCD,vecQC(3),vecQD(3) ! mixing of C and D
!!$!  real(kind=8) :: lambda,vecPQ(3),alphaT,T  ! R_NLM for two electron integral
!!$
!!$  real(kind=dp) :: vecPQ(3)  ! vec{D}
!!$  complex(kind=dp) :: A,B,alphaT,T
!!$  complex(kind=dp), allocatable, dimension(:,:,:,:,:,:) :: NLMthetaJENLM 
!!$  
!!$  real(kind=8), allocatable, dimension(:) :: list_R000j
!!$!  real(kind=8), allocatable, dimension(:,:,:) :: R
!!$  real(kind=8), allocatable, dimension(:,:,:) :: Rj1  ! R_{NLM1}
!!$  integer :: n1_sum(3)  ! upper bounds of N,L,M
!!$  integer :: n2_sum(3)  ! upper bounds of N',L',M'
!!$  integer :: n_sum(3) ! upper bounds of N+N',L+L',M+M'
!!$  integer :: j_max ! upper bound for j when construction R_NLMj
!!$  real(kind=8) :: integral, twoele_tmp
!!$
!!$  if(alpha.eq.0._dp) then  
!!$     thetaJE(:) = 0._dp
!!$  else
!!$
!!$     PI = atan(1.d0)*4.d0
!!$     
!!$     do i=1,3
!!$        n1_sum(i) = n1(i)+nbar1(i)
!!$        n2_sum(i) = n2(i)+nbar2(i)
!!$        n_sum(i)  = n1_sum(i)+n2_sum(i)
!!$     end do
!!$     
!!$     call calc_gaussP(posA,alphaA,posB,alphaB,c_PAB,posP,alphaP,vecPA,vecPB)  ! construct P 
!!$     call calc_gaussP(posC,alphaC,posD,alphaD,c_QCD,posQ,alphaQ,vecQC,vecQD)  ! construct Q
!!$     
!!$     ! construct d,e,f,d',e',f'
!!$     allocate(d1(0:n1(1),0:nbar1(1),0:n1_sum(1)))
!!$     allocate(e1(0:n1(2),0:nbar1(2),0:n1_sum(2)))
!!$     allocate(f1(0:n1(3),0:nbar1(3),0:n1_sum(3)))
!!$     allocate(d2(0:n2(1),0:nbar2(1),0:n2_sum(1)))
!!$     allocate(e2(0:n2(2),0:nbar2(2),0:n2_sum(2)))
!!$     allocate(f2(0:n2(3),0:nbar2(3),0:n2_sum(3)))
!!$     
!!$     call calc_d(alphaP,vecPA(1),vecPB(1),n1(1),nbar1(1),d1)
!!$     call calc_d(alphaP,vecPA(2),vecPB(2),n1(2),nbar1(2),e1)
!!$     call calc_d(alphaP,vecPA(3),vecPB(3),n1(3),nbar1(3),f1)
!!$     call calc_d(alphaQ,vecQC(1),vecQD(1),n2(1),nbar2(1),d2)
!!$     call calc_d(alphaQ,vecQC(2),vecQD(2),n2(2),nbar2(2),e2)
!!$     call calc_d(alphaQ,vecQC(3),vecQD(3),n2(3),nbar2(3),f2)
!!$
!!$     ! compute [NLM|theta_jE|N'L'M']
!!$
!!$     allocate( NLMthetaJENLM(0:n1_sum(1),0:n1_sum(2),0:n1_sum(3),0:n2_sum(1),0:n2_sum(2),0:n2_sum(3)) ) 
!!$
!!$     !---copied from gauss_int_thetaJJ---
!!$     A = IU*alpha/CCC**2
!!$     B = A*(alphaP+alphaQ) +alphaP*alphaQ
!!$     do i=1,3
!!$        vecPQ(i) = posP(i)-posQ(i)  ! vecD
!!$     end do
!!$     alphaT = 1._dp/(1._dp/alphaP +1._dp/alphaQ +1._dp/A)
!!$     !---copy end---
!!$
!!$     ! construct R_{N+N',L+L',M+M'}  (3.34)
!!$!     lambda = 2.d0*PI**(5.d0/2.d0)/(alphaP*alphaQ*sqrt(alphaP+alphaQ))  !M&D eq.(3.31)
!!$     T = 0.d0
!!$     do i=1,3  ! |vec{D}|^2
!!$        vecPQ(i) = posP(i)-posQ(i)
!!$        T = T +vecPQ(i)**2
!!$     end do 
!!$!     alphaT = alphaP*alphaQ/(alphaP+alphaQ)
!!$     T = alphaT*T  ! M&D eq.(3.32)   ( T = alpha_T |vec{D}|^2 )
!!$     
!!$     allocate(R(0:n_sum(1),0:n_sum(2),0:n_sum(3))) ! for two-electron integral  (need R_{N+N',L,M} , R_{N,L+L',M} and R_{N,L,M+M'})
!!$     j_max = n_sum(1)+n_sum(2)+n_sum(3) ! for two-electron integral
!!$     allocate(list_R000j(0:j_max))
!!$     call calc_R000j(alphaT,T,j_max,list_R000j)
!!$     call calc_R(vecPQ,j_max,list_R000j,n_sum(1),n_sum(2),n_sum(3),R) 
!!$     
!!$     twoele = 0.d0
!!$     do i1=0,n1_sum(1)  ! N
!!$        do j1=0,n1_sum(2)  ! L
!!$           do k1=0,n1_sum(3)  ! M
!!$              d1ijk = d1(n1(1),nbar1(1),i1) * e1(n1(2),nbar1(2),j1) * f1(n1(3),nbar1(3),k1)
!!$              twoele_tmp = 0.d0 
!!$              do i2=0,n2_sum(1)  ! N'
!!$                 do j2=0,n2_sum(2) ! L'
!!$                    do k2=0,n2_sum(3) ! M'
!!$                       integral = (-1.d0)**(i2+j2+k2)*R(i1+i2,j1+j2,k1+k2)  ! eq.(3.34) (-1)^{N'+L'+M'}*R_{N+N',L+L',M+M'} 
!!$                       d2ijk = d2(n2(1),nbar2(1),i2) * e2(n2(2),nbar2(2),j2) * f2(n2(3),nbar2(3),k2)
!!$                       twoele_tmp = twoele_tmp +d2ijk*integral
!!$                    end do
!!$                 end do
!!$              end do
!!$              twoele = twoele + twoele_tmp * d1ijk
!!$           end do
!!$        end do
!!$     end do
!!$     twoele = c_PAB*c_QCD*twoele * lambda
!!$   
!!$  end if
!!$end subroutine gauss_int_thetaJE



!====================================================================================
subroutine calc_Rj1(vecPC,n_sum,list_R000j,nx_sum,ny_sum,nz_sum,table_R)
! calculate R_NMLj with j=1
! only depend on nx+nbarx, ny+nbary, nz+nbarz 
!
! the case of "nx_sum = ny_sum = nz_sum = 0" (4 gaussians are all s-type) is not computed.
! 
! 140109
!====================================================================================
  implicit none
  real(kind=8),intent(in) :: vecPC(3)
  integer,intent(in) :: n_sum,nx_sum,ny_sum,nz_sum
  real(kind=8),intent(in) :: list_R000j(0:n_sum)
  real(kind=8),intent(out) :: table_R(0:nx_sum,0:ny_sum,0:nz_sum)
  integer :: i,j,k,l
  real(kind=8) :: table_Rj(0:nx_sum,0:ny_sum,0:nz_sum,0:n_sum)
  real(kind=8) :: R(0:nx_sum,0:ny_sum,0:nz_sum)
  real(kind=8) :: a,b,c

  a = vecPC(1); b = vecPC(2); c = vecPC(3)
  table_Rj = 0.d0

  ! N=0,L=0,M=0  (Eq.(4.4))
  do i=0,n_sum
     table_Rj(0,0,0,i) = list_R000j(i)
  end do
  
  ! use loop variable i for N, j for L, k for M, l for j

  ! Generate nonzero M for N=0,L=0 using Eq.(4.6)
  if(nz_sum.ge.1) then  ! if nz=0,nbarz=0, we do not need to calculate nonzero M
     do k=0,nz_sum-1 ! loop for M
        do l=0,n_sum-(k+1) ! loop for j
           if(k.eq.0) then ! when M=0, there is no M=-1
              table_Rj(0,0,k+1,l) = c*table_Rj(0,0,k,l+1)
           else
              table_Rj(0,0,k+1,l) = c*table_Rj(0,0,k,l+1) +k*table_Rj(0,0,k-1,l+1)
           end if
        end do
     end do
  end if

  ! Generate nonzero L for N=0 and all M using Eq.(4.7)
  do k=0,nz_sum ! loop for all M
     if(ny_sum.ge.1) then ! if ny=0,nbary=0, we do not need to calculate nonzero L
        do j=0,ny_sum-1 ! loop for L
           do l=0,n_sum-(j+1+k) ! loop for j
              if(j.eq.0) then ! L=0
                 table_Rj(0,j+1,k,l) = b*table_Rj(0,j,k,l+1) 
              else
                 table_Rj(0,j+1,k,l) = b*table_Rj(0,j,k,l+1) +j*table_Rj(0,j-1,k,l+1)
              end if
           end do
        end do
     end if
  end do
  
  ! Generate nonzero N for all L&M using Eq.(4.8)
  do k=0,nz_sum ! loop for all M
     do j=0,ny_sum ! loop for all L
        if(nx_sum.ge.1) then ! if nx=0,nbarx=0, we do not need to calculate nonzero N
           do i=0,nx_sum-1 ! loop for N
              do l=0,n_sum-(i+1+j+k) ! loop for j
                 if(i.eq.0) then
                    table_Rj(i+1,j,k,l) = a*table_Rj(i,j,k,l+1)
                 else
                    table_Rj(i+1,j,k,l) = a*table_Rj(i,j,k,l+1) +i*table_Rj(i-1,j,k,l+1)
                 end if
              end do
           end do
        end if
     end do
  end do

  ! R_NLM1
  do i=0,nx_sum ! loop for all N
     do j=0,ny_sum ! loop for all L
        do k=0,nz_sum ! loop for all M
!           R(i,j,k) = table_Rj(i,j,k,0)  ! R_NLM  (j=0 case needed for two-electron integral)
           R(i,j,k) = table_Rj(i,j,k,1) ! R_NLM1 (j=1 case needed for jE integral)
        end do
     end do
  end do

  table_R(:,:,:) = R(:,:,:)

!!$  do i=0,nx_sum ! loop for all N
!!$     do j=0,ny_sum ! loop for all L
!!$        do k=0,nz_sum ! loop for all M
!!$ !          write(*,'(3i5,2es16.6)') i,j,k,R(i,j,k),table_R(i,j,k)
!!$        end do
!!$     end do
!!$  end do

  return
end subroutine calc_Rj1






!=====================================================================
!
! subroutines and functions for calculating I_JJ 
!
!=====================================================================


!======================================================
subroutine setQmat_intKJJ_nz(intKJJ_Qmat_nz)
! K_{JJ, NMPQ} = int (d alpha) I_{JJ, NMPQ}(alpha) 
!
! I_{JJ, NMPQ}(alpha) = sum_{k=1}^3 int dr ds j^k_NM(r) j^k_PQ(s) exp(-I alpha (r-s)^2/c^2)
!
! 140108
!
!======================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use IntegralStorage
  implicit none
  
  type(nonzeromat4legs),allocatable,intent(out) :: intKJJ_Qmat_nz(:)
  
  real(kind=dp) :: alphaJJ
  complex(kind=dp) :: intIJJ(-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2) 
  complex(kind=dp) :: intKJJ(-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2) 
  
  integer :: nn,mm,pp,qq  ! index for 1~4*NBS
  integer :: n,m,p,q ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b,c,d ! + or -
  integer :: aa,bb,cc,dd ! 1 -> +, 2 -> -
  
  character(LEN=80) :: temp
  real(kind=dp) :: int_real,int_complex
  real(kind=dp) :: th
  character(LEN=300) :: readfile
  
  integer :: i
  type(nonzeromat4legs),allocatable :: nonzero_intKJJ(:) ! need large memory but only used when KJJ is computed.

  real(kind=dp) :: d_alphaJJ

  
  if(there_is_intKJJ) then
     readfile = trim(FILESFOLDER)//"/"//file_intKJJ
     write(*,*) " There is intKJJ file. Read integrals from "
     write(*,"(a)") trim(readfile)
     write(10,*) "# There is intKJJ file. Read integrals from "
     write(10,"(a)") "# "//trim(readfile)
     open(unit=204,file=readfile,status='unknown',form='formatted')

     read(204,*) AlphaJJ_MAX, N_alphaJJ
     read(204,*) Diff_u_t
     read(204,*) N_intKJJ, th
     write(*,*) "AlphaJJ_MAX: ", AlphaJJ_MAX
     write(10,*) "# AlphaJJ_MAX: ", AlphaJJ_MAX
     write(*,*) "N_alphaJJ: ", N_alphaJJ
     write(10,*) "# N_alphaJJ: ", N_alphaJJ
     write(*,*) "Number of intKJJ we use:  ", N_intKJJ
     write(10,*) "# Number of intKJJ we use:  ", N_intKJJ
     write(*,*) "threshold of intKJJ:  ", th
     write(10,*) "# threshold of intKJJ:  ", th
     allocate(intKJJ_Qmat_nz(N_intKJJ))
     do i=1,N_intKJJ
        read(204,*) intKJJ_Qmat_nz(i)%a,intKJJ_Qmat_nz(i)%b,intKJJ_Qmat_nz(i)%c, &
             & intKJJ_Qmat_nz(i)%d,int_real,int_complex 
        intKJJ_Qmat_nz(i)%val = cmplx(int_real,int_complex,dp)
     end do
     
     write(*,*) " Reading done."
     
  else ! if there isn't precomputed integrals, compute integrals.
     write(*,*) "There is no intKJJ file. Compute integrals and store at fort.15."
     write(10,*) "# There is no intKJJ file. Compute integrals and store at fort.15."
     write(*,*) "AlphaJJ_MAX: ", AlphaJJ_MAX
     write(10,*) "# AlphaJJ_MAX: ", AlphaJJ_MAX
     write(*,*) "N_alphaJJ: ", N_alphaJJ
     write(10,*) "# N_alphaJJ: ", N_alphaJJ
     write(*,*) "Diff_u_t: ", Diff_u_t
     write(10,*) "# Diff_u_t: ", Diff_u_t
!     stop
     !------------------------------------------------
     ! Compute IJJ(alpha) and integration over alpha
     !------------------------------------------------
     d_alphaJJ = AlphaJJ_MAX/N_alphaJJ
     intKJJ = (0._dp,0._dp)

     do i=1,N_alphaJJ
        alphaJJ = d_alphaJJ*(i-1)
        call calc_intIJJ_mat(alphaJJ,intIJJ) ! get every component of intIJJ (I_jj)
        do q = -NBS,NBS  !q=0 is not used.
           do dd = 1,2
              do p = -NBS,NBS  !p=0 is not used.
                 do cc = 1,2
                    do m = -NBS,NBS  !m=0 is not used.
                       do bb = 1,2
                          do n = -NBS,NBS  !n=0 is not used.
                             do aa = 1,2
                                
                                intKJJ(n,aa,m,bb,p,cc,q,dd) = intKJJ(n,aa,m,bb,p,cc,q,dd)  & 
                                     & + exp( IU*alphaJJ*Diff_u_t**2)      *intIJJ(n,aa,m,bb,p,cc,q,dd) &
                                     & + exp(-IU*alphaJJ*Diff_u_t**2)*conjg(intIJJ(m,bb,n,aa,q,dd,p,cc))
                                
                             end do
                          end do
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
     intKJJ = intKJJ * d_alphaJJ

     !--------------------------------------------------------------------
     ! store components of intKJJ which are larger than TH_intKJJ
     !--------------------------------------------------------------------
     allocate(nonzero_intKJJ((4*NBS)**4))
     i=0
     do nn=1,4*NBS
        do mm=1,4*NBS
           do pp=1,4*NBS
              do qq=1,4*NBS
                 
                 call index_from_Qmat(nn,n,a)
                 call index_from_Qmat(mm,m,b)
                 call index_from_Qmat(pp,p,c)
                 call index_from_Qmat(qq,q,d)

                 ! "+" --> 1, "-" --> 2 to be used in modified routines 
                 if(a.eq."+") aa = 1
                 if(a.eq."-") aa = 2
                 if(b.eq."+") bb = 1
                 if(b.eq."-") bb = 2
                 if(c.eq."+") cc = 1
                 if(c.eq."-") cc = 2
                 if(d.eq."+") dd = 1
                 if(d.eq."-") dd = 2

                 if(abs(intKJJ(n,aa,m,bb,p,cc,q,dd)).gt.TH_intKJJ) then
                    i = i+1
                    nonzero_intKJJ(i)%a = nn
                    nonzero_intKJJ(i)%b = mm
                    nonzero_intKJJ(i)%c = pp
                    nonzero_intKJJ(i)%d = qq
                    nonzero_intKJJ(i)%val = intKJJ(n,aa,m,bb,p,cc,q,dd)
                 end if

              end do
           end do
        end do
     end do
     write(*,*) " Integrals K_JJ done."
     write(10,*) "# Integrals K_JJ done."
     N_intKJJ = i

     write(*,*) "Number of intKJJ we use:  ", N_intKJJ
     write(10,*) "# Number of intKJJ we use:  ", N_intKJJ
     write(*,*) "TH_intKJJ:  ", TH_intKJJ
     write(10,*) "# TH_intKJJ:  ", TH_intKJJ

     write(15,"(1es16.6,1i16)") AlphaJJ_MAX, N_alphaJJ
     write(15,"(1es16.6)") Diff_u_t
     write(15,"(1i16,1es14.4)") N_intKJJ,TH_intKJJ
     allocate(intKJJ_Qmat_nz(N_intKJJ))
     do i=1,N_intKJJ
        write(15,"(4i6,2es16.6)") nonzero_intKJJ(i)%a,nonzero_intKJJ(i)%b,nonzero_intKJJ(i)%c,&
             & nonzero_intKJJ(i)%d,nonzero_intKJJ(i)%val
        intKJJ_Qmat_nz(i) = nonzero_intKJJ(i)
     end do

     deallocate(nonzero_intKJJ)
  end if
     
  return
end subroutine setQmat_intKJJ_nz



!======================================================
subroutine setQmat_intIJJ_nz(intIJJ_Qmat_nz)
! I_{JJ, NMPQ} = sum_{k=1}^3 int dr ds j^k_NM(r) j^k_PQ(s) exp(-I alpha (r-s)^2/c^2)
!
! 140108
!======================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use IntegralStorage
  implicit none
  
  type(nonzeromat4legs),allocatable,intent(out) :: intIJJ_Qmat_nz(:)

  real(kind=dp) :: alphaJJ
  complex(kind=dp) :: intIJJ(-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2) 
  
  integer :: nn,mm,pp,qq  ! index for 1~4*NBS
  integer :: n,m,p,q ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b,c,d ! + or -
  integer :: aa,bb,cc,dd ! 1 -> +, 2 -> -

  character(LEN=80) :: temp
  real(kind=dp) :: int_real,int_complex
  real(kind=dp) :: th
  character(LEN=300) :: readfile

  integer :: i
  type(nonzeromat4legs),allocatable :: nonzero_intIJJ(:) ! need large memory but only used when IJJ is computed.

  if(there_is_intIJJ) then
     readfile = trim(FILESFOLDER)//"/"//file_intIJJ
     write(*,*) " There is intIJJ file. Read integrals from "
     write(*,"(a)") trim(readfile)
     write(10,*) "# There is intIJJ file. Read integrals from "
     write(10,"(a)") "# "//trim(readfile)
     open(unit=204,file=readfile,status='unknown',form='formatted')

     read(204,*) alphaJJ,N_intIJJ, th
     write(*,*) "Number of intIJJ we use:  ", N_intIJJ
     write(10,*) "# Number of intIJJ we use:  ", N_intIJJ
     write(*,*) "threshold of intIJJ:  ", th
     write(10,*) "# threshold of intIJJ:  ", th
     allocate(intIJJ_Qmat_nz(N_intIJJ))
     do i=1,N_intIJJ
        read(204,*) intIJJ_Qmat_nz(i)%a,intIJJ_Qmat_nz(i)%b,intIJJ_Qmat_nz(i)%c, &
             & intIJJ_Qmat_nz(i)%d,int_real,int_complex 
        intIJJ_Qmat_nz(i)%val = cmplx(int_real,int_complex,dp)
     end do
     
     write(*,*) " Reading done."
     
  else ! if there isn't precomputed integrals, compute integrals.
     write(*,*) " There is no intIJJ file. Compute integrals and store at fort.15."
     
!     alphaJJ = 0._dp
     alphaJJ = 1._dp
     call calc_intIJJ_mat(alphaJJ,intIJJ) 

     allocate(nonzero_intIJJ((4*NBS)**4))
     i=0
     do nn=1,4*NBS
        do mm=1,4*NBS
           do pp=1,4*NBS
              do qq=1,4*NBS
                 
                 call index_from_Qmat(nn,n,a)
                 call index_from_Qmat(mm,m,b)
                 call index_from_Qmat(pp,p,c)
                 call index_from_Qmat(qq,q,d)

                 ! "+" --> 1, "-" --> 2 to be used in modified routines 
                 if(a.eq."+") aa = 1
                 if(a.eq."-") aa = 2
                 if(b.eq."+") bb = 1
                 if(b.eq."-") bb = 2
                 if(c.eq."+") cc = 1
                 if(c.eq."-") cc = 2
                 if(d.eq."+") dd = 1
                 if(d.eq."-") dd = 2

                 if(abs(intIJJ(n,aa,m,bb,p,cc,q,dd)).gt.TH_intIJJ) then
                    i = i+1
                    nonzero_intIJJ(i)%a = nn
                    nonzero_intIJJ(i)%b = mm
                    nonzero_intIJJ(i)%c = pp
                    nonzero_intIJJ(i)%d = qq
                    nonzero_intIJJ(i)%val = intIJJ(n,aa,m,bb,p,cc,q,dd)
                 end if

              end do
           end do
        end do
     end do
     write(*,*) " Integrals I_JJ done."
     N_intIJJ = i

     write(*,*) "Number of intIJJ we use:  ", N_intIJJ
     write(10,*) "# Number of intIJJ we use:  ", N_intIJJ
     write(*,*) "TH_intIJJ:  ", TH_intIJJ
     write(10,*) "# TH_intIJJ:  ", TH_intIJJ

     write(15,"(1es16.6,1i16,1es14.4)") alphaJJ,N_intIJJ,TH_intIJJ
!     write(15,"(1i16,1es14.4)") N_intIJJ,TH_intIJJ
     allocate(intIJJ_Qmat_nz(N_intIJJ))
     do i=1,N_intIJJ
        write(15,"(4i6,2es16.6)") nonzero_intIJJ(i)%a,nonzero_intIJJ(i)%b,nonzero_intIJJ(i)%c,&
             & nonzero_intIJJ(i)%d,nonzero_intIJJ(i)%val
        intIJJ_Qmat_nz(i) = nonzero_intIJJ(i)
     end do

     deallocate(nonzero_intIJJ)
  end if
     
  return
end subroutine setQmat_intIJJ_nz


!======================================================
!subroutine setQmat_intIJJ(intIJJ_Qmat)
subroutine setQmat_intIJJ(alphaJJ,intIJJ_Qmat)
! I_{JJ, NMPQ} = sum_{k=1}^3 int dr ds j^k_NM(r) j^k_PQ(s) exp(-I alpha (r-s)^2/c^2)
!
! 131227 ! 
!======================================================
  use Precision
  use DiracOutput
  use IntegralStorage
  implicit none
  
  complex(kind=dp),intent(out) :: intIJJ_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)

  real(kind=dp),intent(in) :: alphaJJ
  complex(kind=dp) :: intIJJ(-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2) 
  
  integer :: nn,mm,pp,qq  ! index for 1~4*NBS
  integer :: n,m,p,q ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b,c,d ! + or -
  integer :: aa,bb,cc,dd ! 1 -> +, 2 -> -

  character(LEN=80) :: temp
  real(kind=dp) :: int_real,int_complex

!!$  if(there_is_twoele) then
!!$     write(*,*) " There is twoele file. Read integrals from ", file_twoele
!!$     open(unit=200,file=file_twoele,status='unknown',form='formatted')
!!$
!!$     do nn=1,4*NBS
!!$        do mm=1,4*NBS
!!$           do pp=1,4*NBS
!!$              do qq=1,4*NBS
!!$
!!$!                 read(200,*) temp,temp,temp,temp, temp,temp,temp,temp, int_real,int_complex            
!!$                 ! use same format as storage (see below)
!!$                 read(200,*) int_real,int_complex            
!!$                 twoele_Qmat(nn,mm,pp,qq) = cmplx(int_real,int_complex,dp)
!!$
!!$              end do
!!$           end do
!!$        end do
!!$     end do
!!$     write(*,*) " Reading done."
!!$     
!!$  else ! if there isn't precomputed integrals, compute integrals.
     write(*,*) " There is no intIJJ file. Compute integrals and store at fort.15."

!     alphaJJ = 0._dp
     call calc_intIJJ_mat(alphaJJ,intIJJ) 

     do nn=1,4*NBS
        do mm=1,4*NBS
           do pp=1,4*NBS
              do qq=1,4*NBS
                 
                 call index_from_Qmat(nn,n,a)
                 call index_from_Qmat(mm,m,b)
                 call index_from_Qmat(pp,p,c)
                 call index_from_Qmat(qq,q,d)

                 ! "+" --> 1, "-" --> 2 to be used in modified routines 
                 if(a.eq."+") aa = 1
                 if(a.eq."-") aa = 2
                 if(b.eq."+") bb = 1
                 if(b.eq."-") bb = 2
                 if(c.eq."+") cc = 1
                 if(c.eq."-") cc = 2
                 if(d.eq."+") dd = 1
                 if(d.eq."-") dd = 2

                 write(15,"(2es16.6)") intIJJ(n,aa,m,bb,p,cc,q,dd)  ! for storage
                 intIJJ_Qmat(nn,mm,pp,qq) = intIJJ(n,aa,m,bb,p,cc,q,dd) 

              end do
           end do
        end do
     end do
     write(*,*) " Integrals I_JJ done."
     
!!$  end if
     
  return
end subroutine setQmat_intIJJ


!=========================================================================
subroutine calc_intIJJ_mat(alphajj,intIJJ)
! I_JJ = Sum_k int dr ds j^k(r)_{na mb} j^k(s)_{pc qd} exp(-i alpha (r-s)^2/c^2)
! modified from calc_inttwoele_mat in sub_int.f90
! 130302
!=========================================================================  
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none
  
  real(kind=dp),intent(in) :: alphaJJ
  complex(kind=dp),intent(out) :: intIJJ(-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2) 
  ! (n^a m^b|p^c q^d)
  ! negative values account for Kramers partners. 0 is not used.

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_pg(4,NMAX_PG)

!  complex(kind=dp) :: intthetaJJ_pg_mat(NBS_S,NBS_S,NBS_S,NBS_S,4,4,4,4) !(ij |theta_jj| kl)^{alpha,beta,gamma,delta}
  complex(kind=dp) :: intthetaJJ_pg_mat(NBS_S,NBS_S,NBS_S,NBS_S,Ngamjj) !(ij |theta_jj| kl)^{alpha,beta,gamma,delta}

  integer :: i,j,k,l ! p.g. indice
  integer :: n,m,p,q,a,b,c,d ! MO indice 
  integer :: alpha,beta,gamma,delta ! spinor indice
  integer :: igj  ! index of GamJJ

  complex(kind=dp) :: tmp1(NBS_S,NBS_S,NBS_S,-NBS:NBS,2,Ngamjj) !(i j |theta_jj| k q^d)^{alpha beta gamma delta}
  complex(kind=dp) :: tmp2(NBS_S,NBS_S,-NBS:NBS,2,-NBS:NBS,2,Ngamjj) !(i j |theta_jj| p^c q^d)^{alpha beta gamma delta}
  complex(kind=dp) :: tmp3(NBS_S,-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2,Ngamjj) !(i m^b |theta_jj| p^c q^d)^{alpha beta gamma delta}
  complex(kind=dp) :: tmp4(-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2,Ngamjj) !(n^a m^b |theta_jj| p^c q^d)^{alpha beta gamma delta}
  complex(kind=dp) :: sum

  call copy_DiracOutput_pg(pg) ! get pg

!  intthetaJJ_pg_mat = (0._dp,0._dp)
!  tmp4 = (0._dp,0._dp)
!  tmp3 = (0._dp,0._dp)
!  tmp2 = (0._dp,0._dp)
!  tmp1 = (0._dp,0._dp)
!  intIJJ = (0._dp,0._dp)

  !(i j |theta_jj| k l)^{alpha beta gamma delta}
  call calc_thetaJJ_pg_mat(alphaJJ,NBS_L,NBS_S,pg,intthetaJJ_pg_mat)

  !(i j |theta_jj| k q^d)^{alpha beta gamma delta}
  do igj = 1,Ngamjj
     delta = GamJJ(igj)%d
     do i=1,NBS_S
        do j=1,NBS_S
           do k=1,NBS_S
              do q = -NBS,NBS  !q=0 is not used.
                 do d = 1,2
                    sum = (0._dp,0._dp)
                    call copy_DiracOutput_cp(q,d,c_pg)
                    do l=1,NBS_S
                       sum = sum +c_pg(delta,l)*intthetaJJ_pg_mat(i,j,k,l,igj)
                       ! c^delta_{q^d l} x (ij |theta_jj| kl)^{alpha beta gamma delta}
                    end do
                    tmp1(i,j,k,q,d,igj) = sum
                 end do
              end do
           end do
        end do
     end do
  end do

  !(i j |theta_jj| p^c q^d)^{alpha beta gamma delta}
  do igj = 1,Ngamjj
     gamma = GamJJ(igj)%c
     do i=1,NBS_S
        do j=1,NBS_S
           do q = -NBS,NBS  !q=0 is not used.
              do d = 1,2
                 do p = -NBS,NBS  !p=0 is not used.
                    do c = 1,2
                       sum = (0._dp,0._dp)
                       call copy_DiracOutput_cp(p,c,c_pg)
                       do k=1,NBS_S
                          sum = sum +conjg(c_pg(gamma,k))*tmp1(i,j,k,q,d,igj)
                          ! (c^gamma_{p^c k})^* x (ij |theta_jj| k q^d)^{alpha beta gamma delta}
                       end do
                       tmp2(i,j,p,c,q,d,igj) = sum
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

  !(i m^b |theta_jj| p^c q^d)^{alpha beta gamma delta}
  do igj = 1,Ngamjj
     beta = GamJJ(igj)%b
     do i=1,NBS_S
        do q = -NBS,NBS  !q=0 is not used.
           do d = 1,2
              do p = -NBS,NBS  !p=0 is not used.
                 do c = 1,2
                    do m = -NBS,NBS  !m=0 is not used.
                       do b = 1,2
                          sum = (0._dp,0._dp)
                          call copy_DiracOutput_cp(m,b,c_pg)
                          do j=1,NBS_S
                             sum = sum +c_pg(beta,j)*tmp2(i,j,p,c,q,d,igj)
                             ! c^beta_{m^b j} x (i j |theta_jj| p^c q^d)^{alpha beta gamma delta}
                          end do
                          tmp3(i,m,b,p,c,q,d,igj) = sum
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

  !(n^a m^b |theta_jj| p^c q^d)^{alpha beta gamma delta}
  do igj = 1,Ngamjj
     alpha = GamJJ(igj)%a
     do q = -NBS,NBS  !q=0 is not used.
        do d = 1,2
           do p = -NBS,NBS  !p=0 is not used.
              do c = 1,2
                 do m = -NBS,NBS  !m=0 is not used.
                    do b = 1,2
                       do n = -NBS,NBS  !n=0 is not used.
                          do a = 1,2
                             sum = (0._dp,0._dp)
                             call copy_DiracOutput_cp(n,a,c_pg)
                             do i=1,NBS_S
                                sum = sum +conjg(c_pg(alpha,i))*tmp3(i,m,b,p,c,q,d,igj)
                                ! (c^alpha_{n^a i})^* x (i m^b |theta_jj| p^c q^d)^{alpha beta gamma delta}
                             end do
                             tmp4(n,a,m,b,p,c,q,d,igj) = sum
                          end do
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

  ! (j_{n^a m^b} |theta_jj| j_{p^c q^d})
  do q = -NBS,NBS  !q=0 is not used.
     do d = 1,2
        do p = -NBS,NBS  !p=0 is not used.
           do c = 1,2
              do m = -NBS,NBS  !m=0 is not used.
                 do b = 1,2
                    do n = -NBS,NBS  !n=0 is not used.
                       do a = 1,2
                          
                          sum = (0._dp,0._dp)
                          do igj = 1,Ngamjj
                             sum = sum  +GamJJ(igj)%val *tmp4(n,a,m,b,p,c,q,d,igj)    
                             ! GammaJJ_{alpha beta gamma delta} x (n^a m^b |theta_jj| p^c q^d)^{alpha beta gamma delta}
                          end do
                          intIJJ(n,a,m,b,p,c,q,d) = sum * CCC**2
                          
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do
  
  return
end subroutine calc_intIJJ_mat


!=========================================================================
subroutine calc_thetaJJ_pg_mat(Alpha,NL,NS,pg,intthetaJJ_pg_mat)
! Every component of 
! (ij | theta_jj | kl)^{alpha beta gamma delta}
! Int dr ds (g^si(r)^+)_i (g^sj(r))_j (g^sk(s)^+)_k (g^sl(s))_l exp(-i alpha (r-s)^2/c^2)
! Note that g depends on large or small components
! g^a = g^L if a=1,2
!     = g^S if a=3,4
!130305
!=========================================================================  
  use Precision
  use DefineTypes
  use Constants
  implicit none
  
  real(kind=dp),intent(in) :: alpha
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
!  complex(kind=dp),intent(out) :: intthetaJJ_pg_mat(NS,NS,NS,NS,4,4,4,4) !(ij |theta_jj| kl)^abcd
  complex(kind=dp),intent(out) :: intthetaJJ_pg_mat(NS,NS,NS,NS,Ngamjj) !(ij |theta_jj| kl)^abcd
  
  integer :: i,j,k,l ! p.g. indice
  integer :: a,b,c,d ! spinor indice
  integer :: n ! index of GamJJ
  integer :: NLS
  complex(kind=dp) :: intthetaJJ_pg
  
  intthetaJJ_pg_mat = (0._dp,0._dp)

  do n=1,Ngamjj
     a = GamJJ(n)%a
     b = GamJJ(n)%b
     c = GamJJ(n)%c
     d = GamJJ(n)%d
     do i=1,NLS(a,NL,NS)
        do j=1,NLS(b,NL,NS)
           do k=1,NLS(c,NL,NS)
              do l=1,NLS(d,NL,NS)
                 intthetaJJ_pg_mat(i,j,k,l,n) = intthetaJJ_pg(alpha,i,j,k,l,a,b,c,d,NL,NS,pg)
              end do
           end do
        end do
     end do
  end do

!!$  do a=1,4
!!$     do b=1,4
!!$        do c=1,4
!!$           do d=1,4
!!$              do i=1,NLS(a,NL,NS)
!!$                 do j=1,NLS(b,NL,NS)
!!$                    do k=1,NLS(c,NL,NS)
!!$                       do l=1,NLS(d,NL,NS)
!!$                          intthetaJJ_pg_mat(i,j,k,l,a,b,c,d) = intthetaJJ_pg(alpha,i,j,k,l,a,b,c,d,NL,NS,pg)
!!$                       end do
!!$                    end do
!!$                 end do
!!$              end do
!!$           end do
!!$        end do
!!$     end do
!!$  end do

  return
end subroutine calc_thetaJJ_pg_mat

!=========================================================================
function intthetaJJ_pg(alpha,i,j,k,l,s_i,s_j,s_k,s_l,NL,NS,pg)
! (Ij | theta_jj | kl)^{alpha beta gamma delta}
! Int dr ds (g^si(r)^+)_i (g^sj(r))_j (g^sk(s)^+)_k (g^sl(s))_l exp(-i alpha (r-s)^2/c^2)
! thetaJJ for primitive gaussians.
! Note that g depends on large or small components
! g^si = g^L if si=1,2
!      = g^S if si=3,4
! 130302
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none

  complex(kind=dp) :: intthetaJJ_pg

  real(kind=dp),intent(in) :: alpha
  integer,intent(in) :: i,j,k,l ! p.g. indice
  integer,intent(in) :: s_i,s_j,s_k,s_l ! spinor indice
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  
  integer :: NLS
  complex(kind=dp) :: thetaJJ
  real(kind=dp) :: posA(3), posB(3), posC(3), posD(3) ! position of center of gaussian 
  real(kind=dp) :: alphaA, alphaB, alphaC, alphaD ! exponents
  integer :: n1(3),nbar1(3),n2(3),nbar2(3)  
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg

  call set_pg(s_i,i,NL,NS,pg,posA,alphaA,n1(1),n1(2),n1(3))
  call set_pg(s_j,j,NL,NS,pg,posB,alphaB,nbar1(1),nbar1(2),nbar1(3))
  call set_pg(s_k,k,NL,NS,pg,posC,alphaC,n2(1),n2(2),n2(3))
  call set_pg(s_l,l,NL,NS,pg,posD,alphaD,nbar2(1),nbar2(2),nbar2(3))
  call gauss_int_thetaJJ(alpha,posA,alphaA,n1,posB,alphaB,nbar1,posC,alphaC,n2,posD,alphaD,nbar2,thetaJJ) 

  f_norm_pg = norm_pg(alphaA,n1(1),n1(2),n1(3))*norm_pg(alphaB,nbar1(1),nbar1(2),nbar1(3)) &
       *norm_pg(alphaC,n2(1),n2(2),n2(3))*norm_pg(alphaD,nbar2(1),nbar2(2),nbar2(3))

  intthetaJJ_pg = thetaJJ*f_norm_pg
  
  return
end function intthetaJJ_pg


!=====================================================================================================
subroutine gauss_int_thetaJJ(alpha,posA,alphaA,n1,posB,alphaB,nbar1,posC,alphaC,n2,posD,alphaD,nbar2,thetaJJ)
! Calculate two I_jj integral for primitive gaussians (AB|theta_jj|CD)
! 130302
! 
! alpha=0  -> use overlap integral
! alpha/=0 -> use formula for [NLM|theta_jj|N'L'M']
!=====================================================================================================
  use Precision
  use Constants  ! for IU and CCC
  implicit none
  
  real(kind=dp),intent(in) :: alpha
  real(kind=dp),intent(in) :: posA(3),posB(3),posC(3),posD(3)
  real(kind=dp),intent(in) :: alphaA,alphaB,alphaC,alphaD
  integer,intent(in) :: n1(3),nbar1(3),n2(3),nbar2(3)  
  complex(kind=dp),intent(out) :: thetaJJ

  integer :: i,j,k
  integer :: i1,j1,k1,i2,j2,k2
  
  real(kind=dp), allocatable, dimension(:,:,:) :: d1,e1,f1,d2,e2,f2
  real(kind=dp) :: d1ijk,d2ijk

  real(kind=dp) :: posP(3),alphaP,c_PAB,vecPA(3),vecPB(3) ! mixing of A and B
  real(kind=dp) :: posQ(3),alphaQ,c_QCD,vecQC(3),vecQD(3) ! mixing of C and D
  real(kind=dp) :: vecPQ(3)
  complex(kind=dp) :: A,B,alphaT
  complex(kind=dp), allocatable, dimension(:,:,:,:,:,:) :: NLMthetaJJNLM 
  integer :: n1_sum(3)  ! upper bounds of N,L,M
  integer :: n2_sum(3)  ! upper bounds of N',L',M'
  complex(kind=dp) :: integral,thetaJJ_tmp

  real(kind=dp) :: overlap1,overlap2

  if(alpha.eq.0._dp) then  ! use overlap integral
     call gauss_int_overlap(posA,alphaA,n1(1),n1(2),n1(3),posB,alphaB,nbar1(1),nbar1(2),nbar1(3),overlap1)
     call gauss_int_overlap(posC,alphaC,n2(1),n2(2),n2(3),posD,alphaD,nbar2(1),nbar2(2),nbar2(3),overlap2)
     thetaJJ = overlap1*overlap2
!     write(*,"(4es16.6,3i6,4es16.6,3i6,1es16.6)") posA(1),posA(2),posA(3),alphaA,n1(1),n1(2),n1(3),posB(1),posB(2),posB(3),alphaB,nbar1(1),nbar1(2),nbar1(3),overlap1
!     write(100,"(2i6,1es16.6)") n1(1)+n1(2)+n1(3),nbar1(1)+nbar1(2)+nbar1(3),overlap1
!     if(thetaJJ.ne.(0._dp,0._dp)) then
!        write(*,"(4es16.6)") overlap1,overlap2,thetaJJ
!     end if
  else
     
     do i=1,3
        n1_sum(i) = n1(i)+nbar1(i)
        n2_sum(i) = n2(i)+nbar2(i)
     end do
     
     call calc_gaussP(posA,alphaA,posB,alphaB,c_PAB,posP,alphaP,vecPA,vecPB)  ! construct P 
     call calc_gaussP(posC,alphaC,posD,alphaD,c_QCD,posQ,alphaQ,vecQC,vecQD)  ! construct Q
     
     ! construct d,e,f,d',e',f'
     allocate(d1(0:n1(1),0:nbar1(1),0:n1_sum(1)))
     allocate(e1(0:n1(2),0:nbar1(2),0:n1_sum(2)))
     allocate(f1(0:n1(3),0:nbar1(3),0:n1_sum(3)))
     allocate(d2(0:n2(1),0:nbar2(1),0:n2_sum(1)))
     allocate(e2(0:n2(2),0:nbar2(2),0:n2_sum(2)))
     allocate(f2(0:n2(3),0:nbar2(3),0:n2_sum(3)))
     
     call calc_d(alphaP,vecPA(1),vecPB(1),n1(1),nbar1(1),d1)
     call calc_d(alphaP,vecPA(2),vecPB(2),n1(2),nbar1(2),e1)
     call calc_d(alphaP,vecPA(3),vecPB(3),n1(3),nbar1(3),f1)
     call calc_d(alphaQ,vecQC(1),vecQD(1),n2(1),nbar2(1),d2)
     call calc_d(alphaQ,vecQC(2),vecQD(2),n2(2),nbar2(2),e2)
     call calc_d(alphaQ,vecQC(3),vecQD(3),n2(3),nbar2(3),f2)
     
     ! compute [NLM|theta_jj|N'L'M']
     A = IU*alpha/CCC**2
     B = A*(alphaP+alphaQ) +alphaP*alphaQ
     do i=1,3
        vecPQ(i) = posP(i)-posQ(i)  ! vecD
     end do
     alphaT = 1._dp/(1._dp/alphaP +1._dp/alphaQ +1._dp/A)
     
     allocate( NLMthetaJJNLM(0:n1_sum(1),0:n1_sum(2),0:n1_sum(3),0:n2_sum(1),0:n2_sum(2),0:n2_sum(3)) ) 
     
     do i1=0,n1_sum(1)  ! N
        do j1=0,n1_sum(2)  ! L
           do k1=0,n1_sum(3)  ! M
              do i2=0,n2_sum(1)  ! N'
                 do j2=0,n2_sum(2) ! L'
                    do k2=0,n2_sum(3) ! M'
                       call calc_thetaJJ(alphaT,B,vecPQ,i1,j1,k1,i2,j2,k2,integral)
                       NLMthetaJJNLM(i1,j1,k1,i2,j2,k2) = integral
                    end do
                 end do
              end do
           end do
        end do
     end do
     
     thetaJJ = (0._dp,0._dp)
     do i1=0,n1_sum(1)  ! N
        do j1=0,n1_sum(2)  ! L
           do k1=0,n1_sum(3)  ! M
              d1ijk = d1(n1(1),nbar1(1),i1) * e1(n1(2),nbar1(2),j1) * f1(n1(3),nbar1(3),k1)
              thetaJJ_tmp = (0._dp,0._dp)
              do i2=0,n2_sum(1)  ! N'
                 do j2=0,n2_sum(2) ! L'
                    do k2=0,n2_sum(3) ! M'
                       d2ijk = d2(n2(1),nbar2(1),i2) * e2(n2(2),nbar2(2),j2) * f2(n2(3),nbar2(3),k2)
                       thetaJJ_tmp = thetaJJ_tmp +d2ijk*NLMthetaJJNLM(i1,j1,k1,i2,j2,k2)
                    end do
                 end do
              end do
              thetaJJ = thetaJJ + thetaJJ_tmp * d1ijk
           end do
        end do
     end do
     
     thetaJJ = c_PAB*c_QCD*thetaJJ
  end if

  return
end subroutine gauss_int_thetaJJ


!========================================================================
subroutine calc_thetaJJ(alphaT,B,vecD,N1,L1,M1,N2,L2,M2,thetaJJ)
! Compute [NLM|theta_jj|N'L'M']
! 130301
!========================================================================
  use Precision
  use Constants  ! for PI
  implicit none

  complex(kind=dp), intent(in) :: alphaT,B
  real(kind=dp),intent(in) :: vecD(3)
  integer,intent(in) :: N1,L1,M1 ! N,L,M
  integer,intent(in) :: N2,L2,M2 ! N',L',M'
  complex(kind=dp), intent(out) :: thetaJJ

  complex(kind=dp) :: hermite,factor
  complex(kind=dp) :: aTDx,aTDy,aTDz
  real(kind=dp) :: D2
  integer :: sum1,sum2

  D2 = vecD(1)**2 +vecD(2)**2 +vecD(3)**2 ! |vecD|^2
  sum1 = N1+L1+M1
  sum2 = N2+L2+M2
  aTDx = sqrt(alphaT)*vecD(1)
  aTDy = sqrt(alphaT)*vecD(2)
  aTDz = sqrt(alphaT)*vecD(3)

  factor = PI**3 *B**(-1.5_dp) *exp(-alphaT*D2) *alphaT**((sum1+sum2)/2._dp) *(-1)**sum1

  thetaJJ = factor *hermite(N1,aTDx)*hermite(L1,aTDy)*hermite(M1,aTDz) &
       & *hermite(N2,aTDx)*hermite(L2,aTDy)*hermite(M2,aTDz) 

  return
end subroutine calc_thetaJJ

!==================================================================
function hermite(n,x)
! Hermite polynomials
! 130301
!==================================================================
  use Precision
  implicit none

  complex(kind=dp) :: hermite
  integer,intent(in) :: n
  complex(kind=dp),intent(in) :: x

  complex(kind=dp) :: y0,y1,y2,y1_sav
  integer :: i

  if(n.lt.0) then
     write(*,*) "n should be 0 or larger in Hermite polynomials."
     stop
  end if

  y0 = (1._dp,0._dp)
  if(n.eq.0) then 
     hermite = y0
  else
     y1 = 2._dp*x
     do i=2,n
        y2 = 2._dp*x*y1 - 2._dp*(i-1)*y0
        y1_sav = y1
        y1 = y2
        y0 = y1_sav
     end do
     hermite = y1
  end if

  return
end function hermite


