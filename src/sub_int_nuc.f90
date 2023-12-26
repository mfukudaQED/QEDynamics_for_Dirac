!==================================================================
subroutine set_NucBasis
! set positions of gaussian functions
!================================================================== 
  use Precision
  use Constants
  use NucBasis
  implicit none

  real(kind=8) :: del_R
  integer :: i,j
  integer :: nx,ny,nz

  NUCTYPE = 1  ! NUCTYPE ! boson(0) of fermion(1)

  NNUC = 1 ! for H

  alpha_N = 1.e3_dp 
  m_N = 1822.89d0  ! real proton mass
!  m_N = 182.289d0 
!  m_N = 1.82289d0
!  m_N = 10000._dp
  
  L_well = 0.2_dp
  L_pg = 0.28_dp
  del_R = L_pg/(N_PGN-1)
  vecR_N = 0._dp

  do i=1,N_PGN
     vecR_N(1,i) = -L_pg/2._dp +del_R*(i-1)
  end do

!  c_nuc(n,m) = c_nuc(n,16-m) if n:odd
!  c_nuc(n,m) = -c_nuc(n,16-m) if n:even

  c_nuc(1,1)  = -0.0245280991468319
  c_nuc(1,2)  = 0.09007664050047938
  c_nuc(1,3)  = -0.1769235561307156
  c_nuc(1,4)  = 0.17591094333388546
  c_nuc(1,5)  = 0.10432249606576617
  c_nuc(1,6)  = 0.18847754568866645
  c_nuc(1,7)  = 0.25291473002415765
  c_nuc(1,8)  = 0.2070511677351923
  c_nuc(1,9)  = 0.2529147300246307
  c_nuc(1,10) = 0.18847754568793418
  c_nuc(1,11) = 0.10432249606658514
  c_nuc(1,12) = 0.1759109433329926
  c_nuc(1,13) = -0.17692355613025004
  c_nuc(1,14) = 0.09007664050016295
  c_nuc(1,15) = -0.024528099146740884

  c_nuc(2,1)  = -0.05362640192898522
  c_nuc(2,2)  = 0.19759233238603274
  c_nuc(2,3)  = -0.3932441350346265
  c_nuc(2,4)  = 0.41637482675045734
  c_nuc(2,5)  = 0.14912348440363207
  c_nuc(2,6)  = 0.32648968433247666
  c_nuc(2,7)  = 0.15119656966039463
  c_nuc(2,8)  = 4.357160010085115d-12
  c_nuc(2,9)  = -0.1511965696688576
  c_nuc(2,10) = -0.326489684324715
  c_nuc(2,11) = -0.1491234844102318
  c_nuc(2,12) = -0.416374826745678
  c_nuc(2,13) = 0.39324413503160716
  c_nuc(2,14) = -0.1975923323845349
  c_nuc(2,15) = 0.05362640192855628

  c_nuc(3,1)  = -0.08165515616907211
  c_nuc(3,2)  = 0.30314311332532246
  c_nuc(3,3)  = -0.6099546317484388
  c_nuc(3,4)  = 0.6573486560923987
  c_nuc(3,5)  = 0.256782395458159
  c_nuc(3,6)  = 0.09999906806383062
  c_nuc(3,7)  = -0.14032228547252457
  c_nuc(3,8)  = -0.5023485826864286
  c_nuc(3,9)  = -0.14032228547151154
  c_nuc(3,10) = 0.09999906806196986
  c_nuc(3,11) = 0.2567823954601683
  c_nuc(3,12) = 0.6573486560906568
  c_nuc(3,13) = -0.6099546317471841
  c_nuc(3,14) = 0.30314311332469596
  c_nuc(3,15) = -0.0816551561688995

  c_nuc(4,1)  = -0.1263824815000815
  c_nuc(4,2)  = 0.4733770499523692
  c_nuc(4,3)  = -0.976968358885807
  c_nuc(4,4)  = 1.1601674382864302
  c_nuc(4,5)  = 0.06921955012595163
  c_nuc(4,6)  = -0.23701338495247773
  c_nuc(4,7)  = -0.612959690314411
  c_nuc(4,8)  = 1.1370723647714727d-11
  c_nuc(4,9)  = 0.6129596902923381
  c_nuc(4,10) = 0.23701338497285201
  c_nuc(4,11) = -0.06921955014308503
  c_nuc(4,12) = -1.1601674382736356
  c_nuc(4,13) = 0.976968358877751
  c_nuc(4,14) = -0.47337704994845325
  c_nuc(4,15) = 0.1263824814989855

  c_nuc(5,1)  = -0.1714772752657301
  c_nuc(5,2)  = 0.6524611836674699
  c_nuc(5,3)  = -1.3861030365789722
  c_nuc(5,4)  = 1.7727998129454685
  c_nuc(5,5)  = -0.22332996820174136
  c_nuc(5,6)  = -1.120155160696469
  c_nuc(5,7)  = 0.20877356525295493
  c_nuc(5,8)  = 0.8029973296663406
  c_nuc(5,9)  = 0.2087735652574163
  c_nuc(5,10) = -1.1201551607043712
  c_nuc(5,11) = -0.2233299681924203
  c_nuc(5,12) = 1.772799812936995
  c_nuc(5,13) = -1.3861030365729665
  c_nuc(5,14) = 0.6524611836643229
  c_nuc(5,15) = -0.17147727526479178

  c_nuc(6,1)  = -0.2588475632364289
  c_nuc(6,2)  = 1.0001792512578729
  c_nuc(6,3)  = -2.2090738222324275
  c_nuc(6,4)  = 3.194577897018155
  c_nuc(6,5)  = -1.7828115450631237
  c_nuc(6,6)  = -1.0081548750480576
  c_nuc(6,7)  = 1.9160855760094742
  c_nuc(6,8)  = 6.117327158851793d-11
  c_nuc(6,9)  = -1.9160855761288627
  c_nuc(6,10) = 1.0081548751578597
  c_nuc(6,11) = 1.7828115449696307
  c_nuc(6,12) = -3.1945778969472487
  c_nuc(6,13) = 2.2090738221869275
  c_nuc(6,14) = -1.000179251235569
  c_nuc(6,15) = 0.2588475632300621

  c_nuc(7,1)  = -0.35701985067488884
  c_nuc(7,2)  = 1.4144467335470534
  c_nuc(7,3)  = -3.2744445591586357
  c_nuc(7,4)  = 5.294113361344026
  c_nuc(7,5)  = -4.845502217537625
  c_nuc(7,6)  = 1.23805766400618
  c_nuc(7,7)  = 3.2416819948523146
  c_nuc(7,8)  = -5.252224106741486
  c_nuc(7,9)  = 3.2416819948903424
  c_nuc(7,10) = 1.2380576639393541
  c_nuc(7,11) = -4.845502217458141
  c_nuc(7,12) = 5.294113361270352
  c_nuc(7,13) = -3.2744445591053
  c_nuc(7,14) = 1.4144467335188429
  c_nuc(7,15) = -0.35701985066639996

  return
end subroutine set_NucBasis

!==================================================================
subroutine set_NucBasis_H2
! set positions of gaussian functions for H2
!================================================================== 
  use Precision
  use Constants
  use NucBasis
  implicit none

  real(kind=8) :: del_R
  integer :: i,j
  integer :: nx,ny,nz

  NUCTYPE = 1  ! NUCTYPE ! boson(0) of fermion(1)

  NNUC = 2 ! for H2

  m_N = 1822.89_dp  ! real proton mass
!  m_N = 182.289_dp 
!  m_N = 1.82289_dp
!  m_N = 10000._dp
  
  alpha_N = 1.e1_dp 
  L_well = 1.6_dp
  L_pg = 2._dp
  del_R = L_pg/(N_PGN-1)
  vecR_N = 0._dp

  do i=1,N_PGN
     vecR_N(3,i) = -L_pg/2._dp +del_R*(i-1)  ! H2 is aligned on z-axis
  end do

!  c_nuc(n,m) = c_nuc(n,16-m) if n:odd
!  c_nuc(n,m) = -c_nuc(n,16-m) if n:even

  c_nuc(1,1)  = -0.11476505289277736
  c_nuc(1,2)  = 0.5583391992199886
  c_nuc(1,3)  = -1.2669012257213546
  c_nuc(1,4)  = 1.5133430340634648
  c_nuc(1,5)  = -0.6712804561916991
  c_nuc(1,6)  = -0.04028118382169025
  c_nuc(1,7)  = 1.2726098347780155
  c_nuc(1,8)  = -1.2253616729352008
  c_nuc(1,9)  = c_nuc(1,7)
  c_nuc(1,10) = c_nuc(1,6)
  c_nuc(1,11) = c_nuc(1,5)
  c_nuc(1,12) = c_nuc(1,4) 
  c_nuc(1,13) = c_nuc(1,3)
  c_nuc(1,14) = c_nuc(1,2)
  c_nuc(1,15) = c_nuc(1,1)

  c_nuc(2,1)  = -0.30099817230061815
  c_nuc(2,2)  = 1.505087859602859
  c_nuc(2,3)  = -3.6498855260976706
  c_nuc(2,4)  = 5.233216184221092
  c_nuc(2,5)  = -4.6384535743233455
  c_nuc(2,6)  =  3.931843412574491
  c_nuc(2,7)  = -1.7989458250164188
  c_nuc(2,8)  = 0.d0
  c_nuc(2,9)  = -c_nuc(2,7)
  c_nuc(2,10) = -c_nuc(2,6)
  c_nuc(2,11) = -c_nuc(2,5)
  c_nuc(2,12) = -c_nuc(2,4)
  c_nuc(2,13) = -c_nuc(2,3)
  c_nuc(2,14) = -c_nuc(2,2)
  c_nuc(2,15) = -c_nuc(2,1)

  c_nuc(3,1)  = -0.39428928965630683
  c_nuc(3,2)  = 1.9473093632456273
  c_nuc(3,3)  = -4.532380503174889
  c_nuc(3,4)  = 5.6849634947477865
  c_nuc(3,5)  = -2.8442535527751556
  c_nuc(3,6)  = -0.15195034897144963
  c_nuc(3,7)  = 2.854269600156633
  c_nuc(3,8)  = -4.708829196095101
  c_nuc(3,9)  = c_nuc(3,7)
  c_nuc(3,10) = c_nuc(3,6)
  c_nuc(3,11) = c_nuc(3,5)
  c_nuc(3,12) = c_nuc(3,4)
  c_nuc(3,13) = c_nuc(3,3)
  c_nuc(3,14) = c_nuc(3,2)
  c_nuc(3,15) = c_nuc(3,1)

  c_nuc(4,1)  = -0.7349313417948198
  c_nuc(4,2)  = 3.741773618823483
  c_nuc(4,3)  = -9.322280537823932
  c_nuc(4,4)  = 13.869871979214302
  c_nuc(4,5)  = -12.411841477706561
  c_nuc(4,6)  = 8.7914553375689
  c_nuc(4,7)  = -5.707545589823882
  c_nuc(4,8)  = 0.d0
  c_nuc(4,9)  = -c_nuc(4,7)
  c_nuc(4,10) = -c_nuc(4,6)
  c_nuc(4,11) = -c_nuc(4,5)
  c_nuc(4,12) = -c_nuc(4,4)
  c_nuc(4,13) = -c_nuc(4,3)
  c_nuc(4,14) = -c_nuc(4,2)
  c_nuc(4,15) = -c_nuc(4,1)

  c_nuc(5,1)  = -0.8993535247825272
  c_nuc(5,2)  = 4.6131008615752975
  c_nuc(5,3)  = -11.459855431777044
  c_nuc(5,4)  = 16.350225247177704
  c_nuc(5,5)  = -11.501421577029728
  c_nuc(5,6)  = 1.4013748840839926
  c_nuc(5,7)  = 3.3328943135752183
  c_nuc(5,8)  = -3.4349102116075523
  c_nuc(5,9)  = c_nuc(5,7)
  c_nuc(5,10) = c_nuc(5,6)
  c_nuc(5,11) = c_nuc(5,5)
  c_nuc(5,12) = c_nuc(5,4)
  c_nuc(5,13) = c_nuc(5,3)
  c_nuc(5,14) = c_nuc(5,2)
  c_nuc(5,15) = c_nuc(5,1)

  c_nuc(6,1)  = -1.6300304531802243
  c_nuc(6,2)  = 8.64733543819561
  c_nuc(6,3)  = -22.98494484691253
  c_nuc(6,4)  = 37.945225884493226
  c_nuc(6,5)  = -38.9492622269624
  c_nuc(6,6)  = 23.59139168503601
  c_nuc(6,7)  = -7.235900206733127
  c_nuc(6,8)  = 0.d0
  c_nuc(6,9)  = -c_nuc(6,7)
  c_nuc(6,10) = -c_nuc(6,6)
  c_nuc(6,11) = -c_nuc(6,5)
  c_nuc(6,12) = -c_nuc(6,4)
  c_nuc(6,13) = -c_nuc(6,3)
  c_nuc(6,14) = -c_nuc(6,2)
  c_nuc(6,15) = -c_nuc(6,1)

  c_nuc(7,1)  = -2.2585613647355376
  c_nuc(7,2)  = 12.464596323287278
  c_nuc(7,3)  = -35.038175009122234
  c_nuc(7,4)  = 63.061023548701506
  c_nuc(7,5)  = -74.65774764413952
  c_nuc(7,6)  = 55.10172097455333
  c_nuc(7,7)  = -19.27987535136487
  c_nuc(7,8)  = 1.3597094443386941
  c_nuc(7,9)  = c_nuc(7,7)
  c_nuc(7,10) = c_nuc(7,6)
  c_nuc(7,11) = c_nuc(7,5)
  c_nuc(7,12) = c_nuc(7,4)
  c_nuc(7,13) = c_nuc(7,3)
  c_nuc(7,14) = c_nuc(7,2)
  c_nuc(7,15) = c_nuc(7,1)

  return
end subroutine set_NucBasis_H2

!========================================================================
!
! routines for nucleus overlap integral
!
!========================================================================

!==================================================================
function intSnuc_mat(i,j)
! overlap integrals for nuclei
! int chi*_i chi_j = c*(i,k) c(j,l) int g_k g_l
! to check orthonormality of chi
!================================================================== 
  use Precision
  use NucBasis
  implicit none

  complex(kind=8) :: intSnuc_mat
  integer,intent(in) :: i,j

  real(kind=8) :: overlap
  real(kind=8) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=8) :: alphaA, alphaB ! exponents
  integer :: nx,nbarx,ny,nbary,nz,nbarz
  real(kind=8) :: norm_pg  ! function
  real(kind=8) :: f_norm_pg

  integer :: k,l

  intSnuc_mat = (0._dp,0._dp)

  !----------------------------
  ! for 1s gauss lobe type
  nx=0;    ny=0;    nz=0
  nbarx=0; nbary=0; nbarz=0
  alphaA = alpha_N
  alphaB = alpha_N
  !----------------------------

  do k=1,N_PGN
     do l=1,N_PGN
        posA(:) = vecR_N(:,k)
        posB(:) = vecR_N(:,l)
        call gauss_int_overlap(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,overlap)
        f_norm_pg = norm_pg(alphaA,nx,ny,nz)*norm_pg(alphaB,nbarx,nbary,nbarz)
        intSnuc_mat = intSnuc_mat +conjg(c_nuc(i,k))*c_nuc(j,l)*overlap*f_norm_pg
     end do
  end do

  return
end function intSnuc_mat

!========================================================================
!
! routines for nucleus kinetic energy integral
!
!========================================================================

!==================================================================
subroutine setmat_Tnuc(Tnuc_mat)
! Compute and set Tnuc_mat = T_kl
!================================================================== 
  use Precision
  use DiracOutput
  use NucBasis 
  implicit none

  complex(kind=8),intent(out) :: Tnuc_mat(NBS_N,NBS_N)

  complex(kind=8) :: intTnuc_mat !function intTnuc_mat(i,j)

  integer :: i,j,ii,jj,s_i,s_j

  do ii=1,NBS_N
     do jj=1,NBS_N
        call convert_index_nuc(NUCTYPE,ii,i,s_i)
        call convert_index_nuc(NUCTYPE,jj,j,s_j)
        if(s_i.eq.s_j) then
           Tnuc_mat(ii,jj) = intTnuc_mat(i,j)
        else
           Tnuc_mat(ii,jj) = (0._dp,0._dp)
        end if
     end do
  end do
  
  return
end subroutine setmat_Tnuc

!==========================================================================================
function intTnuc_mat(i,j)
! kinetic energy integrals for nuclei
! int chi*_i (-1/2 nabla^2) chi_j 
! i,j : labels for atomic nuclei basis functions 
! Note the difference between m_e and m_N
!========================================================================================== 
  use Precision
  use NucBasis
  implicit none

  complex(kind=8) :: intTnuc_mat
  integer,intent(in) :: i,j

  real(kind=8) :: ene_kin
  real(kind=8) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=8) :: alphaA, alphaB ! exponents
  integer :: nx,nbarx,ny,nbary,nz,nbarz
  real(kind=8) :: norm_pg  ! function
  real(kind=8) :: f_norm_pg

  integer :: k,l,ii

  complex(kind=8) :: sum

  sum = (0._dp,0._dp)

  !----------------------------
  ! for 1s gauss lobe type
  nx=0;    ny=0;    nz=0
  nbarx=0; nbary=0; nbarz=0
  alphaA = alpha_N
  alphaB = alpha_N
  !----------------------------

  do k=1,N_PGN
     do l=1,N_PGN
        do ii=1,3
           posA(ii) = vecR_N(ii,k)
           posB(ii) = vecR_N(ii,l)
        end do
        call gauss_int_KE(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,ene_kin) ! nonrelativistic KE (for electron)
        sum =  sum +conjg(c_nuc(i,k))*c_nuc(j,l)*ene_kin
!        write(*,*) k,l,ene_kin
     end do
  end do
  ! normalizatoin does not depende on p.g.
  f_norm_pg = norm_pg(alphaA,nx,ny,nz)*norm_pg(alphaB,nbarx,nbary,nbarz)

  intTnuc_mat = sum*f_norm_pg/m_N  ! multiply by (m_e/m_N) for nucleus
!  write(*,*) i,j,intTnuc_mat

  return
end function intTnuc_mat




!========================================================================
!
! routines for two-nuclei 4c integral
!
!========================================================================

!==================================================================
subroutine setmat_twonuc(twonuc_mat)
! Compute and set twonuc_mat = (i j | k l)
! (NBS_N)^4
! Assume same nuclear species for simplicity.
! computation or reading is only done for boson case
! fermion case is constructed from boson case by symmetry
!================================================================== 
  use Precision
  use DiracOutput ! for FILESFOLDER
  use NucBasis ! for NBS_N
  use IntegralStorage
  implicit none

  complex(kind=8),intent(out) :: twonuc_mat(NBS_N,NBS_N,NBS_N,NBS_N)

  integer :: i,j,k,l ! index for 1~NBS_PHI
  integer :: s_i,s_j,s_k,s_l
  integer :: ii,jj,kk,ll ! index for 1~NBS_N
  integer :: p,q,r,s ! index for 1~N_PGN

  complex(kind=8) :: int
  complex(kind=8):: int_pqrs(N_PGN,N_PGN,N_PGN,N_PGN) ! integral of p.g.

  complex(kind=8) :: sum
  complex(kind=8) :: tmp1(N_PGN,N_PGN,N_PGN,NBS_PHI) ! (pq|rl)
  complex(kind=8) :: tmp2(N_PGN,N_PGN,NBS_PHI,NBS_PHI) ! (pq|kl)
  complex(kind=8) :: tmp3(N_PGN,NBS_PHI,NBS_PHI,NBS_PHI) ! (pj|kl)
  complex(kind=8) :: tmp4(NBS_PHI,NBS_PHI,NBS_PHI,NBS_PHI) ! (ij|kl)

  real(kind=8) :: int_real,int_complex
  character(LEN=300) :: readfile

  if(there_is_twonuc) then
     readfile = trim(FILESFOLDER)//"/"//file_twonuc
     write(*,*) " There is twonuc file. Read integrals from ",file_twonuc
     open(unit=203,file=readfile,status='unknown',form='formatted')

     do i=1,NBS_PHI
        do j=1,NBS_PHI
           do k=1,NBS_PHI
              do l=1,NBS_PHI
                 read(203,*) int_real,int_complex            
                 tmp4(i,j,k,l) = cmplx(int_real,int_complex)
              end do
           end do
        end do
     end do
     write(*,*) " Reading done."
 
  else ! if there isn't precomputed integrals, compute integrals.
     write(*,*) " There is no twonuc file. Compute integrals and store at fort.14."

     ! create (pq|rs)
     do p=1,N_PGN
        do q=1,N_PGN
           do r=1,N_PGN
              do s=1,N_PGN
                 call calc_inttwonuc_pqrs(p,q,r,s,int)
                 int_pqrs(p,q,r,s) = int
!                 if(p.eq.q .and. q.eq.r .and. r.eq.s) then
!                    write(*,*) p,q,r,s,int
!                 end if
              end do
           end do
        end do
     end do

     ! create (pq|rl)
     do p=1,N_PGN
        do q=1,N_PGN
           do r=1,N_PGN
              do l=1,NBS_PHI
                 sum = (0._dp,0._dp)
                 do s=1,N_PGN
                    sum = sum +c_nuc(l,s)*int_pqrs(p,q,r,s)
                 end do
                 tmp1(p,q,r,l) = sum
              end do
           end do
        end do
     end do
     
     ! create (pq|kl)
     do p=1,N_PGN
        do q=1,N_PGN
           do l=1,NBS_PHI
              do k=1,NBS_PHI
                 sum = (0._dp,0._dp)
                 do r=1,N_PGN
                    sum = sum +cmplx(c_nuc(k,r))*tmp1(p,q,r,l)
                 end do
                 tmp2(p,q,k,l) = sum
              end do
           end do
        end do
     end do

     ! create (pj|kl)
     do p=1,N_PGN
        do k=1,NBS_PHI
           do l=1,NBS_PHI
              do j=1,NBS_PHI
                 sum = (0._dp,0._dp)
                 do q=1,N_PGN
                    sum = sum +c_nuc(j,q)*tmp2(p,q,k,l)
                 end do
                 tmp3(p,j,k,l) = sum
              end do
           end do
        end do
     end do
     
     ! create (ij|kl)
     do j=1,NBS_PHI
        do k=1,NBS_PHI
           do l=1,NBS_PHI
              do i=1,NBS_PHI
                 sum = (0._dp,0._dp)
                 do p=1,N_PGN
                    sum = sum +cmplx(c_nuc(i,p))*tmp3(p,j,k,l)
                 end do
                 tmp4(i,j,k,l) = sum
                 write(14,"(2es16.6)") sum   ! for storage                                  
              end do
           end do
        end do
     end do

     write(*,*) " Computation done."

  end if

  ! set (ij|kl) including fermion case
  do ii=1,NBS_N
     do jj=1,NBS_N
        do kk=1,NBS_N
           do ll=1,NBS_N
              call convert_index_nuc(NUCTYPE,ii,i,s_i)
              call convert_index_nuc(NUCTYPE,jj,j,s_j)
              call convert_index_nuc(NUCTYPE,kk,k,s_k)
              call convert_index_nuc(NUCTYPE,ll,l,s_l)
              if((s_i.eq.s_j).and.(s_k.eq.s_l)) then
                 twonuc_mat(ii,jj,kk,ll) = tmp4(i,j,k,l) *(Z_N)**2
              else
                 twonuc_mat(ii,jj,kk,ll) = (0._dp,0._dp)
              end if
           end do
        end do
     end do
  end do

  return
end subroutine setmat_twonuc

!=========================================================================
subroutine calc_inttwonuc_pqrs(p,q,r,s,int_pqrs)
! Routine for 4-center integrals (nucleus-nucleus)
! 4c integration of primitive gaussians
! int dr ds (g^*(r))_p (g(r))_q (1/|r-s|) (g^*(s))_r (g(s))_s
! used in setmat_twonuc
!=========================================================================  
  use Precision
  use NucBasis
  implicit none

  integer,intent(in) :: p,q,r,s  ! labels for p.g.
  complex(kind=8),intent(out) :: int_pqrs
  
  integer :: kk
  real(kind=8) :: twoele
  real(kind=8) :: posA(3), posB(3), posC(3), posD(3) ! position of center of gaussian 
  real(kind=8) :: alphaA, alphaB, alphaC, alphaD ! exponents
  integer :: n1(3),nbar1(3),n2(3),nbar2(3)  
  real(kind=8) :: norm_pg  ! function
  real(kind=8) :: f_norm_pg

  do kk=1,3
     posA(kk) = vecR_N(kk,p)   ! module NucBasis
     posB(kk) = vecR_N(kk,q)   ! module NucBasis
     posC(kk) = vecR_N(kk,r)   ! module NucBasis
     posD(kk) = vecR_N(kk,s)   ! module NucBasis
  end do
  !----------------------------
  ! for 1s gauss lobe type
  n1(1)=0;    n1(2)=0;    n1(3)=0
  nbar1(1)=0; nbar1(2)=0; nbar1(3)=0
  n2(1)=0;    n2(2)=0;    n2(3)=0
  nbar2(1)=0; nbar2(2)=0; nbar2(3)=0
  alphaA = alpha_N
  alphaB = alpha_N
  alphaC = alpha_N
  alphaD = alpha_N
  !----------------------------
  call gauss_int_twoele(posA,alphaA,n1,posB,alphaB,nbar1,posC,alphaC,n2,posD,alphaD,nbar2,twoele)  ! (AB|CD)
  f_norm_pg = norm_pg(alphaA,n1(1),n1(2),n1(3))*norm_pg(alphaB,nbar1(1),nbar1(2),nbar1(3)) &
       *norm_pg(alphaC,n2(1),n2(2),n2(3))*norm_pg(alphaD,nbar2(1),nbar2(2),nbar2(3))
  int_pqrs = twoele*f_norm_pg

  return
end subroutine calc_inttwonuc_pqrs


!========================================================================
!
! routines for nuclei-electron 4c integral
!
!========================================================================

!==================================================================
subroutine setQmat_nucele(nucele_Qmat)
! set nucele_Qmat = (i j | P Q)
! modified 2012.6.6 (to use faster routines)
!================================================================== 
  use Precision
  use DiracOutput
  use NucBasis 
  use IntegralStorage
  use Constants
  implicit none

  complex(kind=8),intent(out) :: nucele_Qmat(NBS_N,NBS_N,4*NBS,4*NBS)

  complex(kind=8) :: intnucele_mat ! --> very slow function intnucele_mat(i,j,p,c,q,d)
  complex(kind=8) :: intnucele(NBS_PHI,NBS_PHI,-NBS:NBS,2,-NBS:NBS,2) ! (i j | p^c q^d)

  integer :: pp,qq  ! index for 1~4*NBS
  integer :: p,q ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: c,d ! + or -
  integer :: cc,dd ! 1 -> +, 2 -> -

  integer :: i,j ! index for 1~NBS_PHI
  integer :: ii,jj ! index for 1~NBS_N
  integer :: s_i,s_j

  complex(kind=8) :: sum
  complex(kind=8) :: tmp4(NBS_PHI,NBS_PHI,4*NBS,4*NBS)

  real(kind=8) :: int_real,int_complex
  character(LEN=300) :: readfile

  ! read or compute nucele integral
  if(there_is_nucele) then
     readfile = trim(FILESFOLDER)//"/"//file_nucele
     write(*,*) " There is nuc-ele file. Read integrals from ",file_nucele
     open(unit=202,file=readfile,status='unknown',form='formatted')

     do i=1,NBS_PHI
        do j=1,NBS_PHI
           do pp=1,4*NBS
              do qq=1,4*NBS

                 ! use same format as storage
                 read(202,*) int_real,int_complex            
                 tmp4(i,j,pp,qq) = cmplx(int_real,int_complex)

              end do
           end do
        end do
     end do
     write(*,*) " Reading done."

  else ! if there isn't precomputed integrals, compute integrals.
     write(*,*) " There is no nuc-ele file. Compute integrals and store at fort.13."

!!$     !--- faster routine but need large memory ---
!!$     call calc_intnucele_mat(intnucele) 
!!$     !------------------------------------------     
     
     do i=1,NBS_PHI
        do j=1,NBS_PHI
           do pp=1,4*NBS
              do qq=1,4*NBS
                 ! convert indices from pp to p^c etc. See routines in sub_makeQmat.f90 for usage.
                 call index_from_Qmat(pp,p,c)
                 call index_from_Qmat(qq,q,d)

                 !--- faster routine but need large memory ---
                 ! "+" --> 1, "-" --> 2 to be used in modified routines 
                 if(c.eq."+") cc = 1
                 if(c.eq."-") cc = 2
                 if(d.eq."+") dd = 1
                 if(d.eq."-") dd = 2
                 !------------------------------------------

                 !----------- very slow way previously used --------------
                 sum = intnucele_mat(i,j,p,c,q,d)
                 tmp4(i,j,pp,qq) = sum
!                 nucele_Qmat(i,j,pp,qq) = sum   ! we assume that we do not compute time evolution for the first time.
                 write(13,"(2es16.6)") sum   ! for storage                 
                 !--------------------------------------------------------

!!$                 !--- faster routine but need large memory ---
!!$                 write(13,"(2es16.6)") intnucele(i,j,p,cc,q,dd)  ! for storage
!!$                 tmp4(i,j,pp,qq) = intnucele(i,j,p,cc,q,dd)
!!$                 !------------------------------------------
              end do
           end do
        end do
     end do

     write(*,*) " Computation done."
     
  end if

  ! set nucele_Qmat
  do ii=1,NBS_N
     do jj=1,NBS_N
        do pp=1,4*NBS
           do qq=1,4*NBS
              call convert_index_nuc(NUCTYPE,ii,i,s_i)
              call convert_index_nuc(NUCTYPE,jj,j,s_j)
              if(s_i.eq.s_j) then  ! for boson case, it is set as s_i=s_j=0 so this condition is true.
                 nucele_Qmat(ii,jj,pp,qq) = tmp4(i,j,pp,qq)*Ze*Z_N
              else
                 nucele_Qmat(ii,jj,pp,qq) = (0._dp,0._dp)
              end if
           end do
        end do
     end do
  end do

  return
end subroutine setQmat_nucele

!=========================================================================
subroutine calc_intnucele_mat(intnucele)
! sum_alpha
! int dr ds (chi(r)^+)_i (chi(r))_j (1/|r-s|)(psi^alpha(s)^+)_pc (psi^alpha(s))_qd 
! (i j | p^c q^d)
! 2012.6.6
!=========================================================================  
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use NucBasis 
  implicit none
  
  complex(kind=8),intent(out) :: intnucele(NBS_PHI,NBS_PHI,-NBS:NBS,2,-NBS:NBS,2) 
  ! (i j |p^c q^d)
  ! negative values account for Kramers partners. 0 is not used.

  type(primitive_gaussian) :: pg
  complex(kind=8) :: c_pg(4,NMAX_PG)

  complex(kind=8) :: intnucele_pg_mat(N_PGN,N_PGN,NBS_S,NBS_S,4) !(kl|rs)^a

  integer :: k,l,r,s ! p.g. indice
  integer :: i,j,p,q,c,d ! MO indice 
  integer :: alpha ! spinor indice

  complex(kind=8) :: tmp1(N_PGN,  N_PGN,   NBS_S    ,-NBS:NBS,2,4) ! (k l | r q^d)^alpha
  complex(kind=8) :: tmp2(N_PGN,  N_PGN,  -NBS:NBS,2,-NBS:NBS,2,4) ! (k l | p^c q^d)^alpha
  complex(kind=8) :: tmp3(N_PGN,  NBS_PHI,-NBS:NBS,2,-NBS:NBS,2,4) ! (k j | p^c q^d)^alpha
  complex(kind=8) :: tmp4(NBS_PHI,NBS_PHI,-NBS:NBS,2,-NBS:NBS,2,4) ! (i j | p^c q^d)^alpha
  complex(kind=8) :: sum

  call copy_DiracOutput_pg(pg) ! get pg

  ! (k l | r s)^alpha
  call calc_intnucele_pg_mat(N_PGN,NBS_L,NBS_S,pg,intnucele_pg_mat)

  ! (k l | r q^d)^alpha   s-->q^d
  do alpha=1,4
     do k=1,N_PGN
        do l=1,N_PGN
           do r=1,NBS_S
              do q = -NBS,NBS  !q=0 is not used.
                 do d = 1,2
                    sum = (0._dp,0._dp)
                    call copy_DiracOutput_cp(q,d,c_pg)
                    do s=1,NBS_S
                       sum = sum +c_pg(alpha,s)*intnucele_pg_mat(k,l,r,s,alpha)
                       ! c^alpha_{q^d s} x (kl|rs)^alpha
                    end do
                    tmp1(k,l,r,q,d,alpha) = sum
                 end do
              end do
           end do
        end do
     end do
  end do

  ! (k l | p^c q^d)^alpha   r-->p^c
  do alpha=1,4
     do k=1,N_PGN
        do l=1,N_PGN
           do q = -NBS,NBS  !q=0 is not used.
              do d = 1,2
                 do p = -NBS,NBS  !p=0 is not used.
                    do c = 1,2
                       sum = (0._dp,0._dp)
                       call copy_DiracOutput_cp(p,c,c_pg)
                       do r=1,NBS_S
                          sum = sum +conjg(c_pg(alpha,r))*tmp1(k,l,r,q,d,alpha)
                          ! (c^alpha_{p^c r})^* x (kl|r q^d)^alpha
                       end do
                       tmp2(k,l,p,c,q,d,alpha) = sum
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

  ! (k j | p^c q^d)^alpha   l-->j
  do alpha=1,4
     do k=1,N_PGN
        do q = -NBS,NBS  !q=0 is not used.
           do d = 1,2
              do p = -NBS,NBS  !p=0 is not used.
                 do c = 1,2
                    do j=1,NBS_PHI
                       sum = (0._dp,0._dp)
                       do l=1,N_PGN
                          sum = sum +c_nuc(j,l)*tmp2(k,l,p,c,q,d,alpha)
                       end do
                       tmp3(k,j,p,c,q,d,alpha) = sum
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

  ! (i j | p^c q^d)^alpha   k-->i
  do alpha=1,4
     do q = -NBS,NBS  !q=0 is not used.
        do d = 1,2
           do p = -NBS,NBS  !p=0 is not used.
              do c = 1,2
                 do j=1,NBS_PHI
                    do i=1,NBS_PHI
                       sum = (0._dp,0._dp)
                       do k=1,N_PGN
                          sum = sum +cmplx(c_nuc(i,k))*tmp3(k,j,p,c,q,d,alpha)
                       end do
                       tmp4(i,j,p,c,q,d,alpha) = sum
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

  ! (i j | p^c q^d)  (sum over alpha)
  do q = -NBS,NBS  !q=0 is not used.
     do d = 1,2
        do p = -NBS,NBS  !p=0 is not used.
           do c = 1,2
              do j=1,NBS_PHI
                 do i=1,NBS_PHI
                    sum = (0._dp,0._dp)
                    do alpha=1,4
                       sum = sum +tmp4(i,j,p,c,q,d,alpha)
                    end do
                    intnucele(i,j,p,c,q,d) = sum
                 end do
              end do
           end do
        end do
     end do
  end do

  return
end subroutine calc_intnucele_mat

!=========================================================================
subroutine calc_intnucele_pg_mat(NPGN,NL,NS,pg,intnucele_pg_mat)
! int dr ds (g_nuc(r)^+)_k (g_nuc(r))_l (1/|r-s|)(g^a(s)^+)_r (g^a(s))_s
! Note that g depends on large or small components
! g^a = g^L if a=1,2
!     = g^S if a=3,4
! compute only necessary spinor combinations
! 4c integration for primitive gaussians.
! 2012.6.6
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none
  
  integer,intent(in) :: NPGN,NL,NS
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(out) :: intnucele_pg_mat(NPGN,NPGN,NS,NS,4) !(kl|rs)^a
  
  integer :: k,l,r,s ! p.g. indice
  integer :: a ! spinor indice
  integer :: NLS
  complex(kind=8) :: intnucele_pg

  intnucele_pg_mat = (0._dp,0._dp)

  do a=1,4
     do k=1,NPGN
        do l=1,NPGN
           do r=1,NLS(a,NL,NS)
              do s=1,NLS(a,NL,NS)
                 intnucele_pg_mat(k,l,r,s,a) = intnucele_pg(k,l,r,s,a,a,NL,NS,pg)
              end do
           end do
        end do
     end do
  end do
  
  return
end subroutine calc_intnucele_pg_mat

!=========================================================================
function intnucele_pg(k,l,r,s,s_r,s_s,NL,NS,pg)
! (k l | r s)^{s_r s_s}
! int dr ds (g_nuc(r)^+)_k (g_nuc(r))_l (1/|r-s|)(g^sr(s)^+)_r (g^ss(s))_s 
! 4c integration for primitive gaussians.
! Note that g depends on large or small components
! g^si = g^L if si=1,2
!      = g^S if si=3,4
! 2012.6.6
!=========================================================================  
  use Precision
  use NucBasis
  use DefineTypes
  implicit none

  complex(kind=8) :: intnucele_pg

  integer,intent(in) :: k,l,r,s ! p.g. indice
  integer,intent(in) :: s_r,s_s ! spinor indice
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  
  integer :: NLS
  real(kind=8) :: twoele
  real(kind=8) :: posA(3), posB(3), posC(3), posD(3) ! position of center of gaussian 
  real(kind=8) :: alphaA, alphaB, alphaC, alphaD ! exponents
  integer :: n1(3),nbar1(3),n2(3),nbar2(3)  
  real(kind=8) :: norm_pg  ! function
  real(kind=8) :: f_norm_pg
  integer :: kk

  do kk=1,3
     posA(kk) = vecR_N(kk,k)   ! module NucBasis
     posB(kk) = vecR_N(kk,l)   ! module NucBasis
  end do
  !----------------------------
  ! for 1s gauss lobe type
  n1(1)=0;    n1(2)=0;    n1(3)=0
  nbar1(1)=0; nbar1(2)=0; nbar1(3)=0
  alphaA = alpha_N
  alphaB = alpha_N
  !----------------------------

  call set_pg(s_r,r,NL,NS,pg,posC,alphaC,n2(1),n2(2),n2(3))
  call set_pg(s_s,s,NL,NS,pg,posD,alphaD,nbar2(1),nbar2(2),nbar2(3))

  call gauss_int_twoele(posA,alphaA,n1,posB,alphaB,nbar1,posC,alphaC,n2,posD,alphaD,nbar2,twoele)  ! (AB|CD)

  f_norm_pg = norm_pg(alphaA,n1(1),n1(2),n1(3))*norm_pg(alphaB,nbar1(1),nbar1(2),nbar1(3)) &
       *norm_pg(alphaC,n2(1),n2(2),n2(3))*norm_pg(alphaD,nbar2(1),nbar2(2),nbar2(3))

  intnucele_pg = twoele*f_norm_pg

  return
end function intnucele_pg



!====================================================================
function intnucele_mat(i,j,p,c,q,d)
! 4c integration of atomic nuclei and electrons
! i,j : labels for atomic nuclei basis functions (-> not orthogonal)
! p,q : labels for molecular orbitals.
! c,d : labels for electron/positron ("+"/"-")
!
!---------------------------
! for slow routines
!---------------------------
!====================================================================
  use Precision
  use NucBasis
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none
 
  complex(kind=8) :: intnucele_mat
  integer,intent(in) :: i,j
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: c,d

  type(primitive_gaussian) :: pg
  complex(kind=8) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  complex(kind=8) :: int_ijpq
  integer :: k

!  call copy_DiracOutput(n,a,m,b,pg,c_n,c_m)  ! set pg,c_n,c_m
  call copy_DiracOutput(p,c,q,d,pg,c_p,c_q)  ! set pg,c_p,c_q

  intnucele_mat = (0._dp,0._dp)
  do k=1,4
     call calc_intnucele_ijpq(i,j,k,k,NBS_L,NBS_S,pg,c_p,c_q,int_ijpq) 
     intnucele_mat = intnucele_mat +int_ijpq
  end do

!  intnucele_mat = intnucele_mat*Ze*Z_N
  intnucele_mat = intnucele_mat  ! 2012.6.6

  return
end function intnucele_mat

!=========================================================================
subroutine calc_intnucele_ijpq(i,j,ip,iq,NL,NS,pg,c_p,c_q,int_ijpq)
! int dr ds (chi^*(r))_i (chi(r))_j (1/|r-s|)(psi(s)^+)_ip (psi(s))_iq 
!
!---------------------------
! for slow routines
!---------------------------
!=========================================================================  
  use Precision
  use NucBasis
  use DefineTypes
  implicit none

  integer,intent(in) :: i,j  ! labels for atomic nuclei basis functions
  integer,intent(in) :: ip,iq ! spinor indice
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
!  complex(kind=8),intent(in) :: c_n(4,NMAX_PG),c_m(4,NMAX_PG)
  complex(kind=8),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  
  complex(kind=8),intent(out) :: int_ijpq
  
  integer :: NLS
  integer :: k,l,kk,m,n
  real(kind=8) :: twoele
  real(kind=8) :: posA(3), posB(3), posC(3), posD(3) ! position of center of gaussian 
  real(kind=8) :: alphaA, alphaB, alphaC, alphaD ! exponents
  integer :: n1(3),nbar1(3),n2(3),nbar2(3)  
  real(kind=8) :: norm_pg  ! function
  real(kind=8) :: f_norm_pg

  int_ijpq = (0._dp,0._dp)

  !----------------------------
  ! for 1s gauss lobe type
  n1(1)=0;    n1(2)=0;    n1(3)=0
  nbar1(1)=0; nbar1(2)=0; nbar1(3)=0
  alphaA = alpha_N
  alphaB = alpha_N
  !----------------------------

  do m=1,N_PGN
     do n=1,N_PGN
        do k=1,NLS(ip,NL,NS)
           do l=1,NLS(iq,NL,NS)
              ! set A and B
              do kk=1,3
                 posA(kk) = vecR_N(kk,m)   ! module NucBasis
                 posB(kk) = vecR_N(kk,n)   ! module NucBasis
              end do
              call set_pg(ip,k,NL,NS,pg,posC,alphaC,n2(1),n2(2),n2(3))
              call set_pg(iq,l,NL,NS,pg,posD,alphaD,nbar2(1),nbar2(2),nbar2(3))
              call gauss_int_twoele(posA,alphaA,n1,posB,alphaB,nbar1,posC,alphaC,n2,posD,alphaD,nbar2,twoele)  ! (AB|CD)
              f_norm_pg = norm_pg(alphaA,n1(1),n1(2),n1(3))*norm_pg(alphaB,nbar1(1),nbar1(2),nbar1(3)) &
                   *norm_pg(alphaC,n2(1),n2(2),n2(3))*norm_pg(alphaD,nbar2(1),nbar2(2),nbar2(3))
              int_ijpq = int_ijpq + conjg(c_nuc(i,m))*c_nuc(j,n)*conjg(c_p(ip,k))*c_q(iq,l)*twoele*f_norm_pg
           end do
        end do
     end do
  end do

  return
end subroutine calc_intnucele_ijpq




!========================================================================
!
! utility routines
!
!========================================================================


!!$!==================================================================
!!$subroutine convert_label_NB(i,nx,ny,nz)
!!$!
!!$! i = (nz-1)*L^2 +(ny-1)*L +nx
!!$!================================================================== 
!!$  use NucBasis
!!$
!!$  implicit none
!!$
!!$  integer,intent(in) :: i
!!$  integer,intent(out) :: nx,ny,nz
!!$  integer :: L
!!$
!!$  L = L_N
!!$
!!$  nx = mod(i,L)
!!$  if(nx.eq.0) then
!!$     nx = L
!!$  end if
!!$  ny = mod((i-nx)/L,L)+1
!!$  nz = (i-(ny-1)*L-nx)/L**2 +1
!!$
!!$  return
!!$end subroutine convert_label_NB







!======================================================
subroutine convert_index_nuc(NT,nn,n,s_n)
! for fermion, we label as 1,1bar,2,2bar,3,3bar,...
! no bar = up --> s_n=1
! bar = down  --> s_n=2
!
! s_n=0 if boson
!======================================================
  use Precision
  use DiracOutput
  implicit none
  
  integer,intent(in) :: NT,nn
  integer,intent(out) :: n,s_n

  if(NT.eq.0) then  ! boson
     n = nn
     s_n = 0
  else ! fermion
     if(mod(nn,2).eq.0) then ! even -> down
        n = nn/2
        s_n = 2
     else ! odd -> up
        n = (nn+1)/2
        s_n = 1
     end if
  end if

  return
end subroutine convert_index_nuc
  

!=====================================
function rombint(f,a,b,tol)
! e.g.)
!  real(kind=8) :: rombint
!  A2 = norm*rombint(integrand_area,ybr,ymax,1.d-6)
!=====================================
!  Rombint returns the integral from a to b of using Romberg integration.
!  The method converges provided that f(x) is continuous in (a,b).
!  f must be double precision and must be declared external in the calling
!  routine.  tol indicates the desired relative accuracy in the integral.
!
  parameter (MAXITER=40,MAXJ=5)
  implicit double precision (a-h,o-z)
  dimension g(MAXJ+1)
  double precision f
  external f
  !
  h=0.5d0*(b-a)
  gmax=h*(f(a)+f(b))
  g(1)=gmax
  nint=1
  error=1.0d20
  i=0
10 i=i+1
!  write(*,*) i,abs(error)
  if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) go to 40
  !  Calculate next trapezoidal rule approximation to integral.
  g0=0.0d0
  do k=1,nint
     g0=g0+f(a+(k+k-1)*h)
  end do
  g0=0.5d0*g(1)+h*g0
  h=0.5d0*h
  nint=nint+nint
  jmax=min(i,MAXJ)
  fourj=1.0d0
  do j=1,jmax
     !  Use Richardson extrapolation.
     fourj=4.0d0*fourj
     g1=g0+(g0-g(j))/(fourj-1.0d0)
     g(j)=g0
     g0=g1
  end do
  if (abs(g0).gt.tol) then
     error=1.0d0-gmax/g0
  else
     error=gmax
  end if
  gmax=g0
  g(jmax+1)=g0
!  if(i>5) stop 
  go to 10
40 rombint=g0
  if (i.gt.MAXITER.and.abs(error).gt.tol) then 
     write(*,*) 'Rombint failed to converge; integral, error=',rombint,error
  end if
  return
end function rombint


