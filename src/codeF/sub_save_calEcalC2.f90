! Last Change:08-Jun-2014.
!====================================================================================
subroutine calc_phys
!====================================================================================
  implicit none
  
  character(LEN=10) :: str1,str2,str3,str4
!  integer,parameter :: nprint = 119 ! interval of print files
!!$  integer,parameter :: nprint = 600 ! interval of print files
!  integer,parameter :: nprint = 20 ! interval of print files
  integer,parameter :: nprint = 100 ! interval of print files
!  integer,parameter :: nprint = 1 ! interval of print files
!  integer,parameter :: nprint = 10 ! interval of print files
!!$  integer,parameter :: filenum = 8! number of print file
!  integer,parameter :: filenum = 24! number of print file
  integer,parameter :: filenum = 50! number of print file
!  integer,parameter :: filenum = 100! number of print file
!  integer,parameter :: filenum = 1! number of print file
!  integer,parameter :: filenum = 10000! number of print file
  integer :: i
  ! location where we carry out calculation
  ! electron
  real(kind=8),parameter :: x0 =  0.d0
  real(kind=8),parameter :: y0 =  0.d0
  real(kind=8),parameter :: z0 =  0.d0

!  call getarg(1,str1) !read "nprint"
!  read(str1,'(i10.10)') nprint  ! 

   write(*,*)'# (x0,y0,z0) =',x0,y0,z0
!     call calc_property_BO(nprint,filenum)
!     call calc_property_BO_point(nprint,filenum,x0,y0,z0)
     call calc_property2


end subroutine calc_phys

!====================================================================================
subroutine save_calE(it,time,NEL,NTIME,calE)
!====================================================================================
  use DiracOutput
  use Constants
  implicit none

  integer,intent(in) :: it  ! it-th time step
  real(kind=8),intent(in) :: time ! it*DeltaT
  complex(kind=8),intent(in) :: calE(4*NBS,4*NBS)  ! calE_NM
  integer,intent(in) :: NEL
  integer(kind=dp),intent(in) :: NTIME

  character(LEN=80) :: filename_E
  character(LEN=80) :: FMT


  if (it.eq.0) then
  write(filename_E,'(a)') 'calE.dat'
  open(unit=1200,file=filename_E)

  write(1200,"(1a10,1es16.6)")'# DeltaT=',DeltaT
  write(1200,"(1a10,1i16)")'# NBS=', NBS 
  write(1200,"(1a10,1i6)") "# NEL : ",NEL
  end if
     ! save electron
     write(FMT,'("(1es25.15,"i0"(2es25.15))")') (4*NBS)**2
     write(1200,FMT) time,calE 

  if (it.eq.NTIME) then
  close(unit=1200)
  end if

end subroutine save_calE

!====================================================================================
subroutine save_calC(it,time,NTIME,calC)
!====================================================================================
  use DiracOutput
  use Constants
  use NucBasis
  implicit none

  integer,intent(in) :: it  ! it-th time step
  real(kind=8),intent(in) :: time ! it*DeltaT
  complex(kind=8),intent(in) :: calC(NBS_N,NBS_N)  ! calE_NM
  integer(kind=dp),intent(in) :: NTIME

  character(LEN=80) :: filename_C
  character(LEN=80) :: FMT


  if (it.eq.0) then
  write(filename_C,'(a)') 'calC.dat'
  open(unit=1201,file=filename_C)

  write(1201,"(1a10,1es16.6)")'# DeltaT=',DeltaT
  write(1201,"(1a10,1i16)")'# NBS_N=', NBS_N 
  end if
     ! save electron
     write(FMT,'("(1es25.15,"i0"(2es25.15e3))")') (NBS_N)**2
     write(1201,FMT) time,calC 

  if (it.eq.NTIME) then
  close(unit=1201)
  end if

end subroutine save_calC

!====================================================================================
subroutine read_calE(it,time,NEL,NTIME,calE)
!====================================================================================
  use DiracOutput
  use Constants
  implicit none

  integer,intent(in) :: it  ! it-th time step
  real(kind=8),intent(in) :: time ! it*DeltaT
  real(kind=8) :: time2 ! it*DeltaT
  complex(kind=8),intent(out) :: calE(4*NBS,4*NBS)  ! calE_NM
  integer,intent(in) :: NEL
  integer(kind=dp),intent(in) :: NTIME

  character(LEN=80) :: filename_E
  character(LEN=80) :: FMT


  if (it.eq.0) then
     write(filename_E,'(a)') 'calE.dat'
     open(unit=1200,file=filename_E,status='old')
   
     read(1200,*)
     read(1200,*)
     read(1200,*)
  end if

     ! read electron
     write(FMT,'("(1es25.15,"i0"(2es25.15))")') (4*NBS)**2
     read(1200,FMT) time2,calE 
!     write(1201,FMT) time,calE 
     if(time.ne.time2) then
        write(*,*)'# error read calE'
     end if

  if (it.eq.NTIME) then
  close(unit=1200)
  end if

end subroutine read_calE

!====================================================================================
subroutine calc_property_BO(nprint,filenum)
!====================================================================================
  use Precision
  Use DiracOutput
  use Constants
  use IntegralStorage
  use NucBasis

  implicit none

  character(LEN=10) :: str1,str2,str3,str4
  character(LEN=80) :: filedens, filetorque, filespin, filespind, filezeta, filechiral
  character(LEN=80) :: filej
  character(LEN=80) :: fileA, filetA
  character(LEN=80) :: fileele,filenuc
  character(LEN=80) :: filename_E,filename_C
  character(LEN=10) :: moji
  character(LEN=16) :: moji2

  integer :: itmp
  real(kind=8) :: tmp
  real(kind=8) :: time
  integer :: icount

  real(kind=8) :: x,y,z,dx,dy,dz
  integer :: it  ! it-th time step
  complex(kind=8),allocatable :: calE(:,:)  ! calE_NM
  complex(kind=8),allocatable :: calEtr(:,:)  ! transpose(calE)
  complex(kind=8),allocatable :: calE0(:,:)  ! calE_NM
  real(kind=8),allocatable :: tmpcalE(:,:,:)
  integer :: i,j,k,l
  integer :: rr,pp,qq
  integer :: nn,mm

  complex(kind=8) :: density_Qmat, rho_Qmat
  complex(kind=8) :: t_Qmat, t_Arad_Qmat, s_Qmat, zeta_Qmat, chiral_Qmat
  complex(kind=8) :: j_Qmat
  complex(kind=8) :: Arad_vec

  complex(kind=8) :: sum
  complex(kind=8) :: sumd, sum_d, sumc
  complex(kind=8) :: sumj(3)
  complex(kind=8) :: sumt(3), sums(3), sumz(3)
  complex(kind=8) :: sumtA(3),sumA(3)
  complex(kind=8),allocatable :: sum0d(:,:,:), sum0c(:,:,:), sum0t(:,:,:,:), sum0s(:,:,:,:), sum0z(:,:,:,:)
  complex(kind=8),allocatable :: sum0j(:,:,:,:)
  complex(kind=8),allocatable :: sum0tA(:,:,:,:), sum0A(:,:,:,:)
  complex(kind=8),allocatable :: tmp_d(:,:,:,:,:), tmp_c(:,:,:,:,:), tmp_t(:,:,:,:,:,:), tmp_s(:,:,:,:,:,:), tmp_z(:,:,:,:,:,:)
  complex(kind=8),allocatable :: tmp_j(:,:,:,:,:,:)
  complex(kind=8) :: pos_e,pos_etot
  complex(kind=8) :: pos_e0

  complex(kind=8) :: tmp_Arad_vec(3)
  complex(kind=8),allocatable :: tmp2_j_Qmat(:,:,:)
  complex(kind=8),allocatable :: tmp_t_Arad_Qmat(:,:,:)

  complex(kind=8) :: totA(3),tott(3),totz(3),totAtz(3),totds(3)

  integer :: meshx, meshy, meshz
  real(kind=8) :: xori, xend, yori, yend, zori, zend

  !----------------------------------------------------
  !   parameters to be read from inifile
  !----------------------------------------------------
  character(LEN=80) :: outfile 
  integer :: NTIME
  integer :: nprint ! print out every nprint step (it)
  integer :: NEL
  !----------------------------------------------------

  character(LEN=80) :: inifile  ! input file which is given in the command line.

  logical :: use_BO_approx    ! .true. --> BO_approximation, .false. --> include nucleus degree of freedom
  logical :: use_exchange    ! .true. --> include exchange terms in diff. eq. , .false. --> without exchange terms
  logical :: use_coherent     ! .true. --> use coherent state for photon, .false. --> no radiation in the initial state.
  integer :: n_coherent ! number of coherent states which are not zero.
  integer :: i_coh ! i_coh-th coherent state is not zero.
  real(kind=dp) :: p0,th,phi ! photon momentum vector in spherical coord.
  character(LEN=1) :: sig ! polarization (+ or -)  
  integer :: ip0,ith,iphi
  real(kind=dp) :: dp0,dth,dphi ! mesh width
  complex(kind=dp),allocatable :: alpha(:)   ! <coherent|hat(a)|coherent>

  integer :: filenum
  character(LEN=80) :: FMT

  call set_GammaMatrix

  !================================================================
  ! set computation parameters
  ! set expansion functions for electrons & nucleus
  ! set photon parameters
  !================================================================
  !---------------------------------------------------------------------------------------------
  ! read computation parameters from the file which is specified in the command line.
  !---------------------------------------------------------------------------------------------
  call getarg(1,inifile)
  write(*,*) "# Read ", trim(inifile)
  open(unit=1000,file=trim(inifile),status='unknown',form='formatted')
  read(1000,*) outfile
  open(unit=10,file=outfile,status='unknown',form='formatted')  ! output file
  read(1000,*) NTIME
  read(1000,*) DeltaT   ! defined as global parameter
  read(1000,*)  !nprint
  read(1000,*) NEL  !  NEL = 1 ! H,H2+
                    !  NEL = 2 ! He,H2
                    !  NEL = 6 ! Li2
  read(1000,*) 
  read(1000,*) there_is_twoele
  !  there_is_twoele = .True. ! true : there is already two electron integrals and have a name defined in module IntegralStorage
  !  there_is_twoele = .false. ! false : compute two electron integrals and store as fort.12
  read(1000,*) there_is_nucele
  !  there_is_nucele = .True. ! true : there is already nucleus-electron 4c integrals and have a name defined in module IntegralStorage
  !  there_is_nucele = .false. ! false : compute nucleus-electron 4c integrals and store as fort.13 (--> nucele.dat)
  read(1000,*) there_is_twonuc
  !  there_is_twonuc = .True. ! true : there is already two nucleus integrals and have a name defined in module IntegralStorage
  !  there_is_twonuc = .false. ! false : compute two nucleus integrals and store as fort.14 (--> twonuc.dat)
  read(1000,*) 
  read(1000,*) use_BO_approx
  ! .true. --> BO_approximation, .false. --> include nucleus degree of freedom
  read(1000,*) 
  read(1000,*) use_exchange
  ! .true. --> include exchange terms in diff. eq. , .false. --> without exchange terms
  read(1000,*) 
  read(1000,*) use_coherent
  ! .true. --> use coherent state for photon, .false. --> no radiation in the initial state.
  if((use_coherent .eq. .true.).and.(use_BO_approx .eq. .false.)) then
     write(*,*) "Coherent state is only supported for BO approximation."
     stop
  end if
  if(use_coherent) then
     read(1000,*) P0MAX  ! global
     read(1000,*) NP0    ! global 
     read(1000,*) NPTH   ! global
     read(1000,*) NPPHI  ! global
     Nph = NP0*NPTH*NPPHI*2 ! j=1~Nph, global
     allocate(alpha(Nph))
     alpha(:) = (0._dp,0._dp)
     read(1000,*) n_coherent
     if(n_coherent.gt.Nph) then
        write(*,*) "n_coherent should be smaller than Nph."
        stop
     end if
     do i=1,n_coherent
        read(1000,*) i_coh
        read(1000,*) alpha(i_coh)
     end do
  end if

 
  write(*,*) "#############################################"
  if(use_BO_approx) then
     write(*,*) "# Use BO approximation."
  else  
     write(*,*) "# No BO approximation."
  end if
  if(use_exchange) then
     write(*,*) "# Use exchange terms."
  else  
     write(*,*) "# Without exchange terms."
  end if
  write(*,*) "#############################################"

  write(*,*) "# Write outputs to ", trim(outfile)
  write(*,"(1a15,1i16)") "# NTIME : ", NTIME
  write(*,"(1a15,1es16.6)") "# DeltaT : ", DeltaT
  write(*,"(1a15,1es16.6)") "# t_end : ", NTIME*DeltaT

  !---------------------------------------------------------------------------------------------
  ! set electron related parameters
  !---------------------------------------------------------------------------------------------
  write(*,*) "#############################################"
  write(*,*) "Set expansion functions for electrons."
  write(*,*) "#############################################"

  call read_qedinp ! read qed.inp
  call read_DiracOutput  ! read and set global variables -> defined in sub_readDiracOutput.f90

  if(mod(NEL,2).eq.0) then
     NOCC = NEL/2  ! closed shell (KR)
  else
     NOCC = (NEL+1)/2 ! open shell (KR)
     ! The occupation number of the NOCCth orbital is 1.
  end if

  write(*,"(1a15,1i6)") "# NEL : ",NEL
  write(*,"(1a15,1i6)") "# NBS : ",NBS
  write(*,"(1a15,1i6)") "# NOCC : ",NOCC

  write(*,"(1a15,1i6)") "# NAT : ",NAT
  do i=1,NAT
     write(*,"(1a15,1i6,3es16.6)") "# cn & xyzc: ",cn(i),xc(i),yc(i),zc(i)
  end do
!  stop
  if(use_coherent) then
     !---------------------------------------------------------------------------------------------
     ! set photon parameters & specify coherent states
     !---------------------------------------------------------------------------------------------
     write(*,*) "#############################################"
     write(*,*) "Set photon parameters."
     write(10,*) "# Set photon parameters."
     write(*,*) "#############################################"
     allocate(ALPHA_COH(Nph))
     ALPHA_COH(:) = alpha(:) 

     write(*,"(1a15,1es16.6)") "# P0MAX: ",P0MAX
     write(10,"(1a15,1es16.6)") "# P0MAX: ",P0MAX
     write(*,"(1a15,1i6)") "# NP0: ",NP0
     write(10,"(1a15,1i6)") "# NP0: ",NP0
     write(*,"(1a15,1i6)") "# NPTH: ",NPTH
     write(10,"(1a15,1i6)") "# NPTH: ",NPTH
     write(*,"(1a15,1i6)") "# NPPHI: ",NPPHI
     write(10,"(1a15,1i6)") "# NPPHI: ",NPPHI
     write(*,"(1a15,1i6)") "# Nph: ",Nph
     write(10,"(1a15,1i6)") "# Nph: ",Nph
     
     write(*,*) "# Photon momentum discretization"
     write(10,*) "# Photon momentum discretization"
     write(*,*) "# Polar coordinates: th measured from z-axis, include both 0 and PI."
     write(*,*) "#                  : phi measured from x-axis, include include 0 but not 2PI."
     write(10,*) "# Polar coordinates: th measured from z-axis, include both 0 and PI."
     write(10,*) "#                  : phi measured from x-axis, include include 0 but not 2PI."

     dp0  = P0MAX/NP0 ! P0MAX specified in GammaMatrix. We do not use p0=0
     dth  = PI/(NPTH-1) ! NPTH specified in GammaMatrix. th include both 0 and PI
     dphi = 2.d0*PI/NPPHI ! NPPHI specified in GammaMatrix. phi include 0 but not 2PI

     write(*,"(1a3,1a6,1a5,1a16,2a16)") "#  "," j ","pol"," p0      "," th/PI ", " phi/PI "
     write(10,"(1a3,1a6,1a5,1a16,2a16)") "#  "," j ","pol"," p0      "," th/PI ", " phi/PI "
     do j=1,Nph
        call convert_label_2(j,NP0,NPTH,NPPHI,ip0,ith,iphi,sig)  ! defined in sub_int.f90
        p0 = dp0*ip0
        th = dth*(ith-1)  ! start from 0
        phi= dphi*(iphi-1) ! start from 0
        write(*,"(1a3,1i6,1a5,1es16.6,2f16.4)") "# ", j, sig, p0, th/PI, phi/PI
        write(10,"(1a3,1i6,1a5,1es16.6,2f16.4)") "# ", j, sig, p0, th/PI, phi/PI
     end do
     
     write(*,*) "# Print non-zero alpha"
     write(10,*) "# Print non-zero alpha"
     j=0
     do i=1,Nph
        if(ALPHA_COH(i).ne.(0._dp,0._dp)) then
           write(*,"(1a8,1i4,2es16.6)") "# alpha",i,ALPHA_COH(i)
           write(10,"(1a8,1i4,2es16.6)") "# alpha",i,ALPHA_COH(i)
           j=j+1
        end if
     end do
     if(j.ne.n_coherent) then
        write(*,*) "n_coherent does not match. Check input file."
        stop
     end if
  end if
!!$  call set_ALPHA_COH  ! set eigenvalue of coherent state


!!$  write(filename_E,'(a)') 'calE.dat'
!!$  write(*,'(a)') filename_E

!!$  open(unit=1200,file=filename_E)
!!$
!!$  !---check parameters--------------------------------------
!!$  write(*,*)'start check parameters'
!!$  read(1200,"(1a10,1es16.6)") moji,DeltaT
!!$  read(1200,"(1a10,1i16)") moji, itmp
!!$  if(NBS.ne.itmp) then
!!$     write(*,*)'error NBS'
!!$     stop
!!$  end if
!!$  read(1200,"(1a16,2es16.6)") moji2, ALPHA_COH(3) 
!!$  read(1200,"(1a16,2es16.6)") moji2, ALPHA_COH(27) 
!!$  read(1200,"(1a16,1es16.6)") moji2,tmp
!!$  if(P0MAX.ne.tmp) then
!!$     write(*,*)'error P0MAX'
!!$     stop
!!$  end if
!!$  read(1200,"(1a10,1i6)") moji,NEL
!!$  
!!$  write(*,"(1a10,1es16.6)")'# DeltaT=',DeltaT
!!$  write(*,"(1a10,1i16)")'# NBS=', NBS 
!!$  write(*,"(1a16,2es16.6)")'# ALPHA_COH(3)=', ALPHA_COH(3) 
!!$  write(*,"(1a16,2es16.6)")'# ALPHA_COH(27)=', ALPHA_COH(27) 
!!$  write(*,"(1a16,1es16.6)") "# P0MAX: ",P0MAX
!!$  write(*,"(1a10,1i6)") "# NEL : ",NEL
!!$  if(mod(NEL,2).eq.0) then
!!$     NOCC = NEL/2  ! closed shell (KR)
!!$  else
!!$     NOCC = (NEL+1)/2 ! open shell (KR)
!!$     ! The occupation number of the NOCCth orbital is 1.
!!$  end if
!!$  write(*,"(1a10,1i6)") "# NOCC : ",NOCC
  !-end check parameters---------------------------------------------------------

  allocate(calE(4*NBS,4*NBS))
!!$  allocate(calEtr(4*NBS,4*NBS))
  allocate(calE0(4*NBS,4*NBS))
!!$  allocate(tmpcalE(2,4*NBS,4*NBS))
  !------------------------------------------------
  ! Initial condition for calE
  ! Kronecker delta for electron occupied.
  calE0(:,:) = (0.d0,0.d0)
  do i=1,NEL
     calE0(i,i) = (1.d0,0.d0)
  end do
  !------------------------------------------------
  !------------------------------------------------
  ! Set mesh
!  meshx = 1
!  meshy = 51
!  meshz = 51
!  meshx = 29
  meshx = 21
  meshy = 21
  meshz = 21
!  xori = -4.d0 ; xend = 10.d0
  xori = -4.d0 ; xend = 6.d0
  yori = -4.d0 ; yend = 6.d0
  zori = -4.d0 ; zend = 6.d0
!  xori = -4.d0 ; xend = 4.d0
!  yori = -4.d0 ; yend = 4.d0
!  zori = -4.d0 ; zend = 4.d0
  if(meshx.ne.1) then
     dx = (xend - xori)/(meshx-1)
  else
     dx = 0.d0
  end if
  if(meshy.ne.1) then
     dy = (yend - yori)/(meshy-1)
  else
     dy = 0.d0
  end if
  if(meshz.ne.1) then
     dz = (zend - zori)/(meshz-1)
  else
     dz = 0.d0
  end if

  allocate(sum0d(meshx,meshy,meshz),sum0t(3,meshx,meshy,meshz), sum0s(3,meshx,meshy,meshz), sum0z(3,meshx,meshy,meshz))
  allocate(sum0A(3,meshx,meshy,meshz), sum0tA(3,meshx,meshy,meshz))
  allocate(sum0c(meshx,meshy,meshz))
  allocate(sum0j(3,meshx,meshy,meshz))
  allocate(tmp_d(meshx,meshy,meshz,4*NBS,4*NBS),tmp_c(meshx,meshy,meshz,4*NBS,4*NBS))
  allocate(tmp_s(3,meshx,meshy,meshz,4*NBS,4*NBS),tmp_t(3,meshx,meshy,meshz,4*NBS,4*NBS),tmp_z(3,meshx,meshy,meshz,4*NBS,4*NBS))
  allocate(tmp_j(3,meshx,meshy,meshz,4*NBS,4*NBS))
  allocate(tmp2_j_Qmat(3,4*NBS,4*NBS),tmp_t_Arad_Qmat(3,4*NBS,4*NBS))

  write(*,*)'meshx =',meshx
  write(*,*)'meshy =',meshy
  write(*,*)'meshz =',meshz
  write(*,*)'dx dy dz =',dx,dy,dz
  !------------------------------------------------

  !------set property matrices---------------------
  write(*,*)'set property matrices'
     do i=1,meshx
        if(meshx.ne.1) then
                x = xori + dx*(i-1)
        else
                x = 0.d0
        end if

        do j=1,meshy
           if(meshy.ne.1) then
                y = yori + dy*(j-1)
           else
                y = 0.d0
           end if
           
           do k=1,meshz
              if(meshz.ne.1) then
                z = zori + dz*(k-1)
              else
                z = 0.d0
              end if

              !$omp parallel do
              do nn=1,4*NBS
                 do mm=1,4*NBS
                    tmp_d(i,j,k,nn,mm) = density_Qmat(x,y,z,nn,mm)
                    tmp_c(i,j,k,nn,mm) = chiral_Qmat(x,y,z,nn,mm)
                 end do
              end do
              !$omp end parallel do

              do l=1,3
                 !$omp parallel do
                 do nn=1,4*NBS
                    do mm=1,4*NBS
                       tmp_j(l,i,j,k,nn,mm) = j_Qmat(l,x,y,z,nn,mm)
                       tmp_t(l,i,j,k,nn,mm) = t_Qmat(l,x,y,z,nn,mm)
                       tmp_s(l,i,j,k,nn,mm) = s_Qmat(l,x,y,z,nn,mm)
                       tmp_z(l,i,j,k,nn,mm) = zeta_Qmat(l,x,y,z,nn,mm)
!                       dsums(l) = dsums(l) +s_Qmat(l,x,y,z,nn,mm)*dcalE(nn,mm)/DeltaT
                    end do
                 end do
                 !$omp end parallel do
              end do
           end do !k=1,meshz
        end do !j=1,meshy
        write(*,*) i,'/',meshx
     end do !i=1,meshx
   !----end set property matrices---------------------


  !------------------------------------------------
  write(filename_E,'(a)') '../calE.dat'
  write(*,'(a)') filename_E

  open(unit=1200,file=filename_E)
  write(FMT,'("(1es25.15,"i0"(2es25.15))")') (4*NBS)**2
     read(1200,*)
     read(1200,*)
     read(1200,*)

  do icount = 0,filenum*nprint
     read(1200,FMT) time,calE 
!!$  !----------------------------------------------------------
!!$  ! set initial condition for density matrices
!!$  !----------------------------------------------------------
!!$
!!$  ! Initial condition for calE
!!$  ! Kronecker delta for electron occupied.
!!$  calE(:,:) = (0._dp,0._dp)
!!$  do i=1,NEL
!!$     calE(i,i) = (1._dp,0._dp)
!!$  end do
!!$  !----------------------------------------------------------
!     write(*,*) time
!     read(1200,*) tmp, calE(:,:)
!!$     read(1200,*) tmp, tmpcalE(:,:,:)
     if((icount==0).or.(mod(icount,nprint)==0)) then
!     write(*,*)'icount = ',icount
!!$        time = tmp
        it = nint(time/DeltaT)
!!$        calEtr(:,:) = tmpcalE(1,:,:) + IU*tmpcalE(2,:,:)
!!$        calE = transpose(calEtr)
  
  !-----------------------for BO approximation--------------------------------------------
  ! calculate initial properties
     if(icount==0) then
            write(filedens,'(a)') 'dens.dat'
            write(filetorque,'(a)') 'torq.dat'
            write(filetA,'(a)') 'tA.dat'
            write(fileA,'(a)') 'A.dat'
            write(filespin,'(a)') 'spin.dat'
            write(filezeta,'(a)') 'zeta.dat'
            write(filechiral,'(a)') 'chir.dat'
            write(filej,'(a)') 'j.dat'
            open(unit=17,file=filedens)
            open(unit=31,file=filezeta)
            open(unit=32,file=filetorque)
            open(unit=320,file=filetA)
            open(unit=321,file=fileA)
            open(unit=33,file=filespin)
            open(unit=37,file=filechiral)
            open(unit=38,file=filej)

            sum0d(:,:,:) = (0.d0,0.d0)
            sum0c(:,:,:) = (0.d0,0.d0)
            sum0t(:,:,:,:) = (0.d0,0.d0)
            sum0tA(:,:,:,:) = (0.d0,0.d0)
            sum0A(:,:,:,:) = (0.d0,0.d0)
            sum0s(:,:,:,:) = (0.d0,0.d0)
            sum0z(:,:,:,:) = (0.d0,0.d0)
            sum0j(:,:,:,:) = (0.d0,0.d0)
     do i=1,meshx
        if(meshx.ne.1) then
                x = xori + dx*(i-1)
        else
                x = 0.d0
        end if

        do j=1,meshy
           if(meshy.ne.1) then
                y = yori + dy*(j-1)
           else
                y = 0.d0
           end if
           
           do k=1,meshz
              if(meshz.ne.1) then
                z = zori + dz*(k-1)
              else
                z = 0.d0
              end if

              do nn=1,4*NBS
                 do mm=1,4*NBS
                    sum0d(i,j,k) = sum0d(i,j,k) +tmp_d(i,j,k,nn,mm)*calE(nn,mm)
                    sum0c(i,j,k) = sum0c(i,j,k) +tmp_c(i,j,k,nn,mm)*calE(nn,mm)
!                    sum0d(i,j,k) = sum0d(i,j,k) +density_Qmat(x,y,z,nn,mm)*calE(nn,mm)
!                    sum0c(i,j,k) = sum0c(i,j,k) +chiral_Qmat(x,y,z,nn,mm)*calE(nn,mm)
                 end do
              end do

              do l=1,3
                 do nn=1,4*NBS
                    do mm=1,4*NBS
                       sum0s(l,i,j,k) = sum0s(l,i,j,k) +tmp_s(l,i,j,k,nn,mm)*calE(nn,mm)
                       sum0z(l,i,j,k) = sum0z(l,i,j,k) +tmp_z(l,i,j,k,nn,mm)*calE(nn,mm)
                       sum0j(l,i,j,k) = sum0j(l,i,j,k) +j_Qmat(l,x,y,z,nn,nn)*calE(nn,mm)
!                       sum0s(l,i,j,k) = sum0s(l,i,j,k) +s_Qmat(l,x,y,z,nn,mm)*calE(nn,mm)
!                       sum0z(l,i,j,k) = sum0z(l,i,j,k) +zeta_Qmat(l,x,y,z,nn,mm)*calE(nn,mm)
                    end do
                 end do
              end do

              !$omp parallel do
              do l=1,3
                tmp_Arad_vec(l) = Arad_vec(0.d0,l,x,y,z)
                do nn=1,4*NBS
                  do mm=1,4*NBS
                    tmp2_j_Qmat(l,nn,mm) = tmp_j(l,i,j,k,nn,mm)
                  end do
                end do
              end do
              !$omp end parallel do

              call calc_t_Arad_Qmat(tmp2_j_Qmat,tmp_Arad_vec,tmp_t_Arad_Qmat)

              !$omp parallel do
              do l=1,3
                 do nn=1,4*NBS
                    do mm=1,4*NBS
                       sum0A(l,i,j,k) = sum0A(l,i,j,k) +tmp_t_Arad_Qmat(l,nn,mm)*calE(nn,mm) !spintorque with Arad effect
                       sum0tA(l,i,j,k) = sum0tA(l,i,j,k) +(tmp_t(l,i,j,k,nn,mm)+tmp_t_Arad_Qmat(l,nn,mm))*calE(nn,mm) !spintorque with Arad effect
                       sum0t(l,i,j,k) = sum0t(l,i,j,k) +tmp_t(l,i,j,k,nn,mm)*calE(nn,mm) !spintorque with Arad effect
                    end do
                 end do
              end do
              !$omp end parallel do
!!$!-calc KRM---------------------------------------------------------------------------
!!$                 sum0d(i,j,k) = sum0d(i,j,k) +density_Qmat(x,y,z,2,2)*calE(2,2)
!!$
!!$              do l=1,3
!!$                 sum0t(l,i,j,k) = sum0t(l,i,j,k) +t_Arad_Qmat(0.d0,l,x,y,z,2,2)*calE(2,2) !spintorque with Arad effect
!!$                 sum0s(l,i,j,k) = sum0s(l,i,j,k) +s_Qmat(l,x,y,z,2,2)*calE(2,2)
!!$                 sum0z(l,i,j,k) = sum0z(l,i,j,k) +zeta_Qmat(l,x,y,z,2,2)*calE(2,2)
!!$              end do
!!$!----------------------------------------------------------------------------

                 write(17,'(5es24.14)') x, y, z, sum0d(i,j,k)
                 write(31,'(9es24.14)') x, y, z, sum0z(1,i,j,k), sum0z(2,i,j,k), sum0z(3,i,j,k)
                 write(32,'(9es24.14)') x, y, z, sum0t(1,i,j,k), sum0t(2,i,j,k), sum0t(3,i,j,k)
                 write(320,'(9es24.14)') x, y, z, sum0tA(1,i,j,k), sum0tA(2,i,j,k), sum0tA(3,i,j,k)
                 write(321,'(9es24.14)') x, y, z, sum0A(1,i,j,k), sum0A(2,i,j,k), sum0A(3,i,j,k)
                 write(33,'(9es24.14)') x, y, z, sum0s(1,i,j,k), sum0s(2,i,j,k), sum0s(3,i,j,k)
                 write(37,'(5es24.14)') x, y, z, sum0c(i,j,k)
                 write(38,'(9es24.14)') x, y, z, sum0j(1,i,j,k), sum0j(2,i,j,k), sum0j(3,i,j,k)
           end do
              write(17,*)
              write(31,*)
              write(32,*)
              write(320,*)
              write(321,*)
              write(33,*)
              write(37,*)
              write(38,*)
        end do
            write(17,*)
            write(31,*)
            write(32,*)
            write(320,*)
            write(321,*)
            write(33,*)
            write(37,*)
            write(38,*)
     end do
        close(unit=17)
        close(unit=31)
        close(unit=32)
        close(unit=320)
        close(unit=321)
        close(unit=33)
        close(unit=37)
        close(unit=38)
!  stop
      end if !icount==0
  ! end calculate initial properties
  !----------------------------------------------------------------

            write(filedens,'(a,i10.10,a)') 'denstime',it,'.dat'
            write(filetorque,'(a,i10.10,a)') 'torqtime',it,'.dat'
            write(filetA,'(a,i10.10,a)') 'tAtime',it,'.dat'
            write(fileA,'(a,i10.10,a)') 'Atime',it,'.dat'
            write(filespin,'(a,i10.10,a)') 'spintime',it,'.dat'
!            write(filespind,'(a,i10.10,a)') 'dspintime',it,'.dat'
            write(filezeta,'(a,i10.10,a)') 'zetatime',it,'.dat'
            write(filechiral,'(a,i10.10,a)') 'chirtime',it,'.dat'
            write(filej,'(a,i10.10,a)') 'jtime',it,'.dat'
            open(unit=17,file=filedens)
            open(unit=31,file=filezeta)
            open(unit=32,file=filetorque)
            open(unit=320,file=filetA)
            open(unit=321,file=fileA)
            open(unit=33,file=filespin)
!            open(unit=34,file=filespind)
            open(unit=37,file=filechiral)
            open(unit=38,file=filej)

            totA(:) = (0.d0,0.d0)
            tott(:) = (0.d0,0.d0)
            totz(:) = (0.d0,0.d0)
            totAtz(:) = (0.d0,0.d0)
            totds(:) = (0.d0,0.d0)

     do i=1,meshx
        if(meshx.ne.1) then
                x = xori + dx*(i-1)
        else
                x = 0.d0
        end if

        do j=1,meshy
           if(meshy.ne.1) then
                y = yori + dy*(j-1)
           else
                y = 0.d0
           end if
           
           do k=1,meshz
              if(meshz.ne.1) then
                z = zori + dz*(k-1)
              else
                z = 0.d0
              end if

            sumd = (0.d0,0.d0)
            sumc = (0.d0,0.d0)
            sumt(:) = (0.d0,0.d0)
            sumtA(:) = (0.d0,0.d0)
            sumA(:) = (0.d0,0.d0)
            sums(:) = (0.d0,0.d0)
            sumz(:) = (0.d0,0.d0)
            sumj(:) = (0.d0,0.d0)
!            dsums(:) = (0.d0,0.d0)

              do nn=1,4*NBS
                 do mm=1,4*NBS
                    sumd = sumd +tmp_d(i,j,k,nn,mm)*calE(nn,mm)
                    sumc = sumc +tmp_c(i,j,k,nn,mm)*calE(nn,mm)
!                    sumd = sumd +density_Qmat(x,y,z,nn,mm)*calE(nn,mm)
!                    sumc = sumc +chiral_Qmat(x,y,z,nn,mm)*calE(nn,mm)
                 end do
              end do

              do l=1,3
                 do nn=1,4*NBS
                    do mm=1,4*NBS
                       sums(l) = sums(l) +tmp_s(l,i,j,k,nn,mm)*calE(nn,mm)
                       sumz(l) = sumz(l) +tmp_z(l,i,j,k,nn,mm)*calE(nn,mm)
                       sumj(l) = sumj(l) +tmp_j(l,i,j,k,nn,mm)*calE(nn,mm)
!                       sums(l) = sums(l) +s_Qmat(l,x,y,z,nn,mm)*calE(nn,mm)
!                       sumz(l) = sumz(l) +zeta_Qmat(l,x,y,z,nn,mm)*calE(nn,mm)
!                       dsums(l) = dsums(l) +tmp_s(l,i,j,k,nn,mm)*dcalE(nn,mm)/DeltaT
                    end do
                 end do
              end do

              do l=1,3
                tmp_Arad_vec(l) = Arad_vec(time,l,x,y,z)
                do nn=1,4*NBS
                  do mm=1,4*NBS
                    tmp2_j_Qmat(l,nn,mm) = tmp_j(l,i,j,k,nn,mm)
                  end do
                end do
              end do

              call calc_t_Arad_Qmat(tmp2_j_Qmat,tmp_Arad_vec,tmp_t_Arad_Qmat)

              do l=1,3
                 do nn=1,4*NBS
                    do mm=1,4*NBS
                       sumA(l) = sumA(l) +tmp_t_Arad_Qmat(l,nn,mm)*calE(nn,mm)
                       sumt(l) = sumt(l) +tmp_t(l,i,j,k,nn,mm)*calE(nn,mm) !spintorque with Arad effect
                       sumtA(l) = sumtA(l) +(tmp_t(l,i,j,k,nn,mm)+tmp_t_Arad_Qmat(l,nn,mm))*calE(nn,mm) !spintorque with Arad effect
                    end do
                 end do
              end do
    
            write(17,'(10es24.14)') it*DeltaT, x, y, z, sumd, sumd-sum0d(i,j,k), sum0d(i,j,k)
            write(31,'(22es24.14)') it*DeltaT, x, y, z, sumz(1), sumz(1)-sum0z(1,i,j,k), sumz(2), sumz(2)-sum0z(2,i,j,k)&
                                                  &, sumz(3), sumz(3)-sum0z(3,i,j,k)&
                                                  &, sum0z(1,i,j,k), sum0z(2,i,j,k), sum0z(3,i,j,k)
            write(32,'(22es24.14)') it*DeltaT, x, y, z, sumt(1), sumt(1)-sum0t(1,i,j,k), sumt(2), sumt(2)-sum0t(2,i,j,k)&
                                                  &, sumt(3), sumt(3)-sum0t(3,i,j,k)&
                                                  &, sum0t(1,i,j,k), sum0t(2,i,j,k), sum0t(3,i,j,k)
            write(320,'(22es24.14)') it*DeltaT, x, y, z, sumtA(1), sumtA(1)-sum0tA(1,i,j,k), sumtA(2), sumtA(2)-sum0tA(2,i,j,k)&
                                                  &, sumtA(3), sumtA(3)-sum0tA(3,i,j,k)&
                                                  &, sum0tA(1,i,j,k), sum0tA(2,i,j,k), sum0tA(3,i,j,k)
            write(321,'(22es24.14)') it*DeltaT, x, y, z, sumA(1), sumA(1)-sum0A(1,i,j,k), sumA(2), sumA(2)-sum0A(2,i,j,k)&
                                                  &, sumA(3), sumA(3)-sum0A(3,i,j,k)&
                                                  &, sum0A(1,i,j,k), sum0A(2,i,j,k), sum0A(3,i,j,k)
            write(33,'(22es24.14)') it*DeltaT, x, y, z, sums(1), sums(1)-sum0s(1,i,j,k), sums(2), sums(2)-sum0s(2,i,j,k)&
                                                  &, sums(3), sums(3)-sum0s(3,i,j,k)&
                                                  &, sum0s(1,i,j,k), sum0s(2,i,j,k), sum0s(3,i,j,k)
!            write(34,'(10es24.14)') it*DeltaT, x, y, z, dsums(1), dsums(2), dsums(3)
            write(37,'(10es24.14)') it*DeltaT, x, y, z, sumc, sumc-sum0c(i,j,k), sum0c(i,j,k)
            write(38,'(22es24.14)') it*DeltaT, x, y, z, sumj(1), sumj(1)-sum0j(1,i,j,k), sumj(2), sumj(2)-sum0j(2,i,j,k)&
                                                  &, sumj(3), sumj(3)-sum0j(3,i,j,k)&
                                                  &, sum0j(1,i,j,k), sum0j(2,i,j,k), sum0j(3,i,j,k)

           do l=1,3
              totA(l) = totA(l) + sumA(l)*dx*dy*dz
              tott(l) = tott(l) + sumt(l)*dx*dy*dz
              totz(l) = totz(l) + sumz(l)*dx*dy*dz
              totAtz(l) = totAtz(l) + abs(sumA(l)+sumt(l)+sumz(l))*dx*dy*dz
              totds(l) = totds(l) + (sums(l) - sum0s(l,i,j,k))*dx*dy*dz
           end do
    

         end do
         write(17,*)
         write(31,*)
         write(32,*)
         write(320,*)
         write(321,*)
         write(33,*)
!         write(34,*)
         write(37,*)
         write(38,*)
       end do
       write(17,*)
       write(31,*)
       write(32,*)
       write(320,*)
       write(321,*)
       write(33,*)
!       write(34,*)
       write(37,*)
       write(38,*)
     end do

     write(1111,'(a4,7es24.14)') 'totA', time, totA(1), totA(2), totA(3)
     write(1112,'(a4,7es24.14)') 'tott', time, tott(1), tott(2), tott(3)
     write(1113,'(a4,7es24.14)') 'totz', time, totz(1), totz(2), totz(3)
     write(1114,'(a6,7es24.14)') 'tottzA', time, tott(1)+totz(1)+totA(1), tott(2)+totz(2)+totA(2), tott(3)+totz(3)+totA(3)
     write(1115,'(a9,7es24.14)') 'abstotAtz', time, totAtz(1), totAtz(2), totAtz(3)
     write(1116,'(a5,7es24.14)') 'totds', time, totds(1), totds(2), totds(3)
        close(unit=17)
        close(unit=31)
        close(unit=32)
        close(unit=320)
        close(unit=321)
        close(unit=33)
!        close(unit=34)
        close(unit=37)
        close(unit=38)
  !-----------------------end 'for BO approximation'--------------------------------------------


!!$  !----------------------- for only density --------------------------------------------
!!$  !------initial state---------------------------
!!$     write(fileele,'(a,i10.10,a)') 'eledenstime',it,'.dat'
!!$
!!$     open(unit=35,file=fileele)
!!$
!!$     do i=1,meshx
!!$        if(meshx.ne.1) then
!!$                x = xori + dx*(i-1)
!!$        else
!!$                x = 0.d0
!!$        end if
!!$
!!$        do j=1,meshy
!!$           if(meshy.ne.1) then
!!$                y = yori + dy*(j-1)
!!$           else
!!$                y = 0.d0
!!$           end if
!!$           
!!$           do k=1,meshz
!!$              if(meshz.ne.1) then
!!$                z = zori + dz*(k-1)
!!$              else
!!$                z = 0.d0
!!$              end if
!!$
!!$                sum = (0.d0,0.d0)
!!$                do nn=1,4*NBS
!!$                   do mm=1,4*NBS
!!$                      sum = sum +density_Qmat(x,y,z,nn,mm)*calE0(nn,mm)
!!$                   end do
!!$                end do
!!$                pos_e0 = sum
!!$
!!$                sum = (0.d0,0.d0)
!!$                do nn=1,4*NBS
!!$                   do mm=1,4*NBS
!!$                      sum = sum +density_Qmat(x,y,z,nn,mm)*calE(nn,mm)
!!$                   end do
!!$                end do
!!$                pos_e = sum
!!$           
!!$            write(35,'(10es16.6)') time, x, y, z, pos_e, pos_e-pos_e0
!!$
!!$         end do
!!$         write(35,*)
!!$       end do
!!$       write(35,*)
!!$     end do
!!$        close(unit=35)

     end if !end if(mod(icount,nprint)==0)
  end do

  deallocate(cn,xc,yc,zc)  ! global
  deallocate(aa_L,xx_L,yy_L,zz_L,nx_L,ny_L,nz_L)
  deallocate(aa_S,xx_S,yy_S,zz_S,nx_S,ny_S,nz_S)
  deallocate(c_La,c_Lb,c_Sa,c_Sb)
  deallocate(d_La,d_Lb,d_Sa,d_Sb)

  deallocate(calE,calEtr,calE0,tmpcalE)

  deallocate(sum0d,sum0t,sum0s, sum0z)
  deallocate(sum0c)
  deallocate(tmp_d,tmp_c)
  deallocate(tmp_s,tmp_t,tmp_z)
  deallocate(tmp_j)

end subroutine calc_property_BO

!====================================================================================
subroutine calc_property_BO_point(nprint,filenum,x0,y0,z0)
!====================================================================================
  use Precision
  Use DiracOutput
  use Constants
  use IntegralStorage
  use NucBasis

  implicit none

  character(LEN=10) :: str1,str2,str3,str4
  character(LEN=80) :: filedens, filetorque, filespin, filespind, filezeta, filechiral
  character(LEN=80) :: fileA, filetA
  character(LEN=80) :: fileele,filenuc
  character(LEN=80) :: filename_E,filename_C
  character(LEN=10) :: moji
  character(LEN=16) :: moji2

  integer :: itmp
  real(kind=8) :: tmp
  real(kind=8) :: time
  integer :: icount

  real(kind=8) :: x,y,z,dx,dy,dz
  integer :: it  ! it-th time step
  complex(kind=8),allocatable :: calE(:,:)  ! calE_NM
  complex(kind=8),allocatable :: calEtr(:,:)  ! transpose(calE)
  complex(kind=8),allocatable :: calE0(:,:)  ! calE_NM
  real(kind=8),allocatable :: tmpcalE(:,:,:)
  integer :: i,j,k,l
  integer :: rr,pp,qq
  integer :: nn,mm

  complex(kind=8) :: density_Qmat, rho_Qmat
  complex(kind=8) :: t_Qmat, t_Arad_Qmat, s_Qmat, zeta_Qmat, chiral_Qmat
  complex(kind=8) :: j_Qmat
  complex(kind=8) :: Arad_vec

  complex(kind=8) :: sum
  complex(kind=8) :: sumd, sum_d, sumc
  complex(kind=8) :: sumt(3), sums(3), sumz(3)
  complex(kind=8) :: sumtA(3),sumA(3)
  complex(kind=8),allocatable :: sum0d(:,:,:), sum0c(:,:,:), sum0t(:,:,:,:), sum0s(:,:,:,:), sum0z(:,:,:,:)
  complex(kind=8),allocatable :: sum0tA(:,:,:,:), sum0A(:,:,:,:)
  complex(kind=8),allocatable :: tmp_d(:,:,:,:,:), tmp_c(:,:,:,:,:), tmp_t(:,:,:,:,:,:), tmp_s(:,:,:,:,:,:), tmp_z(:,:,:,:,:,:)
  complex(kind=8),allocatable :: tmp_j(:,:,:,:,:,:)
  complex(kind=8) :: pos_e,pos_etot
  complex(kind=8) :: pos_e0

  complex(kind=8) :: tmp_Arad_vec(3)
  complex(kind=8),allocatable :: tmp2_j_Qmat(:,:,:)
  complex(kind=8),allocatable :: tmp_t_Arad_Qmat(:,:,:)

  integer :: meshx, meshy, meshz
  real(kind=8) :: xori, xend, yori, yend, zori, zend
  real(kind=8),intent(in) :: x0
  real(kind=8),intent(in) :: y0
  real(kind=8),intent(in) :: z0

  !----------------------------------------------------
  !   parameters to be read from inifile
  !----------------------------------------------------
  character(LEN=80) :: outfile 
  integer :: NTIME
  integer :: nprint ! print out every nprint step (it)
  integer :: NEL
  !----------------------------------------------------

  character(LEN=80) :: inifile  ! input file which is given in the command line.

  logical :: use_BO_approx    ! .true. --> BO_approximation, .false. --> include nucleus degree of freedom
  logical :: use_exchange    ! .true. --> include exchange terms in diff. eq. , .false. --> without exchange terms
  logical :: use_coherent     ! .true. --> use coherent state for photon, .false. --> no radiation in the initial state.
  integer :: n_coherent ! number of coherent states which are not zero.
  integer :: i_coh ! i_coh-th coherent state is not zero.
  real(kind=dp) :: p0,th,phi ! photon momentum vector in spherical coord.
  character(LEN=1) :: sig ! polarization (+ or -)  
  integer :: ip0,ith,iphi
  real(kind=dp) :: dp0,dth,dphi ! mesh width
  complex(kind=dp),allocatable :: alpha(:)   ! <coherent|hat(a)|coherent>

  integer :: filenum
  character(LEN=80) :: FMT

  call set_GammaMatrix

  !================================================================
  ! set computation parameters
  ! set expansion functions for electrons & nucleus
  ! set photon parameters
  !================================================================
  !---------------------------------------------------------------------------------------------
  ! read computation parameters from the file which is specified in the command line.
  !---------------------------------------------------------------------------------------------
  call getarg(1,inifile)
  write(*,*) "# Read ", trim(inifile)
  open(unit=1000,file=trim(inifile),status='unknown',form='formatted')
  read(1000,*) outfile
  open(unit=10,file=outfile,status='unknown',form='formatted')  ! output file
  read(1000,*) NTIME
  read(1000,*) DeltaT   ! defined as global parameter
  read(1000,*)  !nprint
  read(1000,*) NEL  !  NEL = 1 ! H,H2+
                    !  NEL = 2 ! He,H2
                    !  NEL = 6 ! Li2
  read(1000,*) 
  read(1000,*) there_is_twoele
  !  there_is_twoele = .True. ! true : there is already two electron integrals and have a name defined in module IntegralStorage
  !  there_is_twoele = .false. ! false : compute two electron integrals and store as fort.12
  read(1000,*) there_is_nucele
  !  there_is_nucele = .True. ! true : there is already nucleus-electron 4c integrals and have a name defined in module IntegralStorage
  !  there_is_nucele = .false. ! false : compute nucleus-electron 4c integrals and store as fort.13 (--> nucele.dat)
  read(1000,*) there_is_twonuc
  !  there_is_twonuc = .True. ! true : there is already two nucleus integrals and have a name defined in module IntegralStorage
  !  there_is_twonuc = .false. ! false : compute two nucleus integrals and store as fort.14 (--> twonuc.dat)
  read(1000,*) 
  read(1000,*) use_BO_approx
  ! .true. --> BO_approximation, .false. --> include nucleus degree of freedom
  read(1000,*) 
  read(1000,*) use_exchange
  ! .true. --> include exchange terms in diff. eq. , .false. --> without exchange terms
  read(1000,*) 
  read(1000,*) use_coherent
  ! .true. --> use coherent state for photon, .false. --> no radiation in the initial state.
  if((use_coherent .eq. .true.).and.(use_BO_approx .eq. .false.)) then
     write(*,*) "Coherent state is only supported for BO approximation."
     stop
  end if
  if(use_coherent) then
     read(1000,*) P0MAX  ! global
     read(1000,*) NP0    ! global 
     read(1000,*) NPTH   ! global
     read(1000,*) NPPHI  ! global
     Nph = NP0*NPTH*NPPHI*2 ! j=1~Nph, global
     allocate(alpha(Nph))
     alpha(:) = (0._dp,0._dp)
     read(1000,*) n_coherent
     if(n_coherent.gt.Nph) then
        write(*,*) "n_coherent should be smaller than Nph."
        stop
     end if
     do i=1,n_coherent
        read(1000,*) i_coh
        read(1000,*) alpha(i_coh)
     end do
  end if

 
  write(*,*) "#############################################"
  if(use_BO_approx) then
     write(*,*) "# Use BO approximation."
  else  
     write(*,*) "# No BO approximation."
  end if
  if(use_exchange) then
     write(*,*) "# Use exchange terms."
  else  
     write(*,*) "# Without exchange terms."
  end if
  write(*,*) "#############################################"

  write(*,*) "# Write outputs to ", trim(outfile)
  write(*,"(1a15,1i16)") "# NTIME : ", NTIME
  write(*,"(1a15,1es16.6)") "# DeltaT : ", DeltaT
  write(*,"(1a15,1es16.6)") "# t_end : ", NTIME*DeltaT

  !---------------------------------------------------------------------------------------------
  ! set electron related parameters
  !---------------------------------------------------------------------------------------------
  write(*,*) "#############################################"
  write(*,*) "Set expansion functions for electrons."
  write(*,*) "#############################################"

  call read_qedinp ! read qed.inp
  call read_DiracOutput  ! read and set global variables -> defined in sub_readDiracOutput.f90

  if(mod(NEL,2).eq.0) then
     NOCC = NEL/2  ! closed shell (KR)
  else
     NOCC = (NEL+1)/2 ! open shell (KR)
     ! The occupation number of the NOCCth orbital is 1.
  end if

  write(*,"(1a15,1i6)") "# NEL : ",NEL
  write(*,"(1a15,1i6)") "# NBS : ",NBS
  write(*,"(1a15,1i6)") "# NOCC : ",NOCC

  write(*,"(1a15,1i6)") "# NAT : ",NAT
  do i=1,NAT
     write(*,"(1a15,1i6,3es16.6)") "# cn & xyzc: ",cn(i),xc(i),yc(i),zc(i)
  end do
!  stop
  if(use_coherent) then
     !---------------------------------------------------------------------------------------------
     ! set photon parameters & specify coherent states
     !---------------------------------------------------------------------------------------------
     write(*,*) "#############################################"
     write(*,*) "Set photon parameters."
     write(10,*) "# Set photon parameters."
     write(*,*) "#############################################"
     allocate(ALPHA_COH(Nph))
     ALPHA_COH(:) = alpha(:) 

     write(*,"(1a15,1es16.6)") "# P0MAX: ",P0MAX
     write(10,"(1a15,1es16.6)") "# P0MAX: ",P0MAX
     write(*,"(1a15,1i6)") "# NP0: ",NP0
     write(10,"(1a15,1i6)") "# NP0: ",NP0
     write(*,"(1a15,1i6)") "# NPTH: ",NPTH
     write(10,"(1a15,1i6)") "# NPTH: ",NPTH
     write(*,"(1a15,1i6)") "# NPPHI: ",NPPHI
     write(10,"(1a15,1i6)") "# NPPHI: ",NPPHI
     write(*,"(1a15,1i6)") "# Nph: ",Nph
     write(10,"(1a15,1i6)") "# Nph: ",Nph
     
     write(*,*) "# Photon momentum discretization"
     write(10,*) "# Photon momentum discretization"
     write(*,*) "# Polar coordinates: th measured from z-axis, include both 0 and PI."
     write(*,*) "#                  : phi measured from x-axis, include include 0 but not 2PI."
     write(10,*) "# Polar coordinates: th measured from z-axis, include both 0 and PI."
     write(10,*) "#                  : phi measured from x-axis, include include 0 but not 2PI."

     dp0  = P0MAX/NP0 ! P0MAX specified in GammaMatrix. We do not use p0=0
     dth  = PI/(NPTH-1) ! NPTH specified in GammaMatrix. th include both 0 and PI
     dphi = 2.d0*PI/NPPHI ! NPPHI specified in GammaMatrix. phi include 0 but not 2PI

     write(*,"(1a3,1a6,1a5,1a16,2a16)") "#  "," j ","pol"," p0      "," th/PI ", " phi/PI "
     write(10,"(1a3,1a6,1a5,1a16,2a16)") "#  "," j ","pol"," p0      "," th/PI ", " phi/PI "
     do j=1,Nph
        call convert_label_2(j,NP0,NPTH,NPPHI,ip0,ith,iphi,sig)  ! defined in sub_int.f90
        p0 = dp0*ip0
        th = dth*(ith-1)  ! start from 0
        phi= dphi*(iphi-1) ! start from 0
        write(*,"(1a3,1i6,1a5,1es16.6,2f16.4)") "# ", j, sig, p0, th/PI, phi/PI
        write(10,"(1a3,1i6,1a5,1es16.6,2f16.4)") "# ", j, sig, p0, th/PI, phi/PI
     end do
     
     write(*,*) "# Print non-zero alpha"
     write(10,*) "# Print non-zero alpha"
     j=0
     do i=1,Nph
        if(ALPHA_COH(i).ne.(0._dp,0._dp)) then
           write(*,"(1a8,1i4,2es16.6)") "# alpha",i,ALPHA_COH(i)
           write(10,"(1a8,1i4,2es16.6)") "# alpha",i,ALPHA_COH(i)
           j=j+1
        end if
     end do
     if(j.ne.n_coherent) then
        write(*,*) "n_coherent does not match. Check input file."
        stop
     end if
  end if
!!$  call set_ALPHA_COH  ! set eigenvalue of coherent state


!!$  write(filename_E,'(a)') 'calE.dat'
!!$  write(*,'(a)') filename_E

!!$  open(unit=1200,file=filename_E)
!!$
!!$  !---check parameters--------------------------------------
!!$  write(*,*)'start check parameters'
!!$  read(1200,"(1a10,1es16.6)") moji,DeltaT
!!$  read(1200,"(1a10,1i16)") moji, itmp
!!$  if(NBS.ne.itmp) then
!!$     write(*,*)'error NBS'
!!$     stop
!!$  end if
!!$  read(1200,"(1a16,2es16.6)") moji2, ALPHA_COH(3) 
!!$  read(1200,"(1a16,2es16.6)") moji2, ALPHA_COH(27) 
!!$  read(1200,"(1a16,1es16.6)") moji2,tmp
!!$  if(P0MAX.ne.tmp) then
!!$     write(*,*)'error P0MAX'
!!$     stop
!!$  end if
!!$  read(1200,"(1a10,1i6)") moji,NEL
!!$  
!!$  write(*,"(1a10,1es16.6)")'# DeltaT=',DeltaT
!!$  write(*,"(1a10,1i16)")'# NBS=', NBS 
!!$  write(*,"(1a16,2es16.6)")'# ALPHA_COH(3)=', ALPHA_COH(3) 
!!$  write(*,"(1a16,2es16.6)")'# ALPHA_COH(27)=', ALPHA_COH(27) 
!!$  write(*,"(1a16,1es16.6)") "# P0MAX: ",P0MAX
!!$  write(*,"(1a10,1i6)") "# NEL : ",NEL
!!$  if(mod(NEL,2).eq.0) then
!!$     NOCC = NEL/2  ! closed shell (KR)
!!$  else
!!$     NOCC = (NEL+1)/2 ! open shell (KR)
!!$     ! The occupation number of the NOCCth orbital is 1.
!!$  end if
!!$  write(*,"(1a10,1i6)") "# NOCC : ",NOCC
  !-end check parameters---------------------------------------------------------

  allocate(calE(4*NBS,4*NBS))
!!$  allocate(calEtr(4*NBS,4*NBS))
  allocate(calE0(4*NBS,4*NBS))
!!$  allocate(tmpcalE(2,4*NBS,4*NBS))
  !------------------------------------------------
  ! Initial condition for calE
  ! Kronecker delta for electron occupied.
  calE0(:,:) = (0.d0,0.d0)
  do i=1,NEL
     calE0(i,i) = (1.d0,0.d0)
  end do
  !------------------------------------------------
  !------------------------------------------------
  ! Set mesh
!  meshx = 1
!  meshy = 51
!  meshz = 51
  meshx = 1
  meshy = 1
  meshz = 1
  xori = -4.d0 ; xend = 6.d0
  yori = -4.d0 ; yend = 6.d0
  zori = -4.d0 ; zend = 6.d0
!  xori = -4.d0 ; xend = 4.d0
!  yori = -4.d0 ; yend = 4.d0
!  zori = -4.d0 ; zend = 4.d0
  if(meshx.ne.1) then
     dx = (xend - xori)/(meshx-1)
  else
     dx = 0.d0
  end if
  if(meshy.ne.1) then
     dy = (yend - yori)/(meshy-1)
  else
     dy = 0.d0
  end if
  if(meshz.ne.1) then
     dz = (zend - zori)/(meshz-1)
  else
     dz = 0.d0
  end if

  allocate(sum0d(meshx,meshy,meshz),sum0t(3,meshx,meshy,meshz), sum0s(3,meshx,meshy,meshz), sum0z(3,meshx,meshy,meshz))
  allocate(sum0A(3,meshx,meshy,meshz), sum0tA(3,meshx,meshy,meshz))
  allocate(sum0c(meshx,meshy,meshz))
  allocate(tmp_d(meshx,meshy,meshz,4*NBS,4*NBS),tmp_c(meshx,meshy,meshz,4*NBS,4*NBS))
  allocate(tmp_s(3,meshx,meshy,meshz,4*NBS,4*NBS),tmp_t(3,meshx,meshy,meshz,4*NBS,4*NBS),tmp_z(3,meshx,meshy,meshz,4*NBS,4*NBS))
  allocate(tmp_j(3,meshx,meshy,meshz,4*NBS,4*NBS))
  allocate(tmp2_j_Qmat(3,4*NBS,4*NBS),tmp_t_Arad_Qmat(3,4*NBS,4*NBS))

  write(*,*)'meshx =',meshx
  write(*,*)'meshy =',meshy
  write(*,*)'meshz =',meshz
  write(*,*)'dx dy dz =',dx,dy,dz
  !------------------------------------------------

  !------set property matrices---------------------
  write(*,*)'set property matrices'
     do i=1,meshx
        if(meshx.ne.1) then
                x = xori + dx*(i-1)
        else
                x = x0
        end if

        do j=1,meshy
           if(meshy.ne.1) then
                y = yori + dy*(j-1)
           else
                y = y0
           end if
           
           do k=1,meshz
              if(meshz.ne.1) then
                z = zori + dz*(k-1)
              else
                z = z0
              end if

              !$omp parallel do
              do nn=1,4*NBS
                 do mm=1,4*NBS
                    tmp_d(i,j,k,nn,mm) = density_Qmat(x,y,z,nn,mm)
                    tmp_c(i,j,k,nn,mm) = chiral_Qmat(x,y,z,nn,mm)
                 end do
              end do
              !$omp end parallel do

              do l=1,3
                 !$omp parallel do
                 do nn=1,4*NBS
                    do mm=1,4*NBS
                       tmp_j(l,i,j,k,nn,mm) = j_Qmat(l,x,y,z,nn,mm)
                       tmp_t(l,i,j,k,nn,mm) = t_Qmat(l,x,y,z,nn,mm)
                       tmp_s(l,i,j,k,nn,mm) = s_Qmat(l,x,y,z,nn,mm)
                       tmp_z(l,i,j,k,nn,mm) = zeta_Qmat(l,x,y,z,nn,mm)
!                       dsums(l) = dsums(l) +s_Qmat(l,x,y,z,nn,mm)*dcalE(nn,mm)/DeltaT
                    end do
                 end do
                 !$omp end parallel do
              end do
           end do !k=1,meshz
        end do !j=1,meshy
        write(*,*) i,'/',meshx
     end do !i=1,meshx
   !----end set property matrices---------------------


  !------------------------------------------------
  write(filename_E,'(a)') '../calE.dat'
  write(*,'(a)') filename_E

  open(unit=1200,file=filename_E)
  write(FMT,'("(1es25.15,"i0"(2es25.15))")') (4*NBS)**2
     read(1200,*)
     read(1200,*)
     read(1200,*)

  do icount = 0,filenum*nprint
     read(1200,FMT) time,calE 
!!$  !----------------------------------------------------------
!!$  ! set initial condition for density matrices
!!$  !----------------------------------------------------------
!!$
!!$  ! Initial condition for calE
!!$  ! Kronecker delta for electron occupied.
!!$  calE(:,:) = (0._dp,0._dp)
!!$  do i=1,NEL
!!$     calE(i,i) = (1._dp,0._dp)
!!$  end do
!!$  !----------------------------------------------------------
!     read(1200,*) tmp, calE(:,:)
!!$     read(1200,*) tmp, tmpcalE(:,:,:)
!!$     if((icount==0).or.(mod(icount,nprint)==0)) then
!!$!!$        time = tmp
!!$!!$        it = nint(time/DeltaT)
!!$!!$        calEtr(:,:) = tmpcalE(1,:,:) + IU*tmpcalE(2,:,:)
!!$!!$        calE = transpose(calEtr)
!!$  
  !-----------------------for BO approximation--------------------------------------------
  ! calculate initial properties
     if(icount==0) then
!!$            write(filedens,'(a)') 'dens.dat'
!!$            write(filetorque,'(a)') 'torq.dat'
!!$            write(filetA,'(a)') 'tA.dat'
!!$            write(fileA,'(a)') 'A.dat'
!!$            write(filespin,'(a)') 'spin.dat'
!!$            write(filezeta,'(a)') 'zeta.dat'
!!$            write(filechiral,'(a)') 'chir.dat'
!!$            open(unit=17,file=filedens)
!!$            open(unit=31,file=filezeta)
!!$            open(unit=32,file=filetorque)
!!$            open(unit=320,file=filetA)
!!$            open(unit=321,file=fileA)
!!$            open(unit=33,file=filespin)
!!$            open(unit=37,file=filechiral)

            sum0d(:,:,:) = (0.d0,0.d0)
            sum0c(:,:,:) = (0.d0,0.d0)
            sum0t(:,:,:,:) = (0.d0,0.d0)
            sum0tA(:,:,:,:) = (0.d0,0.d0)
            sum0A(:,:,:,:) = (0.d0,0.d0)
            sum0s(:,:,:,:) = (0.d0,0.d0)
            sum0z(:,:,:,:) = (0.d0,0.d0)
     do i=1,meshx
        if(meshx.ne.1) then
                x = xori + dx*(i-1)
        else
                x = x0
        end if

        do j=1,meshy
           if(meshy.ne.1) then
                y = yori + dy*(j-1)
           else
                y = y0
           end if
           
           do k=1,meshz
              if(meshz.ne.1) then
                z = zori + dz*(k-1)
              else
                z = z0
              end if

              do nn=1,4*NBS
                 do mm=1,4*NBS
                    sum0d(i,j,k) = sum0d(i,j,k) +tmp_d(i,j,k,nn,mm)*calE(nn,mm)
                    sum0c(i,j,k) = sum0c(i,j,k) +tmp_c(i,j,k,nn,mm)*calE(nn,mm)
!                    sum0d(i,j,k) = sum0d(i,j,k) +density_Qmat(x,y,z,nn,mm)*calE(nn,mm)
!                    sum0c(i,j,k) = sum0c(i,j,k) +chiral_Qmat(x,y,z,nn,mm)*calE(nn,mm)
                 end do
              end do

              do l=1,3
                 do nn=1,4*NBS
                    do mm=1,4*NBS
                       sum0s(l,i,j,k) = sum0s(l,i,j,k) +tmp_s(l,i,j,k,nn,mm)*calE(nn,mm)
                       sum0z(l,i,j,k) = sum0z(l,i,j,k) +tmp_z(l,i,j,k,nn,mm)*calE(nn,mm)
!                       sum0s(l,i,j,k) = sum0s(l,i,j,k) +s_Qmat(l,x,y,z,nn,mm)*calE(nn,mm)
!                       sum0z(l,i,j,k) = sum0z(l,i,j,k) +zeta_Qmat(l,x,y,z,nn,mm)*calE(nn,mm)
                    end do
                 end do
              end do

              !$omp parallel do
              do l=1,3
                tmp_Arad_vec(l) = Arad_vec(0.d0,l,x,y,z)
                do nn=1,4*NBS
                  do mm=1,4*NBS
                    tmp2_j_Qmat(l,nn,mm) = tmp_j(l,i,j,k,nn,mm)
                  end do
                end do
              end do
              !$omp end parallel do

              call calc_t_Arad_Qmat(tmp2_j_Qmat,tmp_Arad_vec,tmp_t_Arad_Qmat)

              !$omp parallel do
              do l=1,3
                 do nn=1,4*NBS
                    do mm=1,4*NBS
                       sum0A(l,i,j,k) = sum0A(l,i,j,k) +tmp_t_Arad_Qmat(l,nn,mm)*calE(nn,mm) !spintorque with Arad effect
                       sum0tA(l,i,j,k) = sum0tA(l,i,j,k) +(tmp_t(l,i,j,k,nn,mm)+tmp_t_Arad_Qmat(l,nn,mm))*calE(nn,mm) !spintorque with Arad effect
                       sum0t(l,i,j,k) = sum0t(l,i,j,k) +tmp_t(l,i,j,k,nn,mm)*calE(nn,mm) !spintorque with Arad effect
                    end do
                 end do
              end do
              !$omp end parallel do
!!$!!$!-calc KRM---------------------------------------------------------------------------
!!$!!$                 sum0d(i,j,k) = sum0d(i,j,k) +density_Qmat(x,y,z,2,2)*calE(2,2)
!!$!!$
!!$!!$              do l=1,3
!!$!!$                 sum0t(l,i,j,k) = sum0t(l,i,j,k) +t_Arad_Qmat(0.d0,l,x,y,z,2,2)*calE(2,2) !spintorque with Arad effect
!!$!!$                 sum0s(l,i,j,k) = sum0s(l,i,j,k) +s_Qmat(l,x,y,z,2,2)*calE(2,2)
!!$!!$                 sum0z(l,i,j,k) = sum0z(l,i,j,k) +zeta_Qmat(l,x,y,z,2,2)*calE(2,2)
!!$!!$              end do
!!$!!$!----------------------------------------------------------------------------
!!$
!!$!!$                 write(17,'(5es26.16)') x, y, z, sum0d(i,j,k)
!!$!!$                 write(31,'(9es26.16)') x, y, z, sum0z(1,i,j,k), sum0z(2,i,j,k), sum0z(3,i,j,k)
!!$                 write(32,'(9es26.16)') x, y, z, sum0t(1,i,j,k), sum0t(2,i,j,k), sum0t(3,i,j,k)
!!$                 write(320,'(9es26.16)') x, y, z, sum0tA(1,i,j,k), sum0tA(2,i,j,k), sum0tA(3,i,j,k)
!!$                 write(321,'(9es26.16)') x, y, z, sum0A(1,i,j,k), sum0A(2,i,j,k), sum0A(3,i,j,k)
!!$!!$                 write(33,'(9es26.16)') x, y, z, sum0s(1,i,j,k), sum0s(2,i,j,k), sum0s(3,i,j,k)
!!$!!$                 write(37,'(5es26.16)') x, y, z, sum0c(i,j,k)
           end do
!!$              write(17,*)
!!$              write(31,*)
!!$              write(32,*)
!!$              write(320,*)
!!$              write(321,*)
!!$              write(33,*)
!!$              write(37,*)
        end do
!!$            write(17,*)
!!$            write(31,*)
!!$            write(32,*)
!!$            write(320,*)
!!$            write(321,*)
!!$            write(33,*)
!!$            write(37,*)
     end do
!!$        close(unit=17)
!!$        close(unit=31)
!!$        close(unit=32)
!!$        close(unit=320)
!!$        close(unit=321)
!!$        close(unit=33)
!!$        close(unit=37)
!!$!  stop
      end if !icount==0
!!$  ! end calculate initial properties
  !----------------------------------------------------------------

     if(icount==0) then
            write(filedens,'(a,i10.10,a)') 'denstime.dat'
            write(filetorque,'(a,i10.10,a)') 'torqtime.dat'
            write(filetA,'(a,i10.10,a)') 'tAtime.dat'
            write(fileA,'(a,i10.10,a)') 'Atime.dat'
            write(filespin,'(a,i10.10,a)') 'spintime.dat'
!            write(filespind,'(a,i10.10,a)') 'dspintime.dat'
            write(filezeta,'(a,i10.10,a)') 'zetatime.dat'
            write(filechiral,'(a,i10.10,a)') 'chirtime.dat'
            open(unit=17,file=filedens)
            open(unit=31,file=filezeta)
            open(unit=32,file=filetorque)
            open(unit=320,file=filetA)
            open(unit=321,file=fileA)
            open(unit=33,file=filespin)
!            open(unit=34,file=filespind)
            open(unit=37,file=filechiral)
     end if

     do i=1,meshx
        if(meshx.ne.1) then
                x = xori + dx*(i-1)
        else
                x = x0
        end if

        do j=1,meshy
           if(meshy.ne.1) then
                y = yori + dy*(j-1)
           else
                y = y0
           end if
           
           do k=1,meshz
              if(meshz.ne.1) then
                z = zori + dz*(k-1)
              else
                z = z0
              end if

            sumd = (0.d0,0.d0)
            sumc = (0.d0,0.d0)
            sumt(:) = (0.d0,0.d0)
            sumtA(:) = (0.d0,0.d0)
            sumA(:) = (0.d0,0.d0)
            sums(:) = (0.d0,0.d0)
            sumz(:) = (0.d0,0.d0)
!            dsums(:) = (0.d0,0.d0)

              do nn=1,4*NBS
                 do mm=1,4*NBS
                    sumd = sumd +tmp_d(i,j,k,nn,mm)*calE(nn,mm)
                    sumc = sumc +tmp_c(i,j,k,nn,mm)*calE(nn,mm)
!                    sumd = sumd +density_Qmat(x,y,z,nn,mm)*calE(nn,mm)
!                    sumc = sumc +chiral_Qmat(x,y,z,nn,mm)*calE(nn,mm)
                 end do
              end do

              do l=1,3
                 do nn=1,4*NBS
                    do mm=1,4*NBS
                       sums(l) = sums(l) +tmp_s(l,i,j,k,nn,mm)*calE(nn,mm)
                       sumz(l) = sumz(l) +tmp_z(l,i,j,k,nn,mm)*calE(nn,mm)
!                       sums(l) = sums(l) +s_Qmat(l,x,y,z,nn,mm)*calE(nn,mm)
!                       sumz(l) = sumz(l) +zeta_Qmat(l,x,y,z,nn,mm)*calE(nn,mm)
!                       dsums(l) = dsums(l) +tmp_s(l,i,j,k,nn,mm)*dcalE(nn,mm)/DeltaT
                    end do
                 end do
              end do

              do l=1,3
                tmp_Arad_vec(l) = Arad_vec(time,l,x,y,z)
                do nn=1,4*NBS
                  do mm=1,4*NBS
                    tmp2_j_Qmat(l,nn,mm) = tmp_j(l,i,j,k,nn,mm)
                  end do
                end do
              end do

              call calc_t_Arad_Qmat(tmp2_j_Qmat,tmp_Arad_vec,tmp_t_Arad_Qmat)

              do l=1,3
                 do nn=1,4*NBS
                    do mm=1,4*NBS
                       sumA(l) = sumA(l) +tmp_t_Arad_Qmat(l,nn,mm)*calE(nn,mm)
                       sumt(l) = sumt(l) +tmp_t(l,i,j,k,nn,mm)*calE(nn,mm) !spintorque with Arad effect
                       sumtA(l) = sumtA(l) +(tmp_t(l,i,j,k,nn,mm)+tmp_t_Arad_Qmat(l,nn,mm))*calE(nn,mm) !spintorque with Arad effect
                    end do
                 end do
              end do
    
            write(17,'(10es26.16)') time, x, y, z, sumd, sumd-sum0d(i,j,k), sum0d(i,j,k)
            write(31,'(22es26.16)') time, x, y, z, sumz(1), sumz(1)-sum0z(1,i,j,k), sumz(2), sumz(2)-sum0z(2,i,j,k)&
                                                  &, sumz(3), sumz(3)-sum0z(3,i,j,k)&
                                                  &, sum0z(1,i,j,k), sum0z(2,i,j,k), sum0z(3,i,j,k)
            write(32,'(22es26.16)') time, x, y, z, sumt(1), sumt(1)-sum0t(1,i,j,k), sumt(2), sumt(2)-sum0t(2,i,j,k)&
                                                  &, sumt(3), sumt(3)-sum0t(3,i,j,k)&
                                                  &, sum0t(1,i,j,k), sum0t(2,i,j,k), sum0t(3,i,j,k)
            write(320,'(22es26.16)') time, x, y, z, sumtA(1), sumtA(1)-sum0tA(1,i,j,k), sumtA(2), sumtA(2)-sum0tA(2,i,j,k)&
                                                  &, sumtA(3), sumtA(3)-sum0tA(3,i,j,k)&
                                                  &, sum0tA(1,i,j,k), sum0tA(2,i,j,k), sum0tA(3,i,j,k)
            write(321,'(22es26.16)') time, x, y, z, sumA(1), sumA(1)-sum0A(1,i,j,k), sumA(2), sumA(2)-sum0A(2,i,j,k)&
                                                  &, sumA(3), sumA(3)-sum0A(3,i,j,k)&
                                                  &, sum0A(1,i,j,k), sum0A(2,i,j,k), sum0A(3,i,j,k)
            write(33,'(22es26.16)') time, x, y, z, sums(1), sums(1)-sum0s(1,i,j,k), sums(2), sums(2)-sum0s(2,i,j,k)&
                                                  &, sums(3), sums(3)-sum0s(3,i,j,k)&
                                                  &, sum0s(1,i,j,k), sum0s(2,i,j,k), sum0s(3,i,j,k)
!            write(34,'(10es26.16)') time, x, y, z, dsums(1), dsums(2), dsums(3)
            write(37,'(10es26.16)') time, x, y, z, sumc, sumc-sum0c(i,j,k), sum0c(i,j,k)
    

         end do
       end do
     end do

  !-----------------------end 'for BO approximation'--------------------------------------------


!!$  !----------------------- for only density --------------------------------------------
!!$  !------initial state---------------------------
!!$     write(fileele,'(a,i10.10,a)') 'eledenstime',it,'.dat'
!!$
!!$     open(unit=35,file=fileele)
!!$
!!$     do i=1,meshx
!!$        if(meshx.ne.1) then
!!$                x = xori + dx*(i-1)
!!$        else
!!$                x = 0.d0
!!$        end if
!!$
!!$        do j=1,meshy
!!$           if(meshy.ne.1) then
!!$                y = yori + dy*(j-1)
!!$           else
!!$                y = 0.d0
!!$           end if
!!$           
!!$           do k=1,meshz
!!$              if(meshz.ne.1) then
!!$                z = zori + dz*(k-1)
!!$              else
!!$                z = 0.d0
!!$              end if
!!$
!!$                sum = (0.d0,0.d0)
!!$                do nn=1,4*NBS
!!$                   do mm=1,4*NBS
!!$                      sum = sum +density_Qmat(x,y,z,nn,mm)*calE0(nn,mm)
!!$                   end do
!!$                end do
!!$                pos_e0 = sum
!!$
!!$                sum = (0.d0,0.d0)
!!$                do nn=1,4*NBS
!!$                   do mm=1,4*NBS
!!$                      sum = sum +density_Qmat(x,y,z,nn,mm)*calE(nn,mm)
!!$                   end do
!!$                end do
!!$                pos_e = sum
!!$           
!!$            write(35,'(10es16.6)') time, x, y, z, pos_e, pos_e-pos_e0
!!$
!!$         end do
!!$         write(35,*)
!!$       end do
!!$       write(35,*)
!!$     end do
!!$        close(unit=35)

!!$     end if !end if(mod(icount,nprint)==0)
  end do

        close(unit=17)
        close(unit=31)
        close(unit=32)
        close(unit=320)
        close(unit=321)
        close(unit=33)
!        close(unit=34)
        close(unit=37)

  deallocate(cn,xc,yc,zc)  ! global
  deallocate(aa_L,xx_L,yy_L,zz_L,nx_L,ny_L,nz_L)
  deallocate(aa_S,xx_S,yy_S,zz_S,nx_S,ny_S,nz_S)
  deallocate(c_La,c_Lb,c_Sa,c_Sb)
  deallocate(d_La,d_Lb,d_Sa,d_Sb)

  deallocate(calE,calEtr,calE0,tmpcalE)

  deallocate(sum0d,sum0t,sum0s, sum0z)
  deallocate(sum0c)
  deallocate(tmp_d,tmp_c)
  deallocate(tmp_s,tmp_t,tmp_z)
  deallocate(tmp_j)

end subroutine calc_property_BO_point

!====================================================================================
subroutine calc_property2
!====================================================================================
  Use DiracOutput
  use Constants
  use IntegralStorage

  implicit none

  character(LEN=10) :: str1,str2,str3,str4
  character(LEN=80) :: filedens, filetorque, filespin, filespind, filezeta, filechiral, filetz
  character(LEN=80) :: filej
  character(LEN=80) :: fileele,filenuc
  character(LEN=80) :: filename_E,filename_C
  character(LEN=10) :: moji
  character(LEN=16) :: moji2

  integer :: itmp
  real(kind=8) :: tmp
  real(kind=8) :: time
  integer :: icount

  real(kind=8) :: x,y,z,dx,dy,dz
  integer :: it  ! it-th time step
  complex(kind=8),allocatable :: calE(:,:)  ! calE_NM
  complex(kind=8),allocatable :: calEtr(:,:)  ! transpose(calE)
  complex(kind=8),allocatable :: calE0(:,:)  ! calE_NM
  real(kind=8),allocatable :: tmpcalE(:,:,:)
  integer :: i,j,k,l
  integer :: rr,pp,qq
  integer :: nn,mm
  integer :: NEL

  complex(kind=8) :: density_Qmat, rho_Qmat
  complex(kind=8) :: t_Qmat, t_Arad_Qmat, s_Qmat, zeta_Qmat, chiral_Qmat
  complex(kind=8) :: j_Qmat
  complex(kind=8) :: Arad_vec

  complex(kind=8) :: sum
  complex(kind=8) :: sumd, sum_d, sumc
  complex(kind=8) :: sumt(3), sums(3), sumz(3)
  complex(kind=8),allocatable :: sum0d(:,:,:), sum0c(:,:,:), sum0t(:,:,:,:), sum0s(:,:,:,:), sum0z(:,:,:,:)
  complex(kind=8),allocatable :: sum0j(:,:,:,:)
  complex(kind=8),allocatable :: tmp_d(:,:,:,:,:), tmp_c(:,:,:,:,:), tmp_t(:,:,:,:,:,:), tmp_s(:,:,:,:,:,:), tmp_z(:,:,:,:,:,:)
  complex(kind=8),allocatable :: tmp_j(:,:,:,:,:,:)
  complex(kind=8) :: pos_e,pos_etot
  complex(kind=8) :: pos_e0

  complex(kind=8) :: tmp_Arad_vec(3)
  complex(kind=8),allocatable :: tmp2_j_Qmat(:,:,:)
  complex(kind=8),allocatable :: tmp_t_Arad_Qmat(:,:,:)

  integer :: meshx, meshy, meshz
  real(kind=8) :: xori, xend, yori, yend, zori, zend

  integer :: nprint
  integer :: filenum

!  call getarg(1,str1) !read "nprint"
!  read(str1,'(i10.10)') nprint  ! 

  call read_qedinp ! read qed.inp
  call read_DiracOutput  ! read and set global variables -> defined in sub_readDiracOutput.f90
  call set_GammaMatrix

  NEL = 1 !H, H2+, H4+3

  if(mod(NEL,2).eq.0) then
     NOCC = NEL/2  ! closed shell (KR)
  else
     NOCC = (NEL+1)/2 ! open shell (KR)
     ! The occupation number of the NOCCth orbital is 1.
  end if
  write(*,"(1a10,1i6)") "# NOCC : ",NOCC
  !-end check parameters---------------------------------------------------------

!  allocate(calE(4*NBS,4*NBS))
!  allocate(calEtr(4*NBS,4*NBS))
!  allocate(calE0(4*NBS,4*NBS))
!  allocate(tmpcalE(2,4*NBS,4*NBS))
!  !------------------------------------------------
!  ! Initial condition for calE
!  ! Kronecker delta for electron occupied.
!  calE0(:,:) = (0.d0,0.d0)
!  do i=1,NEL
!     calE0(i,i) = (1.d0,0.d0)
!  end do
!  !------------------------------------------------
  !------------------------------------------------
  ! Set mesh
!  meshx = 51
!  meshy = 51
!  meshz = 51
  meshx = 1
  meshy = 1
  meshz = 101
!  meshx = 9
!  meshy = 9
!  meshz = 13
  xori = -2.5d-2 ; xend = 2.5d-2
  yori = -2.5d-2 ; yend = 2.5d-2
  zori = -2.5d-2 ; zend = 2.5d-2
!  xori = -4.d0 ; xend = 4.d0
!  yori = -4.d0 ; yend = 4.d0
!  zori = -4.d0 ; zend = 4.d0
!  xori = -2.d0 ; xend = 2.d0
!  yori = -2.d0 ; yend = 2.d0
!  zori = -3.d0 ; zend = 3.d0
  if(meshx.ne.1) then
     dx = (xend - xori)/(meshx-1)
  else
     dx = 0.d0
  end if
  if(meshy.ne.1) then
     dy = (yend - yori)/(meshy-1)
  else
     dy = 0.d0
  end if
  if(meshz.ne.1) then
     dz = (zend - zori)/(meshz-1)
  else
     dz = 0.d0
  end if

  allocate(sum0d(meshx,meshy,meshz),sum0t(3,meshx,meshy,meshz), sum0s(3,meshx,meshy,meshz), sum0z(3,meshx,meshy,meshz))
  allocate(sum0c(meshx,meshy,meshz))
!  allocate(sum0j(3,meshx,meshy,meshz))
!  allocate(tmp_d(meshx,meshy,meshz,4*NBS,4*NBS),tmp_c(meshx,meshy,meshz,4*NBS,4*NBS))
!  allocate(tmp_s(3,meshx,meshy,meshz,4*NBS,4*NBS),tmp_t(3,meshx,meshy,meshz,4*NBS,4*NBS),tmp_z(3,meshx,meshy,meshz,4*NBS,4*NBS))
!  allocate(tmp_j(3,meshx,meshy,meshz,4*NBS,4*NBS))
!  allocate(tmp2_j_Qmat(3,4*NBS,4*NBS),tmp_t_Arad_Qmat(3,4*NBS,4*NBS))

  write(*,*)'meshx =',meshx
  write(*,*)'meshy =',meshy
  write(*,*)'meshz =',meshz
  write(*,*)'dx dy dz =',dx,dy,dz
  !------------------------------------------------

  !------------------------------------------------

  !-----------------------for BO approximation--------------------------------------------
  ! calculate initial properties
            write(filedens,'(a)') 'dens.dat'
            write(filetorque,'(a)') 'torq.dat'
            write(filespin,'(a)') 'spin.dat'
            write(filezeta,'(a)') 'zeta.dat'
            write(filechiral,'(a)') 'chir.dat'
            write(filetz,'(a)') 'tz.dat'
            write(filej,'(a)') 'j.dat'
            open(unit=17,file=filedens)
            open(unit=31,file=filezeta)
            open(unit=32,file=filetorque)
            open(unit=33,file=filespin)
            open(unit=35,file=filetz)
            open(unit=37,file=filechiral)
            open(unit=38,file=filej)

            sum0d(:,:,:) = (0.d0,0.d0)
            sum0c(:,:,:) = (0.d0,0.d0)
            sum0t(:,:,:,:) = (0.d0,0.d0)
            sum0s(:,:,:,:) = (0.d0,0.d0)
            sum0z(:,:,:,:) = (0.d0,0.d0)
!            sum0j(:,:,:,:) = (0.d0,0.d0)
     do i=1,meshx
        if(meshx.ne.1) then
                x = xori + dx*(i-1)
        else
                x = 0.d0
        end if

        do j=1,meshy
           if(meshy.ne.1) then
                y = yori + dy*(j-1)
           else
                y = 0.d0
           end if
           
           do k=1,meshz
              if(meshz.ne.1) then
                z = zori + dz*(k-1)
              else
                z = 0.d0
              end if

              do nn=1,NEL
                    sum0d(i,j,k) = sum0d(i,j,k) +density_Qmat(x,y,z,nn,nn)
                    sum0c(i,j,k) = sum0c(i,j,k) +chiral_Qmat(x,y,z,nn,nn)
              end do

              do l=1,3
                 do nn=1,NEL
                       sum0s(l,i,j,k) = sum0s(l,i,j,k) +s_Qmat(l,x,y,z,nn,nn)
                       sum0z(l,i,j,k) = sum0z(l,i,j,k) +zeta_Qmat(l,x,y,z,nn,nn)
                       sum0t(l,i,j,k) = sum0t(l,i,j,k) +t_Qmat(l,x,y,z,nn,nn)
!                       sum0j(l,i,j,k) = sum0j(l,i,j,k) +j_Qmat(l,x,y,z,nn,nn)
                 end do
              end do

!!$!-calc KRM---------------------------------------------------------------------------
!!$                 sum0d(i,j,k) = sum0d(i,j,k) +density_Qmat(x,y,z,2,2)*calE(2,2)
!!$
!!$              do l=1,3
!!$                 sum0t(l,i,j,k) = sum0t(l,i,j,k) +t_Arad_Qmat(0.d0,l,x,y,z,2,2)*calE(2,2) !spintorque with Arad effect
!!$                 sum0s(l,i,j,k) = sum0s(l,i,j,k) +s_Qmat(l,x,y,z,2,2)*calE(2,2)
!!$                 sum0z(l,i,j,k) = sum0z(l,i,j,k) +zeta_Qmat(l,x,y,z,2,2)*calE(2,2)
!!$              end do
!!$!----------------------------------------------------------------------------

                 write(17,'(5es24.14)') x, y, z, sum0d(i,j,k)
                 write(31,'(9es24.14)') x, y, z, sum0z(1,i,j,k), sum0z(2,i,j,k), sum0z(3,i,j,k)
                 write(32,'(9es24.14)') x, y, z, sum0t(1,i,j,k), sum0t(2,i,j,k), sum0t(3,i,j,k)
                 write(33,'(9es24.14)') x, y, z, sum0s(1,i,j,k), sum0s(2,i,j,k), sum0s(3,i,j,k)
                 write(35,'(9es24.14)') x, y, z, sum0z(1,i,j,k)+sum0t(1,i,j,k), sum0z(2,i,j,k)+sum0t(2,i,j,k), sum0z(3,i,j,k)+sum0t(3,i,j,k)
                 write(37,'(5es24.14)') x, y, z, sum0c(i,j,k)
!                 write(38,'(9es24.14)') x, y, z, sum0j(1,i,j,k), sum0j(2,i,j,k), sum0j(3,i,j,k)
           end do
              write(17,*)
              write(31,*)
              write(32,*)
              write(33,*)
              write(35,*)
              write(37,*)
              write(38,*)
        end do
            write(17,*)
            write(31,*)
            write(32,*)
            write(33,*)
            write(35,*)
            write(37,*)
            write(38,*)
        write(*,*) i,'/',meshx
     end do
        close(unit=17)
        close(unit=31)
        close(unit=32)
        close(unit=33)
        close(unit=35)
        close(unit=37)
        close(unit=38)
!  stop
  ! end calculate initial properties
  !----------------------------------------------------------------

end subroutine calc_property2
