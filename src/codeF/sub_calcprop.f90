! Last Change: 28-Aug-2014.
!====================================================================================
! 2012.10.12 
!  subroutine prop_calc_menu     ! set mesh etc
!  subroutine calc_prop3d        ! this is main program to calc prop3d
!  subroutine alloc              ! allocate some valuables
!  subroutine set_filename       !
!  subroutine prop_calc_matrices !calc tmp val at each mesh
!  subroutine prop_calc_ini      ! calc initial property at each mesh
!  subroutine prop_calc          ! calc properties at each step
!  subroutine set_paramini       ! read param.ini, qed.inp
!====================================================================================
!====================================================================================
subroutine prop_calc_menu(filename_E, nprint, filenum)
!====================================================================================
  
  use prop_param
  implicit none

  integer i, itmp
  integer,parameter :: N=20
  integer tmp(N)
  character(LEN=80) tmpf(N)

  integer,intent(out) :: nprint ! interval of print files
  integer,intent(out) :: filenum ! number of print file
  character(LEN=80),intent(out) :: filename_E  ! name of calE.dat

  !----------------------------------------------------------------------
!  nprint = 120 ! interval of print files
!  filenum = 40  ! number of print file
!  nprint = 600 ! interval of print files
!  filenum = 8  ! number of print file
  nprint = 50 ! interval of print files
  filenum = 8  ! number of print file
!  nprint = 1 ! interval of print files
!  filenum = 5000  ! number of print file
!  nprint = 1 ! interval of print files
!  filenum = 10000  ! number of print file
  filename_E = '../calE.dat' ! name of calE.dat
  write(*,'(a)') filename_E
  !----------------------------------------------------------------------

  !---calculation range and mesh ----------------------------------------
    ! If you want to calc point, specify mesh=1.
    meshx = 21!21
    meshy = 21!101
    meshz = 21
!    xori = -4.d0 ; xend = 6.d0
!    yori = -4.d0 ; yend = 6.d0
!    zori = -4.d0 ; zend = 6.d0
!    xori = -3.d0 ; xend = 3.d0
!    yori = -3.d0 ; yend = 3.d0
!    zori = -3.d0 ; zend = 3.d0
    xori = -1.d0 ; xend = 1.d0
    yori = -1.d0 ; yend = 1.d0
    zori = -1.d0 ; zend = 1.d0
!    xori = -0.5d0 ; xend = 0.5d0
!    yori = -0.5d0 ; yend = 0.5d0
!    zori = -0.5d0 ; zend = 0.5d0
    write(*,*)'meshx =',meshx
    write(*,*)'meshy =',meshy
    write(*,*)'meshz =',meshz
!    x0 = 1.d0 ; y0 = 1.d0 ; z0 = 1.d0 ! for calc point 
!    x0 = -1.d0 ; y0 = 0.d0 ; z0 = 0.d0 ! for calc point 
!    x0 = 0.d0 ; y0 = 0.d0 ; z0 = 0.1d0 ! for calc point 
    open(unit=301,file="coord.ini",status='old')
      read(301,*) x0
      read(301,*) y0
      read(301,*) z0
    close(unit=301)
    write(*,*)'(x0,y0,z0)=', x0, y0, z0
  !----------------------------------------------------------------------

  !---select property to calculate---------------------------------------
    tmp(:) = 0
    tmp(1)  = 61;  tmpf(1)  = 'dens'    ! electron density
    tmp(2)  = 62;  tmpf(2)  = 'chir'    ! chiral density
    tmp(3)  = 63;  tmpf(3)  = 'spin'    ! spin angular momentum density
    tmp(4)  = 64;  tmpf(4)  = 'torq'    ! spin torque density without Arad
    tmp(5)  = 65;  tmpf(5)  = 'zeta'    ! zeta force density
!!$    tmp(6)  = 66;  tmpf(6)  = 'Atorq'   !(68 is needed) spin torque Arad component
!!$    tmp(7)  = 67;  tmpf(7)  = 'tA'      !(64, 66 and 68 are needed) spin torque with Arad
    tmp(8)  = 68;  tmpf(8)  = 'j'       ! electric current density
    tmp(9)  = 69;  tmpf(9)  = 'divj'    ! divergence j
    tmp(10) = 70;  tmpf(10) = 'stress'  ! stress tensor
!!$    tmp(11) = 71;  tmpf(11) = 'tauA'    ! A component of stress tensor
    tmp(12) = 72;  tmpf(12) = 'pol'     ! polarization vector
    tmp(13) = 73;  tmpf(13) = 'ef'      ! electric field
!!$    tmp(14) = 74;  tmpf(14) = 'tauAM'   ! AM component of stress tensor
!!$    tmp(15) = 75;  tmpf(15) = 'tAM'     ! torque AM component
!!$    tmp(16) = 76;  tmpf(16) = 'tAMt'    ! spin torque with AM
    tmp(17) = 77;  tmpf(17) = 'pi'         ! kinetic momentum without AM
    tmp(18) = 78;  tmpf(18) = 'orbang'     ! orbital angular momentum
    tmp(19) = 79;  tmpf(19) = 'rots'       ! spin vorticity
  !----------------------------------------------------------------------

  itmp = 0
  do i=1,N
    if(tmp(i).ne.0) then
      itmp = itmp + 1
    end if
  end do

  allocate(unitn(itmp),filename(itmp))

  itmp = 0
  do i=1,N
    if(tmp(i).ne.0) then
      itmp = itmp + 1
      unitn(itmp) = tmp(i)
      filename(itmp) = tmpf(i)
    end if
  end do

  itmp=0
  do i=1,size(unitn)
    if(unitn(i).eq.68) itmp=1
  end do
  do i=1,size(unitn)
    if(((itmp==0).and.(unitn(i)==66)) .or.&
      &((itmp==0).and.(unitn(i)==67))) then
        write(*,*)'You have to specify 68 in calcmenu.'
        stop
    end if
  end do

end subroutine prop_calc_menu
!====================================================================================

!====================================================================================
subroutine calc_prop3d
!====================================================================================
  use Precision
  Use DiracOutput
  use Constants
  use IntegralStorage
  use NucBasis
  use prop_param
  use prop_tmp

  implicit none

  integer :: i,j,k,ii,icount
  integer :: it
  real(kind=dp) :: time
  real(kind=dp) :: x,y,z!,dx,dy,dz
  integer :: nprint ! interval of print files
  integer :: filenum ! number of print file
  character(LEN=80) :: filename_E ! name of calE.dat
  character(LEN=80) :: FMT
  complex(kind=dp),allocatable :: calE(:,:)  ! calE_NM
  complex(kind=dp),allocatable :: calE0(:,:)  ! calE0_NM
  character(LEN=80),allocatable :: filedat(:)

  !set param_AM
  call read_param_AM

  call set_paramini ! read param.ini, qed.inp
  call prop_calc_menu(filename_E, nprint, filenum) ! set mesh etc
  call alloc(NBS) ! allocate some valuables
  allocate(calE(4*NBS,4*NBS))
  allocate(calE0(4*NBS,4*NBS))

  dx = 0.d0
  if(meshx.ne.1) dx = (xend - xori)/(meshx-1)
  dy = 0.d0
  if(meshy.ne.1) dy = (yend - yori)/(meshy-1)
  dz = 0.d0
  if(meshz.ne.1) dz = (zend - zori)/(meshz-1)

  write(*,*)'dx dy dz =',dx,dy,dz

  !------set property matrices---------------------
  write(*,*)'set property matrices'
         do i=1,meshx
            x = xori + dx*(i-1)
!            if(meshx.eq.1) x = 2.8345892009d0
            if(meshx.eq.1) x = 0.d0
    
            do j=1,meshy
               y = yori + dy*(j-1)
               if(meshy.eq.1) y = 0.d0
               
               do k=1,meshz
                  z = zori + dz*(k-1)
                  if(meshz.eq.1) z = 0.d0
    
              !calc tmp val at each mesh
              call prop_calc_matrices(i,j,k,x,y,z)

           end do !k=1,meshz
        end do !j=1,meshy
        write(*,*) i,'/',meshx
     end do !i=1,meshx
   write(*,*)'finish set property matrices'
   !----end set property matrices---------------------
  !----------------------------------------------------------
  ! set calE0 = <0|calE(t=0)|0>
  !----------------------------------------------------------
  calE0(:,:) = (0._dp,0._dp)
  do i=2*NBS+1,4*NBS
     calE0(i,i) = (1._dp,0._dp)
  end do

  !read filename_E(calE.dat)
  open(unit=1200,file=filename_E,status='old')
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
  !----------------------------------------------------------
  ! set <phi|:calE:|phi> = <phi|calE|phi> - <0|calE(t=0)|0>
  !----------------------------------------------------------
  calE(:,:) = calE(:,:) - calE0(:,:)

     if((icount==0).or.(mod(icount,nprint)==0)) then
        it = nint(time/DeltaT)

  !----------------------------------------------------------------
     ! calc initial property
     if(icount==0) then
!            write(filedat,'(a)') 'dens.dat'
!            open(unit=unitn,file=filedat)
         write(*,*)'start calc initial properties'

         do i=1,meshx
            x = xori + dx*(i-1)
            if(meshx.eq.1) x = 0.d0
    
            do j=1,meshy
               y = yori + dy*(j-1)
               if(meshy.eq.1) y = 0.d0
               
               do k=1,meshz
                  z = zori + dz*(k-1)
                  if(meshz.eq.1) z = 0.d0

                  ! calc initial property at each mesh
                  call prop_calc_ini(i,j,k,x,y,z,calE)
    
               end do !k=1,meshz
!                  write(unitn,*)
            end do !j=1,meshy
!                write(unitn,*)
         end do !i=1,meshx
!            close(unit=unitn)
         write(*,*)'finish calc initial properties'
       end if !icount==0
     ! end calculate initial properties
  !----------------------------------------------------------------
     !start calculate properties at each step

      call system('mkdir -p data')

      allocate(filedat(size(filename)))
     ! set filename *.dat
      call set_filename(filedat,it)

      do ii=1,size(unitn)
        open(unit=unitn(ii),file=filedat(ii))
      end do

      do i=1,meshx
         x = xori + dx*(i-1)
         if(meshx.eq.1) x = 0.d0
 
         do j=1,meshy
            y = yori + dy*(j-1)
            if(meshy.eq.1) y = 0.d0
            
            do k=1,meshz
               z = zori + dz*(k-1)
               if(meshz.eq.1) z = 0.d0

               call prop_calc(i,j,k,x,y,z,calE,time)

            end do!k=1,meshz
              do ii=1,size(unitn)
                write(unitn(ii),*)
              end do
         end do!j=1,meshy
             do ii=1,size(unitn)
               write(unitn(ii),*)
             end do
      end do!i=1,meshx

      do ii=1,size(unitn)
        close(unit=unitn(ii))
      end do

      deallocate(filedat)

     end if !end if(mod(icount,nprint)==0)
  end do !icount

end subroutine calc_prop3d
!====================================================================================

!====================================================================================
subroutine calc_prop_point
!====================================================================================
  use Precision
  Use DiracOutput
  use Constants
  use IntegralStorage
  use NucBasis
  use prop_param
  use prop_tmp

  implicit none

  integer :: i,j,k,ii,icount
  integer :: it
  real(kind=dp) :: time
  real(kind=dp) :: x,y,z!,dx,dy,dz
  integer :: nprint ! interval of print files
  integer :: filenum ! number of print file
  character(LEN=80) :: filename_E ! name of calE.dat
  character(LEN=80) :: FMT
  complex(kind=dp),allocatable :: calE(:,:)  ! calE_NM
  complex(kind=dp),allocatable :: calE0(:,:)  ! calE0_NM
  character(LEN=80),allocatable :: filedat(:)


  !set param_AM
  call read_param_AM

  call set_paramini ! read param.ini, qed.inp
  call prop_calc_menu(filename_E, nprint, filenum) ! set mesh etc
  call alloc(NBS) ! allocate some valuables
  allocate(calE(4*NBS,4*NBS))
  allocate(calE0(4*NBS,4*NBS))

  dx = 0.d0
  if(meshx.ne.1) dx = (xend - xori)/(meshx-1)
  dy = 0.d0
  if(meshy.ne.1) dy = (yend - yori)/(meshy-1)
  dz = 0.d0
  if(meshz.ne.1) dz = (zend - zori)/(meshz-1)

  write(*,*)'dx dy dz =',dx,dy,dz

  !------set property matrices---------------------
  write(*,*)'set property matrices'
         do i=1,meshx
            x = xori + dx*(i-1)
            if(meshx.eq.1) x = x0
    
            do j=1,meshy
               y = yori + dy*(j-1)
               if(meshy.eq.1) y = y0
               
               do k=1,meshz
                  z = zori + dz*(k-1)
                  if(meshz.eq.1) z = z0
    
              !calc tmp val at each mesh
              call prop_calc_matrices(i,j,k,x,y,z)

           end do !k=1,meshz
        end do !j=1,meshy
        write(*,*) i,'/',meshx
     end do !i=1,meshx
   write(*,*)'finish set property matrices'
   !----end set property matrices---------------------

      call system('mkdir data')

  allocate(filedat(size(filename)))
 ! set filename *.dat
  call set_filename2(filedat,it)
  do ii=1,size(unitn)
    open(unit=unitn(ii),file=filedat(ii))
  end do

  !----------------------------------------------------------
  ! set calE0 = <0|calE(t=0)|0>
  !----------------------------------------------------------
  calE0(:,:) = (0._dp,0._dp)
  do i=2*NBS+1,4*NBS
     calE0(i,i) = (1._dp,0._dp)
  end do

  !read filename_E(calE.dat)
  open(unit=1200,file=filename_E,status='old')
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
  !----------------------------------------------------------
  ! set <phi|:calE:|phi> = <phi|calE|phi> + <0|calE(t=0)|0>
  !----------------------------------------------------------
  calE(:,:) = calE(:,:) - calE0(:,:)


     if((icount==0).or.(mod(icount,nprint)==0)) then
        it = nint(time/DeltaT)

  !----------------------------------------------------------------
     ! calc initial property
     if(icount==0) then
!            write(filedat,'(a)') 'dens.dat'
!            open(unit=unitn,file=filedat)
         write(*,*)'start calc initial properties'

         do i=1,meshx
            x = xori + dx*(i-1)
            if(meshx.eq.1) x = x0
    
            do j=1,meshy
               y = yori + dy*(j-1)
               if(meshy.eq.1) y = y0
               
               do k=1,meshz
                  z = zori + dz*(k-1)
                  if(meshz.eq.1) z = z0

                  ! calc initial property at each mesh
                  call prop_calc_ini(i,j,k,x,y,z,calE)
    
               end do !k=1,meshz
!                  write(unitn,*)
            end do !j=1,meshy
!                write(unitn,*)
         end do !i=1,meshx
!            close(unit=unitn)
         write(*,*)'finish calc initial properties'
       end if !icount==0
     ! end calculate initial properties
  !----------------------------------------------------------------
     !start calculate properties at each step

      do i=1,meshx
         x = xori + dx*(i-1)
         if(meshx.eq.1) x = x0
 
         do j=1,meshy
            y = yori + dy*(j-1)
            if(meshy.eq.1) y = y0
            
            do k=1,meshz
               z = zori + dz*(k-1)
               if(meshz.eq.1) z = z0

               call prop_calc(i,j,k,x,y,z,calE,time)

            end do!k=1,meshz
!              do ii=1,size(unitn)
!                write(unitn(ii),*)
!              end do
         end do!j=1,meshy
!             do ii=1,size(unitn)
!               write(unitn(ii),*)
!             end do
      end do!i=1,meshx

     end if !end if(mod(icount,nprint)==0)
  end do !icount

      do ii=1,size(unitn)
        close(unit=unitn(ii))
      end do

      deallocate(filedat)

end subroutine calc_prop_point
!====================================================================================

!====================================================================================
subroutine alloc(NBS)
!====================================================================================
  
  use prop_tmp
  use prop_param
  implicit none

  integer,intent(in) :: NBS
  integer :: ii

  do ii=1,size(unitn)

    if(unitn(ii).eq.61) then
      allocate(sum0d(meshx,meshy,meshz))
      allocate(tmp_d(meshx,meshy,meshz,4*NBS,4*NBS))

    else if(unitn(ii).eq.62) then
      allocate(sum0c(meshx,meshy,meshz))
      allocate(tmp_c(meshx,meshy,meshz,4*NBS,4*NBS))

    else if(unitn(ii).eq.63) then
      allocate(sum0s(3,meshx,meshy,meshz))
      allocate(tmp_s(3,meshx,meshy,meshz,4*NBS,4*NBS))

    else if(unitn(ii).eq.64) then
      allocate(sum0t(3,meshx,meshy,meshz))
      allocate(tmp_t(3,meshx,meshy,meshz,4*NBS,4*NBS))

    else if(unitn(ii).eq.65) then
      allocate(sum0z(3,meshx,meshy,meshz))
      allocate(tmp_z(3,meshx,meshy,meshz,4*NBS,4*NBS))

    else if(unitn(ii).eq.66) then
      allocate(sum0A(3,meshx,meshy,meshz))

    else if(unitn(ii).eq.67) then
      allocate(sum0tA(3,meshx,meshy,meshz))

    else if(unitn(ii).eq.68) then
      allocate(sum0j(3,meshx,meshy,meshz))
      allocate(tmp_j(3,meshx,meshy,meshz,4*NBS,4*NBS))
      allocate(tmp2_j_Qmat(3,4*NBS,4*NBS))
      allocate(tmp_t_Arad_Qmat(3,4*NBS,4*NBS))

    else if(unitn(ii).eq.69) then
      allocate(sum0divj(meshx,meshy,meshz))
      allocate(tmp_divj(meshx,meshy,meshz,4*NBS,4*NBS))

    else if(unitn(ii).eq.70) then
      allocate(sum0tau(3,3,meshx,meshy,meshz))
      allocate(tmp_tau(3,3,meshx,meshy,meshz,4*NBS,4*NBS))

    else if(unitn(ii).eq.71) then
      allocate(sum0tauA(3,3,meshx,meshy,meshz))

    else if(unitn(ii).eq.72) then
      allocate(sum0pol(3,meshx,meshy,meshz))
      allocate(tmp_E_Qmat(3,meshx,meshy,meshz,4*NBS,4*NBS)) 

    else if(unitn(ii).eq.73) then
      allocate(sum0ef(3,meshx,meshy,meshz))

    else if(unitn(ii).eq.74) then
      allocate(sum0tauAM(3,3,meshx,meshy,meshz))
      allocate(tmp_tauAM(3,3,meshx,meshy,meshz,4*NBS,4*NBS))

    else if(unitn(ii).eq.75) then
      allocate(sum0t_AM(3,meshx,meshy,meshz))
      allocate(tmp_t_AM(3,meshx,meshy,meshz,4*NBS,4*NBS))

    else if(unitn(ii).eq.77) then
      allocate(sum0pi(3,meshx,meshy,meshz))
      allocate(tmp_pi(3,meshx,meshy,meshz,4*NBS,4*NBS))

    else if(unitn(ii).eq.78) then
      allocate(sum0orbang(3,meshx,meshy,meshz))
      allocate(tmp_orbang(3,meshx,meshy,meshz,4*NBS,4*NBS))

    else if(unitn(ii).eq.79) then
      allocate(sum0rots(3,meshx,meshy,meshz))
      allocate(tmp_rots(3,meshx,meshy,meshz,4*NBS,4*NBS))

    end if

  end do !ii

end subroutine alloc
!====================================================================================

!====================================================================================
subroutine set_filename(filedat,it)
!====================================================================================
  use prop_param
  use prop_tmp
  implicit none

  integer, intent(in) :: it
  character(LEN=80) :: filedat(size(filename))
  integer :: i

  do i=1,size(filename)
    write(filedat(i),'(a,a,a4,i10.10,a)') 'data/',trim(filename(i)),'time',it,'.dat'
  end do

end subroutine set_filename
!====================================================================================

!====================================================================================
subroutine set_filename2(filedat)
!====================================================================================
  use prop_param
  use prop_tmp
  implicit none

  character(LEN=80) :: filedat(size(filename))
  integer :: i

  do i=1,size(filename)
    write(filedat(i),'(a,a,a)') 'data/',trim(filename(i)),'point.dat'
  end do

end subroutine set_filename2
!====================================================================================

!====================================================================================
subroutine prop_calc_matrices(i,j,k,x,y,z)
!====================================================================================
  
  use Precision
  Use DiracOutput
  use prop_param
  use prop_tmp
  implicit none

  integer, intent(in) :: i,j,k
  real(kind=dp),intent(in) :: x,y,z
  integer l,m,ii
  integer :: nn,mm
  real(kind=dp) :: vecR(3)

  complex(kind=dp) :: density_Qmat, rho_Qmat
  complex(kind=dp) :: t_Qmat, t_Arad_Qmat, s_Qmat, zeta_Qmat, chiral_Qmat
  complex(kind=dp) :: j_Qmat
  complex(kind=dp) :: Arad_vec
  complex(kind=dp) :: tau_Qmat
  complex(kind=dp) :: divj_Qmat
  complex(kind=dp) :: E_Qmat(3,4*NBS,4*NBS)  ! E^k_NM
  complex(kind=dp) :: tau_AM_Qmat, t_AM_Qmat
  complex(kind=dp) :: pi_Qmat, orbang_Qmat, rots_Qmat

  do ii=1,size(unitn)

    if(unitn(ii).eq.61) then
      !$omp parallel do
      do mm=1,4*NBS
         do nn=1,4*NBS
            tmp_d(i,j,k,nn,mm) = density_Qmat(x,y,z,nn,mm)
         end do
      end do
      !$omp end parallel do
    
    else if(unitn(ii).eq.62) then
      !$omp parallel do
      do mm=1,4*NBS
         do nn=1,4*NBS
            tmp_c(i,j,k,nn,mm) = chiral_Qmat(x,y,z,nn,mm)
         end do
      end do
      !$omp end parallel do
    
    else if(unitn(ii).eq.63) then
      do l=1,3
         !$omp parallel do
         do mm=1,4*NBS
            do nn=1,4*NBS
               tmp_s(l,i,j,k,nn,mm) = s_Qmat(l,x,y,z,nn,mm)
            end do
         end do
         !$omp end parallel do
      end do
    
    else if(unitn(ii).eq.64) then
      do l=1,3
         !$omp parallel do
         do mm=1,4*NBS
            do nn=1,4*NBS
               tmp_t(l,i,j,k,nn,mm) = t_Qmat(l,x,y,z,nn,mm)
            end do
         end do
         !$omp end parallel do
      end do
    
    else if(unitn(ii).eq.65) then
      do l=1,3
         !$omp parallel do
         do mm=1,4*NBS
            do nn=1,4*NBS
               tmp_z(l,i,j,k,nn,mm) = zeta_Qmat(l,x,y,z,nn,mm)
            end do
         end do
         !$omp end parallel do
      end do

    else if((unitn(ii).eq.68)) then
      do l=1,3
         !$omp parallel do
         do mm=1,4*NBS
            do nn=1,4*NBS
               tmp_j(l,i,j,k,nn,mm) = j_Qmat(l,x,y,z,nn,mm)
            end do
         end do
         !$omp end parallel do
      end do

    else if((unitn(ii).eq.69)) then
      !$omp parallel do
      do mm=1,4*NBS
         do nn=1,4*NBS
            tmp_divj(i,j,k,nn,mm) = divj_Qmat(x,y,z,nn,mm)
         end do
      end do
      !$omp end parallel do

    else if((unitn(ii).eq.70)) then
      do l=1,3
      do m=1,3
         !$omp parallel do
         do mm=1,4*NBS
            do nn=1,4*NBS
               tmp_tau(m,l,i,j,k,nn,mm) = tau_Qmat(m,l,x,y,z,nn,mm)
            end do
         end do
         !$omp end parallel do
      end do
      end do

    else if((unitn(ii).eq.72)) then
      vecR(1)=x; vecR(2)=y; vecR(3)=z
      call setQmat_E(vecR,E_Qmat)
      do l=1,3
         do mm=1,4*NBS
            do nn=1,4*NBS
               tmp_E_Qmat(l,i,j,k,nn,mm) = E_Qmat(l,nn,mm)
            end do
         end do
      end do

    else if((unitn(ii).eq.74)) then
      do l=1,3
      do m=1,3
         !$omp parallel do
         do mm=1,4*NBS
            do nn=1,4*NBS
               ! CAUTION!!! this array should be reduced to 7.
               tmp_tauAM(m,l,i,j,k,nn,mm) = tau_AM_Qmat(0.d0,m,l,x,y,z,nn,mm)
            end do
         end do
         !$omp end parallel do
      end do
      end do

    else if((unitn(ii).eq.75)) then
      do l=1,3
         !$omp parallel do
         do mm=1,4*NBS
            do nn=1,4*NBS
               tmp_t_AM(l,i,j,k,nn,mm) = t_AM_Qmat(0.d0,l,x,y,z,nn,mm)
            end do
         end do
         !$omp end parallel do
      end do

    else if((unitn(ii).eq.77)) then
      do l=1,3
         !$omp parallel do
         do mm=1,4*NBS
            do nn=1,4*NBS
               tmp_pi(l,i,j,k,nn,mm) = pi_Qmat(l,x,y,z,nn,mm)
            end do
         end do
         !$omp end parallel do
      end do

    else if((unitn(ii).eq.78)) then
      do l=1,3
         !$omp parallel do
         do mm=1,4*NBS
            do nn=1,4*NBS
               tmp_orbang(l,i,j,k,nn,mm) = orbang_Qmat(l,x,y,z,nn,mm)
            end do
         end do
         !$omp end parallel do
      end do

    else if((unitn(ii).eq.79)) then
      do l=1,3
         !$omp parallel do
         do mm=1,4*NBS
            do nn=1,4*NBS
               tmp_rots(l,i,j,k,nn,mm) = rots_Qmat(l,x,y,z,nn,mm)
            end do
         end do
         !$omp end parallel do
      end do

!    else if(unitn(ii).ne.0) then
!      write(*,*) 'error unitn(',ii,')=',unitn(ii)
    end if !unitn(ii)

  end do
    
end subroutine prop_calc_matrices
!====================================================================================

!====================================================================================
subroutine prop_calc_ini(i,j,k,x,y,z,calE)
!====================================================================================
  
  use Precision
  Use DiracOutput
  use Constants
  use prop_param
  use prop_tmp
  implicit none

  integer, intent(in) :: i,j,k
  real(kind=dp),intent(in) :: x,y,z
  real(kind=dp) :: vecR(3)
  complex(kind=dp), intent(in) :: calE(4*NBS,4*NBS)
  integer :: l,m, ii
  integer :: nn,mm
  complex(kind=dp) :: Arad_vec
  complex(kind=dp) :: tmp_Arad_vec(3)
  complex(kind=dp) :: tmp_tau_Arad_Qmat(3,3,4*NBS,4*NBS)
  complex(kind=dp) :: tmp
  complex(kind=dp) :: E_Qmat(3,4*NBS,4*NBS)  ! E^k_NM
  real(kind=dp) :: tmp_pol(3)
  real(kind=dp) :: dAraddt_4E

  do ii=1,size(unitn)

    if(unitn(ii).eq.61) then
      sum0d(i,j,k) = (0.d0,0.d0)
!      tmp = (0.d0,0.d0)
      do mm=1,4*NBS
         do nn=1,4*NBS
            sum0d(i,j,k) = sum0d(i,j,k) +tmp_d(i,j,k,nn,mm)*calE(nn,mm)
!            tmp = tmp +tmp_d(i,j,k,nn,mm)*calE(nn,mm)
         end do
      end do
!      sum0d(i,j,k) = tmp
    
    else if(unitn(ii).eq.62) then
      sum0c(i,j,k) = (0.d0,0.d0)
      do mm=1,4*NBS
         do nn=1,4*NBS
            sum0c(i,j,k) = sum0c(i,j,k) +tmp_c(i,j,k,nn,mm)*calE(nn,mm)
         end do
      end do
    
    else if(unitn(ii).eq.63) then
      do l=1,3
         sum0s(l,i,j,k) = (0.d0,0.d0)
         do mm=1,4*NBS
            do nn=1,4*NBS
               sum0s(l,i,j,k) = sum0s(l,i,j,k) +tmp_s(l,i,j,k,nn,mm)*calE(nn,mm)
            end do
         end do
      end do
    
    else if(unitn(ii).eq.64) then
      do l=1,3
         sum0t(l,i,j,k) = (0.d0,0.d0)
         do mm=1,4*NBS
            do nn=1,4*NBS
               sum0t(l,i,j,k) = sum0t(l,i,j,k) +tmp_t(l,i,j,k,nn,mm)*calE(nn,mm) !spintorque
            end do
         end do
      end do
    
    else if(unitn(ii).eq.65) then
      do l=1,3
         sum0z(l,i,j,k) = (0.d0,0.d0)
         do mm=1,4*NBS
            do nn=1,4*NBS
               sum0z(l,i,j,k) = sum0z(l,i,j,k) +tmp_z(l,i,j,k,nn,mm)*calE(nn,mm)
            end do
         end do
      end do
    
    else if(unitn(ii).eq.66) then
      !$omp parallel do
      do l=1,3
        tmp_Arad_vec(l) = Arad_vec(0.d0,l,x,y,z)
        do mm=1,4*NBS
          do nn=1,4*NBS
            tmp2_j_Qmat(l,nn,mm) = tmp_j(l,i,j,k,nn,mm)
          end do
        end do
      end do
      !$omp end parallel do
    
      call calc_t_Arad_Qmat(tmp2_j_Qmat,tmp_Arad_vec,tmp_t_Arad_Qmat)
    
      do l=1,3
         sum0A(l,i,j,k) = (0.d0,0.d0)
         do mm=1,4*NBS
            do nn=1,4*NBS
               sum0A(l,i,j,k) = sum0A(l,i,j,k) +tmp_t_Arad_Qmat(l,nn,mm)*calE(nn,mm) !spintorque with Arad effect
            end do
         end do
      end do
    
    else if(unitn(ii).eq.67) then
      do l=1,3
         sum0tA(l,i,j,k) = (0.d0,0.d0)
         do mm=1,4*NBS
            do nn=1,4*NBS
               sum0tA(l,i,j,k) = sum0tA(l,i,j,k) +(tmp_t(l,i,j,k,nn,mm)+tmp_t_Arad_Qmat(l,nn,mm))*calE(nn,mm) !spintorque with Arad effect
            end do
         end do
      end do

    else if(unitn(ii).eq.68) then
      do l=1,3
         sum0j(l,i,j,k) = (0.d0,0.d0)
         do mm=1,4*NBS
            do nn=1,4*NBS
               sum0j(l,i,j,k) = sum0j(l,i,j,k) +tmp_j(l,i,j,k,nn,mm)*calE(nn,mm)
            end do
         end do
      end do

    else if(unitn(ii).eq.69) then
      sum0divj(i,j,k) = (0.d0,0.d0)
      do mm=1,4*NBS
         do nn=1,4*NBS
            sum0divj(i,j,k) = sum0divj(i,j,k) +tmp_divj(i,j,k,nn,mm)*calE(nn,mm)
         end do
      end do

    else if(unitn(ii).eq.70) then
      do l=1,3
      do m=1,3
         sum0tau(m,l,i,j,k) = (0.d0,0.d0)
         do mm=1,4*NBS
            do nn=1,4*NBS
               sum0tau(m,l,i,j,k) = sum0tau(m,l,i,j,k) +tmp_tau(m,l,i,j,k,nn,mm)*calE(nn,mm)
            end do
         end do
      end do
      end do

    else if(unitn(ii).eq.71) then
      !$omp parallel do
      do l=1,3
        tmp_Arad_vec(l) = Arad_vec(0.d0,l,x,y,z)
        do mm=1,4*NBS
           do nn=1,4*NBS
            tmp2_j_Qmat(l,nn,mm) = tmp_j(l,i,j,k,nn,mm)
          end do
        end do
      end do
      !$omp end parallel do
    
      call calc_tau_Arad_Qmat(tmp2_j_Qmat,tmp_Arad_vec,tmp_tau_Arad_Qmat)
    
      do l=1,3
      do m=1,3
        sum0tauA(m,l,i,j,k) = (0.d0,0.d0)
        do mm=1,4*NBS
           do nn=1,4*NBS
               sum0tauA(m,l,i,j,k) = sum0tauA(m,l,i,j,k) +tmp_tau_Arad_Qmat(m,l,nn,mm)*calE(nn,mm)
            end do
         end do
      end do
      end do

    else if(unitn(ii).eq.72) then
      vecR(1)=x; vecR(2)=y; vecR(3)=z
      do l=1,3
         do mm=1,4*NBS
            do nn=1,4*NBS
               E_Qmat(l,nn,mm) = tmp_E_Qmat(l,i,j,k,nn,mm)
            end do
         end do
      end do
      call calc_pol_BO(vecR,calE,E_Qmat,tmp_pol)
      do l=1,3
         sum0pol(l,i,j,k)=tmp_pol(l)
      end do

    else if(unitn(ii).eq.73) then
      vecR(1)=x; vecR(2)=y; vecR(3)=z
      do l=1,3
         sum0ef(l,i,j,k) = dAraddt_4E(l,0.d0,vecR,ALPHA_COH)
      end do

    else if(unitn(ii).eq.74) then
      do l=1,3
      do m=1,3
         sum0tauAM(m,l,i,j,k) = (0.d0,0.d0)
         do mm=1,4*NBS
            do nn=1,4*NBS
               sum0tauAM(m,l,i,j,k) = sum0tauAM(m,l,i,j,k) +tmp_tauAM(m,l,i,j,k,nn,mm)*calE(nn,mm)
            end do
         end do
      end do
      end do

    else if(unitn(ii).eq.75) then
      do l=1,3
         sum0t_AM(l,i,j,k) = (0.d0,0.d0)
         do mm=1,4*NBS
            do nn=1,4*NBS
               sum0t_AM(l,i,j,k) = sum0t_AM(l,i,j,k) +tmp_t_AM(l,i,j,k,nn,mm)*calE(nn,mm)
            end do
         end do
      end do

    else if(unitn(ii).eq.77) then
      do l=1,3
         sum0pi(l,i,j,k) = (0.d0,0.d0)
         do mm=1,4*NBS
            do nn=1,4*NBS
               sum0pi(l,i,j,k) = sum0pi(l,i,j,k) +tmp_pi(l,i,j,k,nn,mm)*calE(nn,mm)
            end do
         end do
      end do

    else if(unitn(ii).eq.78) then
      do l=1,3
         sum0orbang(l,i,j,k) = (0.d0,0.d0)
         do mm=1,4*NBS
            do nn=1,4*NBS
               sum0orbang(l,i,j,k) = sum0orbang(l,i,j,k) +tmp_orbang(l,i,j,k,nn,mm)*calE(nn,mm)
            end do
         end do
      end do

    else if(unitn(ii).eq.79) then
      do l=1,3
         sum0rots(l,i,j,k) = (0.d0,0.d0)
         do mm=1,4*NBS
            do nn=1,4*NBS
               sum0rots(l,i,j,k) = sum0rots(l,i,j,k) +tmp_rots(l,i,j,k,nn,mm)*calE(nn,mm)
            end do
         end do
      end do

    else if(unitn(ii).ne.0) then
      write(*,*) 'error unitn(',ii,')=',unitn(ii)

    end if !unitn(ii)

  end do !ii
    
end subroutine prop_calc_ini
!====================================================================================

!====================================================================================
subroutine prop_calc(i,j,k,x,y,z,calE,time)
!====================================================================================
  
  use Precision
  Use DiracOutput
  use Constants
  use prop_param
  use prop_tmp
  implicit none

  integer, intent(in) :: i,j,k
  complex(kind=dp), intent(in) :: calE(4*NBS,4*NBS)
  integer :: l,m,ii
  integer :: nn,mm
  character(LEN=80) :: FMT1, FMT3, FMT9
  character(LEN=80) :: FMT1R, FMT3R, FMT9R
  complex(kind=dp) :: sumd, sum_d, sumc
  complex(kind=dp) :: sumj(3)
  complex(kind=dp) :: sumt(3), sums(3), sumz(3)
  complex(kind=dp) :: sumtA(3),sumA(3)
  complex(kind=dp) :: Arad_vec
  complex(kind=dp) :: tmp_Arad_vec(3)
  complex(kind=dp) :: sumdivj
  complex(kind=dp) :: sumtau(3,3)
  complex(kind=dp) :: sumtauA(3,3)
  complex(kind=dp) :: sumtauAM(3,3)
  complex(kind=dp) :: sumt_AM(3)
  complex(kind=dp) :: sumpi(3), sumorbang(3),sumrots(3)
  real(kind=dp) :: sumpol(3), sumef(3)

  real(kind=dp),intent(in) :: x,y,z,time
  complex(kind=dp) :: tmp_tau_Arad_Qmat(3,3,4*NBS,4*NBS)
  real(kind=dp) :: vecR(3)
  complex(kind=dp) :: E_Qmat(3,4*NBS,4*NBS)  ! E^k_NM
  real(kind=dp) :: tmp_pol(3)
  real(kind=dp) :: dAraddt_4E
  


  FMT1 = '(10es24.14)'
  FMT3 = '(22es24.14)'
  FMT3R = '(13es24.14)'
  FMT9 = '(58es24.14)'

  do ii=1,size(unitn)

    if(unitn(ii).eq.61) then
      sumd = (0.d0,0.d0)
      do mm=1,4*NBS
         do nn=1,4*NBS
            sumd = sumd +tmp_d(i,j,k,nn,mm)*calE(nn,mm)
         end do
      end do
      write(unitn(ii),FMT1) time, x, y, z, sumd, sumd-sum0d(i,j,k), sum0d(i,j,k)

    else if(unitn(ii).eq.62) then
      sumc = (0.d0,0.d0)
      do mm=1,4*NBS
         do nn=1,4*NBS
            sumc = sumc +tmp_c(i,j,k,nn,mm)*calE(nn,mm)
         end do
      end do
      write(unitn(ii),FMT1) time, x, y, z, sumc, sumc-sum0c(i,j,k), sum0c(i,j,k)
    
    else if(unitn(ii).eq.63) then
      sums(:) = (0.d0,0.d0)
      do l=1,3
         do mm=1,4*NBS
            do nn=1,4*NBS
               sums(l) = sums(l) +tmp_s(l,i,j,k,nn,mm)*calE(nn,mm)
            end do
         end do
      end do
      write(unitn(ii),FMT3) time, x, y, z, sums(1), sums(1)-sum0s(1,i,j,k)&
                                        &, sums(2), sums(2)-sum0s(2,i,j,k)&
                                        &, sums(3), sums(3)-sum0s(3,i,j,k)&
                                        &, sum0s(1,i,j,k), sum0s(2,i,j,k), sum0s(3,i,j,k)
    
    else if(unitn(ii).eq.64) then
      sumt(:) = (0.d0,0.d0)
      do l=1,3
         do mm=1,4*NBS
            do nn=1,4*NBS
               sumt(l) = sumt(l) +tmp_t(l,i,j,k,nn,mm)*calE(nn,mm) !spintorque
            end do
         end do
      end do
      write(unitn(ii),FMT3) time, x, y, z, sumt(1), sumt(1)-sum0t(1,i,j,k)&
                                        &, sumt(2), sumt(2)-sum0t(2,i,j,k)&
                                        &, sumt(3), sumt(3)-sum0t(3,i,j,k)&
                                        &, sum0t(1,i,j,k), sum0t(2,i,j,k), sum0t(3,i,j,k)
    
    else if(unitn(ii).eq.65) then
      sumz(:) = (0.d0,0.d0)
      do l=1,3
         do mm=1,4*NBS
            do nn=1,4*NBS
               sumz(l) = sumz(l) +tmp_z(l,i,j,k,nn,mm)*calE(nn,mm)
            end do
         end do
      end do
      write(unitn(ii),FMT3) time, x, y, z, sumz(1), sumz(1)-sum0z(1,i,j,k)&
                                        &, sumz(2), sumz(2)-sum0z(2,i,j,k)&
                                        &, sumz(3), sumz(3)-sum0z(3,i,j,k)&
                                        &, sum0z(1,i,j,k), sum0z(2,i,j,k), sum0z(3,i,j,k)
    
    else if(unitn(ii).eq.66) then
      !$omp parallel do
      do l=1,3
        tmp_Arad_vec(l) = Arad_vec(time,l,x,y,z)
        do mm=1,4*NBS
           do nn=1,4*NBS
            tmp2_j_Qmat(l,nn,mm) = tmp_j(l,i,j,k,nn,mm)
          end do
        end do
      end do
      !$omp end parallel do
    
      call calc_t_Arad_Qmat(tmp2_j_Qmat,tmp_Arad_vec,tmp_t_Arad_Qmat)
    
      sumA(:) = (0.d0,0.d0)
      do l=1,3
        do mm=1,4*NBS
           do nn=1,4*NBS
               sumA(l) = sumA(l) +tmp_t_Arad_Qmat(l,nn,mm)*calE(nn,mm)
            end do
         end do
      end do
      write(unitn(ii),FMT3) time, x, y, z, sumA(1), sumA(1)-sum0A(1,i,j,k)&
                                        &, sumA(2), sumA(2)-sum0A(2,i,j,k)&
                                        &, sumA(3), sumA(3)-sum0A(3,i,j,k)&
                                        &, sum0A(1,i,j,k), sum0A(2,i,j,k), sum0A(3,i,j,k)
    
    else if(unitn(ii).eq.67) then
      sumtA(:) = (0.d0,0.d0)
      do l=1,3
        do mm=1,4*NBS
           do nn=1,4*NBS
               sumtA(l) = sumtA(l) +(tmp_t(l,i,j,k,nn,mm)+tmp_t_Arad_Qmat(l,nn,mm))*calE(nn,mm) !spintorque with Arad effect
            end do
         end do
      end do
      write(unitn(ii),FMT3) time, x, y, z, sumtA(1), sumtA(1)-sum0tA(1,i,j,k)&
                                        &, sumtA(2), sumtA(2)-sum0tA(2,i,j,k)&
                                        &, sumtA(3), sumtA(3)-sum0tA(3,i,j,k)&
                                        &, sum0tA(1,i,j,k), sum0tA(2,i,j,k), sum0tA(3,i,j,k)
    
    else if(unitn(ii).eq.68) then
      sumj(:) = (0.d0,0.d0)
      do l=1,3
         do mm=1,4*NBS
            do nn=1,4*NBS
               sumj(l) = sumj(l) +tmp_j(l,i,j,k,nn,mm)*calE(nn,mm)
            end do
         end do
      end do
      write(unitn(ii),FMT3) time, x, y, z, sumj(1), sumj(1)-sum0j(1,i,j,k)&
                                        &, sumj(2), sumj(2)-sum0j(2,i,j,k)&
                                        &, sumj(3), sumj(3)-sum0j(3,i,j,k)&
                                        &, sum0j(1,i,j,k), sum0j(2,i,j,k), sum0j(3,i,j,k)
    else if(unitn(ii).eq.69) then
      sumdivj = (0.d0,0.d0)
      do mm=1,4*NBS
         do nn=1,4*NBS
            sumdivj = sumdivj +tmp_divj(i,j,k,nn,mm)*calE(nn,mm)
         end do
      end do
      write(unitn(ii),FMT1) time, x, y, z, sumdivj, sumdivj-sum0divj(i,j,k), sum0divj(i,j,k)

    else if(unitn(ii).eq.70) then
      sumtau(:,:) = (0.d0,0.d0)
      do l=1,3
      do m=1,3
         do mm=1,4*NBS
            do nn=1,4*NBS
               sumtau(m,l) = sumtau(m,l) +tmp_tau(m,l,i,j,k,nn,mm)*calE(nn,mm)
            end do
         end do
      end do
      end do
      write(unitn(ii),FMT9) time, x, y, z, sumtau(1,1), sumtau(1,1)-sum0tau(1,1,i,j,k)&
                                        &, sumtau(1,2), sumtau(1,2)-sum0tau(1,2,i,j,k)&
                                        &, sumtau(1,3), sumtau(1,3)-sum0tau(1,3,i,j,k)&
                                        &, sumtau(2,1), sumtau(2,1)-sum0tau(2,1,i,j,k)&
                                        &, sumtau(2,2), sumtau(2,2)-sum0tau(2,2,i,j,k)&
                                        &, sumtau(2,3), sumtau(2,3)-sum0tau(2,3,i,j,k)&
                                        &, sumtau(3,1), sumtau(3,1)-sum0tau(3,1,i,j,k)&
                                        &, sumtau(3,2), sumtau(3,2)-sum0tau(3,2,i,j,k)&
                                        &, sumtau(3,3), sumtau(3,3)-sum0tau(3,3,i,j,k)&
                                        &, sum0tau(1,1,i,j,k), sum0tau(1,2,i,j,k), sum0tau(1,3,i,j,k)&
                                        &, sum0tau(2,1,i,j,k), sum0tau(2,2,i,j,k), sum0tau(2,3,i,j,k)&
                                        &, sum0tau(3,1,i,j,k), sum0tau(3,2,i,j,k), sum0tau(3,3,i,j,k)

    else if(unitn(ii).eq.71) then
      !$omp parallel do
      do l=1,3
        tmp_Arad_vec(l) = Arad_vec(time,l,x,y,z)
        do mm=1,4*NBS
           do nn=1,4*NBS
            tmp2_j_Qmat(l,nn,mm) = tmp_j(l,i,j,k,nn,mm)
          end do
        end do
      end do
      !$omp end parallel do
    
      call calc_tau_Arad_Qmat(tmp2_j_Qmat,tmp_Arad_vec,tmp_tau_Arad_Qmat)
    
      sumtauA(:,:) = (0.d0,0.d0)
      do l=1,3
      do m=1,3
        do mm=1,4*NBS
           do nn=1,4*NBS
               sumtauA(m,l) = sumtauA(m,l) +tmp_tau_Arad_Qmat(m,l,nn,mm)*calE(nn,mm)
            end do
         end do
      end do
      end do
      write(unitn(ii),FMT9) time, x, y, z, sumtauA(1,1), sumtauA(1,1)-sum0tauA(1,1,i,j,k)&
                                        &, sumtauA(1,2), sumtauA(1,2)-sum0tauA(1,2,i,j,k)&
                                        &, sumtauA(1,3), sumtauA(1,3)-sum0tauA(1,3,i,j,k)&
                                        &, sumtauA(2,1), sumtauA(2,1)-sum0tauA(2,1,i,j,k)&
                                        &, sumtauA(2,2), sumtauA(2,2)-sum0tauA(2,2,i,j,k)&
                                        &, sumtauA(2,3), sumtauA(2,3)-sum0tauA(2,3,i,j,k)&
                                        &, sumtauA(3,1), sumtauA(3,1)-sum0tauA(3,1,i,j,k)&
                                        &, sumtauA(3,2), sumtauA(3,2)-sum0tauA(3,2,i,j,k)&
                                        &, sumtauA(3,3), sumtauA(3,3)-sum0tauA(3,3,i,j,k)&
                                        &, sum0tauA(1,1,i,j,k), sum0tauA(1,2,i,j,k), sum0tauA(1,3,i,j,k)&
                                        &, sum0tauA(2,1,i,j,k), sum0tauA(2,2,i,j,k), sum0tauA(2,3,i,j,k)&
                                        &, sum0tauA(3,1,i,j,k), sum0tauA(3,2,i,j,k), sum0tauA(3,3,i,j,k)

    else if(unitn(ii).eq.72) then
      vecR(1)=x; vecR(2)=y; vecR(3)=z
      do l=1,3
         do mm=1,4*NBS
            do nn=1,4*NBS
               E_Qmat(l,nn,mm) = tmp_E_Qmat(l,i,j,k,nn,mm)
            end do
         end do
      end do
      call calc_pol_BO(vecR,calE,E_Qmat,tmp_pol)
      do l=1,3
         sumpol(l)=tmp_pol(l)
      end do

      write(unitn(ii),FMT3R) time, x, y, z, sumpol(1), sumpol(1)-sum0pol(1,i,j,k)&
                                         &, sumpol(2), sumpol(2)-sum0pol(2,i,j,k)&
                                         &, sumpol(3), sumpol(3)-sum0pol(3,i,j,k)&
                                         &, sum0pol(1,i,j,k), sum0pol(2,i,j,k), sum0pol(3,i,j,k)
    
    else if(unitn(ii).eq.73) then
      vecR(1)=x; vecR(2)=y; vecR(3)=z
      do l=1,3
         sumef(l) = dAraddt_4E(l,time,vecR,ALPHA_COH)
      end do

      write(unitn(ii),FMT3R) time, x, y, z, sumef(1), sumef(1)-sum0ef(1,i,j,k)&
                                         &, sumef(2), sumef(2)-sum0ef(2,i,j,k)&
                                         &, sumef(3), sumef(3)-sum0ef(3,i,j,k)&
                                         &, sum0ef(1,i,j,k), sum0ef(2,i,j,k), sum0ef(3,i,j,k)

    else if(unitn(ii).eq.74) then
      sumtauAM(:,:) = (0.d0,0.d0)
      do l=1,3
      do m=1,3
         do mm=1,4*NBS
            do nn=1,4*NBS
               sumtauAM(m,l) = sumtauAM(m,l) +tmp_tauAM(m,l,i,j,k,nn,mm)*calE(nn,mm)
            end do
         end do
      end do
      end do
      write(unitn(ii),FMT9) time, x, y, z, sumtauAM(1,1), sumtauAM(1,1)-sum0tauAM(1,1,i,j,k)&
                                        &, sumtauAM(1,2), sumtauAM(1,2)-sum0tauAM(1,2,i,j,k)&
                                        &, sumtauAM(1,3), sumtauAM(1,3)-sum0tauAM(1,3,i,j,k)&
                                        &, sumtauAM(2,1), sumtauAM(2,1)-sum0tauAM(2,1,i,j,k)&
                                        &, sumtauAM(2,2), sumtauAM(2,2)-sum0tauAM(2,2,i,j,k)&
                                        &, sumtauAM(2,3), sumtauAM(2,3)-sum0tauAM(2,3,i,j,k)&
                                        &, sumtauAM(3,1), sumtauAM(3,1)-sum0tauAM(3,1,i,j,k)&
                                        &, sumtauAM(3,2), sumtauAM(3,2)-sum0tauAM(3,2,i,j,k)&
                                        &, sumtauAM(3,3), sumtauAM(3,3)-sum0tauAM(3,3,i,j,k)&
                                        &, sum0tauAM(1,1,i,j,k), sum0tauAM(1,2,i,j,k), sum0tauAM(1,3,i,j,k)&
                                        &, sum0tauAM(2,1,i,j,k), sum0tauAM(2,2,i,j,k), sum0tauAM(2,3,i,j,k)&
                                        &, sum0tauAM(3,1,i,j,k), sum0tauAM(3,2,i,j,k), sum0tauAM(3,3,i,j,k)
    
    else if(unitn(ii).eq.75) then
      sumt_AM(:) = (0.d0,0.d0)
      do l=1,3
         do mm=1,4*NBS
            do nn=1,4*NBS
               sumt_AM(l) = sumt_AM(l) +tmp_t_AM(l,i,j,k,nn,mm)*calE(nn,mm)
            end do
         end do
      end do
      write(unitn(ii),FMT3) time, x, y, z, sumt_AM(1), sumt_AM(1)-sum0t_AM(1,i,j,k)&
                                        &, sumt_AM(2), sumt_AM(2)-sum0t_AM(2,i,j,k)&
                                        &, sumt_AM(3), sumt_AM(3)-sum0t_AM(3,i,j,k)&
                                        &, sum0t_AM(1,i,j,k), sum0t_AM(2,i,j,k), sum0t_AM(3,i,j,k)
    
    else if(unitn(ii).eq.77) then
      sumpi(:) = (0.d0,0.d0)
      do l=1,3
         do mm=1,4*NBS
            do nn=1,4*NBS
               sumpi(l) = sumpi(l) +tmp_pi(l,i,j,k,nn,mm)*calE(nn,mm)
            end do
         end do
      end do
      write(unitn(ii),FMT3) time, x, y, z, sumpi(1), sumpi(1)-sum0pi(1,i,j,k)&
                                        &, sumpi(2), sumpi(2)-sum0pi(2,i,j,k)&
                                        &, sumpi(3), sumpi(3)-sum0pi(3,i,j,k)&
                                        &, sum0pi(1,i,j,k), sum0pi(2,i,j,k), sum0pi(3,i,j,k)
    
    else if(unitn(ii).eq.78) then
      sumorbang(:) = (0.d0,0.d0)
      do l=1,3
         do mm=1,4*NBS
            do nn=1,4*NBS
               sumorbang(l) = sumorbang(l) +tmp_orbang(l,i,j,k,nn,mm)*calE(nn,mm)
            end do
         end do
      end do
      write(unitn(ii),FMT3) time, x, y, z, sumorbang(1), sumorbang(1)-sum0orbang(1,i,j,k)&
                                        &, sumorbang(2), sumorbang(2)-sum0orbang(2,i,j,k)&
                                        &, sumorbang(3), sumorbang(3)-sum0orbang(3,i,j,k)&
                                        &, sum0orbang(1,i,j,k), sum0orbang(2,i,j,k), sum0orbang(3,i,j,k)
    
    else if(unitn(ii).eq.79) then
      sumrots(:) = (0.d0,0.d0)
      do l=1,3
         do mm=1,4*NBS
            do nn=1,4*NBS
               sumrots(l) = sumrots(l) +tmp_rots(l,i,j,k,nn,mm)*calE(nn,mm)
            end do
         end do
      end do
      write(unitn(ii),FMT3) time, x, y, z, sumrots(1), sumrots(1)-sum0rots(1,i,j,k)&
                                        &, sumrots(2), sumrots(2)-sum0rots(2,i,j,k)&
                                        &, sumrots(3), sumrots(3)-sum0rots(3,i,j,k)&
                                        &, sum0rots(1,i,j,k), sum0rots(2,i,j,k), sum0rots(3,i,j,k)

    else if(unitn(ii).ne.0) then
      write(*,*) 'error unitn(',ii,')=',unitn(ii)
    end if !unit(ii)
    
    
  end do !ii

end subroutine prop_calc
!====================================================================================

!====================================================================================
subroutine set_paramini
!====================================================================================
  use Precision
  Use DiracOutput
  use Constants
  use IntegralStorage
  use NucBasis

  implicit none

  integer i,j

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

  logical :: use_Sch  ! .true. --> compute in Schwarzschild spacetime, .false. --> flat spacetime
  real(kind=dp) :: R_Sch ! radius at which computation is performed (normalized by Schwarzshild radius)
  real(kind=dp) :: Mass ! mass of the central objects in units of solar mass

  logical :: use_retarded  ! .true. --> include retarded potential, .false. --> ignore retarded potential
  logical :: use_mass_renorm  ! .true. --> perform mass renormalization, .false. --> without mass renormalization

  real(kind=dp) :: alphaJJ


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
  !--------------------------------------------------------------
  ! paramters which are previously given in qed.inp
  read(1000,*)
  read (1000,'(a)') FILESFOLDER
  read(1000,*)
  read (1000,'(a)') FIFI3 !Dirac output file of Dirac10
  read (1000,*) SYMM ! 0=>atom, 1=>molecule using symmmetry calc, 2=>molecule not using symmetry calc
  read (1000,*) KBAL
  !--------------------------------------------------------------
  read(1000,*)
  read(1000,*) outfile
  open(unit=10,file=outfile,status='unknown',form='formatted')  ! output file
  read(1000,*) NTIME
  read(1000,*) DeltaT   ! defined as global parameter
  read(1000,*)  !nprint
  read(1000,*) NEL  !  NEL = 1 ! H,H2+
                    !  NEL = 2 ! He,H2
                    !  NEL = 6 ! Li2
  read(1000,*) 
  read(1000,*) 
  read(1000,*) 
  read(1000,*) 
  read(1000,*) 
  read(1000,*) 

  read(1000,*) 
  read(1000,*) there_is_twoele
  !  there_is_twoele = .True. ! true : there is already two electron integrals and have a name defined in module IntegralStorage
  !  there_is_twoele = .false. ! false : compute two electron integrals and store as fort.12
  if(there_is_twoele) then
     read(1000,*) file_twoele ! defined in module IntegralStorage
  end if
!!$  read(1000,*) there_is_intIJJ
!!$  if(there_is_intIJJ) then
!!$     read(1000,*) file_intIJJ ! defined in module IntegralStorage
!!$  end if
  read(1000,*) there_is_intKJJ
  if(there_is_intKJJ) then
     read(1000,*) file_intKJJ ! defined in module IntegralStorage
  else
     read(1000,*) AlphaJJ_MAX
     read(1000,*) N_alphaJJ
     read(1000,*) Diff_u_t
  end if

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
  read(1000,*) use_retarded
  ! .true. --> include retarded potential, .false. --> ignore retarded potential

  read(1000,*) 
  read(1000,*) use_mass_renorm
  ! .true. --> perform mass renormalization, .false. --> without mass renormalization

  read(1000,*) 
  read(1000,*) use_Sch
  ! .true. --> compute in Schwarzschild spacetime, .false. --> flat spacetime
  if((use_Sch .eq. .true.).and.(use_BO_approx .eq. .false.)) then
     write(*,*) "Computation in curved spacetime is only supported for BO approximation."
     stop
  end if
  read(1000,*) R_Sch
  read(1000,*) Mass
  read(1000,*) NOPT_GRw, NOPT_GR, NOPT_GR_T, NOPT_GR_S, NOPT_GR_M, NOPT_GR_V  ! option for GR calculation 140122, module Constants
!  read(1000,*) NOPT_GR  ! option for GR calculation 140122, module Constants
  if((use_Sch .eq. .true.).and.(R_Sch.le.1._dp)) then
     write(*,*) "R_Sch should be larger than 1."
     stop
  end if

  read(1000,*) 
  read(1000,*) use_coherent
  ! .true. --> use coherent state for photon, .false. --> no radiation in the initial state.
!  if((use_coherent .eq. .true.).and.(use_BO_approx .eq. .false.)) then
!     write(*,*) "Coherent state is only supported for BO approximation."
!     stop
!  end if

  if(use_Sch) then
     if(use_coherent) then
        write(*,*) "Coherent state is only supported for flat spacetime."
        stop
     elseif(use_retarded) then
        write(*,*) "Retarded potential is only supported for flat spacetime."
        stop
     end if
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

     if(n_coherent.gt.Nph) then
        write(*,*) "n_coherent should be smaller than Nph."
        stop
     end if
     do i=1,n_coherent
        read(1000,*) i_coh
        read(1000,*) alpha(i_coh)
     end do
     if(n_coherent.eq.0) then
        write(*,*) "# Self-energy computation."

        if(use_BO_approx .eq. .false.) then
           if(use_mass_renorm) then
              write(*,*) "Mass renormalization is only supported for BO approximation."
              stop
           elseif(use_retarded) then
              write(*,*) "Retarded potential is only supported for BO approximation."
              stop
           end if
        end if

        if(use_mass_renorm) then
           write(*,*) "# With mass renormalization."
        else
           write(*,*) "# Without mass renormalization."
        end if

        if(use_retarded) then
           write(*,*) "Retarded potential is not supported for self-energy computation."
           stop
        end if

     end if
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
  if(use_retarded) then
     write(*,*) "# Use retarded potential."
  else  
     write(*,*) "# Without retarded potential."
  end if
  if(use_Sch) then
     write(*,*) "# Compute in Schwarzschild spacetime."
  else  
     write(*,*) "# Compute in flat spacetime."
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

  if(use_Sch) then
     write(*,*) "#############################################"
     write(*,*) "Parameters for computation in curved spacetime."
     write(*,*) "#############################################"
     write(*,"(1a15,1es16.6)") "# R_Sch : ",R_Sch
     write(*,"(1a15,1es16.6)") "# Mass : ",Mass
     write(*,"(1a15,1i6)") "# NOPT_GR : ", NOPT_GR
     write(*,"(1a15,4i6)") "# NOPT_GR_TSMV : ", NOPT_GR_T, NOPT_GR_S, NOPT_GR_M, NOPT_GR_V
  end if

  if(use_coherent) then
     !---------------------------------------------------------------------------------------------
     ! set photon parameters & specify coherent states
     !---------------------------------------------------------------------------------------------
     write(*,*) "#############################################"
     write(*,*) "Set photon parameters."
     write(*,*) "#############################################"
     allocate(ALPHA_COH(Nph))
     ALPHA_COH(:) = alpha(:) 

     write(*,"(1a15,1es16.6)") "# P0MAX: ",P0MAX
     write(*,"(1a15,1i6)") "# NP0: ",NP0
     write(*,"(1a15,1i6)") "# NPTH: ",NPTH
     write(*,"(1a15,1i6)") "# NPPHI: ",NPPHI
     write(*,"(1a15,1i6)") "# Nph: ",Nph
     
     write(*,*) "# Photon momentum discretization"
     write(*,*) "# Polar coordinates: th measured from z-axis, include both 0 and PI."
     write(*,*) "#                  : phi measured from x-axis, include include 0 but not 2PI."

     dp0  = P0MAX/NP0 ! P0MAX specified in GammaMatrix. We do not use p0=0
     dth  = PI/(NPTH-1) ! NPTH specified in GammaMatrix. th include both 0 and PI
     dphi = 2.d0*PI/NPPHI ! NPPHI specified in GammaMatrix. phi include 0 but not 2PI

     write(*,"(1a3,1a6,1a5,1a16,2a16)") "#  "," j ","pol"," p0      "," th/PI ", " phi/PI "
     do j=1,Nph
        call convert_label_2(j,NP0,NPTH,NPPHI,ip0,ith,iphi,sig)  ! defined in sub_int.f90
        p0 = dp0*ip0
        th = dth*(ith-1)  ! start from 0
        phi= dphi*(iphi-1) ! start from 0
        write(*,"(1a3,1i6,1a5,1es16.6,2f16.4)") "# ", j, sig, p0, th/PI, phi/PI
     end do
     
     write(*,*) "# Print non-zero alpha"
     j=0
     do i=1,Nph
        if(ALPHA_COH(i).ne.(0._dp,0._dp)) then
           write(*,"(1a8,1i4,2es16.6)") "# alpha",i,ALPHA_COH(i)
           j=j+1
        end if
     end do
     if(j.ne.n_coherent) then
        write(*,*) "n_coherent does not match. Check input file."
        stop
     end if
  end if

end subroutine set_paramini
!====================================================================================
!====================================================================================
subroutine set_paramini_tdpt(NEL,NTIME,nprint)
!====================================================================================
  use Precision
  Use DiracOutput
  use Constants
  use IntegralStorage
  use NucBasis

  implicit none

  integer i,j

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
  read(1000,*) nprint
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

end subroutine set_paramini_tdpt
!====================================================================================
