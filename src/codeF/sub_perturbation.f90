! Last Change:08-Jun-2014.
!====================================================================================
! 2013.10.29 
! - subroutine read_ptHFinp
! - subroutine perturbation
! - subroutine renorm_pt
! - subroutine renorm_pt1(c_p0,c_p1,norm)
! - subroutine renorm_pt2(c_p0,c_p1,c_p2,norm)
! - subroutine set_psi_pt0
! - subroutine set_psi_pt1_old
! - subroutine set_psi_pt1
! - subroutine set_psi_pt2
! - function intHrel_pt_Qmat(c_p,c_q)
! - subroutine pm12(a,aa)
! - subroutine energy_sort(eneQ)
!============================================================================
subroutine read_ptHFinp
!read input file (qed.inp)
!============================================================================
  use DiracOutput
  use PTcoef
  use param_AM

  implicit none

  character(len=80) INPUT

!  INPUT = 'qed.inp'
!  open(unit=11,file=INPUT)
     read (*,'(a)') FILESFOLDER
     FIFI1 = trim(FILESFOLDER)//"/"//"basis.txt"
     FIFI2 = trim(FILESFOLDER)//"/"//"vectors.txt"
!     read (*,'(a12)') FIFI1 ! basis.txt of Dirac10
!     read (*,'(a12)') FIFI2 ! vectors.txt of Dirac10
     write(*,'(a)') FIFI1
     write(*,'(a)') FIFI2
     read (*,'(a)') FIFI3 !Dirac output file of Dirac10
     FIFI3 = trim(FILESFOLDER)//"/"//FIFI3
     write(*,'(a)') FIFI3
     read (*,*) SYMM ! 0=>atom, 1=>molecule using symmmetry calc, 2=>molecule not using symmetry calc
     read (*,*) KBAL
     read (*,*) u1,u2,u3,uu,u4,u5,u6
     read (*,*) NEL
     !-- magnetic field for initial state
     read (*,*) Bmag(1)  !|vecBM_x|
     read (*,*) Bmag(2)  !|vecBM_y|
     read (*,*) Bmag(3)  !|vecBM_z|
!     read (*,*) direcB ! 1->x, 2->y, 3->z 
!     read (*,*) Bmag  !|vecBM|
!     !-- magnetic field for additional torq
!     read (*,*) direcB2 ! 1->x, 2->y, 3->z 
!     read (*,*) Bmag2  !|vecBM|
!     stop
!  close(unit=11)

end subroutine read_ptHFinp

!============================================================================
subroutine perturbation
!============================================================================
  use Precision
  Use DiracOutput
  use Constants
  use DefineTypes
  use prop_param
  use PTcoef

  implicit none

  character(LEN=80) :: filedens, filetorque, filespin, filespind, filezeta, filechiral, filetz
  character(LEN=80) :: filej, filetAM, filetAMt

  real(kind=8) :: x,y,z!,dx,dy,dz

  integer :: i,j,k,l
  integer :: rr,pp,qq
  integer :: nn,mm
  integer :: actorb0, actorb

  complex(kind=8) :: density_pt_Qmat, rho_pt_Qmat
  complex(kind=8) :: t_pt_Qmat, t_Arad_pt_Qmat, s_pt_Qmat, zeta_pt_Qmat, chiral_pt_Qmat
  complex(kind=8) :: j_pt_Qmat, t_AM_pt_Qmat
  complex(kind=8) :: intSpin_pt_Qmat

  complex(kind=8),allocatable :: sum0d(:,:,:), sum0c(:,:,:), sum0t(:,:,:,:), sum0s(:,:,:,:), sum0z(:,:,:,:)
  complex(kind=8),allocatable :: sum0j(:,:,:,:), sum0tAM(:,:,:,:)
  complex(kind=8),allocatable :: tmp_d(:,:,:,:,:), tmp_c(:,:,:,:,:), tmp_t(:,:,:,:,:,:), tmp_s(:,:,:,:,:,:), tmp_z(:,:,:,:,:,:)

!  integer :: NEL
  complex(kind=dp) :: intN_pt_Qmat
  double precision,allocatable :: eneQ(:)

  !---------------------------------------------------------------------------------------------
  ! set electron related parameters
  !---------------------------------------------------------------------------------------------
  write(*,*) "#############################################"
  write(*,*) "Set expansion functions for electrons."
  write(*,*) "#############################################"

  call read_ptHFinp ! read qed.inp
  call read_DiracOutput  ! read and set global variables -> defined in sub_readDiracOutput.f90
  call set_GammaMatrix

!  NEL = 1 !H atom

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

  allocate(c_psi(4,NMAX_PG,2*NBS))
  allocate(eneQ0(2*NBS),eneQ1(2*NBS),eneQ2(2*NBS))
  allocate(eneQ(2*NBS),occup(2*NBS))

    actorb0=1
    actorb=NEL
!    actorb0=2
!    actorb=2

!!$    occup(:) = 0.d0
!!$    do i=actorb0,actorb
!!$      occup(i) = 1.d0
!!$    end do

    write(*,*)" # set psi pt0"
    call set_psi_pt0

!!$    write(*,*)" # set psi pt1"
!!$    call set_psi_pt1
!!$
!!$    open(unit=101,file="energy.dat")
!!$    write(101,*) '# orbital, eneQ, eneQ0, eneQ1'
!!$    do i=1,2*NBS
!!$      eneQ(i) = eneQ0(i) + eneQ1(i)
!!$      write(101,"(1i6,3es16.6)") i, eneQ(i), eneQ0(i), eneQ1(i)
!!$    end do
!!$    close(unit=101)
!!$    
!!$    call energy_sort(eneQ)
!!$    write(*,*)" # call energy_sort "
!!$    open(unit=101,file="energy_sort.dat")
!!$    write(101,*) '# orbital, eneQ, eneQ0, eneQ1'
!!$    do i=1,2*NBS
!!$      write(101,"(1i6,1es16.6)") i, eneQ(i)
!!$    end do
!!$    close(unit=101)

!!$    write(*,*)" # set psi pt2"
!!$    call set_psi_pt2
!!$    stop
!!$
!!$    open(unit=101,file="energy.dat")
!!$    write(101,*) '# orbital, eneQ, eneQ0+eneQ1, eneQ0, eneQ1, eneQ2'
!!$    do i=1,2*NBS
!!$      eneQ(i) = eneQ0(i) + eneQ1(i) + eneQ2(i)
!!$      write(101,"(1i6,5es16.6)") i, eneQ(i), eneQ0(i)+eneQ1(i), eneQ0(i), eneQ1(i), eneQ2(i)
!!$    end do
!!$    close(unit=101)

!    write(*,*)" # normalization"
!    call renorm_pt

!!$  do i=1,2*NBS
!!$     do j=1,2*NBS
!!$        write(*,'(2i6,8es16.6)')  i, j,intN_pt_Qmat(i,j)
!!$     end do
!!$  end do
!!$  stop

    open(unit=105,file="cpsi.dat")
    do i=1,2*NBS
      do j=1,NBS_L
        write(105,"(1a,2i6,2es22.14)")"1",i,j,c_psi(1,j,i)
        write(105,"(1a,2i6,2es22.14)")"2",i,j,c_psi(2,j,i)
      end do            
      do j=1,NBS_S      
        write(105,"(1a,2i6,2es22.14)")"3",i,j,c_psi(3,j,i)
        write(105,"(1a,2i6,2es22.14)")"4",i,j,c_psi(4,j,i)
      end do
    end do
    close(unit=105)

!!$    open(unit=105,file="cpsi.dat")
!!$    do i=1,2*NBS
!!$      do j=1,NBS_L
!!$        read(105,"(13x,2es22.14)") c_psi(1,j,i)
!!$        read(105,"(13x,2es22.14)") c_psi(2,j,i)
!!$      end do          
!!$      do j=1,NBS_S    
!!$        read(105,"(13x,2es22.14)") c_psi(3,j,i)
!!$        read(105,"(13x,2es22.14)") c_psi(4,j,i)
!!$      end do
!!$    end do
!!$    close(unit=105)

    open(unit=101,file="intSpin.dat")
    do nn=1,2*NBS
      write(101,"(1i6,6es16.6)") nn, intSpin_pt_Qmat(1,nn,nn), intSpin_pt_Qmat(2,nn,nn), intSpin_pt_Qmat(3,nn,nn)
    end do
    close(unit=101)

!    stop

  !--- start calc physical quantities -----------

!  !--- set magnetic field for additional torq --------
!    Bmag = Bmag2
!    Bdirec = Bdirec2

  !---calculation range and mesh ----------------------------------------
    ! If you want to calc point, specify mesh=1.
    meshx = int(U1/UU) +1
    meshy = int(U2/UU) +1
    meshz = int(U3/UU) +1
    xori = -UU*dble(meshx-1) +U4 ; xend = UU*dble(meshx-1) +U4
    yori = -UU*dble(meshy-1) +U5 ; yend = UU*dble(meshy-1) +U5
    zori = -UU*dble(meshz-1) +U6 ; zend = UU*dble(meshz-1) +U6
    write(*,*)'meshx =',meshx
    write(*,*)'meshy =',meshy
    write(*,*)'meshz =',meshz

    dx = 0.d0
    if(meshx.ne.1) dx = (xend - xori)/(meshx-1)
    dy = 0.d0
    if(meshy.ne.1) dy = (yend - yori)/(meshy-1)
    dz = 0.d0
    if(meshz.ne.1) dz = (zend - zori)/(meshz-1)
  
    write(*,*)'dx dy dz =',dx,dy,dz
  !----------------------------------------------------------------------
  allocate(sum0d(meshx,meshy,meshz),sum0t(3,meshx,meshy,meshz), sum0s(3,meshx,meshy,meshz), sum0z(3,meshx,meshy,meshz))
  allocate(sum0c(meshx,meshy,meshz),sum0j(3,meshx,meshy,meshz),sum0tAM(3,meshx,meshy,meshz))

  ! calculate properties
        write(filedens,'(a)') 'dens.dat'
        write(filetorque,'(a)') 'torq.dat'
        write(filespin,'(a)') 'spin.dat'
        write(filezeta,'(a)') 'zeta.dat'
        write(filechiral,'(a)') 'chir.dat'
        write(filetz,'(a)') 'tztAM.dat'
        write(filej,'(a)') 'j.dat'
        write(filetAM,'(a)') 'tAM.dat'
        write(filetAMt,'(a)') 'tAMt.dat'
        open(unit=17,file=filedens)
        open(unit=31,file=filezeta)
        open(unit=32,file=filetorque)
        open(unit=33,file=filespin)
        open(unit=35,file=filetz)
        open(unit=37,file=filechiral)
        open(unit=38,file=filej)
        open(unit=39,file=filetAM)
        open(unit=40,file=filetAMt)

        sum0d(:,:,:) = (0.d0,0.d0)
        sum0c(:,:,:) = (0.d0,0.d0)
        sum0t(:,:,:,:) = (0.d0,0.d0)
        sum0s(:,:,:,:) = (0.d0,0.d0)
        sum0z(:,:,:,:) = (0.d0,0.d0)
        sum0j(:,:,:,:) = (0.d0,0.d0)
        sum0tAM(:,:,:,:) = (0.d0,0.d0)

         do i=1,meshx
            x = xori + dx*(i-1)
 !           if(meshx.eq.1) x = x0
    
            do j=1,meshy
               y = yori + dy*(j-1)
 !              if(meshy.eq.1) y = y0
               
               do k=1,meshz
                  z = zori + dz*(k-1)
 !                 if(meshz.eq.1) z = z0
    
                  do nn=actorb0,actorb
                        sum0d(i,j,k) = sum0d(i,j,k) +density_pt_Qmat(x,y,z,nn,nn) !*occup(nn)
                        sum0c(i,j,k) = sum0c(i,j,k) +chiral_pt_Qmat(x,y,z,nn,nn) !*occup(nn)
                  end do
    
                  do l=1,3
                     do nn=actorb0,actorb
                           sum0s(l,i,j,k) = sum0s(l,i,j,k) +s_pt_Qmat(l,x,y,z,nn,nn) !*occup(nn)
                           sum0z(l,i,j,k) = sum0z(l,i,j,k) +zeta_pt_Qmat(l,x,y,z,nn,nn) !*occup(nn)
                           sum0t(l,i,j,k) = sum0t(l,i,j,k) +t_pt_Qmat(l,x,y,z,nn,nn) !*occup(nn)
                           sum0j(l,i,j,k) = sum0j(l,i,j,k) +j_pt_Qmat(l,x,y,z,nn,nn) !*occup(nn)
                           sum0tAM(l,i,j,k) = sum0tAM(l,i,j,k) +t_AM_pt_Qmat(0.d0,l,x,y,z,nn,nn) !*occup(nn)
                     end do
                  end do


!                        sum0d(i,j,k) = sum0d(i,j,k) +density_pt_Qmat(x,y,z,4,1) !*occup(1)
!                        sum0c(i,j,k) = sum0c(i,j,k) +chiral_pt_Qmat(x,y,z,4,1) !*occup(1)
!    
!                  do l=1,3
!                           sum0s(l,i,j,k) = sum0s(l,i,j,k) +s_pt_Qmat(l,x,y,z,4,1) !*occup(1)
!                           sum0z(l,i,j,k) = sum0z(l,i,j,k) +zeta_pt_Qmat(l,x,y,z,4,1) !*occup(1)
!                           sum0t(l,i,j,k) = sum0t(l,i,j,k) +t_pt_Qmat(l,x,y,z,4,1) !*occup(1)
!                           sum0j(l,i,j,k) = sum0j(l,i,j,k) +j_pt_Qmat(l,x,y,z,4,1) !*occup(1)
!                           sum0tAM(l,i,j,k) = sum0tAM(l,i,j,k) +t_AM_pt_Qmat(0.d0,l,x,y,z,4,1) !*occup(1)
!                  end do

                 write(17,'(5es24.14)') x, y, z, sum0d(i,j,k)
                 write(31,'(9es24.14)') x, y, z, sum0z(1,i,j,k), sum0z(2,i,j,k), sum0z(3,i,j,k)
                 write(32,'(9es24.14)') x, y, z, sum0t(1,i,j,k), sum0t(2,i,j,k), sum0t(3,i,j,k)
                 write(33,'(9es24.14)') x, y, z, sum0s(1,i,j,k), sum0s(2,i,j,k), sum0s(3,i,j,k)
                 write(35,'(9es24.14)') x, y, z, sum0z(1,i,j,k)+sum0t(1,i,j,k)+sum0tAM(1,i,j,k),&
                                               & sum0z(2,i,j,k)+sum0t(2,i,j,k)+sum0tAM(2,i,j,k), sum0z(3,i,j,k)+sum0t(3,i,j,k)+sum0tAM(3,i,j,k)
                 write(37,'(5es24.14)') x, y, z, sum0c(i,j,k)
                 write(38,'(9es24.14)') x, y, z, sum0j(1,i,j,k), sum0j(2,i,j,k), sum0j(3,i,j,k)
                 write(39,'(9es24.14)') x, y, z, sum0tAM(1,i,j,k), sum0tAM(2,i,j,k), sum0tAM(3,i,j,k)
                 write(40,'(9es24.14)') x, y, z, sum0tAM(1,i,j,k)+sum0t(1,i,j,k), sum0tAM(2,i,j,k)+sum0t(2,i,j,k), sum0tAM(3,i,j,k)+sum0t(3,i,j,k)
    
               end do !k=1,meshz
!                  write(unitn,*)
                  write(17,*)
                  write(31,*)
                  write(32,*)
                  write(33,*)
                  write(35,*)
                  write(37,*)
                  write(38,*)
                  write(39,*)
                  write(40,*)
            end do !j=1,meshy
!                write(unitn,*)
              write(17,*)
              write(31,*)
              write(32,*)
              write(33,*)
              write(35,*)
              write(37,*)
              write(38,*)
              write(39,*)
              write(40,*)
            write(*,*) i,'/',meshx
         end do !i=1,meshx
!            close(unit=unitn)
        close(unit=17)
        close(unit=31)
        close(unit=32)
        close(unit=33)
        close(unit=35)
        close(unit=37)
        close(unit=38)
        close(unit=39)
        close(unit=40)
!  stop
         write(*,*)'finish calc properties'
  !----------------------------------------------------------------

end subroutine perturbation

!============================================================================
subroutine renorm_pt
!============================================================================

  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef

  implicit none

  integer :: i,j,k
  integer :: actorb0, actorb
  complex(kind=dp) :: intN_pt_Qmat
  complex(kind=dp) :: tmp

  actorb0 = 1 ! All orbitals are activated.
  actorb = 2*NBS ! All orbitals are activated.
 
  !$omp parallel do private(tmp,j)
  do i=actorb0,actorb! i runs specified orbitals.
    tmp = 1.d0/dsqrt(dble(intN_pt_Qmat(i,i)))
    do j=1,NBS_L
      c_psi(1,j,i) = c_psi(1,j,i)*tmp
      c_psi(2,j,i) = c_psi(2,j,i)*tmp
    end do 
    do j=1,NBS_S
      c_psi(3,j,i) = c_psi(3,j,i)*tmp
      c_psi(4,j,i) = c_psi(4,j,i)*tmp
    end do 
  end do
  !$omp end parallel do

end subroutine renorm_pt

!============================================================================
subroutine renorm_pt1(c_p0,c_p1,norm)
!============================================================================

  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef

  implicit none

  complex(kind=dp),intent(inout) :: c_p0(4,NMAX_PG),c_p1(4,NMAX_PG)
  complex(kind=dp),intent(out) :: norm
  integer :: i,j,k
  complex(kind=dp) :: intN_pt_Qmat
  complex(kind=dp) :: tmp
  type(primitive_gaussian) :: pg
  complex(kind=dp) :: intN_00, intN_01, intN_10

  call copy_DiracOutput_pg(pg)
  
  intN_pt_Qmat = (0._dp,0._dp)
  do i=1,4
     call calc_intN_pq(i,i,NBS_L,NBS_S,pg,c_p0(:,:),c_p0(:,:),intN_00) 
!     write(*,*) i,intN_00
     intN_pt_Qmat = intN_pt_Qmat +intN_00

     call calc_intN_pq(i,i,NBS_L,NBS_S,pg,c_p0(:,:),c_p1(:,:),intN_01) 
!     write(*,*) i,intN_01
     intN_pt_Qmat = intN_pt_Qmat +intN_01

     call calc_intN_pq(i,i,NBS_L,NBS_S,pg,c_p1(:,:),c_p0(:,:),intN_10) 
!     write(*,*) i,intN_10
     intN_pt_Qmat = intN_pt_Qmat +intN_10
  end do
 
  norm = 1.d0/dsqrt(dble(intN_pt_Qmat))

end subroutine renorm_pt1

!============================================================================
subroutine renorm_pt2(c_p0,c_p1,c_p2,norm)
!============================================================================

  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef

  implicit none

  complex(kind=dp),intent(inout) :: c_p0(4,NMAX_PG),c_p1(4,NMAX_PG),c_p2(4,NMAX_PG)
  complex(kind=dp),intent(out) :: norm
  integer :: i,j,k
  complex(kind=dp) :: intN_pt_Qmat
  complex(kind=dp) :: tmp
  type(primitive_gaussian) :: pg
  complex(kind=dp) :: intN_00, intN_01, intN_10
  complex(kind=dp) :: intN_11, intN_02, intN_20
  
  call copy_DiracOutput_pg(pg)
  
  intN_pt_Qmat = (0._dp,0._dp)
  do i=1,4
     call calc_intN_pq(i,i,NBS_L,NBS_S,pg,c_p0(:,:),c_p0(:,:),intN_00) 
!     write(*,*) i,intN_00
     intN_pt_Qmat = intN_pt_Qmat +intN_00

     call calc_intN_pq(i,i,NBS_L,NBS_S,pg,c_p0(:,:),c_p1(:,:),intN_01) 
!     write(*,*) i,intN_01
     intN_pt_Qmat = intN_pt_Qmat +intN_01

     call calc_intN_pq(i,i,NBS_L,NBS_S,pg,c_p1(:,:),c_p0(:,:),intN_10) 
!     write(*,*) i,intN_10
     intN_pt_Qmat = intN_pt_Qmat +intN_10

     call calc_intN_pq(i,i,NBS_L,NBS_S,pg,c_p0(:,:),c_p2(:,:),intN_02) 
!     write(*,*) i,intN_02
     intN_pt_Qmat = intN_pt_Qmat +intN_02

     call calc_intN_pq(i,i,NBS_L,NBS_S,pg,c_p2(:,:),c_p0(:,:),intN_20) 
!     write(*,*) i,intN_20
     intN_pt_Qmat = intN_pt_Qmat +intN_20

     call calc_intN_pq(i,i,NBS_L,NBS_S,pg,c_p1(:,:),c_p1(:,:),intN_11) 
!     write(*,*) i,intN_11
     intN_pt_Qmat = intN_pt_Qmat +intN_11
  end do
 
  norm = 1.d0/dsqrt(dble(intN_pt_Qmat))

end subroutine renorm_pt2


!============================================================================
subroutine set_psi_pt0
!============================================================================

  use Precision
  use DiracOutput
  use DefineTypes
  use PTcoef

  implicit none

  integer :: i,j,k,l
  integer :: p,n,aa
  character(LEN=1) :: a,b
  integer :: actorb0, actorb
 
  complex(kind=dp) :: c_p(4,NMAX_PG)

  actorb0 = 1 ! All electron orbitals are activated.
  actorb = 2*NBS ! All electron orbitals are activated.
!  actorb = 4*NBS ! All orbitals are activated.

  eneQ0(:) = 0.d0

  do i=1,2*NBS! i runs specified orbitals.
    call index_from_Qmat(i,p,a)
    call pm12(a,aa)
    call copy_DiracOutput_cp(p,aa,c_p)
    if(p.gt.0) then
      eneQ0(i) = e_eig(p)
    elseif(p.lt.0) then ! Kramers pair
      eneQ0(i) = e_eig(-p)
    end if
    do j=1,NBS_L
      c_psi(1,j,i) = c_p(1,j)
      c_psi(2,j,i) = c_p(2,j)
    end do
    do j=1,NBS_S
      c_psi(3,j,i) = c_p(3,j)
      c_psi(4,j,i) = c_p(4,j)
    end do
  end do
  write(*,*)"# finish set psi pt0"

end subroutine set_psi_pt0

!============================================================================
subroutine set_psi_pt1_old
!============================================================================

  use Precision
  use DiracOutput
  use DefineTypes
  use PTcoef

  implicit none

  integer :: i,j,k,l
  integer :: p,n,aa
  character(LEN=1) :: a,b
  integer :: actorb0, actorb
 
  complex(kind=dp) :: c_p(4,NMAX_PG)
  complex(kind=dp) :: c_p0(4,NMAX_PG,2*NBS),c_psi0(4,NMAX_PG,2*NBS)
  complex(kind=dp),allocatable :: c_phi0(:,:,:),c_psi1(:,:,:)
  complex(kind=dp),allocatable :: tmp_intHrel_Qmat(:,:)
  complex(kind=dp) :: intHrel_pt_Qmat
  complex(kind=dp) :: tmp1, tmp2
  integer :: alp, alp2, degA

  complex(kind(0d0)), allocatable :: eigwork(:) ! this may be changed to complex array
  double precision, allocatable :: eigvals(:)
  integer, allocatable :: eigiwork(:)
  double precision,allocatable :: rwork(:)
  integer itmp,itemp,itemp2,itemp3,itemp4

  actorb0 = 1 ! All electron orbitals are activated.
  actorb = 2*NBS ! All electron orbitals are activated.
!  actorb = 4*NBS ! All orbitals are activated.

  c_p0(:,:,:)=(0.d0,0.d0)
  eneQ0(:) = 0.d0
  eneQ1(:) = 0.d0

  do i=1,2*NBS! i runs specified orbitals.
    call index_from_Qmat(i,p,a)
    call pm12(a,aa)
    call copy_DiracOutput_cp(p,aa,c_p)
    if(p.gt.0) then
      eneQ0(i) = e_eig(p)
    elseif(p.lt.0) then ! Kramers pair
      eneQ0(i) = e_eig(-p)
    end if
    do j=1,NBS_L
      c_p0(1,j,i) = c_p(1,j)
      c_p0(2,j,i) = c_p(2,j)
    end do
    do j=1,NBS_S
      c_p0(3,j,i) = c_p(3,j)
      c_p0(4,j,i) = c_p(4,j)
    end do
  end do
  write(*,*)"# finish set psi pt0"
!  do i=1,2*NBS! i runs specified orbitals.
!    write(*,*)i, eneQ0(i)
!  end do

  degA=1
  do i=actorb0,actorb
!      write(*,*)"i=",i,"degA=",degA
    if(degA.ne.1) then
!      write(*,*)"i=",i
      degA = degA -1
      cycle
    end if !degA
!!$ck  if(degA.eq.1) then
    !----decide degenerate number degA----
    do k=1,2*NBS
      if(k.ne.i) then
        if(abs(eneQ0(k)-eneQ0(i))<=1.d-8) then
          degA = degA + 1
        end if
      end if
    end do
!    write(*,*)" i, degA",i,degA
!    stop
    if(degA==1) then
      write(*,*)" Only degenerated calculation is supported."
    end if

    allocate(c_phi0(4,NMAX_PG,degA))
    allocate(c_psi1(4,NMAX_PG,degA))
    c_phi0(:,:,:) = (0.d0,0.d0)
    c_psi1(:,:,:) = (0.d0,0.d0)

    !----  set c_phi0 -------------------------
      do alp=1,degA
        do j=1,NBS_L
          c_phi0(1,j,alp) = c_p0(1,j,i+alp-1)
          c_phi0(2,j,alp) = c_p0(2,j,i+alp-1)
        end do
        do j=1,NBS_S
          c_phi0(3,j,alp) = c_p0(3,j,i+alp-1)
          c_phi0(4,j,alp) = c_p0(4,j,i+alp-1)
        end do
      end do
!      write(*,*)"# finish set c_phi0"
      ! --- calc intHrel_Qmat c_phi0 ----------------
      allocate(tmp_intHrel_Qmat(degA,degA))
      do alp=1,degA
        do alp2=1,degA
          tmp_intHrel_Qmat(alp,alp2) = intHrel_pt_Qmat(c_phi0(:,:,alp),c_phi0(:,:,alp2))
        end do
      end do
!      write(*,*)"# finish set intHrel c_phi0"
      ! check hermite
      do alp=1,degA
        do alp2=1,degA
          if (tmp_intHrel_Qmat(alp,alp2)-dconjg(tmp_intHrel_Qmat(alp,alp2)).ne.0) then
            write(14,*) 'diff',alp,alp2,tmp_intHrel_Qmat(alp,alp2)-dconjg(tmp_intHrel_Qmat(alp2,alp))
          end if
!            write(14,*) 'diff',alp,alp2,tmp_intHrel_Qmat(alp,alp2)-dconjg(tmp_intHrel_Qmat(alp2,alp))
!            write(15,*) alp,alp2,tmp_intHrel_Qmat(alp,alp2),dconjg(tmp_intHrel_Qmat(alp2,alp))
        end do
      end do
!      write(*,*)"# finish check Hermite"
!      stop
      !---- calc eigen problem -----------------

      itemp  = degA*(degA + 2) !lwork
      itemp2 = 5*degA + 3 !liwork
      itemp3 = 2*degA*degA + 5*degA + 1 !lrwork
      allocate(rwork(itemp3))
      allocate(eigwork(1 + 6 * degA + 2 * degA**2) )
      allocate(eigiwork(3 + 5 * degA**2) )
      allocate(eigvals(degA))
!       call dsyevd('V', 'U', degA, dnstmtx, degA, eigvals, &
!                   eigwork, itemp, eigiwork, itemp2, itemp3) 
      call zheevd('V', 'U', degA, tmp_intHrel_Qmat, degA, eigvals, &
                   eigwork, itemp, rwork, itemp3, eigiwork, itemp2, itemp4) 
      
!      write(*,*)"# finish calc eigen values and vectors"
      !--- set energy 1 -----------------
      do alp=1,degA
        eneQ1(i+alp-1) = eigvals(alp)
        write(100,*) alp,eigvals(alp)
      end do

      deallocate(rwork,eigwork,eigiwork,eigvals)

      !--- set c_psi0 -----------------------
        !--- Large component ---------------
        !$omp parallel do private(n,alp,alp2,tmp1,tmp2)
        do j=1,NBS_L
          do n=1,2*NBS
            if(n==i) then
              do alp2=1,degA
                tmp1 = (0.d0,0.d0)
                tmp2 = (0.d0,0.d0)
                do alp=1,degA
                  tmp1 = tmp1 + c_phi0(1,j,alp)*tmp_intHrel_Qmat(alp,alp2)
                  tmp2 = tmp2 + c_phi0(2,j,alp)*tmp_intHrel_Qmat(alp,alp2)
                end do
                c_psi0(1,j,n+alp2-1) = tmp1
                c_psi0(2,j,n+alp2-1) = tmp2
              end do
            else if((n>i).and.(n<i+degA)) then
              cycle
            else
              c_psi0(1,j,n) = c_p0(1,j,n)
              c_psi0(2,j,n) = c_p0(2,j,n)
!              write(*,*)"1",j,n,c_psi0(1,j,n)
!              write(*,*)"2",j,n,c_psi0(2,j,n)
            end if
          end do !n=1,2*NBS
        end do !j=1,NBS_L
        !$omp end parallel do

        !--- Small component ---------------
        !$omp parallel do private(n,alp,alp2,tmp1,tmp2)
        do j=1,NBS_S
          do n=1,2*NBS
            if(n==i) then
              do alp2=1,degA
                tmp1 = (0.d0,0.d0)
                tmp2 = (0.d0,0.d0)
                do alp=1,degA
                  tmp1 = tmp1 + c_phi0(3,j,alp)*tmp_intHrel_Qmat(alp,alp2)
                  tmp2 = tmp2 + c_phi0(4,j,alp)*tmp_intHrel_Qmat(alp,alp2)
                end do
                c_psi0(3,j,n+alp2-1) = tmp1
                c_psi0(4,j,n+alp2-1) = tmp2
              end do
            else if((n>i).and.(n<i+degA)) then
              cycle
            else
              c_psi0(3,j,n) = c_p0(3,j,n)
              c_psi0(4,j,n) = c_p0(4,j,n)
!              write(*,*)"3",j,n,c_psi0(3,j,n)
!              write(*,*)"4",j,n,c_psi0(4,j,n)
            end if
          end do !n=1,2*NBS
        end do !j=1,NBS_S
        !$omp end parallel do
!      write(*,*)"# finish set c_psi0"
!      stop

      deallocate(tmp_intHrel_Qmat)

      !-----  calc intHrel_Qmat c_psi0 --------------
      allocate(tmp_intHrel_Qmat(2*NBS,2*NBS))
      tmp_intHrel_Qmat(:,:)=(0.d0,0.d0)
      itmp=i
      !$omp parallel do
      do l=1,2*NBS
        do k=itmp,itmp+degA-1
          tmp_intHrel_Qmat(k,l) = intHrel_pt_Qmat(c_psi0(:,:,k),c_psi0(:,:,l))
!          write(*,*)k,l,tmp_intHrel_Qmat(k,l)
        end do
      end do
      !$omp end parallel do
!      write(*,*)"# finish calc intHrel for c_psi0"
!      stop
  
      !------------- Large component ---------------
      !$omp parallel do private(n,alp,alp2,tmp1,tmp2)
      do j=1,NBS_L
        do alp=1,degA
          do n=actorb0,actorb
            if((n<i).or.(n>i+degA-1)) then
              tmp1 = (0.d0,0.d0)
              tmp2 = (0.d0,0.d0)
              do alp2=1,degA
                if(alp.ne.alp2) then
                  tmp1 = tmp1 + c_psi0(1,j,i+alp2-1) *tmp_intHrel_Qmat(i+alp2-1,n) /(eneQ1(i+alp2-1)-eneQ1(i+alp-1))
                  tmp2 = tmp2 + c_psi0(2,j,i+alp2-1) *tmp_intHrel_Qmat(i+alp2-1,n) /(eneQ1(i+alp2-1)-eneQ1(i+alp-1))
                end if
              end do! alp2
              tmp1 = tmp1 - c_psi0(1,j,n)
              tmp2 = tmp2 - c_psi0(2,j,n)
              tmp1 = tmp1 * dconjg(tmp_intHrel_Qmat(i+alp-1,n)) /(eneQ0(n)-eneQ0(i))
              tmp2 = tmp2 * dconjg(tmp_intHrel_Qmat(i+alp-1,n)) /(eneQ0(n)-eneQ0(i))
          
              c_psi1(1,j,alp) = c_psi1(1,j,alp) + tmp1
              c_psi1(2,j,alp) = c_psi1(2,j,alp) + tmp2
            end if
          end do !n=actorb0,actorb
!!$          c_psi(1,j,i+alp-1) = c_psi0(1,j,i+alp-1) + lambda*c_psi1(1,j,alp)
!!$          c_psi(2,j,i+alp-1) = c_psi0(2,j,i+alp-1) + lambda*c_psi1(2,j,alp)
          c_psi(1,j,i+alp-1) = c_psi0(1,j,i+alp-1) + c_psi1(1,j,alp)
          c_psi(2,j,i+alp-1) = c_psi0(2,j,i+alp-1) + c_psi1(2,j,alp)
!          write(*,*)"1",j,i+alp-1,c_psi(1,j,i+alp-1)
!          write(*,*)"2",j,i+alp-1,c_psi(2,j,i+alp-1)
        end do !alp
      end do !j=NBS_L
      !$omp end parallel do

      !------------- Small component ---------------
      !$omp parallel do private(n,alp,alp2,tmp1,tmp2)
      do j=1,NBS_S
        do alp=1,degA
          do n=actorb0,actorb
            if((n<i).or.(n>i+degA-1)) then
              tmp1 = (0.d0,0.d0)
              tmp2 = (0.d0,0.d0)
              do alp2=1,degA
                if(alp.ne.alp2) then
                  tmp1 = tmp1 + c_psi0(3,j,i+alp2-1) *tmp_intHrel_Qmat(i+alp2-1,n) /(eneQ1(i+alp2-1)-eneQ1(i+alp-1))
                  tmp2 = tmp2 + c_psi0(4,j,i+alp2-1) *tmp_intHrel_Qmat(i+alp2-1,n) /(eneQ1(i+alp2-1)-eneQ1(i+alp-1))
                end if
              end do
              tmp1 = tmp1 - c_psi0(3,j,n)
              tmp2 = tmp2 - c_psi0(4,j,n)
              tmp1 = tmp1 * dconjg(tmp_intHrel_Qmat(i+alp-1,n)) /(eneQ0(n)-eneQ0(i))
              tmp2 = tmp2 * dconjg(tmp_intHrel_Qmat(i+alp-1,n)) /(eneQ0(n)-eneQ0(i))
          
              c_psi1(3,j,alp) = c_psi1(3,j,alp) + tmp1
              c_psi1(4,j,alp) = c_psi1(4,j,alp) + tmp2
            end if
          end do !n=actorb0,actorb
!!$          c_psi(3,j,i+alp-1) = c_psi0(3,j,i+alp-1) + lambda*c_psi1(3,j,alp)
!!$          c_psi(4,j,i+alp-1) = c_psi0(4,j,i+alp-1) + lambda*c_psi1(4,j,alp)
          c_psi(3,j,i+alp-1) = c_psi0(3,j,i+alp-1) + c_psi1(3,j,alp)
          c_psi(4,j,i+alp-1) = c_psi0(4,j,i+alp-1) + c_psi1(4,j,alp)
!          write(*,*)"3",j,i+alp-1,c_psi(3,j,i+alp-1)
!          write(*,*)"4",j,i+alp-1,c_psi(4,j,i+alp-1)
        end do !alp
      end do !j=NBS_S
      !$omp end parallel do
!      write(*,*)"# finish set c_psi"
!      stop
    write(*,*)"orb =",i,"/",actorb
      deallocate(tmp_intHrel_Qmat)
      deallocate(c_phi0)
      deallocate(c_psi1)
!!$ck  else! degA.ne.1
!!$ck    degA = degA -1
!!$ck  end if
  end do ! i=actorb0,actorb
!  write(*,*) '# end subroutine set psi pt1 old'

end subroutine set_psi_pt1_old

!============================================================================
subroutine set_psi_pt1
!============================================================================

  use Precision
  use DiracOutput
  use DefineTypes
  use PTcoef

  implicit none

  integer :: i,j,k,l
  integer :: p,n,aa
  character(LEN=1) :: a,b
  integer :: actorb0, actorb
 
  complex(kind=dp) :: c_p(4,NMAX_PG)
  complex(kind=dp) :: c_p0(4,NMAX_PG,2*NBS),c_psi0(4,NMAX_PG,2*NBS)
  complex(kind=dp),allocatable :: c_phi0(:,:,:),c_psi1(:,:,:),c_psi2(:,:,:)
  complex(kind=dp),allocatable :: tmp_intHrel_Qmat(:,:)
  complex(kind=dp),allocatable :: int_n0Hir0(:,:),int_n0ir1(:,:),int_ir0ib1(:,:)
  complex(kind=dp) :: intHrel_pt_Qmat
  complex(kind=dp) :: tmp1, tmp2
  complex(kind=dp) :: norm
  integer :: alp, alp2, degA

  complex(kind(0d0)), allocatable :: eigwork(:) ! this may be changed to complex array
  double precision, allocatable :: eigvals(:)
  integer, allocatable :: eigiwork(:)
  double precision,allocatable :: rwork(:)
  integer itmp,itemp,itemp2,itemp3,itemp4

  actorb0 = 1 ! All electron orbitals are activated.
  actorb = 2*NBS ! All electron orbitals are activated.
!  actorb = 4*NBS ! All orbitals are activated.

  c_p0(:,:,:)=(0.d0,0.d0)
  eneQ0(:) = 0.d0
  eneQ1(:) = 0.d0

  do i=1,2*NBS! i runs specified orbitals.
    call index_from_Qmat(i,p,a)
    call pm12(a,aa)
    call copy_DiracOutput_cp(p,aa,c_p)
    if(p.gt.0) then
      eneQ0(i) = e_eig(p)
    elseif(p.lt.0) then ! Kramers pair
      eneQ0(i) = e_eig(-p)
    end if
    do j=1,NBS_L
      c_p0(1,j,i) = c_p(1,j)
      c_p0(2,j,i) = c_p(2,j)
    end do
    do j=1,NBS_S
      c_p0(3,j,i) = c_p(3,j)
      c_p0(4,j,i) = c_p(4,j)
    end do
  end do
  write(*,*)"# finish set psi pt0"
!  do i=1,2*NBS! i runs specified orbitals.
!    write(*,*)i, eneQ0(i)
!  end do
!  stop

  allocate(c_psi1(4,NMAX_PG,2*NBS))
  c_psi1(:,:,:) = (0.d0,0.d0)

  degA=1
  do i=actorb0,actorb
!      write(*,*)"i=",i,"degA=",degA
    if(degA.ne.1) then
!      write(*,*)"i=",i
      degA = degA -1
      cycle
    end if !degA
!!$ck  if(degA.eq.1) then
    !----decide degenerate number degA----
    do k=1,2*NBS
      if(k.ne.i) then
        if(abs(eneQ0(k)-eneQ0(i))<=1.d-8) then
          degA = degA + 1
        end if
      end if
    end do
!    write(*,*)" i, degA",i,degA
!    stop
    if(degA==1) then
      write(*,*)" Only degenerated calculation is supported."
    end if

    allocate(c_phi0(4,NMAX_PG,degA))
    c_phi0(:,:,:) = (0.d0,0.d0)

    !----  set c_phi0 -------------------------
      do alp=1,degA
        do j=1,NBS_L
          c_phi0(1,j,alp) = c_p0(1,j,i+alp-1)
          c_phi0(2,j,alp) = c_p0(2,j,i+alp-1)
        end do
        do j=1,NBS_S
          c_phi0(3,j,alp) = c_p0(3,j,i+alp-1)
          c_phi0(4,j,alp) = c_p0(4,j,i+alp-1)
        end do
      end do
!      write(*,*)"# finish set c_phi0"
      ! --- calc intHrel_Qmat c_phi0 ----------------
      allocate(tmp_intHrel_Qmat(degA,degA))
      do alp=1,degA
        do alp2=1,degA
          tmp_intHrel_Qmat(alp,alp2) = intHrel_pt_Qmat(c_phi0(:,:,alp),c_phi0(:,:,alp2))
        end do
      end do
!      write(*,*)"# finish set intHrel c_phi0"
      ! check hermite
      do alp=1,degA
        do alp2=1,degA
          if (tmp_intHrel_Qmat(alp,alp2)-dconjg(tmp_intHrel_Qmat(alp,alp2)).ne.0) then
            write(14,*) 'diff',alp,alp2,tmp_intHrel_Qmat(alp,alp2)-dconjg(tmp_intHrel_Qmat(alp2,alp))
          end if
!            write(14,*) 'diff',alp,alp2,tmp_intHrel_Qmat(alp,alp2)-dconjg(tmp_intHrel_Qmat(alp2,alp))
!            write(15,*) alp,alp2,tmp_intHrel_Qmat(alp,alp2) !,dconjg(tmp_intHrel_Qmat(alp2,alp))
        end do
      end do
!      write(*,*)"# finish check Hermite"
!      stop
      !---- calc eigen problem -----------------

      itemp  = degA*(degA + 2) !lwork
      itemp2 = 5*degA + 3 !liwork
      itemp3 = 2*degA*degA + 5*degA + 1 !lrwork
      allocate(rwork(itemp3))
      allocate(eigwork(1 + 6 * degA + 2 * degA**2) )
      allocate(eigiwork(3 + 5 * degA**2) )
      allocate(eigvals(degA))
!       call dsyevd('V', 'U', degA, dnstmtx, degA, eigvals, &
!                   eigwork, itemp, eigiwork, itemp2, itemp3) 
      call zheevd('V', 'U', degA, tmp_intHrel_Qmat, degA, eigvals, &
                   eigwork, itemp, rwork, itemp3, eigiwork, itemp2, itemp4) 
      
!      write(*,*)"# finish calc eigen values and vectors"
      !--- set energy 1 -----------------
      do alp=1,degA
        eneQ1(i+alp-1) = eigvals(alp)
        write(100,*) alp,eigvals(alp)
      end do

      deallocate(rwork,eigwork,eigiwork,eigvals)

      !--- set c_psi0 -----------------------
        !--- Large component ---------------
        !$omp parallel do private(n,alp,alp2,tmp1,tmp2)
        do j=1,NBS_L
          do n=actorb0, actorb
            if(n==i) then
              do alp2=1,degA
                tmp1 = (0.d0,0.d0)
                tmp2 = (0.d0,0.d0)
                do alp=1,degA
                  tmp1 = tmp1 + c_phi0(1,j,alp)*tmp_intHrel_Qmat(alp,alp2)
                  tmp2 = tmp2 + c_phi0(2,j,alp)*tmp_intHrel_Qmat(alp,alp2)
                end do
                c_psi0(1,j,n+alp2-1) = tmp1
                c_psi0(2,j,n+alp2-1) = tmp2
              end do
            else if((n>i).and.(n<i+degA)) then
              cycle
            else
              c_psi0(1,j,n) = c_p0(1,j,n)
              c_psi0(2,j,n) = c_p0(2,j,n)
!              write(*,*)"1",j,n,c_psi0(1,j,n)
!              write(*,*)"2",j,n,c_psi0(2,j,n)
            end if
          end do !n=1,2*NBS
        end do !j=1,NBS_L
        !$omp end parallel do

        !--- Small component ---------------
        !$omp parallel do private(n,alp,alp2,tmp1,tmp2)
        do j=1,NBS_S
          do n=actorb0,actorb
            if(n==i) then
              do alp2=1,degA
                tmp1 = (0.d0,0.d0)
                tmp2 = (0.d0,0.d0)
                do alp=1,degA
                  tmp1 = tmp1 + c_phi0(3,j,alp)*tmp_intHrel_Qmat(alp,alp2)
                  tmp2 = tmp2 + c_phi0(4,j,alp)*tmp_intHrel_Qmat(alp,alp2)
                end do
                c_psi0(3,j,n+alp2-1) = tmp1
                c_psi0(4,j,n+alp2-1) = tmp2
              end do
            else if((n>i).and.(n<i+degA)) then
              cycle
            else
              c_psi0(3,j,n) = c_p0(3,j,n)
              c_psi0(4,j,n) = c_p0(4,j,n)
!              write(*,*)"3",j,n,c_psi0(3,j,n)
!              write(*,*)"4",j,n,c_psi0(4,j,n)
            end if
          end do !n=1,2*NBS
        end do !j=1,NBS_S
        !$omp end parallel do
!      write(*,*)"# finish set c_psi0"
!      stop

      deallocate(tmp_intHrel_Qmat)

      !-----  calc int_n0Hir0 = <psi^0_n|H'|psi^0_ir> or <psi^0_n|H'|psi^0_m> --------------
      allocate(int_n0Hir0(2*NBS,2*NBS))
      int_n0Hir0(:,:)=(0.d0,0.d0)
      !$omp parallel do private(n,k)
      do n=actorb0,actorb
        do k=actorb0,actorb
          int_n0Hir0(n,k) = intHrel_pt_Qmat(c_psi0(:,:,n),c_psi0(:,:,k))
!          write(*,*)n,k,int_n0Hir0(n,k)
        end do
      end do
      !$omp end parallel do
!      write(*,*)"# finish calc int_n0Hir0 = <psi^0_n|H'|psi^0_ir>"
!      stop

      !-----  calc int_n0ir1 = <psi^0_n|psi^1_ir> --------------
      allocate(int_n0ir1(2*NBS,2*NBS))
      int_n0ir1(:,:)=(0.d0,0.d0)
!      !$omp parallel do
      open(unit=51,file="n0ib1.dat",position='append')
      do n=actorb0,actorb
        do k=i,i+degA-1
          int_n0ir1(n,k) = - int_n0Hir0(n,k)/(eneQ0(n)-eneQ0(k))
          write(51,"(2i6,2es16.6)")n,k,int_n0ir1(n,k)
        end do
      end do
      close(unit=51)
!      !$omp end parallel do

      !-----  calc int_ir0ib1 = <psi^0_ir|psi^1_ib> --------------
      allocate(int_ir0ib1(2*NBS,2*NBS))
      int_ir0ib1(:,:)=(0.d0,0.d0)
!      !$omp parallel do private(tmp1)
      open(unit=51,file="ir0ib1.dat",position='append')
      do l=i,i+degA-1
        do k=i,i+degA-1
          if(l.ne.k) then
            tmp1=(0.d0,0.d0)
            do n=actorb0,actorb
              if((n<i).or.(n>i+degA-1)) then
                tmp1 = tmp1 - dconjg(int_n0Hir0(n,l))*int_n0ir1(n,k)
              end if
            end do
            int_ir0ib1(l,k) = tmp1/(eneQ1(l)-eneQ1(k))
            write(51,"(2i6,2es16.6)")l,k,int_ir0ib1(l,k)
          end if
        end do
      end do
      close(unit=51)
!      !$omp end parallel do

      !---------------------------
      ! calc psi_pt1
      !---------------------------
      !------------- Large component ---------------
      !$omp parallel do private(n,alp,alp2,tmp1,tmp2)
      do j=1,NBS_L
        do alp=1,degA
          tmp1 = (0.d0,0.d0)
          tmp2 = (0.d0,0.d0)
          do alp2=1,degA
            if(alp.ne.alp2) then
              tmp1 = tmp1 + c_psi0(1,j,i+alp2-1) *int_ir0ib1(i+alp2-1,i+alp-1)
              tmp2 = tmp2 + c_psi0(2,j,i+alp2-1) *int_ir0ib1(i+alp2-1,i+alp-1)
            end if
          end do! alp2
          do n=actorb0,actorb
            if((n<i).or.(n>i+degA-1)) then
              tmp1 = tmp1 + c_psi0(1,j,n) *int_n0ir1(n,i+alp-1)
              tmp2 = tmp2 + c_psi0(2,j,n) *int_n0ir1(n,i+alp-1)
            end if
          end do
          c_psi1(1,j,i+alp-1) = tmp1
          c_psi1(2,j,i+alp-1) = tmp2
!!$          c_psi(1,j,i+alp-1) = c_psi0(1,j,i+alp-1) + lambda*c_psi1(1,j,i+alp-1)
!!$          c_psi(2,j,i+alp-1) = c_psi0(2,j,i+alp-1) + lambda*c_psi1(2,j,i+alp-1)
!!$          c_psi(1,j,i+alp-1) = c_psi0(1,j,i+alp-1) + c_psi1(1,j,i+alp-1)
!!$          c_psi(2,j,i+alp-1) = c_psi0(2,j,i+alp-1) + c_psi1(2,j,i+alp-1)
!!$          write(*,*)"1",j,i+alp-1,c_psi(1,j,i+alp-1)
!!$          write(*,*)"2",j,i+alp-1,c_psi(2,j,i+alp-1)
        end do !alp
      end do !j=NBS_L
      !$omp end parallel do

      !------------- Small component ---------------
      !$omp parallel do private(n,alp,alp2,tmp1,tmp2)
      do j=1,NBS_S
        do alp=1,degA
          tmp1 = (0.d0,0.d0)
          tmp2 = (0.d0,0.d0)
          do alp2=1,degA
            if(alp.ne.alp2) then
              tmp1 = tmp1 + c_psi0(3,j,i+alp2-1) *int_ir0ib1(i+alp2-1,i+alp-1)
              tmp2 = tmp2 + c_psi0(4,j,i+alp2-1) *int_ir0ib1(i+alp2-1,i+alp-1)
            end if
          end do! alp2
          do n=actorb0,actorb
            if((n<i).or.(n>i+degA-1)) then
              tmp1 = tmp1 + c_psi0(3,j,n) *int_n0ir1(n,i+alp-1)
              tmp2 = tmp2 + c_psi0(4,j,n) *int_n0ir1(n,i+alp-1)
            end if
          end do
          c_psi1(3,j,i+alp-1) = tmp1
          c_psi1(4,j,i+alp-1) = tmp2
!!$          c_psi(1,j,i+alp-1) = c_psi0(1,j,i+alp-1) + lambda*c_psi1(1,j,i+alp-1)
!!$          c_psi(2,j,i+alp-1) = c_psi0(2,j,i+alp-1) + lambda*c_psi1(2,j,i+alp-1)
!!$          c_psi(1,j,i+alp-1) = c_psi0(1,j,i+alp-1) + c_psi1(1,j,i+alp-1)
!!$          c_psi(2,j,i+alp-1) = c_psi0(2,j,i+alp-1) + c_psi1(2,j,i+alp-1)
!!$          write(*,*)"1",j,i+alp-1,c_psi(1,j,i+alp-1)
!!$          write(*,*)"2",j,i+alp-1,c_psi(2,j,i+alp-1)
        end do !alp
      end do !j=NBS_S
      !$omp end parallel do
 
      do alp=1,degA
        call renorm_pt1(c_psi0(:,:,i+alp-1),c_psi1(:,:,i+alp-1),norm)
        do j=1,NBS_L
          c_psi(1,j,i+alp-1) = c_psi(1,j,i+alp-1)*norm
          c_psi(2,j,i+alp-1) = c_psi(2,j,i+alp-1)*norm
        end do 
        do j=1,NBS_S
          c_psi(3,j,i+alp-1) = c_psi(3,j,i+alp-1)*norm
          c_psi(4,j,i+alp-1) = c_psi(4,j,i+alp-1)*norm
        end do 
      end do !alp

      !-------- finish calc psi pt1 -------------
    write(*,*)"orb =",i,"/",actorb
      deallocate(int_n0Hir0,int_n0ir1,int_ir0ib1)
      deallocate(c_phi0)
!!$ck  else! degA.ne.1
!!$ck    degA = degA -1
!!$ck  end if
  end do ! i=actorb0,actorb
  deallocate(c_psi1)

end subroutine set_psi_pt1

!============================================================================
subroutine set_psi_pt2
!============================================================================

  use Precision
  use DiracOutput
  use DefineTypes
  use PTcoef

  implicit none

  integer :: i,j,k,l
  integer :: p,n,aa
  character(LEN=1) :: a,b
  integer :: actorb0, actorb
 
  complex(kind=dp) :: c_p(4,NMAX_PG)
  complex(kind=dp) :: c_p0(4,NMAX_PG,2*NBS),c_psi0(4,NMAX_PG,2*NBS)
  complex(kind=dp),allocatable :: c_phi0(:,:,:),c_psi1(:,:,:),c_psi2(:,:,:)
  complex(kind=dp),allocatable :: tmp_intHrel_Qmat(:,:)
  complex(kind=dp),allocatable :: int_n0Hir0(:,:),int_n0ir1(:,:),int_ir0ib1(:,:)
  complex(kind=dp),allocatable :: int_n0Hir1(:,:),int_n0ir2(:,:),int_ir0ib2(:,:)
  complex(kind=dp) :: intHrel_pt_Qmat
  complex(kind=dp) :: tmp1, tmp2
  complex(kind=dp) :: norm
  integer :: alp, alp2, degA

  complex(kind(0d0)), allocatable :: eigwork(:) ! this may be changed to complex array
  double precision, allocatable :: eigvals(:)
  integer, allocatable :: eigiwork(:)
  double precision,allocatable :: rwork(:)
  integer itmp,itemp,itemp2,itemp3,itemp4

  actorb0 = 1 ! All electron orbitals are activated.
  actorb = 2*NBS ! All electron orbitals are activated.
!  actorb = 4*NBS ! All orbitals are activated.

  c_p0(:,:,:)=(0.d0,0.d0)
  eneQ0(:) = 0.d0
  eneQ1(:) = 0.d0
  eneQ2(:) = 0.d0

  do i=1,2*NBS! i runs specified orbitals.
    call index_from_Qmat(i,p,a)
    call pm12(a,aa)
    call copy_DiracOutput_cp(p,aa,c_p)
    if(p.gt.0) then
      eneQ0(i) = e_eig(p)
    elseif(p.lt.0) then ! Kramers pair
      eneQ0(i) = e_eig(-p)
    end if
    do j=1,NBS_L
      c_p0(1,j,i) = c_p(1,j)
      c_p0(2,j,i) = c_p(2,j)
    end do
    do j=1,NBS_S
      c_p0(3,j,i) = c_p(3,j)
      c_p0(4,j,i) = c_p(4,j)
    end do
  end do
  write(*,*)"# finish set psi pt0"
!  do i=1,2*NBS! i runs specified orbitals.
!    write(*,*)i, eneQ0(i)
!  end do
!  stop

  allocate(c_psi1(4,NMAX_PG,2*NBS))
  allocate(c_psi2(4,NMAX_PG,2*NBS))
  c_psi1(:,:,:) = (0.d0,0.d0)
  c_psi2(:,:,:) = (0.d0,0.d0)

  degA=1
  do i=actorb0,actorb
!      write(*,*)"i=",i,"degA=",degA
    if(degA.ne.1) then
!      write(*,*)"i=",i
      degA = degA -1
      cycle
    end if !degA
!!$ck  if(degA.eq.1) then
    !----decide degenerate number degA----
    do k=1,2*NBS
      if(k.ne.i) then
        if(abs(eneQ0(k)-eneQ0(i))<=1.d-8) then
          degA = degA + 1
        end if
      end if
    end do
!    write(*,*)" i, degA",i,degA
!    stop
    if(degA==1) then
      write(*,*)" Only degenerated calculation is supported."
    end if

    allocate(c_phi0(4,NMAX_PG,degA))
!    allocate(c_psi1(4,NMAX_PG,degA))
    c_phi0(:,:,:) = (0.d0,0.d0)

    !----  set c_phi0 -------------------------
      do alp=1,degA
        do j=1,NBS_L
          c_phi0(1,j,alp) = c_p0(1,j,i+alp-1)
          c_phi0(2,j,alp) = c_p0(2,j,i+alp-1)
        end do
        do j=1,NBS_S
          c_phi0(3,j,alp) = c_p0(3,j,i+alp-1)
          c_phi0(4,j,alp) = c_p0(4,j,i+alp-1)
        end do
      end do
!      write(*,*)"# finish set c_phi0"
      ! --- calc intHrel_Qmat c_phi0 ----------------
      allocate(tmp_intHrel_Qmat(degA,degA))
      do alp=1,degA
        do alp2=1,degA
          tmp_intHrel_Qmat(alp,alp2) = intHrel_pt_Qmat(c_phi0(:,:,alp),c_phi0(:,:,alp2))
        end do
      end do
!      write(*,*)"# finish set intHrel c_phi0"
      ! check hermite
      do alp=1,degA
        do alp2=1,degA
          if (tmp_intHrel_Qmat(alp,alp2)-dconjg(tmp_intHrel_Qmat(alp,alp2)).ne.0) then
            write(14,*) 'diff',alp,alp2,tmp_intHrel_Qmat(alp,alp2)-dconjg(tmp_intHrel_Qmat(alp2,alp))
          end if
!            write(14,*) 'diff',alp,alp2,tmp_intHrel_Qmat(alp,alp2)-dconjg(tmp_intHrel_Qmat(alp2,alp))
!            write(15,*) alp,alp2,tmp_intHrel_Qmat(alp,alp2) !,dconjg(tmp_intHrel_Qmat(alp2,alp))
        end do
      end do
!      write(*,*)"# finish check Hermite"
!      stop
      !---- calc eigen problem -----------------

      itemp  = degA*(degA + 2) !lwork
      itemp2 = 5*degA + 3 !liwork
      itemp3 = 2*degA*degA + 5*degA + 1 !lrwork
      allocate(rwork(itemp3))
      allocate(eigwork(1 + 6 * degA + 2 * degA**2) )
      allocate(eigiwork(3 + 5 * degA**2) )
      allocate(eigvals(degA))
!       call dsyevd('V', 'U', degA, dnstmtx, degA, eigvals, &
!                   eigwork, itemp, eigiwork, itemp2, itemp3) 
      call zheevd('V', 'U', degA, tmp_intHrel_Qmat, degA, eigvals, &
                   eigwork, itemp, rwork, itemp3, eigiwork, itemp2, itemp4) 
      
!      write(*,*)"# finish calc eigen values and vectors"
      !--- set energy 1 -----------------
      do alp=1,degA
        eneQ1(i+alp-1) = eigvals(alp)
        write(100,*) alp,eigvals(alp)
      end do

      deallocate(rwork,eigwork,eigiwork,eigvals)

      !--- set c_psi0 -----------------------
        !--- Large component ---------------
        !$omp parallel do private(n,alp,alp2,tmp1,tmp2)
        do j=1,NBS_L
          do n=actorb0, actorb
            if(n==i) then
              do alp2=1,degA
                tmp1 = (0.d0,0.d0)
                tmp2 = (0.d0,0.d0)
                do alp=1,degA
                  tmp1 = tmp1 + c_phi0(1,j,alp)*tmp_intHrel_Qmat(alp,alp2)
                  tmp2 = tmp2 + c_phi0(2,j,alp)*tmp_intHrel_Qmat(alp,alp2)
                end do
                c_psi0(1,j,n+alp2-1) = tmp1
                c_psi0(2,j,n+alp2-1) = tmp2
              end do
            else if((n>i).and.(n<i+degA)) then
              cycle
            else
              c_psi0(1,j,n) = c_p0(1,j,n)
              c_psi0(2,j,n) = c_p0(2,j,n)
!              write(*,*)"1",j,n,c_psi0(1,j,n)
!              write(*,*)"2",j,n,c_psi0(2,j,n)
            end if
          end do !n=1,2*NBS
        end do !j=1,NBS_L
        !$omp end parallel do

        !--- Small component ---------------
        !$omp parallel do private(n,alp,alp2,tmp1,tmp2)
        do j=1,NBS_S
          do n=actorb0,actorb
            if(n==i) then
              do alp2=1,degA
                tmp1 = (0.d0,0.d0)
                tmp2 = (0.d0,0.d0)
                do alp=1,degA
                  tmp1 = tmp1 + c_phi0(3,j,alp)*tmp_intHrel_Qmat(alp,alp2)
                  tmp2 = tmp2 + c_phi0(4,j,alp)*tmp_intHrel_Qmat(alp,alp2)
                end do
                c_psi0(3,j,n+alp2-1) = tmp1
                c_psi0(4,j,n+alp2-1) = tmp2
              end do
            else if((n>i).and.(n<i+degA)) then
              cycle
            else
              c_psi0(3,j,n) = c_p0(3,j,n)
              c_psi0(4,j,n) = c_p0(4,j,n)
!              write(*,*)"3",j,n,c_psi0(3,j,n)
!              write(*,*)"4",j,n,c_psi0(4,j,n)
            end if
          end do !n=1,2*NBS
        end do !j=1,NBS_S
        !$omp end parallel do
!      write(*,*)"# finish set c_psi0"
!      stop

      deallocate(tmp_intHrel_Qmat)

      !-----  calc int_n0Hir0 = <psi^0_n|H'|psi^0_ir> or <psi^0_n|H'|psi^0_m> --------------
      allocate(int_n0Hir0(2*NBS,2*NBS))
      int_n0Hir0(:,:)=(0.d0,0.d0)
      !$omp parallel do private(n,k)
      do n=actorb0,actorb
        do k=actorb0,actorb
          int_n0Hir0(n,k) = intHrel_pt_Qmat(c_psi0(:,:,n),c_psi0(:,:,k))
!          write(*,*)n,k,int_n0Hir0(n,k)
        end do
      end do
      !$omp end parallel do
!      write(*,*)"# finish calc int_n0Hir0 = <psi^0_n|H'|psi^0_ir>"
!      stop

      do alp=1,degA
        tmp1 =(0.d0,0.d0)
        do n=actorb0,actorb
          if((n<i).or.(n>i+degA-1)) then
            tmp1 = tmp1 - dconjg(int_n0Hir0(n,i+alp-1))*int_n0Hir0(n,i+alp-1)/(eneQ0(n)-eneQ0(i+alp-1))
          end if
        end do
        eneQ2(i+alp-1) = dble(tmp1)
        write(*,*)'eneQ2,tmp1',eneQ2(i+alp-1),tmp1
      end do
!      stop

      !-----  calc int_n0ir1 = <psi^0_n|psi^1_ir> --------------
      allocate(int_n0ir1(2*NBS,2*NBS))
      int_n0ir1(:,:)=(0.d0,0.d0)
!      !$omp parallel do
      open(unit=51,file="n0ib1.dat",position='append')
      do n=actorb0,actorb
        do k=i,i+degA-1
          int_n0ir1(n,k) = - int_n0Hir0(n,k)/(eneQ0(n)-eneQ0(k))
          write(51,"(2i6,2es16.6)")n,k,int_n0ir1(n,k)
        end do
      end do
      close(unit=51)
!      !$omp end parallel do

      !-----  calc int_ir0ib1 = <psi^0_ir|psi^1_ib> --------------
      allocate(int_ir0ib1(2*NBS,2*NBS))
      int_ir0ib1(:,:)=(0.d0,0.d0)
!      !$omp parallel do private(tmp1)
      open(unit=51,file="ir0ib1.dat",position='append')
      do l=i,i+degA-1
        do k=i,i+degA-1
          if(l.ne.k) then
            tmp1=(0.d0,0.d0)
            do n=actorb0,actorb
              if((n<i).or.(n>i+degA-1)) then
                tmp1 = tmp1 - dconjg(int_n0Hir0(n,l))*int_n0ir1(n,k)
              end if
            end do
            int_ir0ib1(l,k) = tmp1/(eneQ1(l)-eneQ1(k))
            write(51,"(2i6,2es16.6)")l,k,int_ir0ib1(l,k)
          end if
        end do
      end do
      close(unit=51)
!      !$omp end parallel do

      !---------------------------
      ! calc psi_pt1
      !---------------------------
      !------------- Large component ---------------
      !$omp parallel do private(n,alp,alp2,tmp1,tmp2)
      do j=1,NBS_L
        do alp=1,degA
          tmp1 = (0.d0,0.d0)
          tmp2 = (0.d0,0.d0)
          do alp2=1,degA
            if(alp.ne.alp2) then
              tmp1 = tmp1 + c_psi0(1,j,i+alp2-1) *int_ir0ib1(i+alp2-1,i+alp-1)
              tmp2 = tmp2 + c_psi0(2,j,i+alp2-1) *int_ir0ib1(i+alp2-1,i+alp-1)
            end if
          end do! alp2
          do n=actorb0,actorb
            if((n<i).or.(n>i+degA-1)) then
              tmp1 = tmp1 + c_psi0(1,j,n) *int_n0ir1(n,i+alp-1)
              tmp2 = tmp2 + c_psi0(2,j,n) *int_n0ir1(n,i+alp-1)
            end if
          end do
          c_psi1(1,j,i+alp-1) = tmp1
          c_psi1(2,j,i+alp-1) = tmp2
!!$          c_psi(1,j,i+alp-1) = c_psi0(1,j,i+alp-1) + lambda*c_psi1(1,j,i+alp-1)
!!$          c_psi(2,j,i+alp-1) = c_psi0(2,j,i+alp-1) + lambda*c_psi1(2,j,i+alp-1)
!!$          c_psi(1,j,i+alp-1) = c_psi0(1,j,i+alp-1) + c_psi1(1,j,i+alp-1)
!!$          c_psi(2,j,i+alp-1) = c_psi0(2,j,i+alp-1) + c_psi1(2,j,i+alp-1)
!!$          write(*,*)"1",j,i+alp-1,c_psi(1,j,i+alp-1)
!!$          write(*,*)"2",j,i+alp-1,c_psi(2,j,i+alp-1)
        end do !alp
      end do !j=NBS_L
      !$omp end parallel do

      !------------- Small component ---------------
      !$omp parallel do private(n,alp,alp2,tmp1,tmp2)
      do j=1,NBS_S
        do alp=1,degA
          tmp1 = (0.d0,0.d0)
          tmp2 = (0.d0,0.d0)
          do alp2=1,degA
            if(alp.ne.alp2) then
              tmp1 = tmp1 + c_psi0(3,j,i+alp2-1) *int_ir0ib1(i+alp2-1,i+alp-1)
              tmp2 = tmp2 + c_psi0(4,j,i+alp2-1) *int_ir0ib1(i+alp2-1,i+alp-1)
            end if
          end do! alp2
          do n=actorb0,actorb
            if((n<i).or.(n>i+degA-1)) then
              tmp1 = tmp1 + c_psi0(3,j,n) *int_n0ir1(n,i+alp-1)
              tmp2 = tmp2 + c_psi0(4,j,n) *int_n0ir1(n,i+alp-1)
            end if
          end do
          c_psi1(3,j,i+alp-1) = tmp1
          c_psi1(4,j,i+alp-1) = tmp2
!!$          c_psi(1,j,i+alp-1) = c_psi0(1,j,i+alp-1) + lambda*c_psi1(1,j,i+alp-1)
!!$          c_psi(2,j,i+alp-1) = c_psi0(2,j,i+alp-1) + lambda*c_psi1(2,j,i+alp-1)
!!$          c_psi(1,j,i+alp-1) = c_psi0(1,j,i+alp-1) + c_psi1(1,j,i+alp-1)
!!$          c_psi(2,j,i+alp-1) = c_psi0(2,j,i+alp-1) + c_psi1(2,j,i+alp-1)
!!$          write(*,*)"1",j,i+alp-1,c_psi(1,j,i+alp-1)
!!$          write(*,*)"2",j,i+alp-1,c_psi(2,j,i+alp-1)
        end do !alp
      end do !j=NBS_S
      !$omp end parallel do

      !-------- finish calc psi pt1 -------------
      
      !-----  calc int_n0Hir1 = <psi^0_n|H'|psi^1_ir> --------------
      allocate(int_n0Hir1(2*NBS,2*NBS))
      int_n0Hir1(:,:)=(0.d0,0.d0)
      !$omp parallel do private(tmp1)
      do n=actorb0,actorb
        do k=i,i+degA-1
          tmp1=(0.d0,0.d0)
          do l=i,i+degA-1
            tmp1 = tmp1 + int_n0Hir0(n,l)*int_ir0ib1(l,k)
          end do
          do l=actorb0,actorb
            if((l<i).or.(l>i+degA-1)) then
              tmp1 = tmp1 + int_n0Hir0(n,l)*int_n0ir1(l,k)
            end if
          end do

          int_n0Hir1(n,k) = tmp1
!          write(*,*)n,k,int_n0Hir1(n,k)
        end do
      end do
      !$omp end parallel do
!      write(*,*)"# finish calc int_n0Hir0 = <psi^0_n|H'|psi^1_ir>"
!      stop

      !-----  calc int_n0ir2 = <psi^0_n|psi^2_ir> --------------
      allocate(int_n0ir2(2*NBS,2*NBS))
      int_n0ir2(:,:)=(0.d0,0.d0)
      !$omp parallel do
      do n=actorb0,actorb
        do k=i,i+degA-1
          int_n0ir2(n,k) = (- int_n0Hir1(n,k) + eneQ1(k)*int_n0ir1(n,k))/(eneQ0(n)-eneQ0(k))
!          write(*,*)n,k,int_n0ir2(n,k)
        end do
      end do
      !$omp end parallel do

      !-----  calc int_ir0ib2 = <psi^0_ir|psi^2_ib> --------------
      allocate(int_ir0ib2(2*NBS,2*NBS))
      int_ir0ib2(:,:)=(0.d0,0.d0)
      !$omp parallel do private(tmp1)
      do l=i,i+degA-1
        do k=i,i+degA-1
          if(l.ne.k) then
            tmp1=(0.d0,0.d0)
            do n=actorb0,actorb
              if((n<i).or.(n>i+degA-1)) then
                tmp1 = tmp1 - dconjg(int_n0Hir0(n,l))*int_n0ir2(n,k)
              end if
            end do
            tmp1 = tmp1 + eneQ2(k)*int_ir0ib1(l,k)
            int_ir0ib2(l,k) = tmp1/(eneQ1(l)-eneQ1(k))
!            write(*,*)l,k,int_ir0ib1(l,k)
          end if
        end do
      end do
      !$omp end parallel do


      !---------------------------
      ! calc psi_pt2
      !---------------------------

      !------------- Large component ---------------
      !$omp parallel do private(n,alp,alp2,tmp1,tmp2)
      do j=1,NBS_L
        do alp=1,degA
          tmp1 = (0.d0,0.d0)
          tmp2 = (0.d0,0.d0)
          do alp2=1,degA
            if(alp.ne.alp2) then
              tmp1 = tmp1 + c_psi0(1,j,i+alp2-1) *int_ir0ib2(i+alp2-1,i+alp-1)
              tmp2 = tmp2 + c_psi0(2,j,i+alp2-1) *int_ir0ib2(i+alp2-1,i+alp-1)
            end if
          end do! alp2
          do n=actorb0,actorb
            if((n<i).or.(n>i+degA-1)) then
              tmp1 = tmp1 + c_psi0(1,j,n) *int_n0ir2(n,i+alp-1)
              tmp2 = tmp2 + c_psi0(2,j,n) *int_n0ir2(n,i+alp-1)
            end if
          end do
          c_psi2(1,j,i+alp-1) = tmp1
          c_psi2(2,j,i+alp-1) = tmp2
!!$          c_psi(1,j,i+alp-1) = c_psi0(1,j,i+alp-1) + lambda*c_psi1(1,j,i+alp-1) + lambda*lambda*c_psi2(1,j,i+alp-1)
!!$          c_psi(2,j,i+alp-1) = c_psi0(2,j,i+alp-1) + lambda*c_psi1(2,j,i+alp-1) + lambda*lambda*c_psi2(2,j,i+alp-1)
          c_psi(1,j,i+alp-1) = c_psi0(1,j,i+alp-1) + c_psi1(1,j,i+alp-1) + c_psi2(1,j,i+alp-1)
          c_psi(2,j,i+alp-1) = c_psi0(2,j,i+alp-1) + c_psi1(2,j,i+alp-1) + c_psi2(2,j,i+alp-1)
!          c_psi(1,j,i+alp-1) = c_psi1(1,j,i+alp-1) 
!          c_psi(2,j,i+alp-1) = c_psi1(2,j,i+alp-1) 
!!$          write(*,*)"1",j,i+alp-1,c_psi(1,j,i+alp-1)
!!$          write(*,*)"2",j,i+alp-1,c_psi(2,j,i+alp-1)
        end do !alp
      end do !j=NBS_L
      !$omp end parallel do

      !------------- Small component ---------------
      !$omp parallel do private(n,alp,alp2,tmp1,tmp2)
      do j=1,NBS_S
        do alp=1,degA
          tmp1 = (0.d0,0.d0)
          tmp2 = (0.d0,0.d0)
          do alp2=1,degA
            if(alp.ne.alp2) then
              tmp1 = tmp1 + c_psi0(3,j,i+alp2-1) *int_ir0ib2(i+alp2-1,i+alp-1)
              tmp2 = tmp2 + c_psi0(4,j,i+alp2-1) *int_ir0ib2(i+alp2-1,i+alp-1)
            end if
          end do! alp2
          do n=actorb0,actorb
            if((n<i).or.(n>i+degA-1)) then
              tmp1 = tmp1 + c_psi0(3,j,n) *int_n0ir2(n,i+alp-1)
              tmp2 = tmp2 + c_psi0(4,j,n) *int_n0ir2(n,i+alp-1)
            end if
          end do
          c_psi2(3,j,i+alp-1) = tmp1
          c_psi2(4,j,i+alp-1) = tmp2
!!$          c_psi(3,j,i+alp-1) = c_psi0(3,j,i+alp-1) + lambda*c_psi1(3,j,i+alp-1) + lambda*lambda*c_psi2(3,j,i+alp-1)
!!$          c_psi(4,j,i+alp-1) = c_psi0(4,j,i+alp-1) + lambda*c_psi1(4,j,i+alp-1) + lambda*lambda*c_psi2(4,j,i+alp-1)
          c_psi(3,j,i+alp-1) = c_psi0(3,j,i+alp-1) + c_psi1(3,j,i+alp-1) + c_psi2(3,j,i+alp-1)
          c_psi(4,j,i+alp-1) = c_psi0(4,j,i+alp-1) + c_psi1(4,j,i+alp-1) + c_psi2(4,j,i+alp-1)
!          c_psi(3,j,i+alp-1) = c_psi1(3,j,i+alp-1) 
!          c_psi(4,j,i+alp-1) = c_psi1(4,j,i+alp-1) 
!!$          write(*,*)"3",j,i+alp-1,c_psi(3,j,i+alp-1)
!!$          write(*,*)"4",j,i+alp-1,c_psi(4,j,i+alp-1)
        end do !alp
      end do !j=NBS_S
      !$omp end parallel do

      !-------- finish calc psi pt2 -------------
 
      do alp=1,degA
        call renorm_pt2(c_psi0(:,:,i+alp-1),c_psi1(:,:,i+alp-1),c_psi2(:,:,i+alp-1),norm)
        do j=1,NBS_L
          c_psi(1,j,i+alp-1) = c_psi(1,j,i+alp-1)*norm
          c_psi(2,j,i+alp-1) = c_psi(2,j,i+alp-1)*norm
        end do 
        do j=1,NBS_S
          c_psi(3,j,i+alp-1) = c_psi(3,j,i+alp-1)*norm
          c_psi(4,j,i+alp-1) = c_psi(4,j,i+alp-1)*norm
        end do 
      end do !alp

!      write(*,*)"# finish set c_psi"
!      stop
    write(*,*)"orb =",i,"/",actorb
      deallocate(int_n0Hir0,int_n0ir1,int_ir0ib1)
      deallocate(int_n0Hir1,int_n0ir2,int_ir0ib2)
      deallocate(c_phi0)
!!$ck  else! degA.ne.1
!!$ck    degA = degA -1
!!$ck  end if
  end do ! i=actorb0,actorb
  deallocate(c_psi1)
  deallocate(c_psi2)
!  write(*,*)"# end subroutine  set_psi_pt2"

end subroutine set_psi_pt2

!============================================================================
function intHrel_pt_Qmat(c_p,c_q)
!============================================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
!  use PTcoef
  implicit none
 
  complex(kind=dp),intent(in) :: c_p(4,NMAX_PG), c_q(4,NMAX_PG)
  complex(kind=dp) :: intHrel_pt_Qmat
  type(primitive_gaussian) :: pg

  complex(kind=dp) :: intHrel_pq_14(3),intHrel_pq_23(3),intHrel_pq_32(3),intHrel_pq_41(3)
  complex(kind=dp) :: intHrel_pq_13(3),intHrel_pq_24(3),intHrel_pq_31(3),intHrel_pq_42(3)

  complex(kind=dp) :: func_AMvec

  call set_GammaMatrix

  call copy_DiracOutput_pg(pg)

  call calc_intHrel_pq(1,4,NBS_L,NBS_S,pg,c_p,c_q,intHrel_pq_14) 
  call calc_intHrel_pq(2,3,NBS_L,NBS_S,pg,c_p,c_q,intHrel_pq_23) 
  call calc_intHrel_pq(3,2,NBS_L,NBS_S,pg,c_p,c_q,intHrel_pq_32) 
  call calc_intHrel_pq(4,1,NBS_L,NBS_S,pg,c_p,c_q,intHrel_pq_41) 
  call calc_intHrel_pq(1,3,NBS_L,NBS_S,pg,c_p,c_q,intHrel_pq_13) 
  call calc_intHrel_pq(2,4,NBS_L,NBS_S,pg,c_p,c_q,intHrel_pq_24) 
  call calc_intHrel_pq(3,1,NBS_L,NBS_S,pg,c_p,c_q,intHrel_pq_31) 
  call calc_intHrel_pq(4,2,NBS_L,NBS_S,pg,c_p,c_q,intHrel_pq_42) 
  
  intHrel_pt_Qmat = (0._dp,0._dp)
  
!!$  ! You have to change func_AMvec in sub_AM.f90 too.
!!$  if(direcB==1) then
!!$    !-------------------------------
!!$    ! BM = rot(AM) = (Bmag,0,0)
!!$    ! AM = (0,1/2*Bmag*z,-1/2*Bmag*y)
!!$!    intHrel_pt_Qmat = intHrel_pt_Qmat +func_AMvec(0.d0,2,0.d0,0.d0,1.d0)*(intHrel_pq_13(3)+intHrel_pq_31(3)) -(intHrel_pq_24(3)+intHrel_pq_42(3))
!!$!    intHrel_pt_Qmat = intHrel_pt_Qmat +func_AMvec(0.d0,3,0.d0,1.d0,0.d0)*(intHrel_pq_14(2)+intHrel_pq_23(2)+intHrel_pq_32(2)+intHrel_pq_41(2))
!!$    intHrel_pt_Qmat = intHrel_pt_Qmat -func_AMvec(0.d0,2,0.d0,0.d0,1.d0)*(-IU*(intHrel_pq_14(3)+intHrel_pq_32(3)) +IU*(intHrel_pq_23(3)+intHrel_pq_41(3)))
!!$    intHrel_pt_Qmat = intHrel_pt_Qmat -func_AMvec(0.d0,3,0.d0,1.d0,0.d0)*((intHrel_pq_13(2)+intHrel_pq_31(2)) -(intHrel_pq_24(2)+intHrel_pq_42(2)))
!!$    intHrel_pt_Qmat = intHrel_pt_Qmat *Ze
!!$    !-------------------------------
!!$
!!$  else if(direcB==2) then
!!$    !-------------------------------
!!$    ! BM = rot(AM) = (0,Bmag,0)
!!$    ! AM = (1/2*Bmag*z,0,-1/2*Bmag*x)
!!$!    intHrel_pt_Qmat = intHrel_pt_Qmat +func_AMvec(0.d0,1,0.d0,0.d0,1.d0)*(intHrel_pq_13(3)+intHrel_pq_31(3)) -(intHrel_pq_24(3)+intHrel_pq_42(3))
!!$!    intHrel_pt_Qmat = intHrel_pt_Qmat +func_AMvec(0.d0,3,1.d0,0.d0,0.d0)*(-IU*(intHrel_pq_14(1)+intHrel_pq_32(1)) +IU*(intHrel_pq_23(1)+intHrel_pq_41(1)))
!!$    intHrel_pt_Qmat = intHrel_pt_Qmat -func_AMvec(0.d0,1,0.d0,0.d0,1.d0)*(intHrel_pq_14(3)+intHrel_pq_23(3)+intHrel_pq_32(3)+intHrel_pq_41(3))
!!$    intHrel_pt_Qmat = intHrel_pt_Qmat -func_AMvec(0.d0,3,1.d0,0.d0,0.d0)*((intHrel_pq_13(1)+intHrel_pq_31(1)) -(intHrel_pq_24(1)+intHrel_pq_42(1)))
!!$    intHrel_pt_Qmat = intHrel_pt_Qmat *Ze
!!$    !-------------------------------
!!$
!!$  else if(direcB==3) then
!!$    !-------------------------------
!!$    ! BM = rot(AM) = (0,0,Bmag)
!!$    ! AM = (-1/2*B*y,1/2*B*x,0)
!!$    intHrel_pt_Qmat = intHrel_pt_Qmat -func_AMvec(0.d0,1,0.d0,1.d0,0.d0)*(intHrel_pq_14(2)+intHrel_pq_23(2)+intHrel_pq_32(2)+intHrel_pq_41(2))
!!$    intHrel_pt_Qmat = intHrel_pt_Qmat -func_AMvec(0.d0,2,1.d0,0.d0,0.d0)*(-IU*(intHrel_pq_14(1)+intHrel_pq_32(1)) +IU*(intHrel_pq_23(1)+intHrel_pq_41(1)))
!!$    intHrel_pt_Qmat = intHrel_pt_Qmat *Ze
!!$    !-------------------------------
!!$  else 
!!$    write(*,*)'# direcB should be 1-3.'
!!$  end if

    !-------------------------------
    ! BM = rot(AM) = (Bmag(1),Bmag(2),Bmag(3))
    ! AM = 0.5d0*\vec{BM} \times \vec{r}
    !    = (1/2*Bmagy*z -1/2*Bmagz*y, 1/2*Bmagz*x -1/2*Bmagx*z, 1/2*Bmagx*y -1/2*Bmagy*x)
    intHrel_pt_Qmat = intHrel_pt_Qmat -func_AMvec(0.d0,1,0.d0,0.d0,1.d0)*(intHrel_pq_14(3)+intHrel_pq_23(3)+intHrel_pq_32(3)+intHrel_pq_41(3))
    intHrel_pt_Qmat = intHrel_pt_Qmat -func_AMvec(0.d0,1,0.d0,1.d0,0.d0)*(intHrel_pq_14(2)+intHrel_pq_23(2)+intHrel_pq_32(2)+intHrel_pq_41(2))
    intHrel_pt_Qmat = intHrel_pt_Qmat -func_AMvec(0.d0,2,1.d0,0.d0,0.d0)*(-IU*(intHrel_pq_14(1)+intHrel_pq_32(1)) +IU*(intHrel_pq_23(1)+intHrel_pq_41(1)))
    intHrel_pt_Qmat = intHrel_pt_Qmat -func_AMvec(0.d0,2,0.d0,0.d0,1.d0)*(-IU*(intHrel_pq_14(3)+intHrel_pq_32(3)) +IU*(intHrel_pq_23(3)+intHrel_pq_41(3)))
    intHrel_pt_Qmat = intHrel_pt_Qmat -func_AMvec(0.d0,3,0.d0,1.d0,0.d0)*((intHrel_pq_13(2)+intHrel_pq_31(2)) -(intHrel_pq_24(2)+intHrel_pq_42(2)))
    intHrel_pt_Qmat = intHrel_pt_Qmat -func_AMvec(0.d0,3,1.d0,0.d0,0.d0)*((intHrel_pq_13(1)+intHrel_pq_31(1)) -(intHrel_pq_24(1)+intHrel_pq_42(1)))
    intHrel_pt_Qmat = intHrel_pt_Qmat *Ze
    !-------------------------------

  return
end function intHrel_pt_Qmat

!============================================================================
subroutine pm12(a,aa)
!============================================================================
  implicit none

  character(LEN=1), intent(in) :: a
  integer, intent(out) :: aa

  if(a.eq."+") then
    aa = 1
  else if(a.eq."-") then
    aa = 2
  else
    write(*,*) "# error a in subroutine pm12"
    stop
  end if

end subroutine pm12
!============================================================================
!============================================================================
subroutine energy_sort(eneQ)
!============================================================================
  Use DiracOutput
  use Constants
  use DefineTypes
  use prop_param
  use PTcoef

  implicit none

  double precision, intent(inout) :: eneQ(2*NBS)
  integer :: orbnew(2*NBS)
  integer :: i,j,kk,k,tmporb
  double precision tmpe
  complex(kind(0d0)) tmpc_La(NBS_L,2*NBS),tmpc_Lb(NBS_L,2*NBS)
  complex(kind(0d0)) tmpc_Sa(NBS_S,2*NBS),tmpc_Sb(NBS_S,2*NBS)

  do kk=1,2*NBS
    orbnew(kk) = kk
   do i=1,NBS_L
    tmpc_La(i,kk)=c_psi(1,i,kk)
    tmpc_Lb(i,kk)=c_psi(2,i,kk)
   end do
   do j=1,NBS_S
    tmpc_Sa(j,kk)=c_psi(3,j,kk)
    tmpc_Sb(j,kk)=c_psi(4,j,kk)
   end do
  end do

  do k=1,2*NBS-1
   do kk=k+1,2*NBS
    if(eneQ(k).gt.eneQ(kk)) then
     tmpe=eneQ(k)
     eneQ(k)=eneQ(kk)
     eneQ(kk)=tmpe
     tmporb=orbnew(k)
     orbnew(k)=orbnew(kk)
     orbnew(kk)=tmporb
    end if
   end do
  end do

  do k=1,2*NBS
   do kk=1,2*NBS
    if(k.eq.orbnew(kk)) then
     do i=1,NBS_L
       c_psi(1,i,kk) =tmpc_La(i,k)
       c_psi(2,i,kk) =tmpc_Lb(i,k)
     end do
     do j=1,NBS_S
       c_psi(3,j,kk) =tmpc_Sa(j,k)
       c_psi(4,j,kk) =tmpc_Sb(j,k)
     end do
    end if
   end do
  end do

end subroutine energy_sort
!============================================================================
