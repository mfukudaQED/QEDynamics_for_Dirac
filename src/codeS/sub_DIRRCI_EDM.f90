! Last Change:15-August-2023.
!====================================================================================
! 2013.10.29
! - subroutine read_spinCIinp
! - subroutine calc_CItorq
! - subroutine set_CI_coef(norb)
! 2015.10.05
! 2023.08.15
!============================================================================
!============================================================================
subroutine read_spinCIinp
!read input file (spinci.inp)
!============================================================================
  use DiracOutput
  use param_AM
  use DIRRCI_coef, only : thresh_occdet

  implicit none

  character(len=80) INPUT

!  INPUT = 'spinci.inp'
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
     read (*,'(a)') FIFI4 ! CIcofs.txt
     FIFI4 = trim(FILESFOLDER)//"/"//FIFI4
     write(*,'(a)') FIFI4
     read (*,*) symm ! 0=>atom, 1=>molecule using symmmetry calc, 2=>molecule not using symmetry calc
     read (*,*) kbal
     read (*,*) NROOT
     read (*,*) u1,u2,u3,uu,u4,u5,u6

     !-- magnetic field for initial state
     read (*,*) Bmag(1)  !|vecBM_x|
     read (*,*) Bmag(2)  !|vecBM_y|
     read (*,*) Bmag(3)  !|vecBM_z|

     read (*,*) thresh_occdet
!     stop
!  close(unit=11)

end subroutine read_spinCIinp

!============================================================================
subroutine read_EDMinp(coreorb,actorb,orb,edmflag,DMrw)
!read input file (edm.inp)
!============================================================================
  use DiracOutput
  use param_AM
  use DIRRCI_coef, only : thresh_occdet

  implicit none

  character(len=80) INPUT
  integer,intent(out) :: edmflag(2)
  integer,intent(out) :: coreorb, actorb, orb(3)
  character(len=1),intent(out) ::  DMrw

     read (*,'(a)') FILESFOLDER
     FIFI1 = trim(FILESFOLDER)//"/basis.txt"
     FIFI2 = trim(FILESFOLDER)//"/vectors.txt"
!     read (*,'(a12)') FIFI1 ! basis.txt of Dirac10
!     read (*,'(a12)') FIFI2 ! vectors.txt of Dirac10
     write(*,'(a)') trim(FIFI1)
     write(*,'(a)') trim(FIFI2)
     read (*,'(a)') FIFI3 ! Dirac output file of Dirac10
     FIFI3 = trim(FILESFOLDER)//"/"//FIFI3
     write(*,'(a)') trim(FIFI3)
     read (*,'(a)') FIFI4 ! CIcofs.txt
     FIFI4 = trim(FILESFOLDER)//"/"//FIFI4
     write(*,'(a)') trim(FIFI4)
     read (*,'(a)') FIFI5 ! densmat.txt
     FIFI5 = trim(FILESFOLDER)//"/"//FIFI5
     write(*,'(a)') trim(FIFI5)
     read (*,*) symm ! 1=>atom,homonuclear, 2=>the others but nonsymmetric, 3=>professional option
     read (*,*) kbal
     read (*,*) NROOT
     read (*,*) u1,u2,u3,uu,u4,u5,u6
     read (*,*) coreorb,actorb
     read (*,*) orb(1),orb(2),orb(3)
     read (*,*) edmflag(1) ! integration value
     read (*,*) edmflag(2) ! local quantity

     read (*,*) thresh_occdet
     read (*,'(a1)') DMrw

     !-- external field for initial state
     read (*,*) Bmag(1)   !|vecBM_x|
     read (*,*) Bmag(2)   !|vecBM_y|
     read (*,*) Bmag(3)   !|vecBM_z|
     read (*,*) Eelec(1)  !|vecEM_x|
     read (*,*) Eelec(2)  !|vecEM_y|
     read (*,*) Eelec(3)  !|vecEM_z|
     read (*,*) magdip(1) !|vecmu_x|
     read (*,*) magdip(2) !|vecmu_y|
     read (*,*) magdip(3) !|vecmu_z|

  !---------------------------------------------------------
  ! Input consistency checks
  !---------------------------------------------------------

    if ( (orb(1).lt.1).or.(orb(2).lt.1) ) then
      write(*,*) "Check 'orb' in input."
      stop
    end if
    if ( (DMrw.eq.'K').and.((orb(1).gt.coreorb+actorb).or.(orb(2).gt.coreorb+actorb)) ) then
      write(*,*) "Check 'orb' in input."
      stop
    end if
    if ( (DMrw.ne.'r').and.(DMrw.ne.'w').and.(DMrw.ne.'u').and.(DMrw.ne.'K') ) then
      write(*,*) "Check 'DMrw' in input."
      stop
    end if

    return
end subroutine read_EDMinp

!============================================================================
subroutine calc_CItorq
!============================================================================
  use Precision
  use DiracOutput
  use Constants
  use prop_param
  use param_AM
  use DIRRCI_coef
  use EDM_calculation
  use param_KRCI
  implicit none

  character(LEN=4) :: MOJI4
  character(LEN=12) :: get_time
  character(LEN=24) :: fdate

  integer(kind=8) :: i,j,k,l,m,n
  integer :: p,aa
  integer :: nrcnt
  character(LEN=1) :: a,b

  real(kind=dp),allocatable :: occ(:) ! occupation
  complex(kind=dp),allocatable :: c_psi(:,:,:)

  integer :: edmflag(2)
  integer :: coreorb,actorb,orb(3)
  character(len=1) :: DMrw

  real(kind=8) timea,timeb

  !---------------------------------------------------------------------------------------------
  ! set electron related parameters
  !---------------------------------------------------------------------------------------------

  write(*,*) "#############################################"
  write(*,*) "Set expansion functions for electrons."
  write(*,*) "#############################################"

  write(*,*) '# Calling read_EDMinp'
  call read_EDMinp(coreorb,actorb,orb,edmflag,DMrw)
  write(*,*) '# Calling read_DiracOutput'
  call read_DiracOutput  ! read and set global variables -> defined in sub_readDiracOutput.f90 for DIRRCI
  call set_GammaMatrix

  !---------------------------------------------------------------------------------------------
  ! set natural orbitals
  !---------------------------------------------------------------------------------------------

  write(*,*)
  write(*,'(2x,3a)') '-- Start setting natural orbitals at : ', fdate(), ' --'
  write(*,*)

  allocate(c_psi(NBS0,2*NBS,4))
  c_psi = (0._dp,0._dp)

  do i=1,coreorb
    call index_from_Qmat(i,p,a)
    call pm12(a,aa)
    call copy_DiracOutput_cp_mod(p,aa,c_psi(1:NBS0,i,1:4))
  end do

  if ( (DMrw.eq.'r').or.(DMrw.eq.'w').or.(DMrw.eq.'u') ) then
    call set_CI_coef(actorb,DMrw) ! set c_natu, occup_natu for DIRRCI

    do i=1,actorb
      do j=1,NBS_L
        c_psi(j,i+coreorb,1) = c_natu(1,j,i)
        c_psi(j,i+coreorb,2) = c_natu(2,j,i)
      end do
      do j=1,NBS_S
        c_psi(j,i+coreorb,3) = c_natu(3,j,i)
        c_psi(j,i+coreorb,4) = c_natu(4,j,i)
      end do
    end do
    allocate(occ(coreorb+actorb))
    occ(:) = 0._dp
    do i=1,coreorb
      occ(i) = 1._dp
    end do
    do i=1,actorb
      occ(i+coreorb) = occup_natu(i)
    end do

  else if(DMrw.eq.'K') then

    allocate(occ(coreorb+actorb))
    occ(:) = 0._dp
    do i=1,coreorb
      occ(i) = 1._dp
    end do

    if (actorb.ge.1) then ! for KRCI

      do i=coreorb+1, coreorb+actorb
        call index_from_Qmat(i,p,a)
        call pm12(a,aa)
        call copy_DiracOutput_cp_mod(p,aa,c_psi(1:NBS0,i,1:4))
      end do
      open(unit=45,file=trim(FILESFOLDER)//"/NOcofs.txt")
      do
        read (45,'(a4)') MOJI4
        if (MOJI4.eq.'ROOT') then
          backspace(45)
          read (45,'(4x,i4)') nrcnt
          if (nrcnt.eq.NROOT) exit
        end if
      end do
      read (45,*) ! skip sentence 'Occupation of KR-CI active orbitals:'
      do i = 1, actorb
        read(45,*) occ(i+coreorb)
      end do
      read(45,'(A4)') MOJI4
      if (MOJI4.ne.'END ') then
       write(*,*) 'actorb is wrong.'
       stop
      end if
      close(45)

    end if

  end if ! DMrw

  write(*,*)
  write(*,'(3x,a,i4)')     '- Total occupation number of  core  orbitals : ', int(sum(occ(1:coreorb)))
  write(*,'(3x,a,f20.15)') '- Total occupation number of active orbitals : ', sum(occ(coreorb+1:coreorb+actorb))

  actorb = actorb + coreorb
  write(*,'(3x,a,f20.15)') '- Total occupation number : ', sum(occ(1:actorb))
  write(*,'(3x,a,i4)')     '- Total charge of atoms   : ', sum(cn(1:NAT))

  write(*,*)
  write(*,'(2x,3a)') '-- End setting natural orbitals at : ', fdate(), ' --'
  write(*,*)

  !---------------------------------------------------------------------------------------------
  ! calculate properties
  !---------------------------------------------------------------------------------------------

  call mkdir_today
  call set_pg_normalization_constant

  if(edmflag(1).eq.1) then ! if edmflag(1).ne.1 => skip calculate integration values

  !--- start calc integration values -----------
    write(*,*)
    write(*,'(2x,3a)') '-- Start calc integration values at : ', fdate(), ' --'

    allocate(intPsi_overlap(0:actorb,4,4),intPsi_Lap(0:actorb,4,4))
    allocate(intPsi_x_moment(0:actorb,4,4),intPsi_y_moment(0:actorb,4,4),intPsi_z_moment(0:actorb,4,4))
    allocate(intPsi_xdy_ydx(0:actorb,4,4),intPsi_ydz_zdy(0:actorb,4,4),intPsi_zdx_xdz(0:actorb,4,4))
    allocate(intPsi_ef_x(0:actorb,4,4),intPsi_ef_y(0:actorb,4,4),intPsi_ef_z(0:actorb,4,4))

    call timing_report(.false.)
    call set_int_psi(actorb,occ,c_psi)
    !$omp parallel sections
     !$omp section
     call calc_int_spin(actorb,occ)
     !$omp section
     call calc_int_zeta_potential(actorb,occ)
     !$omp section
     call calc_int_orbital_angular_momentum(actorb,occ)
     !$omp section
     call calc_int_torqAM(actorb,occ)
     !$omp section
     call calc_int_dipole(actorb,occ)
     !$omp section
     call calc_int_torqEDM_ele(actorb,occ)
     !$omp section
     call calc_int_torqEDM_mag(actorb,occ)
     !$omp section
     call calc_int_EeffEDM_ob(actorb,occ)
     !$omp section
     call calc_int_EeffEDM_nuc(actorb,occ)
     !$omp section
     call calc_int_MagHyper(actorb,occ)
    !$omp end parallel sections

    write(*,'(2x,3a)') '-- End calc integration values at : ', fdate(), ' --'
    call timing_report(.true.)

    deallocate(intPsi_overlap,intPsi_Lap)
    deallocate(intPsi_x_moment,intPsi_y_moment,intPsi_z_moment)
    deallocate(intPsi_xdy_ydx,intPsi_ydz_zdy,intPsi_zdx_xdz)
    deallocate(intPsi_ef_x,intPsi_ef_y,intPsi_ef_z)

  end if ! if edmflag(1)

  if(edmflag(2).ne.1) return

  !--- start calc local physical quantities -----------
    write(*,*)
    write(*,'(2x,3a)') '-- Start calc local quantities at : ', fdate(), ' --'
    write(*,*)

  !--- calculation range and mesh ----------------------------------------
    meshx = int(U1/UU) +1
    meshy = int(U2/UU) +1
    meshz = int(U3/UU) +1
    meshyz = meshy*meshz
    mesh = meshx*meshy*meshz
    xori = -UU*(meshx-1)*0.5_dp +U4 ; xend = UU*(meshx-1)*0.5_dp +U4
    yori = -UU*(meshy-1)*0.5_dp +U5 ; yend = UU*(meshy-1)*0.5_dp +U5
    zori = -UU*(meshz-1)*0.5_dp +U6 ; zend = UU*(meshz-1)*0.5_dp +U6
    dx = 0.d0
    if(meshx.ne.1) dx = (xend -xori)/(meshx -1)
    dy = 0.d0
    if(meshy.ne.1) dy = (yend -yori)/(meshy -1)
    dz = 0.d0
    if(meshz.ne.1) dz = (zend -zori)/(meshz -1)

    allocate(x1(meshx),x2(meshy),x3(meshz))
    x1(1:meshx) = (/ (xori +i *dx, i = 0, meshx -1) /)
    x2(1:meshy) = (/ (yori +i *dy, i = 0, meshy -1) /)
    x3(1:meshz) = (/ (zori +i *dz, i = 0, meshz -1) /)

    write(*,'(2x,a)') '* Calculation range and mesh'
    write(*,'(2x,a)') '----------------------------'
    write(*,'(3x,a,3i5)') '# meshx meshy meshz : ', meshx, meshy, meshz
    write(*,'(3x,a,i10)') '# meshyz : ', meshyz
    write(*,'(3x,a,i10)') '# mesh   : ', mesh
    write(*,'(3x,a,2f13.5)') '# xori xend : ', xori, xend
    write(*,'(3x,a,2f13.5)') '# yori yend : ', yori, yend
    write(*,'(3x,a,2f13.5)') '# zori zend : ', zori, zend
    write(*,'(3x,a,3f13.5)') '# dx dy dz  : ', dx, dy, dz

  !----------------------------------------------------------------------
    write(*,*)
    if (orb(1).gt.orb(2)) then
      write(*,'(2x,a)') &
        & '* orbital number to be calculated : 0'
    else
      write(*,'(2x,a,i4,a,i4,a,i4)') &
        & '* orbital number to be calculated : 0 and ', orb(1), ' ..', orb(2), ', ', orb(3)
    end if
    write(*,*)

    allocate(xA(meshx,-2:nmax(1),NAT),yA(meshy,-2:nmax(2),NAT),zA(meshz,-2:nmax(3),NAT))
    allocate(psi_local(mesh,0:actorb,4,4),dpsi_local(mesh,0:actorb,4,4,3))
    allocate(d2psi_local(mesh,0:actorb,4,4,3,3),ddpsi_local(mesh,0:actorb,4,4,3,3))
    allocate(density(mesh),spin(mesh,3),spin_small(mesh,3))
    allocate(grad_N(mesh,3),je(mesh,3))
    allocate(torq(mesh,3))
    allocate(vecA_M(mesh,3),torq_AM(mesh,3))
    allocate(zeta_force(mesh,3))
    allocate(zeta_potential(mesh))
    allocate(torqEDM_ele(mesh,3),torqEDM_mag(mesh,3))
    allocate(vecE_nuc(mesh,3),EeffEDM_nuc(mesh))
    allocate(chirality_density(mesh))
    allocate(tau(mesh,9))


    call timing_report(.false.)
    call set_pg_orbital_angular_momentum
    call set_psi_local(actorb,occ,c_psi)
    !$omp parallel sections
      !$omp section
      call calc_density(0)
      !$omp section
      call calc_spin(0)
      call calc_EeffEDM_nuc
      !$omp section
      call calc_grad_N(0)
      !$omp section
      call calc_spin_torq(0)
      !$omp section
      call calc_je(0)
      call set_vecA_M
      call calc_spin_torq_AM
      !$omp section
      call calc_zeta_potential(0)
      !$omp section
      call calc_zeta_force(0)

    !$omp end parallel sections
    call timing_report(.true.)
    call print_quantities(today,'(15es24.14)',0)
    do n = orb(1),orb(2),orb(3)
      !$omp parallel sections
        !$omp section
        call calc_density(n)
        !$omp section
        call calc_spin(n)
        call calc_EeffEDM_nuc
        !$omp section
        call calc_grad_N(n)
        !$omp section
        call calc_je(n)
        call calc_spin_torq_AM
        !$omp section
        call calc_spin_torq(n)
        !$omp section
        call calc_zeta_potential(n)
        !$omp section
        call calc_zeta_force(n)
      !$omp end parallel sections
      call print_quantities(today,'(15es24.14)',n)
    end do

    write(*,*)
    write(*,'(2x,3a)') '-- End calc local quantities at : ', fdate(), ' --'

    deallocate(x1,x2,x3,xA,yA,zA)
    deallocate(psi_local,dpsi_local)
    deallocate(d2psi_local,ddpsi_local)
    deallocate(density,spin_small)
    deallocate(grad_N,je)
    deallocate(torq)
    deallocate(vecA_M,torq_AM)
    deallocate(zeta_force)
    deallocate(zeta_potential)
    deallocate(torqEDM_ele,torqEDM_mag)
    deallocate(chirality_density)
    deallocate(tau)

  return
end subroutine calc_CItorq

!============================================================================
subroutine set_CI_coef(norb,DMrw)
! get norb, c_natu, occup_natu
!============================================================================
  use Precision
  use DiracOutput
  use DefineTypes
  use DIRRCI_coef
!$ use omp_lib

  implicit none

  integer, intent(out) :: norb
  character(len=1),intent(in) ::  DMrw

  integer :: i,j,k,l
  complex(kind(0d0)), allocatable :: a(:,:,:),a2(:,:,:)
  integer, allocatable :: ras(:,:),n(:,:)
  double precision, allocatable :: enorb(:),eigvals(:)
  complex(kind(0d0)), allocatable :: dnstmtx(:,:), eigwork(:) ! this may be changed to complex array
  complex(kind(0d0)), allocatable :: dnstmtx2(:,:) ! this may be changed to complex array
  complex(kind(0d0)), allocatable :: tmpdnstmtx(:,:,:)
  integer, allocatable :: eigiwork(:)
  double precision,allocatable :: rwork(:)
  complex(kind(0d0)), allocatable :: occdet(:)
  complex(kind(0d0)) :: cdnsty
  integer neact,detcnt
  integer ndet1,ndet2,ne1,ne2,ne3,ne4,tmpck
  integer itemp,itemp2,itemp3,itemp4
  integer detcntnum
  double precision minus1,relsgn,tmp2
  character(len=14) MOJI
  character(len=4) MOJI4
  character(len=1), allocatable :: KRAM(:)
  character(len=12) get_time

  double precision, allocatable :: aa(:,:)
  double precision, allocatable :: xx(:,:),yy(:,:),zz(:,:)

  integer, allocatable :: activeorb(:,:)
!  Logical, allocatable :: activeorb(:,:)
  integer, allocatable :: sumactiveorb(:)
  double precision :: timea, timeb
  integer itmpk,nrcnt
  integer norbh ! norb/2
  character(len=24) today
  integer :: denssub1(2),denssub2(2)

!      call FDATE(today)

!        write(*,*)omp_get_max_threads()
!       !$omp parallel do private(j)
!       do j=1,100
!          write(*,*) omp_get_thread_num()+1
!       end do
!       !$omp end parallel do
!       stop

      open (unit=41,file=FIFI3) !FIFI3 = 'DIRAC.out'
      do
         read (41,'(A14)') MOJI
         if (MOJI.eq.' Core Energy :') exit
      end do
      read (41,1300) norb
      write(*,*) '# norb : ', norb
      close(unit=41)
 1300 FORMAT(47X,I4)

      allocate(enorb(norb))
      allocate(eigwork(1 + 6 * norb + 2 * norb**2) )
      allocate(eigiwork(3 + 5 * norb**2) )
      allocate(eigvals(norb))
      allocate(occup_natu(norb))
      allocate(KRAM(norb))
      allocate(dnstmtx(norb,norb) )
      allocate(dnstmtx2(norb,norb) )

      allocate(a(NBS0,4,norb),a2(NBS0,4,norb))
      allocate(c_natu(4,NBS0,2*NBS))
      c_natu(:,:,:) = (0.d0,0.d0)

      allocate(aa(NBS0,4),n(NBS0,4))
      allocate(xx(NBS0,4),yy(NBS0,4),zz(NBS0,4))

      open(unit=44, file=FIFI4)
         call CIread_1(norb,NROOT,enorb,KRAM,neact,detcnt,detcntnum)
      close(unit=44)

      allocate(ras(neact,detcntnum))
      allocate(occdet(detcntnum))
!      allocate(tmpdnstmtx(detcntnum,norb,norb) )
      allocate(tmpdnstmtx(omp_get_max_threads(),norb,norb) )
      ras(:,:) =0

      if (DMrw.eq.'w') then
        open(unit=44, file=FIFI4)
!           call CIread_2(norb,NROOT,enorb,ras,neact,detcnt,occdet,detcntnum)
           call CIread_3(norb,NROOT,enorb,ras,neact,detcnt,occdet,detcntnum)
        close(unit=44)
      end if

!      write(*,*) 'NBS0=',NBS0,'NBS=',NBS

      allocate(activeorb(norb,detcnt))
      allocate(sumactiveorb(detcnt))

!      !normalization
!      call convert_pg_QEDtoMRDFT(aa,n,xx,yy,zz)

      !convert QED coef to MRDFT coef
      call convert_coef_MRDFTCI(a,KRAM,norb,enorb)
!       do i = 1, nbs0
!       do j = 1, 4
!       do k = 1, norb
!         write(53,'(2e22.14)') a(i,j,k)!!$w
!       end do
!       end do
!       end do
!      stop


     write(*,'(1x,2a)') '* Start setting density matrix : ', get_time()

     if(DMrw=='w') then
!c
!c  Make density matrix
!c
      write(*,'(4x,a)') 'Make Density Matrix'

      call cpu_time(timea)

       activeorb(:,:) = 0
!       activeorb(:,:) = .false.

       !$omp parallel do private(ndet1,j)
       do ndet1=1, detcnt
         do j=1, neact
           activeorb(ras(j,ndet1),ndet1) = 1
!           activeorb(ras(j,ndet1),ndet1) = .true.
         end do ! j->norb
!         write(*,'(30i1)') activeorb(:,ndet1)
       end do ! ndet1 -> detcnt
       !$omp end parallel do

       sumactiveorb(:) = 0
       !$omp parallel do private(ndet1,ne3)
       do ndet1=1, detcnt
         do ne3=1, neact
              sumactiveorb(ndet1) = sumactiveorb(ndet1) +activeorb(ne3,ndet1)
         end do ! ne3->norb
       end do ! ndet1 -> detcnt
       !$omp end parallel do


       minus1 = -1.0d0
       tmpdnstmtx(:,:,:)  = (0.0d0,0.0d0)
       itmpk=0

       !$omp parallel do private(ndet1,ndet2,ne3,ne4,tmpck,denssub1,denssub2) schedule(dynamic)
       do ndet1 = detcnt,2,-1! 1, detcnt
       do ndet2 = ndet1-1,1,-1!1, ndet1

           if(abs(sumactiveorb(ndet1)-sumactiveorb(ndet2)).gt.1)  goto 106
           tmpck = 1
           ne4 = 1

         do ne3=1,norb
           if(activeorb(ne3,ndet1).eq.activeorb(ne3,ndet2)) then
             if (activeorb(ne3,ndet1).eq.1)  ne4 = ne4 + 1
           else
             if(tmpck.eq.3) goto 106
             denssub1(tmpck) = ne3  !denssub1->p,q for D_{pq} norb_index
             denssub2(tmpck) = ne4  !denssub2->p,q for D_{pq} neact_index
             tmpck = tmpck + 1
           end if
           if((tmpck.eq.3).and.(ne4.eq.neact)) goto 105
         end do !ne3

  105    continue
         itmpk = itmpk + 1

!          if((activeorb(denssub1(1),ndet1).eq.1).and.(activeorb(denssub1(2),ndet2).eq.1)) then
!             if(mod(denssub2(1)-denssub2(2),2).eq.0) then
!                tmpdnstmtx(ndet2,denssub1(1),denssub1(2)) = &
!                     & tmpdnstmtx(ndet2,denssub1(1),denssub1(2)) + dconjg(occdet(ndet1))*occdet(ndet2)
!                tmpdnstmtx(ndet2,denssub1(2),denssub1(1)) = &
!                     & tmpdnstmtx(ndet2,denssub1(2),denssub1(1)) + dconjg(occdet(ndet2))*occdet(ndet1)
!             else
!                tmpdnstmtx(ndet2,denssub1(1),denssub1(2)) = &
!                       tmpdnstmtx(ndet2,denssub1(1),denssub1(2)) - dconjg(occdet(ndet1))*occdet(ndet2)
!                tmpdnstmtx(ndet2,denssub1(2),denssub1(1)) = &
!                     & tmpdnstmtx(ndet2,denssub1(2),denssub1(1)) - dconjg(occdet(ndet2))*occdet(ndet1)
!             end if
!          end if
!
!          if((activeorb(denssub1(2),ndet1).eq.1).and.(activeorb(denssub1(1),ndet2).eq.1)) then
!             if(mod(denssub2(2)-denssub2(1),2).eq.0) then
!                tmpdnstmtx(ndet2,denssub1(2),denssub1(1)) = &
!                     & tmpdnstmtx(ndet2,denssub1(2),denssub1(1)) + dconjg(occdet(ndet1))*occdet(ndet2)
!                tmpdnstmtx(ndet2,denssub1(1),denssub1(2)) = &
!                     & tmpdnstmtx(ndet2,denssub1(1),denssub1(2)) + dconjg(occdet(ndet2))*occdet(ndet1)
!             else
!                tmpdnstmtx(ndet2,denssub1(2),denssub1(1)) = &
!                       tmpdnstmtx(ndet2,denssub1(2),denssub1(1)) - dconjg(occdet(ndet1))*occdet(ndet2)
!                tmpdnstmtx(ndet2,denssub1(1),denssub1(2)) = &
!                     & tmpdnstmtx(ndet2,denssub1(1),denssub1(2)) - dconjg(occdet(ndet2))*occdet(ndet1)
!             end if
!          end if

          if((activeorb(denssub1(1),ndet1).eq.1).and.(activeorb(denssub1(2),ndet2).eq.1)) then
             if(mod(denssub2(1)-denssub2(2),2).eq.0) then
                tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) = &
                     & tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) + dconjg(occdet(ndet1))*occdet(ndet2)
                tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) = &
                     & tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) + dconjg(occdet(ndet2))*occdet(ndet1)
             else
                tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) = &
                       tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) - dconjg(occdet(ndet1))*occdet(ndet2)
                tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) = &
                     & tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) - dconjg(occdet(ndet2))*occdet(ndet1)
             end if
          end if

          if((activeorb(denssub1(2),ndet1).eq.1).and.(activeorb(denssub1(1),ndet2).eq.1)) then
             if(mod(denssub2(2)-denssub2(1),2).eq.0) then
                tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) = &
                     & tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) + dconjg(occdet(ndet1))*occdet(ndet2)
                tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) = &
                     & tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) + dconjg(occdet(ndet2))*occdet(ndet1)
             else
                tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) = &
                       tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) - dconjg(occdet(ndet1))*occdet(ndet2)
                tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) = &
                     & tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) - dconjg(occdet(ndet2))*occdet(ndet1)
             end if
          end if

  106        continue

       end do !ndet2
       end do !ndet1
       !$omp end parallel do

       !$omp parallel do private(ndet1,ne1)!! schedule(dynamic)
       do ndet1 = 1,detcnt
          do ne1=1,neact
             tmpdnstmtx(omp_get_thread_num()+1,ras(ne1,ndet1),ras(ne1,ndet1)) = &
                    tmpdnstmtx(omp_get_thread_num()+1,ras(ne1,ndet1),ras(ne1,ndet1)) + dconjg(occdet(ndet1))*occdet(ndet1)
          end do
       end do
       !$omp end parallel do

!       !$omp parallel do private(ndet1,ne1)!! schedule(dynamic)
!       do ndet1 = 1,detcnt
!          do ne1=1,neact
!             tmpdnstmtx(ndet1,ras(ne1,ndet1),ras(ne1,ndet1)) = &
!                    tmpdnstmtx(ndet1,ras(ne1,ndet1),ras(ne1,ndet1)) + dconjg(occdet(ndet1))*occdet(ndet1)
!          end do
!       end do
!       !$omp end parallel do

       write(*,'(4x,a)') 'Make Density Matrix 4'
       dnstmtx(:,:)  = (0.0d0,0.0d0)
       !$omp parallel do private(k,i,j)!! schedule(dynamic)
       do j=1, norb
         do i=1, norb
           do k = 1, omp_get_max_threads()
                dnstmtx(i,j) = dnstmtx(i,j) + tmpdnstmtx(k,i,j)
           end do
!           do ndet1 = 1, detcnt
!                dnstmtx(i,j) = dnstmtx(i,j) + tmpdnstmtx(ndet1,i,j)
!           end do
         end do
       end do
       !$omp end parallel do


       call cpu_time(timeb)

      write(*,'(4x,a)')       'Finish : Density Matrix'
      write(*,'(4x,a,f15.5)') 'CPU time is ', timeb - timea
      write(*,'(4x,a,i15)')   'Products of determinant : ', itmpk

     else if(DMrw=='u') then
!c
!c  Use density matrix calculated by DIRAC
!c

      write(*,'(4x,a)') 'Use density matrix calculated by DIRAC'
      write(*,'(4x,a)') 'Read from '//trim(FIFI5)

       open(unit=45,file=FIFI5)

       do
          read (45,'(a4)') MOJI4
          if (MOJI4.eq.'ROOT') then
            backspace(45)
            read (45,'(4x,i4)') nrcnt
              if (nrcnt.eq.NROOT) exit
          end if
       end do

      write(16,'(a4,i4)') 'ROOT', nrcnt

       do i = 1, norb
          do j = 1, norb
             read(45,*) nrcnt, nrcnt, dnstmtx(i,j)
             write(16,*) i, j, dnstmtx(i,j) ! read check
          end do
       end do

       close(45)

      write(*,'(4x,a)') 'Finish : Density Matrix'

     end if ! DMrw w,u

     if(DMrw=='w'.or.DMrw=='u') then

      write(*,'(4x,a)') 'Make Natural Orbital'

      cdnsty = (0.0d0,0.0d0)
      do i=1, norb
      do j=1, norb
          dnstmtx2 (i,j) = dnstmtx(i,j)
          if (i.eq.j) then
            write (*,*) i,j,dnstmtx(i,j)
            write (15,*) i,j,dnstmtx(i,j)
            cdnsty = cdnsty + dnstmtx(i,j)
          end if
          write(13,*) i,j,dnstmtx(i,j)
          if (abs(dnstmtx(i,j)-dconjg(dnstmtx(j,i))).gt.1.d-15) then
            write(*,*) 'diff',i,j,dnstmtx(i,j)-dconjg(dnstmtx(j,i))
            write(14,*) i,j,dnstmtx(i,j)-dconjg(dnstmtx(j,i))
          end if
       end do !j
       end do !i
       write(*,*) cdnsty
       write(15,*) 'total',cdnsty

       itemp  = norb*(norb + 2) !lwork
       itemp2 = 5*norb + 3 !liwork
       itemp3 = 2*norb*norb + 5*norb + 1 !lrwork
       allocate(rwork(itemp3))
!       call dsyevd('V', 'U', norb, dnstmtx, norb, eigvals, &
!                   eigwork, itemp, eigiwork, itemp2, itemp3)
       call zheevd('V', 'U', norb, dnstmtx, norb, eigvals, &
                   eigwork, itemp, rwork, itemp3, eigiwork, itemp2, itemp4)

      call cpu_time(timea)

      write(*,'(4x,a)')       'Finish : Natural Orbital'
      write(*,'(4x,a,f15.5)') 'CPU time is : ', timea - timeb

       write(100,*) eigvals
       do j = 1, norb
       do i = 1, norb
         write(101,*) i,j,dnstmtx(i,j)
       end do
       end do

       do i = 1, norb
       do j=  1, norb
         tmp2 = 0.0d0
       do k = 1, norb
       do l=  1, norb
         tmp2 = tmp2 + dble(dconjg(dnstmtx(k,i)) * dnstmtx2(k,l) * dnstmtx(l,j))
       end do !ndet2
       end do !ndet1
         if(i.eq.j) write(100,*) tmp2, eigvals(i) - tmp2, (eigvals(i) - tmp2)/tmp2
         if(abs(tmp2).gt.10.0d-15) write(14,*) i,j,tmp2
       end do !ndet2
       end do !ndet1

       a2(:,:,:) = dcmplx(0.0D0,0.0D0)

       do i = 1, nbs0
       do j = 1, 4
       do l = 1, norb
       do k = 1, norb
         a2(i,j,l) = a2(i,j,l) + dnstmtx(k,l) * a(i,j,k)
       end do !k
       end do !l
       end do !j
       end do !i

     end if ! DMrw w or u

       if(DMrw=='r') then
       end if

       open(unit=50,file='c_natu.dat')
       open(unit=52,file='occup_natu.dat')
       if(DMrw=='w'.or.DMrw=='u') then
          do i = 1, nbs0
          do j = 1, 4
          do k = 1, norb
            write(50,'(2e22.14)') a2(i,j,k)!!$w
!            write(53,'(2e22.14)') a(i,j,k) !!$w
!!$r             read(50,'(2e22.14)') a2(i,j,k)
!!!$r            write(51,'(2e22.14)') a2(i,j,k)
!            a(i,j,k) = a2(i,j,k)
            c_natu(j,i,k) = a2(i,j,k) !i<->j is not bug
          end do !k
          end do !j
          end do !i
          do k = 1, norb
!!$r             read(52,'(e22.14)') eigvals(k)
             write(52,'(e22.14)') eigvals(k)!!$w
             occup_natu(k) = eigvals(k)
          end do

!stop

       else if(DMrw=='r') then
          write(*,'(4x,a)') 'Reading density matrix...'
          do i = 1, nbs0
          do j = 1, 4
          do k = 1, norb
             read(50,'(2e22.14)') a2(i,j,k)
             c_natu(j,i,k) = a2(i,j,k) !i<->j is not bug
          end do !k
          end do !j
          end do !i
          do k = 1, norb
             read(52,'(e22.14)') eigvals(k)
             occup_natu(k) = eigvals(k)
          end do

       else
          write(*,*) 'error DMrw'
          stop
       end if
       close(unit=50)
       close(unit=52)

       write(*,'(1x,2a)') '* End setting density matrix : ', get_time()
!       stop

!       write(*,*)'call calc_calEnpqm :', today
!       call calc_calEnpqm(norb,neact,detcnt,ras,occdet)
!       write(*,*)'end calc_calEnpqm :', today

      deallocate(enorb)
      deallocate(eigwork)
      deallocate(eigiwork)
      deallocate(eigvals)
      deallocate(KRAM)
      deallocate(dnstmtx)
      deallocate(dnstmtx2)
      deallocate(a,a2)
      deallocate(aa,n)
      deallocate(xx,yy,zz)
      deallocate(ras)
      deallocate(occdet)
      deallocate(tmpdnstmtx)

end subroutine set_CI_coef

!============================================================================
subroutine calc_calEnpqm(norb,neact,detcnt,ras,occdet)
!subroutine calc_calEnpqm(norb,neact,detcnt,ras,occdet,calEnpqm)
!============================================================================
  use Precision
  use DiracOutput
  use DIRRCI_coef
  use AA_nucatt_pg, only : calEnpqm

  implicit none

  integer, intent(in) :: norb
  integer, intent(in) :: neact
  integer, intent(in) :: detcnt
  integer, intent(in) :: ras(neact,detcnt)
  complex(kind=dp),intent(in) :: occdet(detcnt)

!  complex(kind=dp),intent(out) :: calEnpqm(norb,norb,norb,norb)


  integer :: i,j,k,l
  complex(kind=dp)  :: tmpcalEnpqm(detcnt,norb,norb,norb,norb)
  integer ndet1,ndet2,ne11,ne12,ne21,ne22,ne3,ne4,tmpck
  double precision minus1,relsgn,tmp2

       minus1 = -1.0d0
       tmpcalEnpqm(:,:,:,:,:)  = (0._dp,0._dp)

       !$omp parallel do private(ndet1,ndet2,ne11,ne12,ne21,ne22,ne3,ne4,tmpck,relsgn)
       do ndet1 = 1, detcnt
       do ndet2 = 1, detcnt

         do ne11 = 1, neact
         do ne12 = 1, neact
           if(ne12.eq.ne11) cycle
         do ne21 = 1, neact
         do ne22 = 1, neact
           if(ne22.eq.ne21) cycle
           tmpck = 0
!        check the other states
           relsgn = 1.0d0
           do ne3=1,neact
           do ne4=1,neact
             if (ne3.eq.ne11) cycle
             if (ne3.eq.ne12) cycle
             if (ne4.eq.ne21) cycle
             if (ne4.eq.ne22) cycle
             if (ras(ne3,ndet1).eq.ras(ne4,ndet2)) then
               tmpck=tmpck+1
               if(mod(ne3-ne4,2).ne.0) relsgn = relsgn * minus1
             end if
           end do !ne4
           end do !ne3
           !!!CAUTION
           if (tmpck.lt.neact-2) cycle
!           if ((tmpck.lt.neact-2).or.(tmpck.gt.neact)) cycle

           tmpcalEnpqm(ndet1,ras(ne11,ndet1),ras(ne12,ndet1),ras(ne21,ndet2),ras(ne22,ndet2)) = &
                & tmpcalEnpqm(ndet1,ras(ne11,ndet1),ras(ne12,ndet1),ras(ne21,ndet2),ras(ne22,ndet2)) &
                & + dconjg(occdet(ndet1))*occdet(ndet2)*relsgn

         end do !ne22
         end do !ne21
         end do !ne12
         end do !ne11

       end do !ndet2
       write(*,*)'ndet1',ndet1,'/',detcnt
       end do !ndet1
       !$omp end parallel do

       allocate(calEnpqm(norb,norb,norb,norb))
       calEnpqm(:,:,:,:)  = (0.0d0,0.0d0)
       do i=1, norb
       do j=1, norb
       do k=1, norb
       do l=1, norb
         do ndet1 = 1, detcnt
            calEnpqm(i,j,k,l) = calEnpqm(i,j,k,l) + tmpcalEnpqm(ndet1,i,j,k,l)
         end do !ndet1
       end do
       end do
       end do
       end do

end subroutine calc_calEnpqm

!============================================================================
subroutine calc_CIphys_norb(actorb0,actorb,occ)
  use Precision
  Use DiracOutput
  use Constants
  use prop_param

  implicit none

  real(kind=8) :: x,y,z

  integer :: i,j,k,l
  integer :: rr,pp,qq
  integer :: nn,mm
  integer :: actorb0, actorb
  character(LEN=80) :: filetorq,filezeta,filespin

  complex(kind=dp) :: density_pt_Qmat, rho_pt_Qmat !function
  complex(kind=dp) :: t_pt_Qmat, t_Arad_pt_Qmat, s_pt_Qmat, zeta_pt_Qmat, chiral_pt_Qmat !function
  complex(kind=dp) :: ssmall_pt_Qmat, spol_pt_Qmat !function
  complex(kind=dp) :: j_pt_Qmat, t_AM_pt_Qmat !function
  complex(kind=dp) :: divj_pt_Qmat, AA_pt_Qmat !function
  complex(kind=dp) :: t_AA_pt_Qmat, Eeff_EDM_ele_pt_Qmat, Eeff_EDM_nuc_pt_Qmat !function
  complex(kind=dp) :: EDMtorqE_pt_Qmat,EDMtorqB_pt_Qmat !function
  complex(kind=dp) :: Hso_nuc_pt_Qmat, Hso_ele_pt_Qmat !function
  complex(kind=dp) :: Eele_pt_Qmat(3) !not function
  complex(kind=dp) :: js_pt_Qmat !function
  complex(kind=dp) :: intSpin_pt_Qmat
  complex(kind=dp) :: intspol_pt_Qmat, intgamk_pt_Qmat !function
  complex(kind=dp) :: intjmoment_pt_Qmat !function
  complex(kind=dp) :: intEeff_EDM_ele_pt_Qmat !function
  complex(kind=dp) :: intEeff_EDM_nuc_pt_Qmat !function

  real(kind=dp) :: occ(actorb) ! occupation of natural orb (norb)

  complex(kind=dp), allocatable :: tmp_t0(:,:,:,:,:), tmp_z0(:,:,:,:,:), tmp_s0(:,:,:,:,:)


  allocate(tmp_z0(3,meshx,meshy,meshz,actorb))
  allocate(tmp_t0(3,meshx,meshy,meshz,actorb))
  allocate(tmp_s0(3,meshx,meshy,meshz,actorb))

    dx = 0.d0
    if(meshx.ne.1) dx = (xend - xori)/(meshx-1)
    dy = 0.d0
    if(meshy.ne.1) dy = (yend - yori)/(meshy-1)
    dz = 0.d0
    if(meshz.ne.1) dz = (zend - zori)/(meshz-1)

    write(*,*)'dx dy dz =',dx,dy,dz
         !$omp parallel do private(x,y,z,nn,l,i,j,k)
         do i=1,meshx
            x = xori + dx*(i-1)
!            if(meshx.eq.1) x = x0
            if(meshx.eq.1) x = 0.d0

            do j=1,meshy
               y = yori + dy*(j-1)
!               if(meshy.eq.1) y = y0
               if(meshy.eq.1) y = 0.d0

               do k=1,meshz
                  z = zori + dz*(k-1)
!                  if(meshz.eq.1) z = z0
                  if(meshz.eq.1) z = 0.d0

                     do l=1,3
                        do nn=actorb0,actorb
                           tmp_z0(l,i,j,k,nn) = zeta_pt_Qmat(l,x,y,z,nn,nn)
                           tmp_t0(l,i,j,k,nn) = t_pt_Qmat(l,x,y,z,nn,nn)
                           tmp_s0(l,i,j,k,nn) = s_pt_Qmat(l,x,y,z,nn,nn)
                        end do
                     end do

!        write(*,*)'z', k,'/',meshz
               end do !k=1,meshz
            end do !j=1,meshy
            write(*,*) i,'/',meshx
         end do !i=1,meshx
         !$omp end parallel do


      do nn=actorb0,actorb
        write(filezeta,'(a,i5.5,a)')'data/zeta_',nn,'.dat'
        write(filetorq,'(a,i5.5,a)')'data/torq_',nn,'.dat'
        write(filespin,'(a,i5.5,a)')'data/spin_',nn,'.dat'
        open(unit=31,file=filezeta)
        open(unit=32,file=filetorq)
        open(unit=33,file=filespin)

         do i=1,meshx
            x = xori + dx*(i-1)
 !           if(meshx.eq.1) x = x0

            do j=1,meshy
               y = yori + dy*(j-1)
 !              if(meshy.eq.1) y = y0

               do k=1,meshz
                  z = zori + dz*(k-1)
 !                 if(meshz.eq.1) z = z0

                 write(31,'(15es24.14)') x, y, z, tmp_z0(1,i,j,k,nn), tmp_z0(2,i,j,k,nn), tmp_z0(3,i,j,k,nn) &
                                             & , tmp_z0(1,i,j,k,nn)*occ(nn), tmp_z0(2,i,j,k,nn)*occ(nn), tmp_z0(3,i,j,k,nn)*occ(nn)
                 write(32,'(15es24.14)') x, y, z, tmp_t0(1,i,j,k,nn), tmp_t0(2,i,j,k,nn), tmp_t0(3,i,j,k,nn) &
                                             & , tmp_t0(1,i,j,k,nn)*occ(nn), tmp_t0(2,i,j,k,nn)*occ(nn), tmp_t0(3,i,j,k,nn)*occ(nn)
                 write(33,'(15es24.14)') x, y, z, tmp_s0(1,i,j,k,nn), tmp_s0(2,i,j,k,nn), tmp_s0(3,i,j,k,nn) &
                                             & , tmp_s0(1,i,j,k,nn)*occ(nn), tmp_s0(2,i,j,k,nn)*occ(nn), tmp_s0(3,i,j,k,nn)*occ(nn)

               end do
             end do
           end do
        end do !nn

end subroutine calc_CIphys_norb
