! Last Change:08-Jun-2014.
!============================================================================
subroutine tdpt
! HF, H'=-ZeeAM^i \gamma^0 \gamma^i
! save calE
!============================================================================
  use Precision
  Use DiracOutput
  use Constants
  use DefineTypes, only: NMAX_PG
  implicit none

  integer :: i,j,k,l,n,p,aa
  character(LEN=1) :: a,b

  integer :: it ! time counter
  real(kind=dp) :: time 

  !----------------------------------------------------
  !   parameters to be read from inifile
  !----------------------------------------------------
  integer(kind=dp) :: NTIME
  integer :: nprint ! print out every nprint step (it)
  integer :: NEL
  !----------------------------------------------------

  complex(kind=dp) :: c_p(4,NMAX_PG)
  complex(kind=dp),allocatable :: calE(:,:)  ! expectation value of (e^dagger e)
  complex(kind=dp),allocatable :: calE0(:,:)  ! expectation value of (e^dagger e)
  complex(kind=dp),allocatable :: c_p0(:,:,:)
  real(kind=dp),allocatable :: eneQ0(:)
  complex(kind=dp),allocatable :: e(:),e_0(:),e_1(:)

  complex(kind=dp),allocatable :: int_k0Hn0(:,:)
  complex(kind=dp) :: intHrel_pt_Qmat

  !set param_AM
  call read_param_AM

  call set_paramini_tdpt(NEL,NTIME,nprint) ! read param.ini, qed.inp, DiracOutput, set gammamatrix
  call set_GammaMatrix

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
  allocate(eneQ0(2*NBS))
  allocate(c_p0(4,NMAX_PG,2*NBS))
  allocate(calE(4*NBS,4*NBS))
  allocate(calE0(4*NBS,4*NBS))
  allocate(e(2*NBS),e_0(2*NBS),e_1(2*NBS))
  allocate(int_k0Hn0(2*NBS,2*NBS))

  c_p0(:,:,:)=(0.d0,0.d0)
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


!!$      !-----  calc int_k0Hn0 = <psi^0_k|H'|psi^0_n> --------------
!!$      ! H'=-ZeeAM^i \gamma^0 \gamma^i
!!$      int_k0Hn0(:,:)=(0.d0,0.d0)
!!$!      !$omp parallel do private(n,k)
!!$      do n=1,2*NBS
!!$        do k=1,2*NBS
!!$          int_k0Hn0(k,n) = intHrel_pt_Qmat(c_p0(:,:,k),c_p0(:,:,n))
!!$!          write(*,*)k,n,int_n0Hir0(k,n)
!!$        end do
!!$      end do
!!$!      !$omp end parallel do

      e_0(:) = (0.d0,0.d0)
      do n=1,NEL
        e_0(n) = (1.d0,0.d0)
      end do

  write(*,*) "#############################################"
  write(*,*) "Start time evolution loop ... "
  write(*,*) "#############################################"

  do it = 0,NTIME  
     time = DeltaT*it  ! DeltaT is defined in GammaMatrix and set in the main routine.

     !print out every nprint steps.
     if(mod(it,nprint).eq.0) then

!!$      e_1(:) = (0.d0,0.d0)
!!$      do k=1,2*NBS
!!$        do n=1,2*NBS
!!$          if (abs(eneQ0(k)-eneQ0(n)).lt.1.d-10) then
!!$            e_1(k) = e_1(k) - IU *int_k0Hn0(k,n) *time
!!$          else
!!$            e_1(k) = e_1(k) - int_k0Hn0(k,n) *(exp(IU*(eneQ0(k)-eneQ0(n))*time)-1.d0)/(eneQ0(k)-eneQ0(n))
!!$          end if
!!$        end do
!!$      end do
!!$
!!$      do k=1,2*NBS
!!$        e(k) = e_0(k) + e_1(k)
!!$      end do
!!$
!!$      calE(:,:) = (0.d0,0.d0)
!!$      do n=2*NBS+1,4*NBS
!!$         calE(n,n) = (1._dp,0._dp)
!!$      end do
!!$      do n=1,2*NBS
!!$        do k=1,2*NBS
!!$          calE(k,n) = dconjg(e(k))*e(n)
!!$        end do
!!$      end do

!----------------------------------------------
      calE0(:,:) = (0.d0,0.d0)
      calE(:,:) = (0.d0,0.d0)
      do n=1,NEL
         calE(n,n) = (1.d0,0.d0)
      end do
      do n=2*NBS+1,4*NBS
         calE0(n,n) = (1._dp,0._dp)
         calE(n,n) = (1._dp,0._dp)
      end do
      do n=1,2*NBS
        do k=1,2*NBS
          calE(k,n) = calE0(k,n)*exp(IU*(eneQ0(k)-eneQ0(n))*time)
        end do
      end do
      call save_calE(it,time,NEL,NTIME,calE)

     end if
  end do

  deallocate(eneQ0)
  deallocate(c_p0)
  deallocate(calE,calE0)
  deallocate(e,e_0,e_1)
  deallocate(int_k0Hn0)

  write(*,*)' # end calc tdpt'

end subroutine tdpt
