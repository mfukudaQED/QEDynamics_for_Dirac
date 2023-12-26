! Last Change:28-Sep-2015.
!====================================================================================
! 2013.11.07!
! - subroutine calc_prop_int
! - function intSpin_mat(k,p,a,q,b)
! - subroutine calc_intdensmom_mat(p,a,q,b,intdensmom_mat)
! - subroutine calc_intOrbang_mat(p,a,q,b,intOrbang_mat)
! - subroutine calc_intOrbang_AM_mat(p,a,q,b,intOrbang_AM_mat)
! - function inttau_mat(k,l,p,a,q,b)
! - function intjmoment_mat(i,j,p,a,q,b)
! - function intMagHyp_mat(posR,p,a,q,b)
! - subroutine calc_inttau_AM(tmp_intjmom_Qmat,inttau_AM_Qmat)
! - subroutine calc_inttorq_AM(inttau_AM_Qmat,inttorq_AM_Qmat)
! - subroutine calc_inttau_Arad(t,tmp_intFj_Qmat,inttau_Arad_Qmat)
! - subroutine calc_inttorq_Arad(inttau_Arad_Qmat,inttorq_Arad_Qmat)
! - function intHrel_mat(p,a,q,b)
! - subroutine calc_intHrel_pq(ip,iq,NL,NS,pg,c_p,c_q,intHrel_pq)
! - subroutine calc_intorbang_pq(ip,iq,NL,NS,pg,c_p,c_q,intorbang_pq)
! - subroutine calc_intmoment_pq(ip,iq,NL,NS,pg,c_p,c_q,intmoment_pq)
! - subroutine calc_intmoment_wpos_pq(posC,ip,iq,NL,NS,pg,c_p,c_q,intmoment_pq)
! - subroutine calc_int_2nd_moment_pq(ip,iq,NL,NS,pg,c_p,c_q,int2moment_pq)
! - function intrhoA0M_mat(p,a,q,b)
!====================================================================================
!====================================================================================
subroutine calc_prop_int
!====================================================================================
  use Precision
  Use DiracOutput
  use Constants
  use IntegralStorage
  use NucBasis
  use prop_param
  use prop_tmp

  implicit none

  integer :: i,j,k,l,ii,icount
  integer :: nn,mm
  integer :: it
  real(kind=dp) :: time
!  real(kind=dp) :: x,y,z,dx,dy,dz
  integer :: filenum ! number of print file
  character(LEN=80) :: FMT
  complex(kind=dp),allocatable :: calE(:,:)  ! calE_NM
  complex(kind=dp),allocatable :: calE0(:,:)  ! calE0_NM
  character(LEN=80),allocatable :: filedat(:)

  integer :: nprint ! interval of print files
  character(LEN=80) :: filename_E  ! name of calE.dat
  complex(kind=dp) :: intspin(3), inttorq(3)
  complex(kind=dp) :: intspin0(3),inttorq0(3)
  complex(kind=dp) :: intdensmom(3), intOrbang(3), intOrbang_AM(3)
  complex(kind=dp) :: intdensmom0(3),intOrbang0(3),intOrbang_AM0(3)
  complex(kind=dp) :: intSpin_Qmat
  complex(kind=dp) :: inttorq_Qmat
  complex(kind=dp) :: intFj_Qmat
  complex(kind=dp) :: intjmoment_Qmat !function
  complex(kind=dp) :: intMagHyp_Qmat
  complex(kind=dp),allocatable :: tmp_intSpin_Qmat(:,:,:)
  complex(kind=dp),allocatable :: tmp_intFj_Qmat(:,:,:,:)
  complex(kind=dp),allocatable :: tmp_inttorq_Qmat(:,:,:)
!  complex(kind=dp),allocatable :: inttau_Arad_Qmat(:,:,:,:)
!  complex(kind=dp),allocatable :: inttorq_Arad_Qmat(:,:,:)
  complex(kind=dp),allocatable :: inttau_AM_Qmat(:,:,:,:)
  complex(kind=dp),allocatable :: inttorq_AM_Qmat(:,:,:)
  complex(kind=dp),allocatable :: tmp_intjmoment_Qmat(:,:,:,:)
  complex(kind=dp),allocatable :: tmp_intdensmom_Qmat(:,:,:)
  complex(kind=dp),allocatable :: tmp_intOrbang_Qmat(:,:,:)
  complex(kind=dp),allocatable :: tmp_intOrbang_AM_Qmat(:,:,:)

  real(kind=dp) :: posR(3)

  call set_paramini ! read param.ini, qed.inp
  call read_param_AM ! read param_AM.ini

  allocate(calE(4*NBS,4*NBS))
  allocate(calE0(4*NBS,4*NBS))
  allocate(tmp_intSpin_Qmat(3,4*NBS,4*NBS))
  allocate(tmp_intFj_Qmat(Nph,3,4*NBS,4*NBS))
!  allocate(inttau_Arad_Qmat(3,3,4*NBS,4*NBS))
  allocate(tmp_inttorq_Qmat(3,4*NBS,4*NBS))
!  allocate(inttorq_Arad_Qmat(3,4*NBS,4*NBS))
  allocate(tmp_intjmoment_Qmat(3,3,4*NBS,4*NBS))
  allocate(tmp_intdensmom_Qmat(3,4*NBS,4*NBS))
  allocate(tmp_intOrbang_Qmat(3,4*NBS,4*NBS))
  allocate(tmp_intOrbang_AM_Qmat(3,4*NBS,4*NBS))
  allocate(inttau_AM_Qmat(3,3,4*NBS,4*NBS))
  allocate(inttorq_AM_Qmat(3,4*NBS,4*NBS))
  
  do k=1,3
    !$omp parallel do
    do nn=1,4*NBS
      do mm=1,4*NBS
        tmp_intSpin_Qmat(k,nn,mm) = intSpin_Qmat(k,nn,mm)
        tmp_inttorq_Qmat(k,nn,mm) = inttorq_Qmat(k,nn,mm)
!        do j=1,Nph
!          tmp_intFj_Qmat(j,k,nn,mm) = intFj_Qmat(j,k,nn,mm)
!        end do
!!$      do l=1,NAT
!!$         posR(1)=xc(l); posR(2)=yc(l); posR(3)=zc(l)
!!$         write(*,*) intMagHyp_Qmat(posR,nn,mm)
!!$      end do
      call calc_intdensmom_Qmat(nn,mm,tmp_intdensmom_Qmat(:,nn,mm))
      call calc_intOrbang_Qmat(nn,mm,tmp_intOrbang_Qmat(:,nn,mm))
      call calc_intOrbang_AM_Qmat(nn,mm,tmp_intOrbang_AM_Qmat(:,nn,mm))
      end do
    end do
    !$omp end parallel do
  end do
!  stop

  do k=1,3
  do l=1,3
    !$omp parallel do
    do nn=1,4*NBS
      do mm=1,4*NBS
        tmp_intjmoment_Qmat(k,l,nn,mm) = intjmoment_Qmat(k,l,nn,mm)
      end do
    end do
    !$omp end parallel do
  end do
  end do

  call calc_inttau_AM(tmp_intjmoment_Qmat,inttau_AM_Qmat)
  call calc_inttorq_AM(inttau_AM_Qmat,inttorq_AM_Qmat)

  nprint = 1 ! interval of print files
  filenum = 10000  ! number of print file
  filename_E = '../calE.dat' ! name of calE.dat
  write(*,'(a)') filename_E
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

        intspin0(:) = (0.d0,0.d0)
        inttorq0(:) = (0.d0,0.d0)
  do icount = 0,filenum*nprint
     read(1200,FMT) time,calE 

  !----------------------------------------------------------
  ! set <phi|:calE:|phi> = <phi|calE|phi> - <0|calE(t=0)|0>
  !----------------------------------------------------------
  calE(:,:) = calE(:,:) - calE0(:,:)

     if((icount==0).or.(mod(icount,nprint)==0)) then
        it = nint(time/DeltaT)

    !------------------------------------------------------
    ! calc integral of spin angular momentum density
    !------------------------------------------------------
        intspin(:) = (0.d0,0.d0)
        do k=1,3
         do nn=1,4*NBS
           do mm=1,4*NBS
              intspin(k) = intspin(k) + tmp_intSpin_Qmat(k,nn,mm)*calE(nn,mm)
            end do
          end do
          if(icount==0) then
            intspin0(k) = intspin(k)
          end if
        end do
        write(75,'(13es24.14)') time, intspin(1), intspin(1)-intspin0(1) &
                                   &, intspin(2), intspin(2)-intspin0(2) &
                                   &, intspin(3), intspin(3)-intspin0(3) 

    !------------------------------------------------------
    ! calc integral of spin torque density
    !------------------------------------------------------

!        call calc_inttau_Arad(time,tmp_intFj_Qmat,inttau_Arad_Qmat)
!        call calc_inttorq_Arad(inttau_Arad_Qmat,inttorq_Arad_Qmat)

        inttorq(:) = (0.d0,0.d0)
        do k=1,3
         do nn=1,4*NBS
           do mm=1,4*NBS
!              inttorq(k) = inttorq(k) + (tmp_inttorq_Qmat(k,nn,mm) + inttorq_Arad_Qmat(k,nn,mm))*calE(nn,mm)
              inttorq(k) = inttorq(k) + (tmp_inttorq_Qmat(k,nn,mm) + inttorq_AM_Qmat(k,nn,mm))*calE(nn,mm)
!              inttorq(k) = inttorq(k) + (inttorq_AM_Qmat(k,nn,mm))*calE(nn,mm)
            end do
          end do
          if(icount==0) then
            inttorq0(k) = inttorq(k)
          end if
        end do
        write(76,'(13es24.14)') time, inttorq(1), inttorq(1)-inttorq0(1) &
                                   &, inttorq(2), inttorq(2)-inttorq0(2) &
                                   &, inttorq(3), inttorq(3)-inttorq0(3) 

    !------------------------------------------------------
    ! calc integral of momentum of electron density
    !------------------------------------------------------
        intdensmom(:) = (0.d0,0.d0)
        do k=1,3
         do nn=1,4*NBS
           do mm=1,4*NBS
              intdensmom(k) = intdensmom(k) + tmp_intdensmom_Qmat(k,nn,mm)*calE(nn,mm)
            end do
          end do
          if(icount==0) then
            intdensmom0(k) = intdensmom(k)
          end if
        end do
        write(77,'(13es24.14)') time, intdensmom(1), intdensmom(1)-intdensmom0(1) &
                                   &, intdensmom(2), intdensmom(2)-intdensmom0(2) &
                                   &, intdensmom(3), intdensmom(3)-intdensmom0(3) 

    !------------------------------------------------------
    ! calc integral of orbital angular momentum density
    !------------------------------------------------------
        intOrbang(:) = (0.d0,0.d0)
        do k=1,3
         do nn=1,4*NBS
           do mm=1,4*NBS
              intOrbang(k) = intOrbang(k) + tmp_intOrbang_Qmat(k,nn,mm)*calE(nn,mm)
            end do
          end do
          if(icount==0) then
            intOrbang0(k) = intOrbang(k)
          end if
        end do !k
        write(78,'(13es24.14)') time, intOrbang(1), intOrbang(1)-intOrbang0(1) &
                                   &, intOrbang(2), intOrbang(2)-intOrbang0(2) &
                                   &, intOrbang(3), intOrbang(3)-intOrbang0(3) 

    !------------------------------------------------------
    ! calc integral of AM part of the orbital angular momentum density
    !------------------------------------------------------
        intOrbang_AM(:) = (0.d0,0.d0)
        do k=1,3
         do nn=1,4*NBS
           do mm=1,4*NBS
              intOrbang_AM(k) = intOrbang_AM(k) + tmp_intOrbang_AM_Qmat(k,nn,mm)*calE(nn,mm)
            end do
          end do
          if(icount==0) then
            intOrbang_AM0(k) = intOrbang_AM(k)
          end if
        end do !k
        write(79,'(13es24.14)') time, intOrbang_AM(1), intOrbang_AM(1)-intOrbang_AM0(1) &
                                   &, intOrbang_AM(2), intOrbang_AM(2)-intOrbang_AM0(2) &
                                   &, intOrbang_AM(3), intOrbang_AM(3)-intOrbang_AM0(3) 



     end if !end if(mod(icount,nprint)==0)
  end do !icount

  write(*,*)'# fort.75 => intspin   : t, intspinx, intspinx-intspinx0, intspiny, intspiny-intspiny0, intspinz, intspinz-intspinz0'
  write(*,*)'# fort.76 => inttorq   : t, inttorqx, inttorqx-inttorqx0, inttorqy, inttorqy-inttorqy0, inttorqz, inttorqz-inttorqz0'
  write(*,*)'# fort.77 => intdensmom: t, intdensmomx, intdensmomx-intdensmomx0, intdensmomy, intdensmomy-intdensmomy0, intdensmomz, intdensmomz-intdensmomz0'
  write(*,*)'# fort.78 => intOrbang : t, intOrbangx, intOrbangx-intOrbangx0, intOrbangy, intOrbangy-intOrbangy0, intOrbangz, intOrbangz-intOrbangz0'
  write(*,*)'# fort.79 => intOrbang_AM : t, intOrbang_AMx, intOrbang_AMx-intOrbang_AMx0, intOrbang_AMy, intOrbang_AMy-intOrbang_AMy0, intOrbang_AMz, intOrbang_AMz-intOrbang_AMz0'

end subroutine calc_prop_int

!============================================================
function intSpin_mat(k,p,a,q,b)
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none
 
  complex(kind=dp) :: intSpin_mat
  integer,intent(in) :: k
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  complex(kind=dp) :: intSpin_pq_11,intSpin_pq_22,intSpin_pq_33,intSpin_pq_44
  complex(kind=dp) :: intSpin_pq_12,intSpin_pq_21,intSpin_pq_34,intSpin_pq_43
  
  call set_GammaMatrix
  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  
  intSpin_mat = (0._dp,0._dp)
  if(k==1) then
     call calc_intN_pq(1,2,NBS_L,NBS_S,pg,c_p,c_q,intSpin_pq_12) 
     call calc_intN_pq(2,1,NBS_L,NBS_S,pg,c_p,c_q,intSpin_pq_21) 
     call calc_intN_pq(3,4,NBS_L,NBS_S,pg,c_p,c_q,intSpin_pq_34) 
     call calc_intN_pq(4,3,NBS_L,NBS_S,pg,c_p,c_q,intSpin_pq_43) 
    
     intSpin_mat = (intSpin_pq_12 +intSpin_pq_21 +intSpin_pq_34 +intSpin_pq_43 )*0.5d0

  else if(k==2) then
     call calc_intN_pq(1,2,NBS_L,NBS_S,pg,c_p,c_q,intSpin_pq_12) 
     call calc_intN_pq(2,1,NBS_L,NBS_S,pg,c_p,c_q,intSpin_pq_21) 
     call calc_intN_pq(3,4,NBS_L,NBS_S,pg,c_p,c_q,intSpin_pq_34) 
     call calc_intN_pq(4,3,NBS_L,NBS_S,pg,c_p,c_q,intSpin_pq_43) 
    
     intSpin_mat = (-IU*intSpin_pq_12 +IU*intSpin_pq_21 -IU*intSpin_pq_34 +IU*intSpin_pq_43 )*0.5d0

  else if(k==3) then
     call calc_intN_pq(1,1,NBS_L,NBS_S,pg,c_p,c_q,intSpin_pq_11) 
     call calc_intN_pq(2,2,NBS_L,NBS_S,pg,c_p,c_q,intSpin_pq_22) 
     call calc_intN_pq(3,3,NBS_L,NBS_S,pg,c_p,c_q,intSpin_pq_33) 
     call calc_intN_pq(4,4,NBS_L,NBS_S,pg,c_p,c_q,intSpin_pq_44) 
    
     intSpin_mat = (intSpin_pq_11 -intSpin_pq_22 +intSpin_pq_33 -intSpin_pq_44 )*0.5d0
   end if

  return
end function intSpin_mat

!============================================================
subroutine calc_intdensmom_mat(p,a,q,b,intdensmom_mat)
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
! intdensmom_mat(i) = \int d^3r \psi_p^a x^i  \psi_q^b
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none
 
  complex(kind=dp),intent(out) :: intdensmom_mat(3)
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  integer :: k
  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  complex(kind=dp) :: intmom_pq_11(3),intmom_pq_22(3),intmom_pq_33(3),intmom_pq_44(3)
  
  call set_GammaMatrix
  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  
     call calc_intmoment_pq(1,1,NBS_L,NBS_S,pg,c_p,c_q,intmom_pq_11) 
     call calc_intmoment_pq(2,2,NBS_L,NBS_S,pg,c_p,c_q,intmom_pq_22) 
     call calc_intmoment_pq(3,3,NBS_L,NBS_S,pg,c_p,c_q,intmom_pq_33) 
     call calc_intmoment_pq(4,4,NBS_L,NBS_S,pg,c_p,c_q,intmom_pq_44) 
    
  intdensmom_mat(:) = (0._dp,0._dp)
     do k=1,3
        intdensmom_mat(k) = intmom_pq_11(k) +intmom_pq_22(k) +intmom_pq_33(k) +intmom_pq_44(k)
     end do

  return
end subroutine calc_intdensmom_mat

!============================================================
subroutine calc_intOrbang_mat(p,a,q,b,intOrbang_mat)
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
! intOrbang_mat = -i \hbar \int d^3r \psi_p^a x \times \nabla \psi_q^b
! The vector potential term is NOT included.
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none
 
  complex(kind=dp),intent(out) :: intOrbang_mat(3)
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  integer :: k
  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  complex(kind=dp) :: intOrbang_pq_11(3),intOrbang_pq_22(3),intOrbang_pq_33(3),intOrbang_pq_44(3)
  
  call set_GammaMatrix
  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  
     call calc_intorbang_pq(1,1,NBS_L,NBS_S,pg,c_p,c_q,intOrbang_pq_11) 
     call calc_intorbang_pq(2,2,NBS_L,NBS_S,pg,c_p,c_q,intOrbang_pq_22) 
     call calc_intorbang_pq(3,3,NBS_L,NBS_S,pg,c_p,c_q,intOrbang_pq_33) 
     call calc_intorbang_pq(4,4,NBS_L,NBS_S,pg,c_p,c_q,intOrbang_pq_44) 
    
  intOrbang_mat(:) = (0._dp,0._dp)
     do k=1,3
        intOrbang_mat(k) = intOrbang_pq_11(k) +intOrbang_pq_22(k) +intOrbang_pq_33(k) +intOrbang_pq_44(k)
     end do

  return
end subroutine calc_intOrbang_mat

!============================================================
subroutine calc_intOrbang_AM_mat(p,a,q,b,intOrbang_AM_mat)
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
! intOrbang_AM_mat = -(Zee/2c)* \int d^3r \psi_p^a [x^j x^j Bmag(i) - x^i x^j Bmag(j)] \psi_q^b
! AM = 1/2 * B \times r
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use param_AM
  implicit none
 
  complex(kind=dp),intent(out) :: intOrbang_AM_mat(3)
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  integer :: k,l
  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  complex(kind=dp) :: int2moment_pq_11(6),int2moment_pq_22(6),int2moment_pq_33(6),int2moment_pq_44(6)
  complex(kind=dp) :: tmp(6), tmpmat(3,3)
  
  call set_GammaMatrix
  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  
     call calc_int_2nd_moment_pq(1,1,NBS_L,NBS_S,pg,c_p,c_q,int2moment_pq_11) 
     call calc_int_2nd_moment_pq(2,2,NBS_L,NBS_S,pg,c_p,c_q,int2moment_pq_22) 
     call calc_int_2nd_moment_pq(3,3,NBS_L,NBS_S,pg,c_p,c_q,int2moment_pq_33) 
     call calc_int_2nd_moment_pq(4,4,NBS_L,NBS_S,pg,c_p,c_q,int2moment_pq_44) 

     do k=1,6
        tmp(k) = int2moment_pq_11(k) +int2moment_pq_22(k) +int2moment_pq_33(k) +int2moment_pq_44(k)
     end do
     tmpmat(1,1) = tmp(1)
     tmpmat(2,2) = tmp(2)
     tmpmat(3,3) = tmp(3)
     tmpmat(1,2) = tmp(4)
     tmpmat(2,1) = tmp(4)
     tmpmat(2,3) = tmp(5)
     tmpmat(3,2) = tmp(5)
     tmpmat(3,1) = tmp(6)
     tmpmat(1,3) = tmp(6)
    
  intOrbang_AM_mat(:) = (0._dp,0._dp)
     do k=1,3
        do l=1,3
          intOrbang_AM_mat(k) = tmpmat(l,l)*Bmag(k) - tmpmat(k,l)*Bmag(l)
        end do
        intOrbang_AM_mat(k) = - Ze*intOrbang_AM_mat(k)/CCC/2.d0
     end do

!     do k=1,3
!        write(*,*) k, intOrbang_AM_mat(k)
!     end do
!     stop

  return
end subroutine calc_intOrbang_AM_mat

!============================================================
function inttau_mat(k,l,p,a,q,b)
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants !IU, CCC
  implicit none
 
  complex(kind=dp) :: inttau_mat
  integer,intent(in) :: k,l
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  complex(kind=dp) :: inttau_pq_14(3),inttau_pq_23(3),inttau_pq_32(3),inttau_pq_41(3)
  complex(kind=dp) :: inttau_pq_13(3),inttau_pq_24(3),inttau_pq_31(3),inttau_pq_42(3)
  complex(kind=dp) :: inttau_mat1,inttau_mat2,inttau_mat3
  integer :: i
  
  call set_GammaMatrix
  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  
  call calc_intT_pq(1,4,NBS_L,NBS_S,pg,c_p,c_q,inttau_pq_14) 
  call calc_intT_pq(2,3,NBS_L,NBS_S,pg,c_p,c_q,inttau_pq_23) 
  call calc_intT_pq(3,2,NBS_L,NBS_S,pg,c_p,c_q,inttau_pq_32) 
  call calc_intT_pq(4,1,NBS_L,NBS_S,pg,c_p,c_q,inttau_pq_41) 
  call calc_intT_pq(1,3,NBS_L,NBS_S,pg,c_p,c_q,inttau_pq_13) 
  call calc_intT_pq(2,4,NBS_L,NBS_S,pg,c_p,c_q,inttau_pq_24) 
  call calc_intT_pq(3,1,NBS_L,NBS_S,pg,c_p,c_q,inttau_pq_31) 
  call calc_intT_pq(4,2,NBS_L,NBS_S,pg,c_p,c_q,inttau_pq_42) 

  if(l==1) then
    inttau_mat = inttau_pq_14(k)+inttau_pq_23(k)+inttau_pq_32(k)+inttau_pq_41(k)
  else if(l==2) then
    inttau_mat = -IU*(inttau_pq_14(k)+inttau_pq_32(k)) +IU*(inttau_pq_23(k)+inttau_pq_41(k))
  else if(l==3) then
    inttau_mat = (inttau_pq_13(k)+inttau_pq_31(k)) -(inttau_pq_24(k)+inttau_pq_42(k))
  end if

  inttau_mat = inttau_mat*IU*CCC

  return
end function inttau_mat

!============================================================
function intjmoment_mat(i,j,p,a,q,b)
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
! \psi^\dagger x^i \gamma^0 \gamma^j \psi
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants  ! use IU
  implicit none

  complex(kind=8) :: intjmoment_mat
  integer,intent(in) :: i,j
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  complex(kind=dp) :: intmom_pq_14(3),intmom_pq_23(3),intmom_pq_32(3),intmom_pq_41(3)
  complex(kind=dp) :: intmom_pq_13(3),intmom_pq_24(3),intmom_pq_31(3),intmom_pq_42(3)

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in intjmoment_mat."
     stop
  end if

  if(j.lt.1 .or. j.gt.3) then
     write(*,*) "j should be 1-3 in intjmoment_mat."
     stop
  end if
  
  call set_GammaMatrix
  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  
  intjmoment_mat = (0._dp,0._dp)
  if(j==1) then
     call calc_intmoment_pq(1,4,NBS_L,NBS_S,pg,c_p,c_q,intmom_pq_14) 
     call calc_intmoment_pq(2,3,NBS_L,NBS_S,pg,c_p,c_q,intmom_pq_23) 
     call calc_intmoment_pq(3,2,NBS_L,NBS_S,pg,c_p,c_q,intmom_pq_32) 
     call calc_intmoment_pq(4,1,NBS_L,NBS_S,pg,c_p,c_q,intmom_pq_41) 
    
     intjmoment_mat = intmom_pq_14(i) +intmom_pq_23(i) +intmom_pq_32(i) +intmom_pq_41(i) 

  else if(j==2) then
     call calc_intmoment_pq(1,4,NBS_L,NBS_S,pg,c_p,c_q,intmom_pq_14) 
     call calc_intmoment_pq(2,3,NBS_L,NBS_S,pg,c_p,c_q,intmom_pq_23) 
     call calc_intmoment_pq(3,2,NBS_L,NBS_S,pg,c_p,c_q,intmom_pq_32) 
     call calc_intmoment_pq(4,1,NBS_L,NBS_S,pg,c_p,c_q,intmom_pq_41) 
    
     intjmoment_mat = -IU*intmom_pq_14(i) +IU*intmom_pq_23(i) -IU*intmom_pq_32(i) +IU*intmom_pq_41(i)

  else if(j==3) then
     call calc_intmoment_pq(1,3,NBS_L,NBS_S,pg,c_p,c_q,intmom_pq_13) 
     call calc_intmoment_pq(2,4,NBS_L,NBS_S,pg,c_p,c_q,intmom_pq_24) 
     call calc_intmoment_pq(3,1,NBS_L,NBS_S,pg,c_p,c_q,intmom_pq_31) 
     call calc_intmoment_pq(4,2,NBS_L,NBS_S,pg,c_p,c_q,intmom_pq_42) 
    
     intjmoment_mat = intmom_pq_13(i) -intmom_pq_24(i) +intmom_pq_31(i) -intmom_pq_42(i)
   end if

  return
end function intjmoment_mat

!============================================================
function intMagHyp_mat(posR,p,a,q,b)
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
! Z_ee \mu^i epsilon_{ijk} \int d^3r \psi^\dagger \gamma^0 \gamma^j r^k/|\vec{r}|^3 \psi
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants  ! use IU
  use param_AM !magdip
  implicit none

  complex(kind=8) :: intMagHyp_mat
  integer,intent(in) :: p,q
  real(kind=dp),intent(in) :: posR(3)
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  complex(kind=dp) :: intE_pq_14(3),intE_pq_23(3),intE_pq_32(3),intE_pq_41(3)
  complex(kind=dp) :: intE_pq_13(3),intE_pq_24(3),intE_pq_31(3),intE_pq_42(3)
  complex(kind=dp) :: tmp_mat(3,3)
  integer :: k

!  magdip(1) = (0._dp, 0._dp)
!  magdip(2) = (0._dp, 0._dp)
!  magdip(3) = (0._dp, 0._dp)
  
  call set_GammaMatrix
  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  
  tmp_mat(:,:) = (0._dp,0._dp)
     call calc_intE_pq(1,4,posR,NBS_L,NBS_S,pg,c_p,c_q,intE_pq_14) 
     call calc_intE_pq(2,3,posR,NBS_L,NBS_S,pg,c_p,c_q,intE_pq_23) 
     call calc_intE_pq(3,2,posR,NBS_L,NBS_S,pg,c_p,c_q,intE_pq_32) 
     call calc_intE_pq(4,1,posR,NBS_L,NBS_S,pg,c_p,c_q,intE_pq_41) 
     call calc_intE_pq(1,3,posR,NBS_L,NBS_S,pg,c_p,c_q,intE_pq_13) 
     call calc_intE_pq(2,4,posR,NBS_L,NBS_S,pg,c_p,c_q,intE_pq_24) 
     call calc_intE_pq(3,1,posR,NBS_L,NBS_S,pg,c_p,c_q,intE_pq_31) 
     call calc_intE_pq(4,2,posR,NBS_L,NBS_S,pg,c_p,c_q,intE_pq_42) 
    
     do k=1,3
        tmp_mat(1,k) = intE_pq_14(k) +intE_pq_23(k) +intE_pq_32(k) +intE_pq_41(k) 
        tmp_mat(2,k) = -IU*intE_pq_14(k) +IU*intE_pq_23(k) -IU*intE_pq_32(k) +IU*intE_pq_41(k)
        tmp_mat(3,k) = intE_pq_13(k) -intE_pq_24(k) +intE_pq_31(k) -intE_pq_42(k)
     end do

     intMagHyp_mat = magdip(1)*( tmp_mat(2,3)-tmp_mat(3,2) )&
                   &+magdip(2)*( tmp_mat(3,1)-tmp_mat(1,3) )&
                   &+magdip(3)*( tmp_mat(1,2)-tmp_mat(2,1) )

     intMagHyp_mat = Ze*intMagHyp_mat

  return
end function intMagHyp_mat

!============================================================
subroutine calc_inttau_AM(tmp_intjmom_Qmat,inttau_AM_Qmat)
! t : time = timestep*DeltaT
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
! AM = 1/2 * B \times r
! tau^{kl} = j^l AM^k /CCC
!============================================================

  use Precision
  use DefineTypes
  use DiracOutput
  use Constants  ! use IU
  use param_AM
  implicit none

  integer :: k,l
  integer :: nn,mm
  integer :: j
  character(LEN=1) :: a,b
  complex(kind=dp),intent(out) :: inttau_AM_Qmat(3,3,4*NBS,4*NBS)
  complex(kind=dp),intent(in) :: tmp_intjmom_Qmat(3,3,4*NBS,4*NBS)

  inttau_AM_Qmat(:,:,:,:) = (0.d0,0.d0)

  do l=1,3
    do nn=1,4*NBS
      do mm=1,4*NBS
        inttau_AM_Qmat(1,l,nn,mm) = Bmag(2)*tmp_intjmom_Qmat(3,l,nn,mm) - Bmag(3)*tmp_intjmom_Qmat(2,l,nn,mm)
        inttau_AM_Qmat(2,l,nn,mm) = Bmag(3)*tmp_intjmom_Qmat(1,l,nn,mm) - Bmag(1)*tmp_intjmom_Qmat(3,l,nn,mm)
        inttau_AM_Qmat(3,l,nn,mm) = Bmag(1)*tmp_intjmom_Qmat(2,l,nn,mm) - Bmag(2)*tmp_intjmom_Qmat(1,l,nn,mm)
      end do
    end do
  end do

  inttau_AM_Qmat(:,:,:,:) = inttau_AM_Qmat(:,:,:,:)*0.5_dp*Ze ! intjmom_Qmat doesn't include Ze, CCC

end subroutine calc_inttau_AM

!============================================================
subroutine calc_inttorq_AM(inttau_AM_Qmat,inttorq_AM_Qmat)
! t : time = timestep*DeltaT
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================

  use Precision
  use DefineTypes
  use DiracOutput
  implicit none

  integer :: nn,mm
  character(LEN=1) :: a,b
  complex(kind=dp),intent(in) :: inttau_AM_Qmat(3,3,4*NBS,4*NBS)
  complex(kind=dp),intent(out) :: inttorq_AM_Qmat(3,4*NBS,4*NBS)

  do nn=1,4*NBS
    do mm=1,4*NBS
      inttorq_AM_Qmat(1,nn,mm) = -inttau_AM_Qmat(2,3,nn,mm) + inttau_AM_Qmat(3,2,nn,mm)
      inttorq_AM_Qmat(2,nn,mm) =  inttau_AM_Qmat(1,3,nn,mm) - inttau_AM_Qmat(3,1,nn,mm)
      inttorq_AM_Qmat(3,nn,mm) = -inttau_AM_Qmat(1,2,nn,mm) + inttau_AM_Qmat(2,1,nn,mm)
    end do
  end do

end subroutine calc_inttorq_AM

!============================================================
subroutine calc_inttau_Arad(t,tmp_intFj_Qmat,inttau_Arad_Qmat)
! t : time = timestep*DeltaT
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================

  use Precision
  use DefineTypes
  use DiracOutput
  use Constants  ! use P0MAX,NP0,NPTH,NPPHI,Nph,ALPHA_COH
  implicit none

  real(kind=dp),intent(in) :: t  ! time
  integer :: k,l
  real(kind=dp) :: DeltaPj,P0j
  complex(kind=dp) :: vecPj(3) ! (circular) polarization vector (depends on sig)
  complex(kind=dp) :: vecpolj(3) ! (circular) polarization vector (depends on sig)
  integer :: nn,mm
  integer :: j
  character(LEN=1) :: a,b
  complex(kind=dp) :: inttau_Arad_Qmat(3,3,4*NBS,4*NBS)
  complex(kind=dp),intent(in) :: tmp_intFj_Qmat(Nph,3,4*NBS,4*NBS)

  inttau_Arad_Qmat = (0.d0,0.d0)

  do k=1,3
  do l=1,3
  do nn=1,4*NBS
    do mm=1,4*NBS
      do j=1,Nph
!        call calc_DeltaPj_and_P0j(j,DeltaPj,P0j)
        call calc_mode_pj(j,P0j,vecPj,vecpolj)
!        inttau_Arad_Qmat(k,l,nn,mm) = inttau_Arad_Qmat(k,l,nn,mm) &
!                                    &+ DeltaPj*tmp_intFj_Qmat(j,l,nn,mm) *cdexp(-IU*CCC*P0j*t)*vecpolj(k) *ALPHA_COH(j)
!        inttau_Arad_Qmat(k,l,nn,mm) = inttau_Arad_Qmat(k,l,nn,mm) &
!                                    &+ DeltaPj*dconjg(tmp_intFj_Qmat(j,l,mm,nn)) *cdexp(IU*CCC*P0j*t)*dconjg(vecpolj(k) *ALPHA_COH(j))
        inttau_Arad_Qmat(k,l,nn,mm) = inttau_Arad_Qmat(k,l,nn,mm) &
                                    &+ tmp_intFj_Qmat(j,l,nn,mm) *cdexp(-IU*CCC*P0j*t)*vecpolj(k) *ALPHA_COH(j)/(2.d0*PI*dsqrt(P0j*CCC))
        inttau_Arad_Qmat(k,l,nn,mm) = inttau_Arad_Qmat(k,l,nn,mm) &
                                    &+ dconjg(tmp_intFj_Qmat(j,l,mm,nn)) *cdexp(IU*CCC*P0j*t)*dconjg(vecpolj(k) *ALPHA_COH(j))/(2.d0*PI*dsqrt(P0j*CCC))
      end do
    end do
  end do
  end do
  end do

end subroutine calc_inttau_Arad

!============================================================
subroutine calc_inttorq_Arad(inttau_Arad_Qmat,inttorq_Arad_Qmat)
! t : time = timestep*DeltaT
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================

  use Precision
  use DefineTypes
  use DiracOutput
  use Constants  ! use P0MAX,NP0,NPTH,NPPHI,Nph,ALPHA_COH
  implicit none

  integer :: nn,mm
  character(LEN=1) :: a,b
  complex(kind=dp),intent(in) :: inttau_Arad_Qmat(3,3,4*NBS,4*NBS)
  complex(kind=dp),intent(out) :: inttorq_Arad_Qmat(3,4*NBS,4*NBS)

  do nn=1,4*NBS
    do mm=1,4*NBS
      inttorq_Arad_Qmat(1,nn,mm) = -inttau_Arad_Qmat(2,3,nn,mm) + inttau_Arad_Qmat(3,2,nn,mm)
      inttorq_Arad_Qmat(2,nn,mm) =  inttau_Arad_Qmat(1,3,nn,mm) - inttau_Arad_Qmat(3,1,nn,mm)
      inttorq_Arad_Qmat(3,nn,mm) = -inttau_Arad_Qmat(1,2,nn,mm) + inttau_Arad_Qmat(2,1,nn,mm)
    end do
  end do

end subroutine calc_inttorq_Arad

!============================================================
function intHrel_mat(p,a,q,b)
! <psi(i)|H'|psi(j)> = \int dr^3 psi_i^dagger Ze e AM_l \gamma^0 \gamma^l psi_j
!                    = \int dr^3 psi_i^dagger Ze e (-1)*AM^l \gamma^0 \gamma^l psi_j
!                    = \int dr^3 (-1/c) je_{ij}^l AM^l
! AM = (-1/2*B*y,1/2*B*x,0)
! B = rot(AM) = (0,0,B)
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none
 
  complex(kind=dp) :: intHrel_mat
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=dp) :: intHrel_pq_14(3),intHrel_pq_23(3),intHrel_pq_32(3),intHrel_pq_41(3)
  complex(kind=dp) :: intHrel_pq_13(3),intHrel_pq_24(3),intHrel_pq_31(3),intHrel_pq_42(3)

  complex(kind=dp) :: func_AMvec
  
  call set_GammaMatrix
  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  
  call calc_intHrel_pq(1,4,NBS_L,NBS_S,pg,c_p,c_q,intHrel_pq_14) 
  call calc_intHrel_pq(2,3,NBS_L,NBS_S,pg,c_p,c_q,intHrel_pq_23) 
  call calc_intHrel_pq(3,2,NBS_L,NBS_S,pg,c_p,c_q,intHrel_pq_32) 
  call calc_intHrel_pq(4,1,NBS_L,NBS_S,pg,c_p,c_q,intHrel_pq_41) 
  call calc_intHrel_pq(1,3,NBS_L,NBS_S,pg,c_p,c_q,intHrel_pq_13) 
  call calc_intHrel_pq(2,4,NBS_L,NBS_S,pg,c_p,c_q,intHrel_pq_24) 
  call calc_intHrel_pq(3,1,NBS_L,NBS_S,pg,c_p,c_q,intHrel_pq_31) 
  call calc_intHrel_pq(4,2,NBS_L,NBS_S,pg,c_p,c_q,intHrel_pq_42) 
  
  intHrel_mat = (0._dp,0._dp)

  ! You have to change func_AMvec in sub_AM.f90 too.
!    !-------------------------------
!    ! BM = rot(AM) = (Bmag,0,0)
!    ! AM = (0,1/2*Bmag*z,-1/2*Bmag*y)
!    intHrel_mat = intHrel_mat -func_AMvec(0.d0,2,0.d0,0.d0,1.d0)*(-IU*(intHrel_pq_14(3)+intHrel_pq_32(3)) +IU*(intHrel_pq_23(3)+intHrel_pq_41(3)))
!    intHrel_mat = intHrel_mat -func_AMvec(0.d0,3,0.d0,1.d0,0.d0)*((intHrel_pq_13(2)+intHrel_pq_31(2)) -(intHrel_pq_24(2)+intHrel_pq_42(2)))
!    intHrel_mat = intHrel_mat *Ze
!    !-------------------------------
!
!    !-------------------------------
!    ! BM = rot(AM) = (0,Bmag,0)
!    ! AM = (1/2*Bmag*z,0,-1/2*Bmag*x)
!    intHrel_mat = intHrel_mat -func_AMvec(0.d0,1,0.d0,0.d0,1.d0)*(intHrel_pq_14(3)+intHrel_pq_23(3)+intHrel_pq_32(3)+intHrel_pq_41(3))
!    intHrel_mat = intHrel_mat -func_AMvec(0.d0,3,1.d0,0.d0,0.d0)*((intHrel_pq_13(1)+intHrel_pq_31(1)) -(intHrel_pq_24(1)+intHrel_pq_42(1)))
!    intHrel_mat = intHrel_mat *Ze
!    !-------------------------------
!
!    !-------------------------------
!    ! BM = rot(AM) = (0,0,Bmag)
!    ! AM = (-1/2*B*y,1/2*B*x,0)
!    intHrel_mat = intHrel_mat -func_AMvec(0.d0,1,0.d0,1.d0,0.d0)*(intHrel_pq_14(2)+intHrel_pq_23(2)+intHrel_pq_32(2)+intHrel_pq_41(2))
!    intHrel_mat = intHrel_mat -func_AMvec(0.d0,2,1.d0,0.d0,0.d0)*(-IU*(intHrel_pq_14(1)+intHrel_pq_32(1)) +IU*(intHrel_pq_23(1)+intHrel_pq_41(1)))
!    intHrel_mat = intHrel_mat *Ze
!    !-------------------------------


    !-------------------------------
    ! BM = rot(AM) = (Bmag(1),Bmag(2),Bmag(3))
    ! AM = 0.5d0*\vec{BM} \times \vec{r}
    !    = (1/2*Bmagy*z -1/2*Bmagz*y, 1/2*Bmagz*x -1/2*Bmagx*z, 1/2*Bmagx*y -1/2*Bmagy*x)
    intHrel_mat = intHrel_mat -func_AMvec(0.d0,1,0.d0,0.d0,1.d0)*(intHrel_pq_14(3)+intHrel_pq_23(3)+intHrel_pq_32(3)+intHrel_pq_41(3))
    intHrel_mat = intHrel_mat -func_AMvec(0.d0,1,0.d0,1.d0,0.d0)*(intHrel_pq_14(2)+intHrel_pq_23(2)+intHrel_pq_32(2)+intHrel_pq_41(2))
    intHrel_mat = intHrel_mat -func_AMvec(0.d0,2,1.d0,0.d0,0.d0)*(-IU*(intHrel_pq_14(1)+intHrel_pq_32(1)) +IU*(intHrel_pq_23(1)+intHrel_pq_41(1)))
    intHrel_mat = intHrel_mat -func_AMvec(0.d0,2,0.d0,0.d0,1.d0)*(-IU*(intHrel_pq_14(3)+intHrel_pq_32(3)) +IU*(intHrel_pq_23(3)+intHrel_pq_41(3)))
    intHrel_mat = intHrel_mat -func_AMvec(0.d0,3,0.d0,1.d0,0.d0)*((intHrel_pq_13(2)+intHrel_pq_31(2)) -(intHrel_pq_24(2)+intHrel_pq_42(2)))
    intHrel_mat = intHrel_mat -func_AMvec(0.d0,3,1.d0,0.d0,0.d0)*((intHrel_pq_13(1)+intHrel_pq_31(1)) -(intHrel_pq_24(1)+intHrel_pq_42(1)))
    intHrel_mat = intHrel_mat *Ze
    !-------------------------------

  return
end function intHrel_mat

!=========================================================================
subroutine calc_intHrel_pq(ip,iq,NL,NS,pg,c_p,c_q,intHrel_pq)
! int (psi^+)_ip (psi)_iq 
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none

  integer,intent(in) :: ip,iq ! spinor indice
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=dp),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  
  complex(kind=dp),intent(out) :: intHrel_pq(3)
  
  integer :: NLS
  integer :: i,j,k
  real(kind=dp) :: overlap
  real(kind=dp) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=dp) :: alphaA, alphaB ! exponents
  integer :: nx,nbarx,ny,nbary,nz,nbarz
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg

  real(kind=8) :: posC(3) 
  real(kind=8) :: moment(3)

  posC(:) = 0._dp
  intHrel_pq(:) = (0._dp,0._dp)

  do i=1,NLS(ip,NL,NS)
     do j=1,NLS(iq,NL,NS)
        call set_pg(ip,i,NL,NS,pg,posA,alphaA,nx,ny,nz)
        call set_pg(iq,j,NL,NS,pg,posB,alphaB,nbarx,nbary,nbarz)
        call gauss_int_moment(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,posC,moment)
        f_norm_pg = norm_pg(alphaA,nx,ny,nz)*norm_pg(alphaB,nbarx,nbary,nbarz)
        do k=1,3
          intHrel_pq(k) = intHrel_pq(k) + conjg(c_p(ip,i))*c_q(iq,j)*moment(k)*f_norm_pg
        end do
     end do
  end do
  return
end subroutine calc_intHrel_pq

!=========================================================================
subroutine calc_intorbang_pq(ip,iq,NL,NS,pg,c_p,c_q,intorbang_pq)
! -i \hbar \int d^3r \psi_p^a x \times \nabla \psi_q^b
!=========================================================================  
  use Precision
  use DefineTypes
  use Constants
  implicit none

  integer,intent(in) :: ip,iq ! spinor indice
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=dp),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  
  complex(kind=dp),intent(out) :: intorbang_pq(3)
  
  integer :: NLS
  integer :: i,j,k
  real(kind=dp) :: overlap
  real(kind=dp) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=dp) :: alphaA, alphaB ! exponents
  integer :: nx,nbarx,ny,nbary,nz,nbarz
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg

  real(kind=8) :: posC(3) 
  real(kind=8) :: xigradj(3,3)

  posC(:) = 0._dp
  intorbang_pq(:) = (0._dp,0._dp)

  do i=1,NLS(ip,NL,NS)
     do j=1,NLS(iq,NL,NS)
        call set_pg(ip,i,NL,NS,pg,posA,alphaA,nx,ny,nz)
        call set_pg(iq,j,NL,NS,pg,posB,alphaB,nbarx,nbary,nbarz)
        call gauss_int_xigradj(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,xigradj)
        f_norm_pg = norm_pg(alphaA,nx,ny,nz)*norm_pg(alphaB,nbarx,nbary,nbarz)
        do k=1,3
           if(k.eq.1) then
              intorbang_pq(k) = intorbang_pq(k) + conjg(c_p(ip,i))*c_q(iq,j)*f_norm_pg*(xigradj(2,3)-xigradj(3,2))
           elseif(k.eq.2) then
              intorbang_pq(k) = intorbang_pq(k) + conjg(c_p(ip,i))*c_q(iq,j)*f_norm_pg*(xigradj(3,1)-xigradj(1,3))
           elseif(k.eq.3) then
              intorbang_pq(k) = intorbang_pq(k) + conjg(c_p(ip,i))*c_q(iq,j)*f_norm_pg*(xigradj(1,2)-xigradj(2,1))
           end if
        end do
     end do
  end do

  do k=1,3
    intorbang_pq(k) = - intorbang_pq(k) *IU
  end do
  return
end subroutine calc_intorbang_pq

!=========================================================================
subroutine calc_intmoment_pq(ip,iq,NL,NS,pg,c_p,c_q,intmoment_pq)
! int (psi^+)_ip vec{r} (psi)_iq 
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none

  integer,intent(in) :: ip,iq ! spinor indice
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=dp),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  
  complex(kind=dp),intent(out) :: intmoment_pq(3)
  
  integer :: NLS
  integer :: i,j,k
  real(kind=dp) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=dp) :: alphaA, alphaB ! exponents
  integer :: nx,nbarx,ny,nbary,nz,nbarz
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg

  real(kind=8) :: posC(3) 
  real(kind=8) :: moment(3)

  intmoment_pq(:) = (0._dp,0._dp)
  posC(:) = 0._dp

  do i=1,NLS(ip,NL,NS)
     do j=1,NLS(iq,NL,NS)
        call set_pg(ip,i,NL,NS,pg,posA,alphaA,nx,ny,nz)
        call set_pg(iq,j,NL,NS,pg,posB,alphaB,nbarx,nbary,nbarz)
        call gauss_int_moment(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,posC,moment)
        f_norm_pg = norm_pg(alphaA,nx,ny,nz)*norm_pg(alphaB,nbarx,nbary,nbarz)
        do k=1,3
          intmoment_pq(k) = intmoment_pq(k) + conjg(c_p(ip,i))*c_q(iq,j)*moment(k)*f_norm_pg
        end do
     end do
  end do
  return
end subroutine calc_intmoment_pq

!=========================================================================
subroutine calc_intmoment_wpos_pq(posC,ip,iq,NL,NS,pg,c_p,c_q,intmoment_pq)
! int (psi^+)_ip vec{r} (psi)_iq 
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none

  integer,intent(in) :: ip,iq ! spinor indice
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=dp),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  real(kind=8),intent(in) :: posC(3) 
  
  complex(kind=dp),intent(out) :: intmoment_pq(3)
  
  integer :: NLS
  integer :: i,j,k
  real(kind=dp) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=dp) :: alphaA, alphaB ! exponents
  integer :: nx,nbarx,ny,nbary,nz,nbarz
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg

  real(kind=8) :: moment(3)

  intmoment_pq(:) = (0._dp,0._dp)
!  posC(:) = 0._dp

  do i=1,NLS(ip,NL,NS)
     do j=1,NLS(iq,NL,NS)
        call set_pg(ip,i,NL,NS,pg,posA,alphaA,nx,ny,nz)
        call set_pg(iq,j,NL,NS,pg,posB,alphaB,nbarx,nbary,nbarz)
        call gauss_int_moment(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,posC,moment)
        f_norm_pg = norm_pg(alphaA,nx,ny,nz)*norm_pg(alphaB,nbarx,nbary,nbarz)
        do k=1,3
          intmoment_pq(k) = intmoment_pq(k) + conjg(c_p(ip,i))*c_q(iq,j)*moment(k)*f_norm_pg
        end do
     end do
  end do
  return
end subroutine calc_intmoment_wpos_pq

!=========================================================================
subroutine calc_int_2nd_moment_pq(ip,iq,NL,NS,pg,c_p,c_q,int2moment_pq)
! int (psi^+)_ip x^i x^j (psi)_iq 
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none

  integer,intent(in) :: ip,iq ! spinor indice
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=dp),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  
  complex(kind=dp),intent(out) :: int2moment_pq(6)
  
  integer :: NLS
  integer :: i,j,k
  real(kind=dp) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=dp) :: alphaA, alphaB ! exponents
  integer :: nx,nbarx,ny,nbary,nz,nbarz
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg

  real(kind=8) :: posC(3) 
  real(kind=8) :: moment(6)

  int2moment_pq(:) = (0._dp,0._dp)
  posC(:) = 0._dp

  do i=1,NLS(ip,NL,NS)
     do j=1,NLS(iq,NL,NS)
        call set_pg(ip,i,NL,NS,pg,posA,alphaA,nx,ny,nz)
        call set_pg(iq,j,NL,NS,pg,posB,alphaB,nbarx,nbary,nbarz)
        call gauss_int_2nd_moment(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,posC,moment)
        f_norm_pg = norm_pg(alphaA,nx,ny,nz)*norm_pg(alphaB,nbarx,nbary,nbarz)
        do k=1,6
          int2moment_pq(k) = int2moment_pq(k) + conjg(c_p(ip,i))*c_q(iq,j)*moment(k)*f_norm_pg
        end do
     end do
  end do
  return
end subroutine calc_int_2nd_moment_pq

!!$!=========================================================================
!!$subroutine calc_A_pq(ip,iq,posRA,NL,NS,pg,c_p,c_q,AA_pq,intE_pq)
!!$! Int (psi^+)_ip (psi)_iq /|r-R|
!!$!=========================================================================  
!!$  use Precision
!!$  use DefineTypes
!!$  implicit none
!!$
!!$  integer,intent(in) :: ip,iq ! spinor indice
!!$  real(kind=dp),intent(in) :: posRA(3)  ! position of nucleus
!!$  integer,intent(in) :: NL,NS
!!$  type(primitive_gaussian),intent(in) :: pg
!!$  complex(kind=dp),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
!!$  
!!$  complex(kind=dp),intent(out) :: AA_pq
!!$  complex(kind=dp),intent(out) :: intE_pq(3)
!!$  
!!$  integer :: NLS
!!$  integer :: i,j,k
!!$  real(kind=dp) :: overlap,ef(3) ! --> dummy variables here 
!!$  real(kind=dp) :: nucatt
!!$  real(kind=dp) :: posA(3), posB(3) ! position of center of gaussian 
!!$  real(kind=dp) :: alphaA, alphaB ! exponents
!!$  integer :: nx,nbarx,ny,nbary,nz,nbarz
!!$  real(kind=dp) :: norm_pg  ! function
!!$  real(kind=dp) :: f_norm_pg
!!$
!!$
!!$  AA_pq = (0._dp,0._dp)
!!$  intE_pq(1:3) = (0._dp,0._dp)
!!$
!!$  do i=1,NLS(ip,NL,NS)
!!$     do j=1,NLS(iq,NL,NS)
!!$        call set_pg(ip,i,NL,NS,pg,posA,alphaA,nx,ny,nz)
!!$        call set_pg(iq,j,NL,NS,pg,posB,alphaB,nbarx,nbary,nbarz)
!!$        call gauss_int(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,posRA,overlap,nucatt,ef)
!!$        f_norm_pg = norm_pg(alphaA,nx,ny,nz)*norm_pg(alphaB,nbarx,nbary,nbarz)
!!$        AA_pq = AA_pq + conjg(c_p(ip,i))*c_q(iq,j)*nucatt*f_norm_pg
!!$        do k=1,3
!!$           intE_pq(k) = intE_pq(k) + conjg(c_p(ip,i))*c_q(iq,j)*ef(k)*f_norm_pg
!!$        end do
!!$     end do
!!$  end do
!!$  return
!!$end subroutine calc_AA_pq
!!$

!============================================================
function intrhoA0M_mat(p,a,q,b)
! intrhoA0M = \int rho A0M d^3r
! A0M = \vec{E}_M \cdot \vec{r}
! \vec{E}_M = - grad(A0M)
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use param_AM !Eelec
  implicit none
 
  complex(kind=dp) :: intrhoA0M_mat
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=dp) :: intmom_pq_11(3),intmom_pq_22(3),intmom_pq_33(3),intmom_pq_44(3)
  integer :: k

  
  call set_GammaMatrix
  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  
  call calc_intmoment_pq(1,1,NBS_L,NBS_S,pg,c_p,c_q,intmom_pq_11) 
  call calc_intmoment_pq(2,2,NBS_L,NBS_S,pg,c_p,c_q,intmom_pq_22) 
  call calc_intmoment_pq(3,3,NBS_L,NBS_S,pg,c_p,c_q,intmom_pq_33) 
  call calc_intmoment_pq(4,4,NBS_L,NBS_S,pg,c_p,c_q,intmom_pq_44) 

  intrhoA0M_mat = (0._dp, 0._dp)
  do k=1,3
     intrhoA0M_mat = intrhoA0M_mat + Eelec(k)*(intmom_pq_11(k)+intmom_pq_22(k)+intmom_pq_33(k)+intmom_pq_44(k))
  end do
  
  return
end function intrhoA0M_mat
