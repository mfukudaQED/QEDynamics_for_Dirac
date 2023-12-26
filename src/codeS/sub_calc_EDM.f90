! Last Change:28-Oct-2015.03
!============================================================
! 2015.07.06
!============================================================
! SUBROUTINE
!============================================================
SUBROUTINE SET_PG_NORMALIZATION_CONSTANT
!============================================================
  use Precision
  use Constants
  use DiracOutput
  implicit none

  real(kind=dp) :: c1(0:7)
  integer :: i

  allocate(c_Lnorm(NBS_L),c_Snorm(NBS_S))

  c1(0) = sqrt(2._dp/PI)
  do i = 1,7
    c1(i) = c1(i-1)/(2*i-1)
  end do

  do i = 1,NBS_L
    c_Lnorm(i) = sqrt(sqrt(aa_L(i)**3) *(4._dp *aa_L(i))**(nx_L(i) +ny_L(i) +nz_L(i))&
                   &*c1(nx_L(i)) *c1(ny_L(i)) *c1(nz_L(i)))
  end do
  do i = 1,NBS_S
    c_Snorm(i) = sqrt(sqrt(aa_S(i)**3) *(4._dp *aa_S(i))**(nx_S(i) +ny_S(i) +nz_S(i))&
                   &*c1(nx_S(i)) *c1(ny_S(i)) *c1(nz_S(i)))
  end do

  return
END SUBROUTINE SET_PG_NORMALIZATION_CONSTANT

!============================================================
SUBROUTINE CALC_COEF_XX(coef,coefLL,coefSS,coefLS)
!============================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=dp),intent(in) :: coef(NBS0,4)
  complex(kind=dp),intent(out) :: coefLL(NBS_LL,1:2,1:2)
  complex(kind=dp),intent(out) :: coefSS(NBS_SS,3:4,3:4)
  complex(kind=dp),intent(out) :: coefLS(NBS_LS,1:2,3:4)

  integer :: i,j,a,b,m,n

  !-- LL -----------------------------------------------------------------
  do b = 1,2
    do a = 1,2

      do m = 1,NBS_L
        coefLL(m,a,b) = conjg(coef(m,a)) *coef(m,b)
      end do
      m = NBS_L + 1
      n = NBS_LT + 1
      do j = 2,NBS_L
        do i = 1,j-1
          coefLL(m,a,b) = conjg(coef(i,a)) *coef(j,b)
          m = m + 1
        end do
        do i = 1,j-1
          coefLL(n,a,b) = conjg(coef(j,a)) *coef(i,b)
          n = n + 1
        end do
      end do

    end do
  end do

  !-- SS -----------------------------------------------------------------
  do b = 3,4
    do a = 3,4

      do m = 1,NBS_S
        coefSS(m,a,b) = conjg(coef(m,a)) *coef(m,b)
      end do
      m = NBS_S + 1
      n = NBS_ST + 1
      do j = 2,NBS_S
        do i = 1,j-1
          coefSS(m,a,b) = conjg(coef(i,a)) *coef(j,b)
          m = m + 1
        end do
        do i = 1,j-1
          coefSS(n,a,b) = conjg(coef(j,a)) *coef(i,b)
          n = n + 1
        end do
      end do

    end do
  end do

  !-- LS -----------------------------------------------------------------
  do b = 3,4
    do a = 1,2

      m = 1
      do j = 1,NBS_S
        do i = 1,NBS_L
          coefLS(m,a,b) = conjg(coef(i,a)) *coef(j,b)
          m = m + 1
        end do
      end do

    end do
  end do

  return
END SUBROUTINE CALC_COEF_XX

!=======================================================
SUBROUTINE COPY_DIRACOUTPUT_CP_MOD(p,a,c_p)
! copy coefficients according to the value of p and a
! 2010.12.2 modified to include Kramers pair
!    If negative values are in p or q, they denote Kramers pair of -p or -q MO.
!    We do similar transformation for positron coefficients (just formally).
! 2012.5.19 modified from copy_DiracOutput
!    a is specified by 1(+,electron) or 2(-,positron)
!    to be used in modified twoele integration routines.
!=======================================================
  use DefineTypes
  use DiracOutput
  implicit none

  integer,intent(in) :: p,a
  complex(kind=8),intent(out) :: c_p(NBS0,4)

  complex(kind=8) :: cLa_p(NBS_L),cLb_p(NBS_L),cSa_p(NBS_S),cSb_p(NBS_S)
  complex(kind=8) :: dLa_p(NBS_L),dLb_p(NBS_L),dSa_p(NBS_S),dSb_p(NBS_S)
  integer :: j

  ! c_La,c_Lb,c_Sa,c_Sb for electron, d_La,d_Lb,d_Sa,d_Sb for positron
  ! La->1, Lb->2, Sa->3, Sb->4, no matter components are really large or small

  ! Kramers pair (MO with bar) is obtained by transforming
  ! bar(c_La) = -(c_Lb)*
  ! bar(c_Lb) =  (c_La)*
  ! bar(c_Sa) = -(c_Sb)*
  ! bar(c_Sb) =  (c_Sa)*

  if(p.gt.0) then
     do j=1,NBS_L
        cLa_p(j) = c_La(j,p); cLb_p(j) = c_Lb(j,p)
        dLa_p(j) = d_La(j,p); dLb_p(j) = d_Lb(j,p)
     end do
     do j=1,NBS_S
        cSa_p(j) = c_Sa(j,p); cSb_p(j) = c_Sb(j,p)
        dSa_p(j) = d_Sa(j,p); dSb_p(j) = d_Sb(j,p)
     end do
  elseif(p.lt.0) then ! Kramers pair
     do j=1,NBS_L
        cLa_p(j) = -conjg(c_Lb(j,-p)); cLb_p(j) = conjg(c_La(j,-p))
        dLa_p(j) = -conjg(d_Lb(j,-p)); dLb_p(j) = conjg(d_La(j,-p))
     end do
     do j=1,NBS_S
        cSa_p(j) = -conjg(c_Sb(j,-p)); cSb_p(j) = conjg(c_Sa(j,-p))
        dSa_p(j) = -conjg(d_Sb(j,-p)); dSb_p(j) = conjg(d_Sa(j,-p))
     end do
  elseif(p.eq.0) then
     c_p = (0.d0,0.d0)
     !stop
     return
  end if

  c_p=(0.d0,0.d0)

  if(a.eq.1) then
     do j=1,NBS_L
        c_p(j,1) = cLa_p(j); c_p(j,2) = cLb_p(j)
     end do
     do j=1,NBS_S
        c_p(j,3) = cSa_p(j); c_p(j,4) = cSb_p(j)
     end do
  elseif(a.eq.2) then
     do j=1,NBS_L
        c_p(j,1) = dLa_p(j); c_p(j,2) = dLb_p(j)
     end do
     do j=1,NBS_S
        c_p(j,3) = dSa_p(j); c_p(j,4) = dSb_p(j)
     end do
  else
     write(*,*) "Use + or - in function ...."
     stop
  end if

  return
END SUBROUTINE COPY_DIRACOUTPUT_CP_MOD

!=======================================================
SUBROUTINE TIMING_REPORT(report)
!=======================================================
  use month_to_day
  implicit none

  logical,intent(in) :: report
  integer,save :: t1(8)
  integer :: t2(8)
  real(kind=8),save :: timea
  real(kind=8) :: timeb,walltime

  if(.not.report) then
    call cpu_time(timea)
    call date_and_time(values=t1)
  else
    call cpu_time(timeb)
    call date_and_time(values=t2)
    if(mod(t2(1),4).eq.0) days(2) = 29
    ! Assume that the calculation time does not exceed one year.
    walltime = (((sum(days(1:t2(2))) -sum(days(1:t1(2))) +t2(3) -t1(3))*24.d0 +t2(5) -t1(5))*60.d0 +t2(6) -t1(6))*60.d0 +t2(7) -t1(7) +(t2(8) -t1(8))*1.d-3
    write(*,'(3x,a,f13.3)') '-  CPU time is :', timeb -timea
    write(*,'(3x,a,f13.3)') '- Wall time is :', walltime
  end if

  return
END SUBROUTINE TIMING_REPORT

!============================================================
SUBROUTINE MKDIR_TODAY
!============================================================
  use EDM_calculation, only : today
  implicit none

  call date_and_time(date=today)
  open(unit=111,file='today.sh')
  write(111,'(a)') &
   "#!/bin/bash",&
   "if [ -e "//today//" ]; then",&
   "  i=0; for f in "//today//"* ; do j=0${i}; mv ${f} tmpdir.${j:(-2)}; i=$((i + 1)); done",&
   "  i=0; for f in tmpdir* ; do j=0${i}; mv ${f} "//today//".${j:(-2)}; i=$((i + 1)); done",&
   "fi",&
   "mkdir "//today
  close(111)
  call system('sh ./today.sh')
  call system('rm -f ./today.sh')

  return
END SUBROUTINE MKDIR_TODAY

!============================================================
! FUNCTION
!============================================================
FUNCTION GET_TIME
!============================================================
  implicit none

  character(len=12) :: get_time
  character(len=10) :: t

  call date_and_time(time=t)
  write(get_time,'(a12)') t(1:2)//':'//t(3:4)//':'//t(5:10)

  return
END FUNCTION GET_TIME
