!============================================================
! SUBROUTINE
!============================================================

!============================================================
SUBROUTINE SET_INT_PSI(actorb,occ,coef)
!============================================================
  use Precision
  use DiracOutput
  !$ use omp_lib
  implicit none

  integer,intent(in) :: actorb
  real(kind=dp),intent(in) :: occ(actorb)
  complex(kind=dp) :: coef(NBS0,2*NBS,4)

  integer :: flag, id, Nthreads
  complex(kind=dp),allocatable,dimension(:,:,:) :: coefLL
  complex(kind=dp),allocatable,dimension(:,:,:) :: coefSS
  complex(kind=dp),allocatable,dimension(:,:,:) :: coefLS

  real(kind=dp) :: overlapLL (NBS_LT), overlapSS (NBS_ST), overlapLS (NBS_LS)
  real(kind=dp) :: x_momentLL(NBS_LT), x_momentSS(NBS_ST), x_momentLS(NBS_LS)
  real(kind=dp) :: y_momentLL(NBS_LT), y_momentSS(NBS_ST), y_momentLS(NBS_LS)
  real(kind=dp) :: z_momentLL(NBS_LT), z_momentSS(NBS_ST), z_momentLS(NBS_LS)
  real(kind=dp) :: xdy_ydxLL(NBS_LL), xdy_ydxSS(NBS_SS)
  real(kind=dp) :: ydz_zdyLL(NBS_LL), ydz_zdySS(NBS_SS)
  real(kind=dp) :: zdx_xdzLL(NBS_LL), zdx_xdzSS(NBS_SS)
  real(kind=dp) :: LapLS(NBS_LS), LapSL(NBS_LS)
  real(kind=dp) :: ef_xSS(NBS_ST),ef_ySS(NBS_ST),ef_zSS(NBS_ST)
  real(kind=dp) :: ef_xLS(NBS_LS),ef_yLS(NBS_LS),ef_zLS(NBS_LS)
  integer :: n

  call calc_int_gaussian( overlapLL,  overlapSS,  overlapLS,  &
                        & x_momentLL, x_momentSS, x_momentLS, &
                        & y_momentLL, y_momentSS, y_momentLS, &
                        & z_momentLL, z_momentSS, z_momentLS, &
                        & xdy_ydxLL, xdy_ydxSS,               &
                        & ydz_zdyLL, ydz_zdySS,               &
                        & zdx_xdzLL, zdx_xdzSS,               &
                        & LapLS, LapSL,                       &
                        & ef_xSS, ef_ySS, ef_zSS,             &
                        & ef_xLS, ef_yLS, ef_zLS )

   !Nthreads = omp_get_num_threads()

  flag = 0
  !$omp parallel do private(n,coefLL,coefSS,coefLS,flag)
  do n = 1,actorb
    if (flag==0) then
      allocate(coefLL(NBS_LL,1:2,1:2))
      allocate(coefSS(NBS_SS,3:4,3:4))
      allocate(coefLS(NBS_LS,1:2,3:4))
      flag = 1
    end if

    call calc_coef_XX(coef(1:NBS0,n,1:4),coefLL(:,:,:),coefSS(:,:,:),coefLS(:,:,:))
    call calc_int_psi( n, coefLL(:,:,:),coefSS(:,:,:),coefLS(:,:,:),          &
                     &  overlapLL,  overlapSS,  overlapLS, &
                     & x_momentLL, x_momentSS, x_momentLS, &
                     & y_momentLL, y_momentSS, y_momentLS, &
                     & z_momentLL, z_momentSS, z_momentLS, &
                     & xdy_ydxLL, xdy_ydxSS,               &
                     & ydz_zdyLL, ydz_zdySS,               &
                     & zdx_xdzLL, zdx_xdzSS,               &
                     & LapLS, LapSL,                       &
                     & ef_xSS, ef_ySS, ef_zSS,             &
                     & ef_xLS, ef_yLS, ef_zLS )
  end do
  !$omp end parallel do
  call contract_occupation(actorb,occ)

  return
END SUBROUTINE SET_INT_PSI

!============================================================
SUBROUTINE CALC_GAUSS_P(posA,posB,alphaA,alphaB,posP,vecPA,vecPB,alphaP,cAB,cP)
!============================================================
  use Precision
  use Constants
  implicit none

  real(kind=dp),intent(in) :: posA(3),posB(3)
  real(kind=dp),intent(in) :: alphaA,alphaB
  real(kind=dp),intent(out) :: posP(3),vecPA(3),vecPB(3),alphaP,cAB,cP

  alphaP = alphaA +alphaB
  posP(1:3) = posA(1:3) +alphaB/alphaP *(posB(1:3) -posA(1:3))
  vecPA(1:3) = posP(1:3) -posA(1:3)
  vecPB(1:3) = posP(1:3) -posB(1:3)
  cAB = exp(-alphaA*alphaB/alphaP *sum((posA(1:3) -posB(1:3))**2))
  cP = sqrt((PI/alphaP)**3)

  return
END SUBROUTINE CALC_GAUSS_P

!============================================================
SUBROUTINE GAUSSIAN_INTEGRAL_1DIM(vecPAi,vecPBi,posPi,alphaB,alphaP,nA,nB, &
                                & overlap,moment,derivative_1st,derivative_2nd)
!============================================================
  use Precision
  use Binomial
  implicit none

  real(kind=dp),intent(in) :: vecPAi,vecPBi,posPi,alphaB,alphaP
  integer,intent(in) :: nA,nB
  real(kind=dp),intent(out) :: overlap,moment,derivative_1st,derivative_2nd

  real(kind=dp) :: workA(0:nA),workB(0:nB+2)
  real(kind=dp) :: workP(0:nA+nB+2)
  real(kind=dp) :: d1,d2
  integer :: s,t

  do s = 0,nA
    workA(s) = nCr(s,nA) *vecPAi**(nA-s)
  end do
  do s = 0,nB
    workB(s) = nCr(s,nB) *vecPBi**(nB-s)
  end do
  do s = 0,nA+nB+2,2
    workP(s) = d_factorial(s-1)/sqrt((2*alphaP)**s)
  end do

  !--Overlap----------------------------------------------------------
  overlap = 0._dp
  do t = 0,nB,2
    do s = 0,nA,2
      overlap = overlap +workA(s) *workB(t) *workP(s+t)
    end do
  end do
  do t = 1,nB,2
    do s = 1,nA,2
      overlap = overlap +workA(s) *workB(t) *workP(s+t)
    end do
  end do

  !--Moment----------------------------------------------------------
  moment = 0._dp
  do s = 1,nA,2
    do t = 0,nB,2
      moment = moment +workA(s) *workB(t) *workP(s+t+1)
    end do
  end do
  do t = 1,nB,2
    do s = 0,nA,2
      moment = moment +workA(s) *workB(t) *workP(s+t+1)
    end do
  end do

  moment = moment +posPi *overlap

  !--Derivative_1st----------------------------------------------------------
  do t = 0,nB+1
    workB(t) = nCr(t,nB+1) *vecPBi**(nB+1-t)
  end do
  d1 = 0._dp
  do s = 0,nA,2
    do t = 0,nB+1,2
      d1 = d1 +workA(s) *workB(t) *workP(s+t)
    end do
  end do
  do s = 1,nA,2
    do t = 1,nB+1,2
      d1 = d1 +workA(s) *workB(t) *workP(s+t)
    end do
  end do

  do t = 0,nB-1
    workB(t) = nCr(t,nB-1) *vecPBi**(nB-1-t)
  end do
  d2 = 0._dp
  do t = 0,nB-1,2
    do s = 0,nA,2
      d2 = d2 +workA(s) *workB(t) *workP(s+t)
    end do
  end do
  do t = 1,nB-1,2
    do s = 1,nA,2
      d2 = d2 +workA(s) *workB(t) *workP(s+t)
    end do
  end do

  derivative_1st = -2._dp*alphaB*d1 +nB*d2

  !--Derivative_2nd----------------------------------------------------------
  do t = 0,nB+2
    workB(t) = nCr(t,nB+2) *vecPBi**(nB+2-t)
  end do
  d1 = 0._dp
  do s = 0,nA,2
    do t = 0,nB+2,2
      d1 = d1 +workA(s) *workB(t) *workP(s+t)
    end do
  end do
  do s = 1,nA,2
    do t = 1,nB+2,2
      d1 = d1 +workA(s) *workB(t) *workP(s+t)
    end do
  end do

  do t = 0,nB-2
    workB(t) = nCr(t,nB-2) *vecPBi**(nB-2-t)
  end do
  d2 = 0._dp
  do t = 0,nB-2,2
    do s = 0,nA,2
      d2 = d2 +workA(s) *workB(t) *workP(s+t)
    end do
  end do
  do t = 1,nB-2,2
    do s = 1,nA,2
      d2 = d2 +workA(s) *workB(t) *workP(s+t)
    end do
  end do

  derivative_2nd = -2._dp*alphaB*(2*nB+1)*overlap +4._dp*alphaB**2*d1 +nB*(nB-1)*d2

  return
END SUBROUTINE GAUSSIAN_INTEGRAL_1DIM

!============================================================
SUBROUTINE GAUSSIAN_INTEGRAL_EF(vecPA,nx,ny,nz,vecPB,nbarx,nbary,nbarz,alphaP,posP,posC,ef)
! [input]
! nx,ny,nz : type of p.g. A 
! nbarx,nbary,nbarz : type of p.g. B
! posC : position of nucleus (or where you want to evaluate electric field)
!
! [output]
! ef(1~3) : x,y,z components of electric field
!============================================================
  use Precision
  use Constants
  implicit none

  real(kind=dp),intent(in) :: vecPA(3),vecPB(3),posP(3)
  real(kind=dp),intent(in) :: alphaP
  real(kind=dp),intent(in) :: posC(3) ! position of nucleus
  integer,intent(in) :: nx,nbarx,ny,nbary,nz,nbarz

  real(kind=dp),intent(out) :: ef(3)

  real(kind=dp) :: vecPC(3) ! P-C
  real(kind=dp) :: d(0:nx+nbarx,0:nbarx,0:nx),e(0:ny+nbary,0:nbary,0:ny),f(0:nz+nbarz,0:nbarz,0:nz)
  real(kind=dp) :: R(0:nx+nbarx+1,0:ny+nbary+1,0:nz+nbarz+1)
  real(kind=dp) :: list_R000j(0:nx+nbarx+ny+nbary+nz+nbarz+3)
  real(kind=dp) :: T,dijk,const
  integer :: i,j,k
  integer :: nx_sum,ny_sum,nz_sum,j_max

  call calculate_d(alphaP,vecPA(1),vecPB(1),nx,nbarx,d)
  call calculate_d(alphaP,vecPA(2),vecPB(2),ny,nbary,e)
  call calculate_d(alphaP,vecPA(3),vecPB(3),nz,nbarz,f)

  nx_sum = nx+nbarx
  ny_sum = ny+nbary
  nz_sum = nz+nbarz
  j_max = nx_sum+ny_sum+nz_sum+3

  vecPC(1:3) = posP(1:3) -posC(1:3)
  T = alphaP*sum(vecPC(1:3)**2)

  call calc_R000j(alphaP,T,j_max,list_R000j)
  call calc_R(vecPC,j_max,list_R000j,nx_sum+1,ny_sum+1,nz_sum+1,R) ! for Electric field

  const = 2._dp*PI/alphaP

  ef(:) = 0.d0
  do k=0,nz_sum
    do j=0,ny_sum
      do i=0,nx_sum
        dijk = d(i,nbarx,nx)*e(j,nbary,ny)*f(k,nbarz,nz)
        ef(1) = ef(1) -dijk*const*R(i+1,j,k) !  [NLM | xc rc^-3] (3.21)
        ef(2) = ef(2) -dijk*const*R(i,j+1,k) !  [NLM | yc rc^-3] (3.21)
        ef(3) = ef(3) -dijk*const*R(i,j,k+1) !  [NLM | zc rc^-3] (3.21)
      end do
    end do
  end do

END SUBROUTINE GAUSSIAN_INTEGRAL_EF

!============================================================
SUBROUTINE CALCULATE_D(alphaP,PAx,PBx,n,nbar,d)
! calculate coefficients for Lambda polynomial expansion
! (p.220 of Eq.(2.17)~(2.22) McMurchie & Davidson)
!========================================================================
  use Precision
  implicit none
  
  real(kind=dp), intent(in) :: alphaP,PAx,PBx
  integer, intent(in) :: n,nbar
  real(kind=dp), intent(out) :: d(0:n+nbar,0:nbar,0:n)

  real(kind=dp) :: const
  integer :: i,j,k

  const = 1._dp/(2._dp*alphaP)

  d = 0._dp
  d(0,0,0) = 1._dp ! (2.22) McMurchie & Davidson

  ! Make table for n=0 (loop over nbar & N) 
  do j=1,nbar
    d(0,j,0) = PBx*d(0,j-1,0) +d(1,j-1,0)
    do k=1,j-1  ! loop for N 
      d(k,j,0) = d(k-1,j-1,0)*const +PBx*d(k,j-1,0) +(k+1)*d(k+1,j-1,0)
    end do
    d(j,j,0) = d(j-1,j-1,0)*const
  end do
  
  ! loop over n & N starting from n=0 table (for every nbar)
  do i=1,n
    do j=0,nbar
      d(0,j,i) = PAx*d(0,j,i-1) +d(1,j,i-1)
      do k=1,j+i-1
        d(k,j,i) = d(k-1,j,i-1)*const +PAx*d(k,j,i-1) +(k+1)*d(k+1,j,i-1)
      end do
      d(j+i,j,i) = d(j+i-1,j,i-1)*const
    end do
  end do

  return
END SUBROUTINE CALCULATE_D

!============================================================
SUBROUTINE CALC_INT_GAUSSIAN( overlapLL,  overlapSS,  overlapLS,  &
                            & x_momentLL, x_momentSS, x_momentLS, &
                            & y_momentLL, y_momentSS, y_momentLS, &
                            & z_momentLL, z_momentSS, z_momentLS, &
                            & xdy_ydxLL, xdy_ydxSS,               &
                            & ydz_zdyLL, ydz_zdySS,               &
                            & zdx_xdzLL, zdx_xdzSS,               &
                            & LapLS, LapSL,                       &
                            & ef_xSS, ef_ySS, ef_zSS,             &
                            & ef_xLS, ef_yLS, ef_zLS )
!============================================================
  use Precision
  use DiracOutput
  use Constants
  implicit none

  real(kind=dp),intent(out) :: overlapLL (NBS_LT), overlapSS (NBS_ST), overlapLS (NBS_LS)
  real(kind=dp),intent(out) :: x_momentLL(NBS_LT), x_momentSS(NBS_ST), x_momentLS(NBS_LS)
  real(kind=dp),intent(out) :: y_momentLL(NBS_LT), y_momentSS(NBS_ST), y_momentLS(NBS_LS)
  real(kind=dp),intent(out) :: z_momentLL(NBS_LT), z_momentSS(NBS_ST), z_momentLS(NBS_LS)
  real(kind=dp),intent(out) :: xdy_ydxLL(NBS_LL), xdy_ydxSS(NBS_SS)
  real(kind=dp),intent(out) :: ydz_zdyLL(NBS_LL), ydz_zdySS(NBS_SS)
  real(kind=dp),intent(out) :: zdx_xdzLL(NBS_LL), zdx_xdzSS(NBS_SS)
  real(kind=dp),intent(out) :: LapLS(NBS_LS), LapSL(NBS_LS)
  real(kind=dp),intent(out) :: ef_xSS(NBS_ST),ef_ySS(NBS_ST),ef_zSS(NBS_ST)
  real(kind=dp),intent(out) :: ef_xLS(NBS_LS),ef_yLS(NBS_LS),ef_zLS(NBS_LS)

  real(kind=dp) :: overlap_x,overlap_y,overlap_z
  real(kind=dp) :: moment_x,moment_y,moment_z
  real(kind=dp) :: d_x,d_y,d_z
  real(kind=dp) :: d_xx,d_yy,d_zz,ef(3),ef_sum(3)
  real(kind=dp) :: posP(3),vecPA(3),vecPB(3),alphaP,cAB,cP
  real(kind=dp) :: c_ef,c_integral
  integer :: i,j,k,m,n

  overlapLL(1:NBS_L) = 1._dp
  overlapSS(1:NBS_S) = 1._dp

  !$omp parallel do &
  !$omp private(i,j,k,m,n,posP,vecPA,vecPB,alphaP,cAB,cP,c_ef,c_integral) &
  !$omp private(overlap_x,overlap_y,overlap_z,moment_x,moment_y,moment_z) &
  !$omp private(d_x,d_y,d_z,d_xx,d_yy,d_zz,ef,ef_sum)
  do j = 1,NBS0
    if (j.le.NBS_L) then
    !-- LL diagonal ------------------------------------------------------------------------------------------------------------
      alphaP = aa_L(j) +aa_L(j)
      call gaussian_integral_1dim(0._dp,0._dp,xx_L(j),aa_L(j),alphaP,nx_L(j),nx_L(j),overlap_x,moment_x,d_x,d_xx)
      call gaussian_integral_1dim(0._dp,0._dp,yy_L(j),aa_L(j),alphaP,ny_L(j),ny_L(j),overlap_y,moment_y,d_y,d_yy)
      call gaussian_integral_1dim(0._dp,0._dp,zz_L(j),aa_L(j),alphaP,nz_L(j),nz_L(j),overlap_z,moment_z,d_z,d_zz)
      c_ef = c_Lnorm(j) *c_Lnorm(j)
      c_integral = c_ef *sqrt((PI/alphaP)**3)
      x_momentLL(j) = c_integral *moment_x*overlap_y*overlap_z
      y_momentLL(j) = c_integral *overlap_x*moment_y*overlap_z
      z_momentLL(j) = c_integral *overlap_x*overlap_y*moment_z
      xdy_ydxLL(j)  = c_integral *(moment_x*d_y -moment_y*d_x)*overlap_z
      ydz_zdyLL(j)  = c_integral *(moment_y*d_z -moment_z*d_y)*overlap_x
      zdx_xdzLL(j)  = c_integral *(moment_z*d_x -moment_x*d_z)*overlap_y
    end if
    if (j.le.NBS_S) then
    !-- SS diagonal ------------------------------------------------------------------------------------------------------------
      alphaP = aa_S(j) +aa_S(j)
      call gaussian_integral_1dim(0._dp,0._dp,xx_S(j),aa_S(j),alphaP,nx_S(j),nx_S(j),overlap_x,moment_x,d_x,d_xx)
      call gaussian_integral_1dim(0._dp,0._dp,yy_S(j),aa_S(j),alphaP,ny_S(j),ny_S(j),overlap_y,moment_y,d_y,d_yy)
      call gaussian_integral_1dim(0._dp,0._dp,zz_S(j),aa_S(j),alphaP,nz_S(j),nz_S(j),overlap_z,moment_z,d_z,d_zz)
      c_ef = c_Snorm(j) *c_Snorm(j)
      c_integral = c_ef *sqrt((PI/alphaP)**3)
      x_momentSS(j) = c_integral *moment_x*overlap_y*overlap_z
      y_momentSS(j) = c_integral *overlap_x*moment_y*overlap_z
      z_momentSS(j) = c_integral *overlap_x*overlap_y*moment_z
      xdy_ydxSS(j)  = c_integral *(moment_x*d_y -moment_y*d_x)*overlap_z
      ydz_zdySS(j)  = c_integral *(moment_y*d_z -moment_z*d_y)*overlap_x
      zdx_xdzSS(j)  = c_integral *(moment_z*d_x -moment_x*d_z)*overlap_y
      ef_sum(:) = 0._dp
      do k = 1,NAT
        call gaussian_integral_ef((/0._dp,0._dp,0._dp/),nx_S(j),ny_S(j),nz_S(j),(/0._dp,0._dp,0._dp/),nx_S(j),ny_S(j),nz_S(j),&
                                 &alphaP,(/xx_S(j),yy_S(j),zz_S(j)/),(/xc(k),yc(k),zc(k)/),ef)
        ef_sum(:) = ef_sum(:) +cn(k)*ef(:)
      end do
      ef_xSS(j) = c_ef *ef_sum(1) ; ef_ySS(j) = c_ef *ef_sum(2) ; ef_zSS(j) = c_ef *ef_sum(3)

      m = (j-1)*NBS_L +1
      do i = 1,NBS_L
    !-- LS ---------------------------------------------------------------------------------------------------------------------
        call calc_gauss_P((/xx_L(i),yy_L(i),zz_L(i)/),(/xx_S(j),yy_S(j),zz_S(j)/),aa_L(i),aa_S(j),posP,vecPA,vecPB,alphaP,cAB,cP)
        call gaussian_integral_1dim(vecPA(1),vecPB(1),posP(1),aa_S(j),alphaP,nx_L(i),nx_S(j),overlap_x,moment_x,d_x,d_xx)
        call gaussian_integral_1dim(vecPA(2),vecPB(2),posP(2),aa_S(j),alphaP,ny_L(i),ny_S(j),overlap_y,moment_y,d_y,d_yy)
        call gaussian_integral_1dim(vecPA(3),vecPB(3),posP(3),aa_S(j),alphaP,nz_L(i),nz_S(j),overlap_z,moment_z,d_z,d_zz)
        c_ef = c_Lnorm(i) *c_Snorm(j) *cAB
        c_integral = c_ef *cP
        overlapLS(m) = c_integral *overlap_x*overlap_y*overlap_z
        x_momentLS(m) = c_integral *moment_x*overlap_y*overlap_z
        y_momentLS(m) = c_integral *overlap_x*moment_y*overlap_z
        z_momentLS(m) = c_integral *overlap_x*overlap_y*moment_z
        LapLS(m) = c_integral *(d_xx*overlap_y*overlap_z +overlap_x*d_yy*overlap_z +overlap_x*overlap_y*d_zz)
        call gaussian_integral_ef(vecPA,nx_L(i),ny_L(i),nz_L(i),vecPB,nx_S(j),ny_S(j),nz_S(j),alphaP,posP,(/xc(1),yc(1),zc(1)/),ef)
        ef_xLS(m) = c_ef *ef(1); ef_yLS(m) = c_ef *ef(2); ef_zLS(m) = c_ef *ef(3)
    !-- SL ----------------------------------------------------------------------------------------------------------------------
        call calc_gauss_P((/xx_S(j),yy_S(j),zz_S(j)/),(/xx_L(i),yy_L(i),zz_L(i)/),aa_S(j),aa_L(i),posP,vecPA,vecPB,alphaP,cAB,cP)
        call gaussian_integral_1dim(vecPA(1),vecPB(1),posP(1),aa_L(i),alphaP,nx_S(j),nx_L(i),overlap_x,moment_x,d_x,d_xx)
        call gaussian_integral_1dim(vecPA(2),vecPB(2),posP(2),aa_L(i),alphaP,ny_S(j),ny_L(i),overlap_y,moment_y,d_y,d_yy)
        call gaussian_integral_1dim(vecPA(3),vecPB(3),posP(3),aa_L(i),alphaP,nz_S(j),nz_L(i),overlap_z,moment_z,d_z,d_zz)
        c_integral = c_Snorm(j) *c_Lnorm(i) *cP *cAB
        LapSL(m) = c_integral *(d_xx*overlap_y*overlap_z +overlap_x*d_yy*overlap_z +overlap_x*overlap_y*d_zz)
        m = m + 1
      end do
    end if

    if (j.eq.1) cycle

    n = nint((j-1)*(j-2)*0.5)

    if (j.le.NBS_L) then
    !-- LL triangular ----------------------------------------------------------------------------------------------------------
      m = n +NBS_L +1
      do i = 1,j-1
        call calc_gauss_P((/xx_L(i),yy_L(i),zz_L(i)/),(/xx_L(j),yy_L(j),zz_L(j)/),aa_L(i),aa_L(j),posP,vecPA,vecPB,alphaP,cAB,cP)
        call gaussian_integral_1dim(vecPA(1),vecPB(1),posP(1),aa_L(j),alphaP,nx_L(i),nx_L(j),overlap_x,moment_x,d_x,d_xx)
        call gaussian_integral_1dim(vecPA(2),vecPB(2),posP(2),aa_L(j),alphaP,ny_L(i),ny_L(j),overlap_y,moment_y,d_y,d_yy)
        call gaussian_integral_1dim(vecPA(3),vecPB(3),posP(3),aa_L(j),alphaP,nz_L(i),nz_L(j),overlap_z,moment_z,d_z,d_zz)
        c_ef = c_Lnorm(i) *c_Lnorm(j) *cAB
        c_integral = c_ef *cP
        overlapLL(m)  = c_integral *overlap_x*overlap_y*overlap_z
        x_momentLL(m) = c_integral *moment_x*overlap_y*overlap_z
        y_momentLL(m) = c_integral *overlap_x*moment_y*overlap_z
        z_momentLL(m) = c_integral *overlap_x*overlap_y*moment_z
        xdy_ydxLL(m)  = c_integral *(moment_x*d_y -moment_y*d_x)*overlap_z
        ydz_zdyLL(m)  = c_integral *(moment_y*d_z -moment_z*d_y)*overlap_x
        zdx_xdzLL(m)  = c_integral *(moment_z*d_x -moment_x*d_z)*overlap_y
        m = m + 1
      end do
      m = n +NBS_LT +1
      do i = 1,j-1
        call calc_gauss_P((/xx_L(j),yy_L(j),zz_L(j)/),(/xx_L(i),yy_L(i),zz_L(i)/),aa_L(j),aa_L(i),posP,vecPA,vecPB,alphaP,cAB,cP)
        call gaussian_integral_1dim(vecPA(1),vecPB(1),posP(1),aa_L(i),alphaP,nx_L(j),nx_L(i),overlap_x,moment_x,d_x,d_xx)
        call gaussian_integral_1dim(vecPA(2),vecPB(2),posP(2),aa_L(i),alphaP,ny_L(j),ny_L(i),overlap_y,moment_y,d_y,d_yy)
        call gaussian_integral_1dim(vecPA(3),vecPB(3),posP(3),aa_L(i),alphaP,nz_L(j),nz_L(i),overlap_z,moment_z,d_z,d_zz)
        c_integral = c_Lnorm(j) *c_Lnorm(i) *cP *cAB
        xdy_ydxLL(m) = c_integral *(moment_x*d_y -moment_y*d_x)*overlap_z
        ydz_zdyLL(m) = c_integral *(moment_y*d_z -moment_z*d_y)*overlap_x
        zdx_xdzLL(m) = c_integral *(moment_z*d_x -moment_x*d_z)*overlap_y
        m = m + 1
      end do
    end if

    if (j.le.NBS_S) then
    !-- SS triangular ----------------------------------------------------------------------------------------------------------
      m = n +NBS_S +1
      do i = 1,j-1
        call calc_gauss_P((/xx_S(i),yy_S(i),zz_S(i)/),(/xx_S(j),yy_S(j),zz_S(j)/),aa_S(i),aa_S(j),posP,vecPA,vecPB,alphaP,cAB,cP)
        call gaussian_integral_1dim(vecPA(1),vecPB(1),posP(1),aa_S(j),alphaP,nx_S(i),nx_S(j),overlap_x,moment_x,d_x,d_xx)
        call gaussian_integral_1dim(vecPA(2),vecPB(2),posP(2),aa_S(j),alphaP,ny_S(i),ny_S(j),overlap_y,moment_y,d_y,d_yy)
        call gaussian_integral_1dim(vecPA(3),vecPB(3),posP(3),aa_S(j),alphaP,nz_S(i),nz_S(j),overlap_z,moment_z,d_z,d_zz)
        c_ef = c_Snorm(i) *c_Snorm(j) *cAB
        c_integral = c_ef *cP
        overlapSS(m)  = c_integral *overlap_x*overlap_y*overlap_z
        x_momentSS(m) = c_integral *moment_x*overlap_y*overlap_z
        y_momentSS(m) = c_integral *overlap_x*moment_y*overlap_z
        z_momentSS(m) = c_integral *overlap_x*overlap_y*moment_z
        xdy_ydxSS(m)  = c_integral *(moment_x*d_y -moment_y*d_x)*overlap_z
        ydz_zdySS(m)  = c_integral *(moment_y*d_z -moment_z*d_y)*overlap_x
        zdx_xdzSS(m)  = c_integral *(moment_z*d_x -moment_x*d_z)*overlap_y
        ef_sum(:) = 0._dp
        do k = 1,NAT
          call gaussian_integral_ef(vecPA,nx_S(i),ny_S(i),nz_S(i),vecPB,nx_S(j),ny_S(j),nz_S(j),alphaP,posP,(/xc(k),yc(k),zc(k)/),ef)
          ef_sum(:) = ef_sum(:) +cn(k)*ef(:)
        end do
        ef_xSS(m) = c_ef *ef_sum(1); ef_ySS(m) = c_ef *ef_sum(2); ef_zSS(m) = c_ef *ef_sum(3)
        m = m + 1
      end do
      m = n +NBS_ST +1
      do i = 1,j-1
        call calc_gauss_P((/xx_S(j),yy_S(j),zz_S(j)/),(/xx_S(i),yy_S(i),zz_S(i)/),aa_S(j),aa_S(i),posP,vecPA,vecPB,alphaP,cAB,cP)
        call gaussian_integral_1dim(vecPA(1),vecPB(1),posP(1),aa_S(i),alphaP,nx_S(j),nx_S(i),overlap_x,moment_x,d_x,d_xx)
        call gaussian_integral_1dim(vecPA(2),vecPB(2),posP(2),aa_S(i),alphaP,ny_S(j),ny_S(i),overlap_y,moment_y,d_y,d_yy)
        call gaussian_integral_1dim(vecPA(3),vecPB(3),posP(3),aa_S(i),alphaP,nz_S(j),nz_S(i),overlap_z,moment_z,d_z,d_zz)
        c_integral = c_Snorm(j) *c_Snorm(i) *cAB *cP
        xdy_ydxSS(m) = c_integral *(moment_x*d_y -moment_y*d_x)*overlap_z
        ydz_zdySS(m) = c_integral *(moment_y*d_z -moment_z*d_y)*overlap_x
        zdx_xdzSS(m) = c_integral *(moment_z*d_x -moment_x*d_z)*overlap_y
        m = m + 1
      end do
    end if

  end do
  !$omp end parallel do

  return
END SUBROUTINE CALC_INT_GAUSSIAN

!============================================================
SUBROUTINE CALC_INT_PSI( n, coefLL, coefSS, coefLS,          &
                       &  overlapLL,  overlapSS,  overlapLS, &
                       & x_momentLL, x_momentSS, x_momentLS, &
                       & y_momentLL, y_momentSS, y_momentLS, &
                       & z_momentLL, z_momentSS, z_momentLS, &
                       & xdy_ydxLL, xdy_ydxSS,               &
                       & ydz_zdyLL, ydz_zdySS,               &
                       & zdx_xdzLL, zdx_xdzSS,               &
                       & LapLS, LapSL,                       &
                       & ef_xSS, ef_ySS, ef_zSS,             &
                       & ef_xLS, ef_yLS, ef_zLS )
!============================================================
  use Precision
  use DiracOutput
  use EDM_calculation
  implicit none

  integer,intent(in) :: n
  complex(kind=dp),intent(in) :: coefLL(NBS_LL,1:2,1:2)
  complex(kind=dp),intent(in) :: coefSS(NBS_SS,3:4,3:4)
  complex(kind=dp),intent(in) :: coefLS(NBS_LS,1:2,3:4)
  real(kind=dp),intent(in) :: overlapLL (NBS_LT), overlapSS (NBS_ST), overlapLS (NBS_LS)
  real(kind=dp),intent(in) :: x_momentLL(NBS_LT), x_momentSS(NBS_ST), x_momentLS(NBS_LS)
  real(kind=dp),intent(in) :: y_momentLL(NBS_LT), y_momentSS(NBS_ST), y_momentLS(NBS_LS)
  real(kind=dp),intent(in) :: z_momentLL(NBS_LT), z_momentSS(NBS_ST), z_momentLS(NBS_LS)
  real(kind=dp),intent(in) :: xdy_ydxLL(NBS_LL), xdy_ydxSS(NBS_SS)
  real(kind=dp),intent(in) :: ydz_zdyLL(NBS_LL), ydz_zdySS(NBS_SS)
  real(kind=dp),intent(in) :: zdx_xdzLL(NBS_LL), zdx_xdzSS(NBS_SS)
  real(kind=dp),intent(in) :: LapLS(NBS_LS), LapSL(NBS_LS)
  real(kind=dp),intent(in) :: ef_xSS(NBS_ST),ef_ySS(NBS_ST),ef_zSS(NBS_ST)
  real(kind=dp),intent(in) :: ef_xLS(NBS_LS),ef_yLS(NBS_LS),ef_zLS(NBS_LS)

  complex(kind=dp) :: work_D (1:4,1:4),work_T (1:4,1:4)
  complex(kind=dp) :: work_D1(1:4,1:4),work_T1(1:4,1:4)
  complex(kind=dp) :: work_D2(1:4,1:4),work_T2(1:4,1:4)
  complex(kind=dp) :: work_D3(1:4,1:4),work_T3(1:4,1:4)
  complex(kind=dp) :: work_D4(3:4,3:4),work_T4(3:4,3:4)
  complex(kind=dp) :: work_D5(3:4,3:4),work_T5(3:4,3:4)
  complex(kind=dp) :: work_D6(3:4,3:4),work_T6(3:4,3:4)
  integer :: a,b

  !--LL--------------------------------------------------------------------------------
  do b = 1,2
    do a = 1,2
      work_D (a,b) = dot_product(overlapLL (1:NBS_L), coefLL(1:NBS_L,a,b))
      work_D1(a,b) = dot_product(x_momentLL(1:NBS_L), coefLL(1:NBS_L,a,b))
      work_D2(a,b) = dot_product(y_momentLL(1:NBS_L), coefLL(1:NBS_L,a,b))
      work_D3(a,b) = dot_product(z_momentLL(1:NBS_L), coefLL(1:NBS_L,a,b))
      work_T (a,b) = dot_product(overlapLL (NBS_L+1:NBS_LT), coefLL(NBS_L+1:NBS_LT,a,b))
      work_T1(a,b) = dot_product(x_momentLL(NBS_L+1:NBS_LT), coefLL(NBS_L+1:NBS_LT,a,b))
      work_T2(a,b) = dot_product(y_momentLL(NBS_L+1:NBS_LT), coefLL(NBS_L+1:NBS_LT,a,b))
      work_T3(a,b) = dot_product(z_momentLL(NBS_L+1:NBS_LT), coefLL(NBS_L+1:NBS_LT,a,b))
    end do
  end do
  do b = 1,2
    do a = 1,2
      intPsi_overlap (n,a,b) = work_D (a,b) + work_T (a,b) + conjg(work_T (b,a))
      intPsi_x_moment(n,a,b) = work_D1(a,b) + work_T1(a,b) + conjg(work_T1(b,a))
      intPsi_y_moment(n,a,b) = work_D2(a,b) + work_T2(a,b) + conjg(work_T2(b,a))
      intPsi_z_moment(n,a,b) = work_D3(a,b) + work_T3(a,b) + conjg(work_T3(b,a))
    end do
  end do

  do a = 1,2
    intPsi_xdy_ydx(n,a,a) = dot_product(xdy_ydxLL(1:NBS_LL), coefLL(1:NBS_LL,a,a))
    intPsi_ydz_zdy(n,a,a) = dot_product(ydz_zdyLL(1:NBS_LL), coefLL(1:NBS_LL,a,a))
    intPsi_zdx_xdz(n,a,a) = dot_product(zdx_xdzLL(1:NBS_LL), coefLL(1:NBS_LL,a,a))
  end do

  !--SS--------------------------------------------------------------------------------
  do b = 3,4
    do a = 3,4
      work_D (a,b) = dot_product(overlapSS (1:NBS_S), coefSS(1:NBS_S,a,b))
      work_D1(a,b) = dot_product(x_momentSS(1:NBS_S), coefSS(1:NBS_S,a,b))
      work_D2(a,b) = dot_product(y_momentSS(1:NBS_S), coefSS(1:NBS_S,a,b))
      work_D3(a,b) = dot_product(z_momentSS(1:NBS_S), coefSS(1:NBS_S,a,b))
      work_D4(a,b) = dot_product(ef_xSS(1:NBS_S), coefSS(1:NBS_S,a,b))
      work_D5(a,b) = dot_product(ef_ySS(1:NBS_S), coefSS(1:NBS_S,a,b))
      work_D6(a,b) = dot_product(ef_zSS(1:NBS_S), coefSS(1:NBS_S,a,b))
      work_T (a,b) = dot_product(overlapSS (NBS_S+1:NBS_ST), coefSS(NBS_S+1:NBS_ST,a,b))
      work_T1(a,b) = dot_product(x_momentSS(NBS_S+1:NBS_ST), coefSS(NBS_S+1:NBS_ST,a,b))
      work_T2(a,b) = dot_product(y_momentSS(NBS_S+1:NBS_ST), coefSS(NBS_S+1:NBS_ST,a,b))
      work_T3(a,b) = dot_product(z_momentSS(NBS_S+1:NBS_ST), coefSS(NBS_S+1:NBS_ST,a,b))
      work_T4(a,b) = dot_product(ef_xSS(NBS_S+1:NBS_ST), coefSS(NBS_S+1:NBS_ST,a,b))
      work_T5(a,b) = dot_product(ef_ySS(NBS_S+1:NBS_ST), coefSS(NBS_S+1:NBS_ST,a,b))
      work_T6(a,b) = dot_product(ef_zSS(NBS_S+1:NBS_ST), coefSS(NBS_S+1:NBS_ST,a,b))
    end do
  end do
  do b = 3,4
    do a = 3,4
      intPsi_overlap (n,a,b) = work_D (a,b) + work_T (a,b) + conjg(work_T (b,a))
      intPsi_x_moment(n,a,b) = work_D1(a,b) + work_T1(a,b) + conjg(work_T1(b,a))
      intPsi_y_moment(n,a,b) = work_D2(a,b) + work_T2(a,b) + conjg(work_T2(b,a))
      intPsi_z_moment(n,a,b) = work_D3(a,b) + work_T3(a,b) + conjg(work_T3(b,a))
      intPsi_ef_x(n,a,b) = work_D4(a,b) + work_T4(a,b) + conjg(work_T4(b,a))
      intPsi_ef_y(n,a,b) = work_D5(a,b) + work_T5(a,b) + conjg(work_T5(b,a))
      intPsi_ef_z(n,a,b) = work_D6(a,b) + work_T6(a,b) + conjg(work_T6(b,a))
    end do
  end do

  do a = 3,4
    intPsi_xdy_ydx(n,a,a) = dot_product(xdy_ydxSS(1:NBS_SS), coefSS(1:NBS_SS,a,a))
    intPsi_ydz_zdy(n,a,a) = dot_product(ydz_zdySS(1:NBS_SS), coefSS(1:NBS_SS,a,a))
    intPsi_zdx_xdz(n,a,a) = dot_product(zdx_xdzSS(1:NBS_SS), coefSS(1:NBS_SS,a,a))
  end do

  !--LS--------------------------------------------------------------------------------
  do b = 3,4
    do a = 1,2
      intPsi_overlap (n,a,b) = dot_product(overlapLS (1:NBS_LS), coefLS(1:NBS_LS,a,b))
      intPsi_x_moment(n,a,b) = dot_product(x_momentLS(1:NBS_LS), coefLS(1:NBS_LS,a,b))
      intPsi_y_moment(n,a,b) = dot_product(y_momentLS(1:NBS_LS), coefLS(1:NBS_LS,a,b))
      intPsi_z_moment(n,a,b) = dot_product(z_momentLS(1:NBS_LS), coefLS(1:NBS_LS,a,b))
      intPsi_ef_x(n,a,b) = dot_product(ef_xLS(1:NBS_LS), coefLS(1:NBS_LS,a,b))
      intPsi_ef_y(n,a,b) = dot_product(ef_yLS(1:NBS_LS), coefLS(1:NBS_LS,a,b))
      intPsi_ef_z(n,a,b) = dot_product(ef_zLS(1:NBS_LS), coefLS(1:NBS_LS,a,b))
      intPsi_Lap(n,a,b) = dot_product(LapLS(1:NBS_LS), coefLS(1:NBS_LS,a,b))
    end do
  end do

  !--SL--------------------------------------------------------------------------------
  do b = 1,2
    do a = 3,4
      intPsi_Lap(n,a,b) = dot_product(LapSL(1:NBS_LS), conjg(coefLS(1:NBS_LS,b,a)))
    end do
  end do

  return
END SUBROUTINE CALC_INT_PSI

!============================================================
SUBROUTINE CONTRACT_OCCUPATION(actorb,occ)
!============================================================
  use Precision
  use DiracOutput
  use EDM_calculation
  implicit none

  integer,intent(in) :: actorb
  real(kind=dp),intent(in) :: occ(actorb)

  integer :: a,b

  !--LL--------------------------------------------------------------------------------
  do b = 1,2
    do a = 1,2
      intPsi_overlap (0,a,b) = dot_product(occ(1:actorb), intPsi_overlap (1:actorb,a,b))
      intPsi_x_moment(0,a,b) = dot_product(occ(1:actorb), intPsi_x_moment(1:actorb,a,b))
      intPsi_y_moment(0,a,b) = dot_product(occ(1:actorb), intPsi_y_moment(1:actorb,a,b))
      intPsi_z_moment(0,a,b) = dot_product(occ(1:actorb), intPsi_z_moment(1:actorb,a,b))
    end do
  end do
  do a = 1,2
    intPsi_xdy_ydx(0,a,a) = dot_product(occ(1:actorb), intPsi_xdy_ydx(1:actorb,a,a))
    intPsi_ydz_zdy(0,a,a) = dot_product(occ(1:actorb), intPsi_ydz_zdy(1:actorb,a,a))
    intPsi_zdx_xdz(0,a,a) = dot_product(occ(1:actorb), intPsi_zdx_xdz(1:actorb,a,a))
  end do

  !--SS--------------------------------------------------------------------------------
  do b = 3,4
    do a = 3,4
      intPsi_overlap (0,a,b) = dot_product(occ(1:actorb), intPsi_overlap (1:actorb,a,b))
      intPsi_x_moment(0,a,b) = dot_product(occ(1:actorb), intPsi_x_moment(1:actorb,a,b))
      intPsi_y_moment(0,a,b) = dot_product(occ(1:actorb), intPsi_y_moment(1:actorb,a,b))
      intPsi_z_moment(0,a,b) = dot_product(occ(1:actorb), intPsi_z_moment(1:actorb,a,b))
      intPsi_ef_x(0,a,b) = dot_product(occ(1:actorb), intPsi_ef_x(1:actorb,a,b))
      intPsi_ef_y(0,a,b) = dot_product(occ(1:actorb), intPsi_ef_y(1:actorb,a,b))
      intPsi_ef_z(0,a,b) = dot_product(occ(1:actorb), intPsi_ef_z(1:actorb,a,b))
    end do
  end do
  do a = 3,4
    intPsi_xdy_ydx(0,a,a) = dot_product(occ(1:actorb), intPsi_xdy_ydx(1:actorb,a,a))
    intPsi_ydz_zdy(0,a,a) = dot_product(occ(1:actorb), intPsi_ydz_zdy(1:actorb,a,a))
    intPsi_zdx_xdz(0,a,a) = dot_product(occ(1:actorb), intPsi_zdx_xdz(1:actorb,a,a))
  end do

  !--LS--------------------------------------------------------------------------------
  do b = 3,4
    do a = 1,2
      intPsi_overlap (0,a,b) = dot_product(occ(1:actorb), intPsi_overlap (1:actorb,a,b))
      intPsi_x_moment(0,a,b) = dot_product(occ(1:actorb), intPsi_x_moment(1:actorb,a,b))
      intPsi_y_moment(0,a,b) = dot_product(occ(1:actorb), intPsi_y_moment(1:actorb,a,b))
      intPsi_z_moment(0,a,b) = dot_product(occ(1:actorb), intPsi_z_moment(1:actorb,a,b))
      intPsi_ef_x(0,a,b) = dot_product(occ(1:actorb), intPsi_ef_x(1:actorb,a,b))
      intPsi_ef_y(0,a,b) = dot_product(occ(1:actorb), intPsi_ef_y(1:actorb,a,b))
      intPsi_ef_z(0,a,b) = dot_product(occ(1:actorb), intPsi_ef_z(1:actorb,a,b))
      intPsi_Lap(0,a,b) = dot_product(occ(1:actorb), intPsi_Lap(1:actorb,a,b))
    end do
  end do

  !--SL--------------------------------------------------------------------------------
  do b = 1,2
    do a = 3,4
      intPsi_overlap (0:actorb,a,b) = conjg(intPsi_overlap (0:actorb,b,a))
      intPsi_x_moment(0:actorb,a,b) = conjg(intPsi_x_moment(0:actorb,b,a))
      intPsi_y_moment(0:actorb,a,b) = conjg(intPsi_y_moment(0:actorb,b,a))
      intPsi_z_moment(0:actorb,a,b) = conjg(intPsi_z_moment(0:actorb,b,a))
      intPsi_ef_x(0:actorb,a,b) = conjg(intPsi_ef_x(0:actorb,b,a))
      intPsi_ef_y(0:actorb,a,b) = conjg(intPsi_ef_y(0:actorb,b,a))
      intPsi_ef_z(0:actorb,a,b) = conjg(intPsi_ef_z(0:actorb,b,a))
      intPsi_Lap(0,a,b) = dot_product(occ(1:actorb), intPsi_Lap(1:actorb,a,b))
    end do
  end do

  return
END SUBROUTINE CONTRACT_OCCUPATION

!============================================================
SUBROUTINE CALC_INT_SPIN(actorb,occ)
!============================================================
  use Precision
  use Constants
  use EDM_calculation
  implicit none

  integer,intent(in) :: actorb
  real(kind=dp),intent(in) :: occ(actorb)

  complex(kind=dp) :: intSpin(0:actorb,3)
  integer :: m,n

  intSpin(0:actorb,1) = intPsi_overlap(0:actorb,1,2) +intPsi_overlap(0:actorb,2,1)&
                      &+intPsi_overlap(0:actorb,3,4) +intPsi_overlap(0:actorb,4,3)
  intSpin(0:actorb,2) = IU*(-intPsi_overlap(0:actorb,1,2) +intPsi_overlap(0:actorb,2,1)&
                           &-intPsi_overlap(0:actorb,3,4) +intPsi_overlap(0:actorb,4,3))
  intSpin(0:actorb,3) = intPsi_overlap(0:actorb,1,1) -intPsi_overlap(0:actorb,2,2)&
                      &+intPsi_overlap(0:actorb,3,3) -intPsi_overlap(0:actorb,4,4)

  intSpin = 0.5_dp *intSpin

  open(unit=31,file=today//'/int_Spin.dat')
  write(31,'(a)') repeat('#',120)
  write(31,'(a)') '  Spin'
  write(31,'(a)') repeat('#',120)
  write(31,*)
  write(31,'(a7,7a16)') 'Orbital', 'x Real', 'x Imag', 'y Real', 'y Imag', 'z Real', 'z Imag', 'Occupation'
  do m = 1,actorb
    write(31,'(i7,7es16.6)') m, (intSpin(m,n), n = 1,3), occ(m)
  end do
  write(31,'(a7,6es16.6)') 'Total :', (intSpin(0,n), n = 1,3)
  close(31)

  return
END SUBROUTINE CALC_INT_SPIN

!============================================================
SUBROUTINE CALC_INT_ZETA_POTENTIAL(actorb,occ)
!============================================================
  use Precision
  use Constants
  use EDM_calculation
  implicit none

  integer,intent(in) :: actorb
  real(kind=dp),intent(in) :: occ(actorb)

  complex(kind=dp) :: intZetaPotential(0:actorb)
  integer :: m

  intZetaPotential(0:actorb) = intPsi_overlap(0:actorb,1,3) +intPsi_overlap(0:actorb,2,4)&
                             &+intPsi_overlap(0:actorb,3,1) +intPsi_overlap(0:actorb,4,2)

  open(unit=320,file=today//'/int_Electron_Chirality.dat')
  write(320,'(a)') repeat('#',120)
  write(320,'(a)') '  Electron Chirality'
  write(320,'(a)') repeat('#',120)
  write(320,*)
  write(320,'(a7,4a16)') 'Orbital', 'Real', 'Imag', 'Occupation'
  do m = 1,actorb
    write(320,'(i7,3es16.6)') m, intZetaPotential(m), occ(m)
  end do
  write(320,'(a7,2es16.6)') 'Total :', intZetaPotential(0)    !integrated electron chirality
  close(320)

  intZetaPotential = 0.5_dp *CCC *intZetaPotential

  open(unit=32,file=today//'/int_ZetaPotential.dat')
  write(32,'(a)') repeat('#',120)
  write(32,'(a)') '  Zeta Potential'
  write(32,'(a)') repeat('#',120)
  write(32,*)
  write(32,'(a7,4a16)') 'Orbital', 'Real', 'Imag', 'Occupation'
  do m = 1,actorb
    write(32,'(i7,3es16.6)') m, intZetaPotential(m), occ(m)
  end do
  write(32,'(a7,2es16.6)') 'Total :', intZetaPotential(0)
  close(32)

  return
END SUBROUTINE CALC_INT_ZETA_POTENTIAL

!============================================================
SUBROUTINE CALC_INT_ORBITAL_ANGULAR_MOMENTUM(actorb,occ)
!============================================================
  use Precision
  use Constants
  use EDM_calculation
  implicit none

  integer,intent(in) :: actorb
  real(kind=dp),intent(in) :: occ(actorb)

  complex(kind=dp) :: vec(0:actorb,3)
  integer :: m,n

  vec(0:actorb,1) = intPsi_ydz_zdy(0:actorb,1,1) +intPsi_ydz_zdy(0:actorb,2,2)&
                  &+intPsi_ydz_zdy(0:actorb,3,3) +intPsi_ydz_zdy(0:actorb,4,4)
  vec(0:actorb,2) = intPsi_zdx_xdz(0:actorb,1,1) +intPsi_zdx_xdz(0:actorb,2,2)&
                  &+intPsi_zdx_xdz(0:actorb,3,3) +intPsi_zdx_xdz(0:actorb,4,4)
  vec(0:actorb,3) = intPsi_xdy_ydx(0:actorb,1,1) +intPsi_xdy_ydx(0:actorb,2,2)&
                  &+intPsi_xdy_ydx(0:actorb,3,3) +intPsi_xdy_ydx(0:actorb,4,4)

  vec = -IU *vec

  open(unit=33,file=today//'/int_OrbAng.dat')
  write(33,'(a)') repeat('#',120)
  write(33,'(a)') '  Orbital Angular Momentum'
  write(33,'(a)') repeat('#',120)
  write(33,*)
  write(33,'(a7,7a16)') 'Orbital', 'x Real', 'x Imag', 'y Real', 'y Imag', 'z Real', 'z Imag', 'Occupation'
  do m = 1,actorb
    write(33,'(i7,7es16.6)') m, (vec(m,n), n = 1,3), occ(m)
  end do
  write(33,'(a7,6es16.6)') 'Total :', (vec(0,n), n = 1,3)
  close(33)

  return
END SUBROUTINE CALC_INT_ORBITAL_ANGULAR_MOMENTUM

!============================================================
SUBROUTINE CALC_INT_TORQAM(actorb,occ)
!============================================================
  use Precision
  use Constants
  use param_AM
  use EDM_calculation
  implicit none

  integer,intent(in) :: actorb
  real(kind=dp),intent(in) :: occ(actorb)

  complex(kind=dp) :: work(0:actorb,3,3)
  complex(kind=dp) :: int_tauAM (0:actorb,3,3)
  complex(kind=dp) :: int_torqAM(0:actorb,3)
  integer :: m,n

  work(0:actorb,1,1) = intPsi_x_moment(0:actorb,1,4) +intPsi_x_moment(0:actorb,2,3)&
                     &+intPsi_x_moment(0:actorb,3,2) +intPsi_x_moment(0:actorb,4,1)
  work(0:actorb,1,2) = IU*(-intPsi_x_moment(0:actorb,1,4) +intPsi_x_moment(0:actorb,2,3)&
                          &-intPsi_x_moment(0:actorb,3,2) +intPsi_x_moment(0:actorb,4,1))
  work(0:actorb,1,3) = intPsi_x_moment(0:actorb,1,3) -intPsi_x_moment(0:actorb,2,4)&
                     &+intPsi_x_moment(0:actorb,3,1) -intPsi_x_moment(0:actorb,4,2)
  work(0:actorb,2,1) = intPsi_y_moment(0:actorb,1,4) +intPsi_y_moment(0:actorb,2,3)&
                     &+intPsi_y_moment(0:actorb,3,2) +intPsi_y_moment(0:actorb,4,1)
  work(0:actorb,2,2) = IU*(-intPsi_y_moment(0:actorb,1,4) +intPsi_y_moment(0:actorb,2,3)&
                          &-intPsi_y_moment(0:actorb,3,2) +intPsi_y_moment(0:actorb,4,1))
  work(0:actorb,2,3) = intPsi_y_moment(0:actorb,1,3) -intPsi_y_moment(0:actorb,2,4)&
                     &+intPsi_y_moment(0:actorb,3,1) -intPsi_y_moment(0:actorb,4,2)
  work(0:actorb,3,1) = intPsi_z_moment(0:actorb,1,4) +intPsi_z_moment(0:actorb,2,3)&
                     &+intPsi_z_moment(0:actorb,3,2) +intPsi_z_moment(0:actorb,4,1)
  work(0:actorb,3,2) = IU*(-intPsi_z_moment(0:actorb,1,4) +intPsi_z_moment(0:actorb,2,3)&
                          &-intPsi_z_moment(0:actorb,3,2) +intPsi_z_moment(0:actorb,4,1))
  work(0:actorb,3,3) = intPsi_z_moment(0:actorb,1,3) -intPsi_z_moment(0:actorb,2,4)&
                     &+intPsi_z_moment(0:actorb,3,1) -intPsi_z_moment(0:actorb,4,2)

  int_tauAM(0:actorb,1,2) = Bmag(2)*work(0:actorb,3,2) -Bmag(3)*work(0:actorb,2,2)
  int_tauAM(0:actorb,1,3) = Bmag(2)*work(0:actorb,3,3) -Bmag(3)*work(0:actorb,2,3)
  int_tauAM(0:actorb,2,1) = Bmag(3)*work(0:actorb,1,1) -Bmag(1)*work(0:actorb,3,1)
  int_tauAM(0:actorb,2,3) = Bmag(3)*work(0:actorb,1,3) -Bmag(1)*work(0:actorb,3,3)
  int_tauAM(0:actorb,3,1) = Bmag(1)*work(0:actorb,2,1) -Bmag(2)*work(0:actorb,1,1)
  int_tauAM(0:actorb,3,2) = Bmag(1)*work(0:actorb,2,2) -Bmag(2)*work(0:actorb,1,2)

  int_tauAM = 0.5_dp*Ze*int_tauAM

  int_torqAM(0:actorb,1) = -int_tauAM(0:actorb,2,3) +int_tauAM(0:actorb,3,2)
  int_torqAM(0:actorb,2) = -int_tauAM(0:actorb,3,1) +int_tauAM(0:actorb,1,3)
  int_torqAM(0:actorb,3) = -int_tauAM(0:actorb,1,2) +int_tauAM(0:actorb,2,1)

  open(unit=34,file=today//'/int_torqAM.dat')
  write(34,'(a)') repeat('#',120)
  write(34,'(a,es13.5,a,es13.5,a,es13.5,a)')&
                 &'  torqAM  vecB_M = (',  Bmag(1), ', ', Bmag(2), ', ', Bmag(3), ')'
  write(34,'(a)') repeat('#',120)
  write(34,*)
  write(34,'(a7,7a16)') 'Orbital', 'x Real', 'x Imag', 'y Real', 'y Imag', 'z Real', 'z Imag', 'Occupation'
  do m = 1,actorb
    write(34,'(i7,7es16.6)') m, (int_torqAM(m,n), n = 1,3), occ(m)
  end do
  write(34,'(a7,6es16.6)') 'Total :', (int_torqAM(0,n), n = 1,3)
  close(34)

  return
END SUBROUTINE CALC_INT_TORQAM

!============================================================
SUBROUTINE CALC_INT_DIPOLE(actorb,occ)
!============================================================
  use Precision
  use EDM_calculation
  use DiracOutput
  use Constants
  implicit none

  integer,intent(in) :: actorb
  real(kind=dp),intent(in) :: occ(actorb)

  complex(kind=dp) :: dipole(0:actorb,3)
  complex(kind=dp) :: dipole_Debye(0:actorb,3)
  complex(kind=dp) :: dipole_nuc(NAT,3)
  complex(kind=dp) :: dipole_nuc_Debye(NAT,3)
  integer :: m,n

  dipole(0:actorb,1) = intPsi_x_moment(0:actorb,1,1) +intPsi_x_moment(0:actorb,2,2)&
                     &+intPsi_x_moment(0:actorb,3,3) +intPsi_x_moment(0:actorb,4,4)
  dipole(0:actorb,2) = intPsi_y_moment(0:actorb,1,1) +intPsi_y_moment(0:actorb,2,2)&
                     &+intPsi_y_moment(0:actorb,3,3) +intPsi_y_moment(0:actorb,4,4)
  dipole(0:actorb,3) = intPsi_z_moment(0:actorb,1,1) +intPsi_z_moment(0:actorb,2,2)&
                     &+intPsi_z_moment(0:actorb,3,3) +intPsi_z_moment(0:actorb,4,4)

  dipole = Ze *dipole
  dipole_Debye = dipole *2.54177000

  do m = 1,NAT
    dipole_nuc(m,1:3) = cn(m) *(/xc(m),yc(m),zc(m)/)
    dipole_nuc_Debye(m,1:3) = dipole_nuc(m,1:3) *2.54177000
  end do

  open(unit=35,file=today//'/int_dipole.dat')
  write(35,'(a)') repeat('#',120)
  write(35,'(a)') '  Dipole Moment'
  write(35,'(a)') repeat('#',120)
  write(35,*)
  write(35,'(a7,7a16)') 'Orbital', 'x Real', 'x Imag', 'y Real', 'y Imag', 'z Real', 'z Imag', 'Occupation'
  do m = 1,actorb
    write(35,'(i7,7es16.6)') m, (dipole(m,n), n = 1,3), occ(m)
  end do
  write(35,'(2x,a)') 'Electronic contribution :'       
  write(35,'(7x,6es16.6)') (dipole(0,n), n = 1,3)      
  write(35,'(2x,a)') 'Electronic contribution[Debye] :'
  write(35,'(7x,6es16.6)') (dipole_Debye(0,n), n = 1,3)
  write(35,*)
  write(35,'(2x,a)') 'Nuclear contribution :'
  do m = 1,NAT
    write(35,'(i5,a,6es16.6)') m, ' :', (dipole_nuc(m,n), n = 1,3)
    dipole(0,1:3) = dipole(0,1:3) +dipole_nuc(m,1:3)
  end do
  write(35,'(2x,a)') 'Nuclear contribution[Debye] :'
  do m = 1,NAT
    write(35,'(i5,a,6es16.6)') m, ' :', (dipole_nuc_Debye(m,n), n = 1,3)
    dipole_Debye(0,1:3) = dipole_Debye(0,1:3) +dipole_nuc_Debye(m,1:3)
  end do
  write(35,*)
  write(35,'(a7,6es16.6)') 'Total :', (dipole(0,n), n = 1,3)
  write(35,'(a14)') 'Total[Debye] :'                   
  write(35,'(7x,6es16.6)') (dipole_Debye(0,n), n = 1,3)
  close(35)

  return
END SUBROUTINE CALC_INT_DIPOLE

!============================================================
SUBROUTINE CALC_INT_TORQEDM_ELE(actorb,occ)
!============================================================
  use Precision
  use Constants
  use param_AM
  use EDM_calculation
  implicit none

  integer,intent(in) :: actorb
  real(kind=dp),intent(in) :: occ(actorb)

  complex(kind=dp) :: work(0:actorb,3)
  complex(kind=dp) :: int_torqEDM_ele(0:actorb,3)
  integer :: m,n

  work(0:actorb,1) = intPsi_overlap(0:actorb,1,2) +intPsi_overlap(0:actorb,2,3)&
                   &-intPsi_overlap(0:actorb,3,4) -intPsi_overlap(0:actorb,4,3)
  work(0:actorb,2) = IU*(-intPsi_overlap(0:actorb,1,2) +intPsi_overlap(0:actorb,2,1)&
                        &+intPsi_overlap(0:actorb,3,4) -intPsi_overlap(0:actorb,4,3))
  work(0:actorb,3) = intPsi_overlap(0:actorb,1,1) -intPsi_overlap(0:actorb,2,2)&
                   &-intPsi_overlap(0:actorb,3,3) +intPsi_overlap(0:actorb,4,4)

  int_torqEDM_ele(0:actorb,1) = work(0:actorb,2) *Eelec(3) -work(0:actorb,3) *Eelec(2)
  int_torqEDM_ele(0:actorb,2) = work(0:actorb,3) *Eelec(1) -work(0:actorb,1) *Eelec(3)
  int_torqEDM_ele(0:actorb,3) = work(0:actorb,1) *Eelec(2) -work(0:actorb,2) *Eelec(1)

  open(unit=36,file=today//'/int_torqEDMele.dat')
  write(36,'(a)') repeat('#',120)
  write(36,'(a,es13.5,a,es13.5,a,es13.5,a)')&
                 &'  torqEDM_ele  vecE_M = (',  Eelec(1), ', ', Eelec(2), ', ', Eelec(3), ')'
  write(36,'(a)') repeat('#',120)
  write(36,*)
  write(36,'(a7,7a16)') 'Orbital', 'x Real', 'x Imag', 'y Real', 'y Imag', 'z Real', 'z Imag', 'Occupation'
  do m = 1,actorb
    write(36,'(i7,7es16.6)') m, (int_torqEDM_ele(m,n), n = 1,3), occ(m)
  end do
  write(36,'(a7,6es16.6)') 'Total :', (int_torqEDM_ele(0,n), n = 1,3)
  close(36)

  return
END SUBROUTINE CALC_INT_TORQEDM_ELE

!============================================================
SUBROUTINE CALC_INT_TORQEDM_MAG(actorb,occ)
!============================================================
  use Precision
  use Constants
  use param_AM
  use EDM_calculation
  implicit none

  integer,intent(in) :: actorb
  real(kind=dp),intent(in) :: occ(actorb)

  complex(kind=dp) :: work(0:actorb,3)
  complex(kind=dp) :: int_torqEDM_mag(0:actorb,3)
  integer :: m,n

  work(0:actorb,1) = intPsi_overlap(0:actorb,1,4) +intPsi_overlap(0:actorb,2,3)&
                   &-intPsi_overlap(0:actorb,3,2) -intPsi_overlap(0:actorb,4,1)
  work(0:actorb,2) = IU*(-intPsi_overlap(0:actorb,1,4) +intPsi_overlap(0:actorb,2,3)&
                        &+intPsi_overlap(0:actorb,3,2) -intPsi_overlap(0:actorb,4,1))
  work(0:actorb,3) = intPsi_overlap(0:actorb,1,3) -intPsi_overlap(0:actorb,2,4)&
                   &-intPsi_overlap(0:actorb,3,1) +intPsi_overlap(0:actorb,4,2)

  int_torqEDM_mag(0:actorb,1) = work(0:actorb,2) *Bmag(3) -work(0:actorb,3) *Bmag(2)
  int_torqEDM_mag(0:actorb,2) = work(0:actorb,3) *Bmag(1) -work(0:actorb,1) *Bmag(3)
  int_torqEDM_mag(0:actorb,3) = work(0:actorb,1) *Bmag(2) -work(0:actorb,2) *Bmag(1)

  int_torqEDM_mag = IU *int_torqEDM_mag

  open(unit=37,file=today//'/int_torqEDMmag.dat')
  write(37,'(a)') repeat('#',120)
  write(37,'(a,es13.5,a,es13.5,a,es13.5,a)')&
                  &'  torqEDM_mag  vecB_M = (',  Bmag(1), ', ', Bmag(2), ', ', Bmag(3), ')'
  write(37,'(a)') repeat('#',120)
  write(37,*)
  write(37,'(a7,7a16)') 'Orbital', 'x Real', 'x Imag', 'y Real', 'y Imag', 'z Real', 'z Imag', 'Occupation'
  do m = 1,actorb
    write(37,'(i7,7es16.6)') m, (int_torqEDM_mag(m,n), n = 1,3), occ(m)
  end do
  write(37,'(a7,6es16.6)') 'Total :', (int_torqEDM_mag(0,n), n = 1,3)
  close(37)

  return
END SUBROUTINE CALC_INT_TORQEDM_MAG

!============================================================
SUBROUTINE CALC_INT_EEFFEDM_NUC(actorb,occ)
!============================================================
  use Precision
  use Constants
  use DiracOutput
  use EDM_calculation
  implicit none

  integer,intent(in) :: actorb
  real(kind=dp),intent(in) :: occ(actorb)

  complex(kind=dp) :: int_EeffEDM_nuc(0:actorb)
  integer :: m

  int_EeffEDM_nuc(0:actorb) = intPsi_ef_x(0:actorb,3,4) +intPsi_ef_x(0:actorb,4,3)&
                       &+IU*(-intPsi_ef_y(0:actorb,3,4) +intPsi_ef_y(0:actorb,4,3))&
                            &+intPsi_ef_z(0:actorb,3,3) -intPsi_ef_z(0:actorb,4,4)

  int_EeffEDM_nuc = -2._dp*int_EeffEDM_nuc

  open(unit=38,file=today//'/int_EeffEDMnuc.dat')
  write(38,'(a)') repeat('#',120)
  write(38,'(a)') '  EeffEDM_nuc'
  write(38,'(a)') repeat('#',120)
  write(38,*)
  write(38,'(a7,3a25)') 'Orbital', 'Real', 'Imag', 'Occupation'
  do m = 1,actorb
    write(38,'(i7,3es25.15)') m, int_EeffEDM_nuc(m), occ(m)
  end do
  write(38,'(a7,2es25.15)') 'Total :', int_EeffEDM_nuc(0)
  write(38,'(a14,2es25.15)') 'Total[GV/cm] :', int_EeffEDM_nuc(0) *5.1422067476
  close(38)

  return
END SUBROUTINE CALC_INT_EEFFEDM_NUC

!============================================================
SUBROUTINE CALC_INT_EEFFEDM_OB(actorb,occ)
!============================================================
  use Precision
  use Constants
  use EDM_calculation
  implicit none

  integer,intent(in) :: actorb
  real(kind=dp),intent(in) :: occ(actorb)

  complex(kind=dp) :: int_EeffEDM_ob(0:actorb)
  integer :: m

  int_EeffEDM_ob(0:actorb) = intPsi_Lap(0:actorb,1,3) +intPsi_Lap(0:actorb,2,4)&
                           &-intPsi_Lap(0:actorb,3,1) -intPsi_Lap(0:actorb,4,2)

  int_EeffEDM_ob = -2._dp*IU*CCC*int_EeffEDM_ob/Ze

  open(unit=39,file=today//'/int_EeffEDMob.dat')
  write(39,'(a)') repeat('#',120)
  write(39,'(a)') '  EeffEDM_ob'
  write(39,'(a)') repeat('#',120)
  write(39,*)
  write(39,'(a7,3a25)') 'Orbital', 'Real', 'Imag', 'Occupation'
  do m = 1,actorb
    write(39,'(i7,3es25.15)') m, int_EeffEDM_ob(m), occ(m)
  end do
  write(39,'(a7,2es25.15)') 'Total :', int_EeffEDM_ob(0)
  write(39,'(a14,2es25.15)') 'Total[GV/cm] :', int_EeffEDM_ob(0) *5.1422067476
  close(39)

  return
END SUBROUTINE CALC_INT_EEFFEDM_OB

!============================================================
SUBROUTINE CALC_INT_MAGHYPER(actorb,occ)
!============================================================
  use Precision
  use Constants
  use param_AM
  use EDM_calculation
  implicit none

  integer,intent(in) :: actorb
  real(kind=dp),intent(in) :: occ(actorb)

  complex(kind=dp) :: int_MagHyper(0:actorb)
  complex(kind=dp) :: work(0:actorb,6)
  integer :: m

  work(0:actorb,1) = intPsi_ef_y(0:actorb,1,4) +intPsi_ef_y(0:actorb,2,3)&
                   &+intPsi_ef_y(0:actorb,3,2) +intPsi_ef_y(0:actorb,4,1)
  work(0:actorb,2) = IU*(-intPsi_ef_x(0:actorb,1,4) +intPsi_ef_x(0:actorb,2,3)&
                        &-intPsi_ef_x(0:actorb,3,2) +intPsi_ef_x(0:actorb,4,1))
  work(0:actorb,3) = IU*(-intPsi_ef_z(0:actorb,1,4) +intPsi_ef_z(0:actorb,2,3)&
                        &-intPsi_ef_z(0:actorb,3,2) +intPsi_ef_z(0:actorb,4,1))
  work(0:actorb,4) = intPsi_ef_y(0:actorb,1,3) -intPsi_ef_y(0:actorb,2,4)&
                   &+intPsi_ef_y(0:actorb,3,1) -intPsi_ef_y(0:actorb,4,2)
  work(0:actorb,5) = intPsi_ef_x(0:actorb,1,3) -intPsi_ef_x(0:actorb,2,4)&
                   &+intPsi_ef_x(0:actorb,3,1) -intPsi_ef_x(0:actorb,4,2)
  work(0:actorb,6) = intPsi_ef_z(0:actorb,1,4) +intPsi_ef_z(0:actorb,2,3)&
                   &+intPsi_ef_z(0:actorb,3,2) +intPsi_ef_z(0:actorb,4,1)

  int_MagHyper(0:actorb) = magdip(3)*(work(0:actorb,1) -work(0:actorb,2))&
                         &+magdip(1)*(work(0:actorb,3) -work(0:actorb,4))&
                         &+magdip(2)*(work(0:actorb,5) -work(0:actorb,6))

  int_MagHyper = Ze*int_MagHyper

  open(unit=40,file=today//'/int_MagHyper.dat')
  write(40,'(a)') repeat('#',120)
  write(40,'(a,es13.5,a,es13.5,a,es13.5,a)')&
                 &'  Magnetic Hyperfine Structure Constant  magdip_nuc = (',  magdip(1), ', ', magdip(2), ', ', magdip(3), ')'
  write(40,'(a)') repeat('#',120)
  write(40,*)
  write(40,'(a7,3a25)') 'Orbital', 'Real', 'Imag', 'Occupation'
  do m = 1,actorb
    write(40,'(i7,3es25.15)') m, int_MagHyper(m), occ(m)
  end do
  write(40,'(a7,2es25.15)') 'Total :', int_MagHyper(0)
  close(40)

  return
END SUBROUTINE CALC_INT_MAGHYPER
