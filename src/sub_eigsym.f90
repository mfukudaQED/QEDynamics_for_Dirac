!==============================================================================
! subroutines for eigenvalues and eigenvectors real symmetric matrix.
!
! 2012.9.13 
!  -copied from H2p_5.f90 and diag.f90 in 02_H2p (other project)
!  -modified so that eig(1)<eig(2)<eig(3) 
! 
!
! f90 version of tqli, tred2, pythag
!==============================================================================



!======================================================================
!subroutine calc_eigen(Ts,d,a)
subroutine calc_eigen(Ts,eig,vec)
!
! For given 3x3 matrix "Ts", calculate eigenvalues "d" and eigenvectors "a".
! d's are sorted in order of d(1)>d(2)>d(3)
! 
! [a(1,i), a(2,i), a(3,i)] is the eigenvector for i-th eigenvalue d(i)
!
!======================================================================
  implicit none

  real(kind=8), intent(in) :: Ts(3,3)
  REAL(kind=8), intent(out) :: eig(3),vec(3,3)
  REAL(kind=8) :: e(3),d(3),a(3,3)
  integer :: i,j,ii

  do i=1,3
     do j=1,3
        a(i,j) = Ts(i,j)
!        write(*,*) a(i,j)
     end do
  end do

  call tred2(a,3,3,d,e)
  call tqli(d,e,3,3,a)
  !do ii=1,3
  !write(*,'(4es16.6)') d(ii),a(1,ii),a(2,ii),a(3,ii) 
  !end do
  
  ! sort w.r.t. eigenvalue (In the order of d(1) > d(2) > d(3) )
  call eigsrt(d,a,3,3)  

 !            do ii=1,3
 !               write(*,'(4es16.6)') d(ii),a(1,ii),a(2,ii),a(3,ii) 
 !            end do
  !do ii=1,3
  !     write(*,'(4es16.6)',advance='no') d(ii),a(1,ii),a(2,ii),a(3,ii) 
  !end do
  !  write(*,*)

  ! re-sort
  eig(3) = d(1)
  eig(2) = d(2)
  eig(1) = d(3)
  do j=1,3
     vec(j,3) = a(j,1)
     vec(j,2) = a(j,2)
     vec(j,1) = a(j,3)
  end do

  return
end subroutine calc_eigen

!=====================================
subroutine tqli(d,e,iend,np,z)
!=====================================
  integer :: i,iend,iter,k,l,m,np
  real(kind=8) :: b,c,dd,f,g,p,r,s,pythag_f90
  real(kind=8) ::  d(np),e(np),z(np,np)

  do i = 2,iend
     e(i-1) = e(i)
  end do
  e(n) = 0.d0
  do l = 1,iend
     iter = 0
10    do m = l,iend-1
        dd = abs(d(m)) + abs(d(m+1))
        if (abs(e(m)) + dd.eq.dd) goto 20
     end do
     m = iend
20    if(m.ne.l)then
 
        if(iter.eq.30) stop 'Many iterations in subroutine tqli (re)'
        iter = iter+1
        g = (d(l+1)-d(l))/(2.*e(l))
        r = pythag_f90(g,1.d0)
        g = d(m)-d(l)+e(l)/(g+sign(r,g))
        s = 1.d0
        c = 1.d0
        p = 0.d0
        do i = m-1,l,-1
           f = s*e(i)
           b = c*e(i)
           r = pythag_f90(f,g)
           e(i+1) = r
           if(r.eq.0.d0)then
              d(i+1) = d(i+1)-p
              e(m) = 0.d0
              goto 10
           endif
           s = f/r
           c = g/r
           g = d(i+1)-p
           r = (d(i)-g)*s+2.d0*c*b
           p = s*r
           d(i+1) = g+p
           g = c*r-b
!  
           do k = 1,iend
              f = z(k,i+1)
              z(k,i+1) = s*z(k,i)+c*f
              z(k,i) = c*z(k,i)-s*f
           end do
! 
        end do
        d(l) = d(l)-p
        e(l) = g
        e(m) = 0.d0
        goto 10
     endif
  end do
  return
end subroutine tqli
!! This subroutine is made by refering numerical recipe

!=====================================
subroutine tred2(a,n,np,d,e)
!=====================================
  integer :: i,j,k,l,n,np
  real(kind=8) :: a(np,np),d(np),e(np)
  real(kind=8) :: f,g,h,hh,scale
!
  do i = n,2,-1
     l = i-1
     h = 0.d0
     scale = 0.d0
     if(l.gt.1)then
        do k = 1,l
           scale = scale+abs(a(i,k))
        end do
        if(scale.eq.0.d0)then
           e(i) = a(i,l)
        else
           do k = 1,l
              a(i,k) = a(i,k) / scale
              h = h + a(i,k)**2
           end do
           f = a(i,l)
           g = -sign(sqrt(h),f)
           e(i) = scale*g
           h = h - f*g
           a(i,l) = f - g
           f = 0.d0
           do j = 1,l
!
              a(j,i) = a(i,j)/h
              g = 0.d0
              do k = 1,j
                 g = g +a(j,k)*a(i,k)
              end do
              do k = j+1,l
                 g = g +a(k,j)*a(i,k)
              end do
              e(j) = g/h
              f = f +e(j)*a(i,j)
           end do
           hh = f/(h+h)
           do j = 1,l
              f = a(i,j)
              g = e(j) -hh*f
              e(j) = g
              do k = 1,j
                 a(j,k) = a(j,k) -f*e(k) -g*a(i,k)
              end do
           end do
        endif
     else
        e(i) = a(i,l)
     endif
     d(i) = h
  end do
!
  d(1) = 0.d0
  e(1) = 0.d0
  do i = 1,n
!
     l = i-1
     if(d(i).ne.0.d0)then
        do j = 1,l
           g = 0.d0
           do k = 1,l
              g = g + a(i,k) * a(k,j)
           end do
           do k = 1,l
              a(k,j) = a(k,j) - g * a(k,i)
           end do
        end do
     endif
!
     d(i) = a(i,i)
!
     a(i,i) = 1.d0
     do j = 1,l
        a(i,j) = 0.d0
        a(j,i) = 0.d0
     end do
!
  end do
  return
end subroutine tred2
!! This subroutine is made by refering numerical recipe

!=====================================
function pythag_f90(a1,a2)
!=====================================
  real(kind=8) :: a1,a2,pythag_f90
  real(kind=8) :: absr1,absr2

  absr1 = abs(a1)
  absr2 = abs(a2)

  if(absr1.gt.absr2)then
     pythag_f90 = absr1 * sqrt( 1.d0 + (absr1/absr2)**2.0 )
  else
     if(absb.eq.0.d0)then
        pythag_f90=0.d0
     else
        pythag_f90 = absr2 * sqrt( 1.d0 + (absr1/absr2)**2.0 )
     endif
  endif

  return

end function pythag_f90

!===============================
subroutine eigsrt(d,v,n,m)
!===============================
  integer n,m  
  real(kind=8) :: d(m),v(m,m)
  integer :: i,j,k
  real(kind=8) :: temp
! m :number of eigenvalue

do i=1, n-1
     k = i
     temp = d(i)
  do j=i+1, n
     if(d(j).ge.p)then
       k = j
       temp = d(j)
     endif
  end do
     if(k.ne.i)then
        d(k) = d(i)
        d(i) = temp
        do j=1,n
           temp = v(j,i)
           v(j,i) = v(j,k)
           v(j,k) = temp
        end do
     endif
end do
return


end subroutine eigsrt
