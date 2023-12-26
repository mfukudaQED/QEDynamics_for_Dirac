! Last Change: 12-Oct-2015.
!============================================================================
!  subroutine CIread_1                    2012/04/16 
!  subroutine CIread_2                    2012/04/16 
!  subroutine convert_coef_MRDFTCI        2012/04/16 
!  subroutine convert_coef_QEDtoMRDFT     2012/04/16 
!  subroutine convert_pg_QEDtoMRDFT       2012/04/16 
!============================================================================

!============================================================================
subroutine CIread_1(norb,NROOT,enorb,KRAM,neact,detcnt,detcntnum)
!============================================================================
       implicit none

       integer norb,i,neact,NROOT,nrcnt,detcnt,j,k,l
       integer cnt
       double precision enorb(norb)
       character(len=1) REPCH(14),KRAM(norb)
       character(len=14) MOJI
       character(len=4) MOJI4
       character(len=10) MOJI10
       character(len=30) MOJI30
       character(len=4) cras
       integer detcntnum
       integer :: countflag = 0

!       open(unit=44, file=FIFI4)

       read (44,4100) MOJI ! number of orbitals
       write(*,*) "Number of active orbitals  :", trim(MOJI)

       do i=1,norb
       KRAM(i) = 'F'
       read(44,1100) (REPCH(j),j=1,14),enorb(i)
!       write(*,1100) (REPCH(j),j=1,14),enorb(i)
         do j=1,14
           if (REPCH(j).eq.'2') KRAM(i) ='T'
         end do
!         write(*,*) KRAM(i)
       end do

       read (44,'(i4,a30)') neact, MOJI30
         if (MOJI30.ne."     Number of active electron") then
            write(*,*) 'error  Number of active electron'
            stop
         end if
       write(*,*) 'Number of active electrons :', neact

!
       read (44,'(i8,a30)') detcntnum, MOJI30
         if (MOJI30.ne.'         Number of determinant') then
           backspace(44)
           countflag = 1
         else
           backspace(44)
         end if

       do
!          read (44,'(a4,i4)') MOJI4, nrcnt
!          if ((MOJI4.eq.'ROOT').and.(nrcnt.eq.NROOT)) exit
          read (44,'(a4)') MOJI4
          if (MOJI4.eq.'ROOT') then
            backspace(44)
            read (44,'(4x,i4)') nrcnt
              if (nrcnt.eq.NROOT) exit
          end if
       end do
       write(*,*) 'Determinant of Root :', nrcnt

       if (countflag.eq.1) then
         detcntnum = 0
         do
            read (44,'(a10)') MOJI10
            if(MOJI10.eq.' Strings: ') then
               detcntnum = detcntnum +1
            else if(MOJI10.eq.'    0     ') then
                exit
            end if
         end do
       end if

       detcnt = detcntnum

       write (*,*) 'Count strings :', detcntnum
!--END--


! 1100 FORMAT(20X,14A,6X,G20.10)
 ! This format have to be specified in the same way of "dirrci/cidirr.F"
 1100 FORMAT(20X,14A,6X,G22.14)
 1200 format(i4)
 1300 FORMAT(15X,6(I4,6X))
 1310 FORMAT(6(I4,6X))
 1350 FORMAT(10X,100(I4,6X))
 1400 FORMAT(9X,D22.15)
 1401 FORMAT(8X,2D23.15)
 4100 format(a14)

end subroutine CIread_1

!============================================================================
subroutine CIread_2(norb,NROOT,enorb,ras,neact,detcnt,occdet,detcntnum)
!============================================================================
       use DIRRCI_coef, only : thresh_occdet
       implicit none

       integer norb,i,neact,NROOT,nrcnt,detcnt,j,k,l
       integer detcntnum
       integer cnt
       integer :: ras(neact,detcntnum)
       integer :: sumkram(detcntnum)
       double precision enorb(norb)
       complex(kind(0d0)) :: occdet(detcntnum)
       character(len=14) MOJI
       character(len=4) MOJI4
       character(len=10) MOJI10
       character(len=30) MOJI30
       character(len=4) cras
!       integer :: I1E,I1H,I2E,I2H,I3E

       do
!          read (44,'(a4,i4)') MOJI4, nrcnt
!          if ((MOJI4.eq.'ROOT').and.(nrcnt.eq.NROOT)) exit
          read (44,'(a4)') MOJI4
          if (MOJI4.eq.'ROOT') then
            backspace(44)
            read (44,'(4x,i4)') nrcnt
              if (nrcnt.eq.NROOT) exit
          end if
       end do

       detcnt = 1
       sumkram(:) = 0

       do j=1,detcntnum!read ras start
          read (44,1401) occdet(detcnt)
!          write(*,*)'occdet', occdet(detcnt)
            if (abs(occdet(detcnt)).gt.thresh_occdet) then
               read (44,1350) (ras(i,detcnt),i=1,neact) !Strings
!               read (44,'(14X,5I4)') I1E,I1H,I2E,I2H,I3E

               detcnt = detcnt + 1
            else
               read (44,*)
!               read (44,*)
            end if

       end do !end read ras

            detcnt = detcnt -1
!            if(detcnt.eq.detcntnum) then
!               write (*,*) 'Number of determinant :', detcnt
!            else
!               write (*,*) 'error Number of determinant :', detcnt, detcntnum
!               stop
!            end if

       write(*,*)'detcnt=',detcnt,'detcntnum=',detcntnum


!       do j=1, detcnt
!          write(*,*) 'Determinant :',j,(ras(i,j),i=1,neact)
!          write(*,*) 'Coefficient :',occdet(j)
!       end do

!       do i=1,norb
!       write(*,1100) (REPCH(j),j=1,14),enorb(i)
!       end do

!      close(unit=44)

 1100 FORMAT(20X,14A,6X,G20.10)
 1200 format(i4)
 1300 FORMAT(15X,6(I4,6X))
 1310 FORMAT(6(I4,6X))
 1350 FORMAT(10X,100(I4,6X))
 1400 FORMAT(9X,D22.15)
 1401 FORMAT(8X,2D23.15)
 4100 format(a14)

end subroutine CIread_2

!============================================================================
subroutine CIread_3(norb,NROOT,enorb,ras,neact,detcnt,occdet,detcntnum)
!============================================================================
       use DIRRCI_coef, only : thresh_occdet
       implicit none

       integer norb,i,neact,NROOT,nrcnt,detcnt,j,k,l
       integer detcntnum
       integer cnt
       integer :: ras(neact,detcntnum)
       double precision enorb(norb)
       complex(kind(0d0)) :: occdet(detcntnum)
       character(len=14) MOJI
       character(len=4) MOJI4
       character(len=10) MOJI10
       character(len=30) MOJI30
       character(len=4) cras
       integer :: I1E,I1H,I2E,I2H,I3E

       do
!          read (44,'(a4,i4)') MOJI4, nrcnt
!          if ((MOJI4.eq.'ROOT').and.(nrcnt.eq.NROOT)) exit
          read (44,'(a4)') MOJI4
          if (MOJI4.eq.'ROOT') then
            backspace(44)
            read (44,'(4x,i4)') nrcnt
              if (nrcnt.eq.NROOT) exit
          end if
       end do

       detcnt = 1

       do j=1,detcntnum!read ras start
          read (44,1401) occdet(detcnt)
!          write(*,*)'occdet', occdet(detcnt)
            if (abs(occdet(detcnt)).gt.thresh_occdet) then
               read (44,1350) (ras(i,detcnt),i=1,neact) !Strings
               read (44,'(14X,5I4)') I1E,I1H,I2E,I2H,I3E

               detcnt = detcnt + 1
            else
               read (44,*)
               read (44,*)
            end if

       end do !end read ras

            detcnt = detcnt -1
!            if(detcnt.eq.detcntnum) then
!               write (*,*) 'Number of determinant :', detcnt
!            else
!               write (*,*) 'error Number of determinant :', detcnt, detcntnum
!               stop
!            end if

       write(*,*) '# detcnt :',detcnt
       write(*,*) '# detcntnum :',detcntnum

!       do j=1, detcnt
!          write(*,*) 'Determinant :',j,(ras(i,j),i=1,neact)
!          write(*,*) 'Coefficient :',occdet(j)
!       end do

!       do i=1,norb
!       write(*,1100) (REPCH(j),j=1,14),enorb(i)
!       end do

!      close(unit=44)

 1100 FORMAT(20X,14A,6X,G20.10)
 1200 format(i4)
 1300 FORMAT(15X,6(I4,6X))
 1310 FORMAT(6(I4,6X))
 1350 FORMAT(10X,100(I4,6X))
 1400 FORMAT(9X,D22.15)
 1401 FORMAT(8X,2D23.15)
 4100 format(a14)

end subroutine CIread_3

!============================================================================
subroutine CIread_KRCI1(coreorb,norb,neact,detcntall)
!============================================================================
  use Precision
  use DIRRCI_coef, only : thresh_occdet
  implicit none

  integer,intent(inout) :: coreorb
  integer,intent(out) :: norb,neact
  integer(kind=dp),intent(out) :: detcntall

  character(len=26) MOJI26


  read (44,1450) coreorb
  coreorb = coreorb *2
  read (44,1450) norb
  norb = norb *2
  read (44,1450) neact
  read (44,1451) detcntall
  write(*,*)'# coreorb  :',coreorb
  write(*,*)'# norb     :',norb
  write(*,*)'# neact    :',neact
  write(*,*)'# detcntall:',detcntall

 1450 FORMAT(27X,I4)
 1451 FORMAT(27X,I12)

end subroutine CIread_KRCI1

!============================================================================
subroutine CIread_KRCI1atom(coreorbg,coreorbu,coreorb,norbg,norbu,norb,neact,detcntall)
!============================================================================
  use Precision
  use DIRRCI_coef, only : thresh_occdet
  implicit none

  integer,intent(out) :: coreorbg,norbg
  integer,intent(out) :: coreorbu,norbu
  integer :: coreorb,norb,neact
  integer(kind=dp),intent(out) :: detcntall

  character(len=26) MOJI26


  read (44,1450) coreorbg,coreorbu
  coreorb = coreorbg + coreorbu
  coreorbg = coreorbg *2
  coreorbu = coreorbu *2
  coreorb = coreorb *2
  read (44,1450) norbg, norbu
  norb = norbg + norbu
  norbg = norbg *2
  norbu = norbu *2
  norb = norb *2
  read (44,1450) neact
  read (44,1451) detcntall
  write(*,*)'# coreorbg :',coreorbg
  write(*,*)'# coreorbu :',coreorbu
  write(*,*)'# coreorb  :',coreorb
  write(*,*)'# norbg :',norbg
  write(*,*)'# norbu :',norbu
  write(*,*)'# norb     :',norb
  write(*,*)'# neact    :',neact
  write(*,*)'# detcntall:',detcntall

 1450 FORMAT(27X,2I4)
 1451 FORMAT(27X,I12)

end subroutine CIread_KRCI1atom

!============================================================================
subroutine CIread_KRCI2(NROOT,neact,detcntall,detcnt,gasA,gasB,nael,nbel,occdet)
!============================================================================
       use Precision
       use DIRRCI_coef, only : thresh_occdet
       implicit none

       integer,intent(in) :: NROOT
       integer,intent(in) :: neact
       integer(kind=dp),intent(in) :: detcntall
       integer(kind=dp),intent(out) :: detcnt
       integer,intent(out) ::  gasA(neact,detcntall),gasB(neact,detcntall)
       integer,intent(out) ::  nael(detcntall), nbel(detcntall)
       complex(kind(0d0)),intent(out) :: occdet(detcntall)
       integer :: i,j
       integer :: nrcnt
       character(len=4) MOJI4
       character(len=6) MOJI6

       do
          read(44,'(A4)') MOJI4
          if(MOJI4.eq.'ROOT') then
             backspace(44)
             read(44,'(4X,I4)') nrcnt
             if(nrcnt.eq.NROOT) exit
          end if
       end do
       write(*,*)'# NROOT =', nrcnt

       occdet(:) = (0.d0,0.d0)
       gasA(:,:) = 0; gasB(:,:) = 0
       nael(:) = 0; nbel(:) = 0
       detcnt = 1

       do i=1,detcntall
          read (44,1500) occdet(detcnt), nael(detcnt), nbel(detcnt)
!          write(*,1500) occdet(detcnt), nael(detcnt), nbel(detcnt)
          if (abs(occdet(detcnt)).gt.thresh_occdet) then
             read (44,1400) (gasA(j,detcnt),MOJI6,j=1,nael(detcnt))
             read (44,1400) (gasB(j,detcnt),MOJI6,j=1,nbel(detcnt))
             do j=1,nael(detcnt)
                gasA(j,detcnt) = gasA(j,detcnt) *2 -1 !odd
             end do
             do j=1,nbel(detcnt)
                gasB(j,detcnt) = gasB(j,detcnt) *2    !even -> Kramers partner
             end do
!             write (*,1400) (gasA(j,detcnt),MOJI6,j=1,nael(detcnt))
!             write (*,1400) (gasB(j,detcnt),MOJI6,j=1,nbel(detcnt))
             detcnt = detcnt + 1
          else
             read(44,*)
             read(44,*)
          end if
       end do

       detcnt = detcnt -1
       write(*,*)'detcnt=',detcnt,'detcntall=',detcntall

 1400 FORMAT(7X,100(I4,A6))
 1500 FORMAT(7X,2E22.15,2I5)


end subroutine CIread_KRCI2

!============================================================================
subroutine convert_coef_MRDFTCI(a,KRAM,norb,enorb)
!============================================================================
  use DiracOutput
  implicit none

  integer i,j,k
  integer ck
  integer,intent(in) :: norb
  double precision,intent(in) :: enorb(norb)
  character(LEN=1),intent(in) :: KRAM(norb)
  character(LEN=1) KRAMF(NBS_E),KRAMT(NBS_E)
  complex(kind(0d0)),intent(inout) :: a(NBS0,4,norb)

  write(*,*) '* Start convert coef'
  KRAMF(:) = 'T'
  KRAMT(:) = 'T'

  do i=1,norb
     ck=0
     do j=1,NBS_E
        if(ck.eq.0) then
!           if(enorb(i).eq.e_eig(j)) then
           if(abs(enorb(i)-e_eig(j)).lt.1.d-13) then
              if(KRAM(i).eq.'F') then
                 if(KRAMF(j).eq.'T') then
                    call convert_coef_QEDtoMRDFT(a,KRAM(i),norb,i,j)
!!$                    do k=1,NBS_L
!!$                       a(k,1,i) = c_La(k,j)
!!$                       a(k,2,i) = c_Lb(k,j)
!!$                    end do
!!$                    do k=1,NBS_S
!!$                       a(k,3,i) = c_Sa(k,j)
!!$                       a(k,4,i) = c_Sb(k,j)
!!$                    end do
!!$                    if(NBS_L.gt.NBS_S) then
!!$                       do k=NBS_S+1,NBS0
!!$                             a(k,3,i) = (0.d0,0.d0)
!!$                             a(k,4,i) = (0.d0,0.d0)
!!$                       end do
!!$                    else
!!$                       do k=NBS_L+1,NBS0
!!$                             a(k,1,i) = (0.d0,0.d0)
!!$                             a(k,2,i) = (0.d0,0.d0)
!!$                       end do
!!$                    end if
                    KRAMF(j)='O'
                    ck=1
                 end if
              else if(KRAM(i).eq.'T') then !Kramers pair
                 if(KRAMT(j).eq.'T') then
                    call convert_coef_QEDtoMRDFT(a,KRAM(i),norb,i,j)
!!$                    do k=1,NBS_L
!!$                       a(k,1,i) = -dconjg(c_Lb(k,j))
!!$                       a(k,2,i) =  dconjg(c_La(k,j))
!!$                    end do
!!$                    do k=1,NBS_S
!!$                       a(k,3,i) = -dconjg(c_Sb(k,j))
!!$                       a(k,4,i) =  dconjg(c_Sa(k,j))
!!$                    end do
!!$                    if(NBS_L.gt.NBS_S) then
!!$                       do k=NBS_S+1,NBS0
!!$                             a(k,3,i) = (0.d0,0.d0)
!!$                             a(k,4,i) = (0.d0,0.d0)
!!$                       end do
!!$                    else
!!$                       do k=NBS_L+1,NBS0
!!$                             a(k,1,i) = (0.d0,0.d0)
!!$                             a(k,2,i) = (0.d0,0.d0)
!!$                       end do
!!$                    end if
                    KRAMT(j)='O'
                    ck=1
                 end if
              else
                 write(*,*)'error KRAM'
                 stop
              end if
           end if
        end if
     end do
  end do

!  do i=1,norb
!     write(*,'(a5,i3,a1,a1)') 'KRAM(',i,')',KRAM(i)
!  end do
!  do j=1,NBS_E
!     write(*,'(a6,i3,a1,a1)') 'KRAMF(',j,')',KRAMF(j)
!     write(*,'(a6,i3,a1,a1)') 'KRAMT(',j,')',KRAMT(j)
!  end do

  write(*,*) '* End convert coef'

end subroutine convert_coef_MRDFTCI

!============================================================================
subroutine convert_coef_MRDFTCI2(a,KRAM,norb,enorb)
!============================================================================
  use DiracOutput
  implicit none

  integer i,j,k
  integer ck
  integer,intent(in) :: norb
  double precision,intent(in) :: enorb(norb)
  character(LEN=1),intent(in) :: KRAM(norb)
  character(LEN=1) KRAMF(NBS_E),KRAMT(NBS_E)
  complex(kind(0d0)),intent(inout) :: a(NBS0,4,norb)

  write(*,*) '* Start convert coef'
  KRAMF(:) = 'T'
  KRAMT(:) = 'T'

  do i=1,norb
     ck=0
     do j=1,NBS_E
        if(ck.eq.0) then
!           if(enorb(i).eq.e_eig(j)) then
           if(abs(enorb(i)-e_eig(j)).lt.1.d-13) then
              if(KRAM(i).eq.'F') then
                 if(KRAMF(j).eq.'T') then
                    call convert_coef_QEDtoMRDFT(a,KRAM(i),norb,i,j)
                    KRAMF(j)='O'
                    ck=1
                 end if
              else if(KRAM(i).eq.'T') then !Kramers pair
                 if(KRAMT(j).eq.'T') then
                    call convert_coef_QEDtoMRDFT(a,KRAM(i),norb,i,j)
                    KRAMT(j)='O'
                    ck=1
                 end if
              else
                 write(*,*)'error KRAM'
                 stop
              end if
           end if
        end if
     end do
  end do

!  do i=1,norb
!     write(*,'(a5,i3,a1,a1)') 'KRAM(',i,')',KRAM(i)
!  end do
!  do j=1,NBS_E
!     write(*,'(a6,i3,a1,a1)') 'KRAMF(',j,')',KRAMF(j)
!     write(*,'(a6,i3,a1,a1)') 'KRAMT(',j,')',KRAMT(j)
!  end do

  write(*,*) '* End convert coef'

end subroutine convert_coef_MRDFTCI2

!============================================================================
subroutine convert_coef_QEDtoMRDFT(a,kramers,norb,inumnorb,jnumnorb)
!============================================================================
  use DiracOutput
  implicit none

  integer i,k
  integer,intent(in) :: norb,inumnorb, jnumnorb
  character(LEN=1),intent(in) :: kramers
  complex(kind(0d0)),intent(inout) :: a(NBS0,4,norb)

!  write(*,*)'norb,inumnorb,jnumnorb',norb,inumnorb,jnumnorb
  if(kramers.eq.'F') then
     do i=1,NBS_L
        a(i,1,inumnorb) = c_La(i,jnumnorb)
        a(i,2,inumnorb) = c_Lb(i,jnumnorb)
!        write(*,*)'F',i,a(i,1,inumnorb),c_La(i,jnumnorb)
     end do
     do i=1,NBS_S
        a(i,3,inumnorb) = c_Sa(i,jnumnorb)
        a(i,4,inumnorb) = c_Sb(i,jnumnorb)
     end do

  else if(kramers.eq.'T') then
     do i=1,NBS_L
        a(i,1,inumnorb) = -dconjg(c_Lb(i,jnumnorb))
        a(i,2,inumnorb) =  dconjg(c_La(i,jnumnorb))
!        write(*,*)'T',i,a(i,1,inumnorb),c_La(i,jnumnorb)
     end do
     do i=1,NBS_S
        a(i,3,inumnorb) = -dconjg(c_Sb(i,jnumnorb))
        a(i,4,inumnorb) =  dconjg(c_Sa(i,jnumnorb))
     end do
  else
     write(*,*)'error kramers in convert_QEDtoMRDFT'
     stop
  end if

  if(NBS_L.gt.NBS_S) then
     do i=NBS_S+1,NBS0
        do k=3,4
           a(i,k,inumnorb) = (0.d0,0.d0)
        end do
     end do
  else
     do i=NBS_L+1,NBS0
        do k=1,2
           a(i,k,inumnorb) = (0.d0,0.d0)
        end do
     end do
  end if

end subroutine convert_coef_QEDtoMRDFT

!============================================================================
subroutine convert_pg_QEDtoMRDFT(aa,n,xx,yy,zz)
!============================================================================
  use DiracOutput
  implicit none

  integer i,j,k
  integer,intent(out) :: n(NBS0,4)
  double precision,intent(out) :: aa(NBS0,4)
  double precision,intent(out) :: xx(NBS0,4),yy(NBS0,4),zz(NBS0,4)
  double precision :: bb(NBS0,4)

  do i=1,NBS_L
     do k=1,2
        aa(i,k) = aa_L(i)
        xx(i,k) = xx_L(i)
        yy(i,k) = yy_L(i)
        zz(i,k) = zz_L(i)
        call translaten2(nx_L(i),ny_L(i),nz_L(i),n(i,k))
     end do
  end do
  do i=1,NBS_S
     do k=3,4
        aa(i,k) = aa_S(i)
        xx(i,k) = xx_S(i)
        yy(i,k) = yy_S(i)
        zz(i,k) = zz_S(i)
        call translaten2(nx_S(i),ny_S(i),nz_S(i),n(i,k))
     end do
  end do

  if(NBS_L.gt.NBS_S) then
     do i=NBS_S+1,NBS0
        do k=3,4
           xx(i,k) = 0.0D0
           yy(i,k) = 0.0D0
           zz(i,k) = 0.0D0
           aa(i,k) = 1.D0
           n(i,k)  = -1.D0
        end do
     end do
  else
     do i=NBS_L+1,NBS0
        do k=1,2
           xx(i,k) = 0.0D0
           yy(i,k) = 0.0D0
           zz(i,k) = 0.0D0
           aa(i,k) = 1.D0
           n(i,k)  = -1.D0
        end do
     end do
  end if

  ! normalization for MRDFT
  bb(:,:)=1.d0
  call normal(aa,bb,n,NBS0)

!      do j=1,NBS0
!         do k=1,4
!            write(*,*)j,k,bb(j,k)
!         end do
!      end do
!      stop

  do j=1,NBS_E
     do k=1,NBS_L
        c_La(k,j) = c_La(k,j)*bb(k,1)
        c_Lb(k,j) = c_Lb(k,j)*bb(k,2)
     end do
     do k=1,NBS_S
        c_Sa(k,j) = c_Sa(k,j)*bb(k,3)
        c_Sb(k,j) = c_Sb(k,j)*bb(k,4)
     end do
  end do

end subroutine convert_pg_QEDtoMRDFT
