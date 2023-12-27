
!============================================================================
subroutine set_KRCI_coef(coreorbg,coreorbu,coreorb,norbg,norbu,norb,DMrw)
! get norb, c_natu, occup_natu
!============================================================================
  use Precision
  use DiracOutput
  use DefineTypes 
  use DIRRCI_coef
  use param_KRCI
!$ use omp_lib

  implicit none

  integer, intent(out) :: norb
  character(len=1),intent(in) ::  DMrw
  integer,intent(inout) :: coreorb
  integer, intent(out) :: coreorbg,coreorbu,norbg,norbu

  integer(kind=dp) :: detcntall,detcnt
  integer :: neact
  integer,allocatable :: gasA(:,:),gasB(:,:)
  integer,allocatable :: nael(:),nbel(:)
  complex(kind=dp),allocatable :: occdet(:)

  integer(kind=dp) :: ndetK,ndetL
  integer :: n_KL
  integer :: ne3, ne4, tmpck, ne1a, ne1b
  integer :: i,j,k,l
  integer :: denssub1(2),denssub2(2)
  integer :: neactA,neactB
  integer, allocatable :: activeorbA(:,:),activeorbB(:,:)

  integer :: itemp,itemp2,itemp3,itemp4
  integer, allocatable :: eigiwork(:)
  real(kind=dp), allocatable :: rwork(:)
  real(kind=dp), allocatable :: eigvals(:)
  complex(kind=dp), allocatable :: dnstmtx(:,:), eigwork(:) ! this may be changed to complex array
  complex(kind=dp), allocatable :: dnstmtx2(:,:) ! this may be changed to complex array
  complex(kind=dp), allocatable :: tmpdnstmtx(:,:,:)
  character(len=24) today
  double precision :: timea, timeb
  double precision :: tmp2
  complex(kind=dp) :: cdnsty
  integer :: p,aa
  character(LEN=1) :: a,b
  complex(kind=dp), allocatable :: tmpc_natu(:,:,:) ! complex coefficients (4, NBS_L or NBS_S,2*NBS)
!  complex(kind=dp), allocatable :: tmpc_natu2(:,:,:) ! complex coefficients (4, NBS_L or NBS_S,2*NBS)
  integer :: itmp
  integer :: iorbold, iorbnew

  if(SYMM.eq.2) then 
     open(unit=44, file=FIFI4)
        call CIread_KRCI1(coreorb,norb,neact,detcntall)
     close(unit=44)
   else !symmetry (D2h atom)
     open(unit=44, file=FIFI4)
        call CIread_KRCI1atom(coreorbg,coreorbu,coreorb,norbg,norbu,norb,neact,detcntall)
     close(unit=44)
   end if

  allocate(occdet(detcntall))
  allocate(gasA(neact,detcntall),gasB(neact,detcntall))
  allocate(nael(detcntall), nbel(detcntall)) 

  open(unit=44, file=FIFI4)
     call CIread_KRCI2(NROOT,neact,detcntall,detcnt,gasA,gasB,nael,nbel,occdet)
  close(unit=44)
  ! detcntall : number of all determinants
  ! detcnt    : number of thresholded determinants
  ! neact = nael + nbel

!  stop

  allocate(occup_natu(norb))
  allocate(c_natu(4,NBS0,2*NBS))
  allocate(tmpc_natu(4,NBS0,2*NBS))
!  allocate(tmpc_natu2(4,NBS0,2*NBS))
  occup_natu(:) = 0._dp
  c_natu(:,:,:) = (0._dp,0._dp)
  tmpc_natu(:,:,:) = (0._dp,0._dp)
!  tmpc_natu2(:,:,:) = (0._dp,0._dp)

  allocate(eigwork(1 + 6 * norb + 2 * norb**2))
  allocate(eigiwork(3 + 5 * norb**2))
  allocate(eigvals(norb))
  allocate(dnstmtx(norb,norb))
  allocate(dnstmtx2(norb,norb))
  allocate(tmpdnstmtx(omp_get_max_threads(),norb,norb) )

  allocate(activeorbA(norb,detcnt))
  allocate(activeorbB(norb,detcnt))

  !for diagonalization
  itemp  = norb*(norb + 2) !lwork
  itemp2 = 5*norb + 3 !liwork
  itemp3 = 2*norb*norb + 5*norb + 1 !lrwork
  allocate(rwork(itemp3))

!c   
!c    Make density matrix
!c   
   if(DMrw=='w') then
      call FDATE(today)
      write(*,*)'Start setting density matrix :', today
     
     write(*,*) ' '
     write(*,*) 'Finish : read orbital infomation '
     write(*,*) ' '
     write(*,*) 'Make Density Matrix '
     write(*,*) ' '

     activeorbA(:,:) = 0
     activeorbB(:,:) = 0

     !$omp parallel do private(ndetL,j)
     do ndetL=1, detcnt
       do j=1, nael(ndetL)
         activeorbA(gasA(j,ndetL),ndetL) = 1  !activeorbA('odd ',detcnt)
       end do 
       do j=1, nbel(ndetL)
         activeorbB(gasB(j,ndetL),ndetL) = 1  !activeorbB('even',detcnt)
       end do 
     end do ! ndetL -> detcnt 
     !$omp end parallel do

     tmpdnstmtx(:,:,:)  = (0._dp,0._dp)

     !$omp parallel do private(ndetK,ndetL,ne1a,ne1b,ne3,ne4,tmpck,neactA,neactB,&
                           !$omp &n_KL,denssub1,denssub2) schedule(dynamic)
     do ndetL=detcnt,2,-1
    lpK:do ndetK=ndetL-1,1,-1
           n_KL = nael(ndetK) - nael(ndetL)
           if(n_KL.eq.0) then !D_{pq} or D_{barp barq}
              !D_{pq}
              do ne1b=nbel(ndetL),1,-1
                 if(gasB(ne1b,ndetK).ne.gasB(ne1b,ndetL)) goto 105
              end do
              ne4=1; tmpck=1
              neactA = nael(ndetL)
              do ne3=1,norb,2 ! odd
                 if(activeorbA(ne3,ndetK).eq.activeorbA(ne3,ndetL)) then
                    if(activeorbA(ne3,ndetK).eq.1)  ne4 = ne4 + 1 
                 else
                    if(tmpck.eq.3) cycle lpK ! <3>
                    denssub1(tmpck) = ne3  !denssub1->p,q for D_{pq} norb_index
                    denssub2(tmpck) = ne4  !denssub2->p,q for D_{pq} neact_index
                    tmpck = tmpck +1
                 end if
                 if((tmpck.eq.3).and.(ne4.eq.neactA)) exit
              end do
    
              if((activeorbA(denssub1(1),ndetK).eq.1).and.(activeorbA(denssub1(2),ndetL).eq.1)) then 
                 if(mod(denssub2(1)+denssub2(2),2).eq.0) then
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) = &
                         & tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) + dconjg(occdet(ndetK))*occdet(ndetL)
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) = &
                         & tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) + dconjg(occdet(ndetL))*occdet(ndetK)
                 else
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) = &
                           tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) - dconjg(occdet(ndetK))*occdet(ndetL)
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) = &
                         & tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) - dconjg(occdet(ndetL))*occdet(ndetK)
                 end if
              end if
     
              if((activeorbA(denssub1(2),ndetK).eq.1).and.(activeorbA(denssub1(1),ndetL).eq.1)) then
                 if(mod(denssub2(2)+denssub2(1),2).eq.0) then
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) = &
                         & tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) + dconjg(occdet(ndetK))*occdet(ndetL)
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) = &
                         & tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) + dconjg(occdet(ndetL))*occdet(ndetK)
                 else
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) = &
                           tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) - dconjg(occdet(ndetK))*occdet(ndetL)
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) = &
                         & tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) - dconjg(occdet(ndetL))*occdet(ndetK)
                 end if
              end if
    
  105         continue
              !D_{barp barq}
              do ne1a=nael(ndetL),1,-1
                 if(gasA(ne1a,ndetK).ne.gasA(ne1a,ndetL)) cycle lpK
              end do
              ne4=1; tmpck=1
              neactB = nbel(ndetL)
              do ne3=2,norb,2 ! even for loop kramers partners
                 if(activeorbB(ne3,ndetK).eq.activeorbB(ne3,ndetL)) then
                    if(activeorbB(ne3,ndetK).eq.1)  ne4 = ne4 + 1 
                 else
                    if(tmpck.eq.3) cycle lpK ! <3>
                    denssub1(tmpck) = ne3  !denssub1->p,q for D_{pq} norb_index
                    denssub2(tmpck) = ne4  !denssub2->p,q for D_{pq} neact_index
                    tmpck = tmpck +1
                 end if
                 if((tmpck.eq.3).and.(ne4.eq.neactB)) exit
              end do
    
              if((activeorbB(denssub1(1),ndetK).eq.1).and.(activeorbB(denssub1(2),ndetL).eq.1)) then 
                 if(mod(denssub2(1)+denssub2(2),2).eq.0) then
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) = &
                         & tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) + dconjg(occdet(ndetK))*occdet(ndetL)
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) = &
                         & tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) + dconjg(occdet(ndetL))*occdet(ndetK)
                 else
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) = &
                           tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) - dconjg(occdet(ndetK))*occdet(ndetL)
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) = &
                         & tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) - dconjg(occdet(ndetL))*occdet(ndetK)
                 end if
              end if
     
              if((activeorbB(denssub1(2),ndetK).eq.1).and.(activeorbB(denssub1(1),ndetL).eq.1)) then
                 if(mod(denssub2(2)+denssub2(1),2).eq.0) then
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) = &
                         & tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) + dconjg(occdet(ndetK))*occdet(ndetL)
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) = &
                         & tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) + dconjg(occdet(ndetL))*occdet(ndetK)
                 else
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) = &
                           tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) - dconjg(occdet(ndetK))*occdet(ndetL)
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) = &
                         & tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) - dconjg(occdet(ndetL))*occdet(ndetK)
                 end if
              end if
                 
     
           else if(n_KL.eq.-1) then ! D_{barp q}, n_K<n_L, barn_K>barn_L
              ne4=1; tmpck=1
              neactA = nael(ndetL)! e_q|L> 
              do ne3=1,norb,2 ! odd
                 if(activeorbA(ne3,ndetK).eq.activeorbA(ne3,ndetL)) then
                    if(activeorbA(ne3,ndetK).eq.1)  ne4 = ne4 + 1 
                 else
                    if(tmpck.eq.2) cycle lpK ! <2>
                    denssub1(2) = ne3  !denssub1(2)->q for D_{barp q} norb_index
                    denssub2(2) = ne4  !denssub2(2)->q for D_{barp q} neact_index
                    tmpck = tmpck +1
                 end if
                 if((tmpck.eq.2).and.(ne4.eq.neactA)) exit
              end do
    
              ne4=1; tmpck=1
              neactB = nbel(ndetK) ! <K|e_barp^\dagger
              do ne3=2,norb,2! even for loop kramers partners
                 if(activeorbB(ne3,ndetK).eq.activeorbB(ne3,ndetL)) then
                    if(activeorbB(ne3,ndetK).eq.1)  ne4 = ne4 + 1 
                 else
                    if(tmpck.eq.2) cycle lpK  ! <2>
                    denssub1(1) = ne3  !denssub1(1)->barp for D_{barp q} norb_index
                    denssub2(1) = ne4  !denssub2(1)->barp for D_{barp q} neact_index
                    tmpck = tmpck +1
                 end if
                 if((tmpck.eq.2).and.(ne4.eq.neactB)) exit
              end do
    
              if((activeorbB(denssub1(1),ndetK).eq.1).and.(activeorbA(denssub1(2),ndetL).eq.1)) then 
                 if(mod(denssub2(1)+denssub2(2),2).eq.0) then
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) = &
                         & tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) + dconjg(occdet(ndetK))*occdet(ndetL)
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) = &
                         & tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) + dconjg(occdet(ndetL))*occdet(ndetK)
                 else
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) = &
                           tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) - dconjg(occdet(ndetK))*occdet(ndetL)
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) = &
                         & tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) - dconjg(occdet(ndetL))*occdet(ndetK)
                 end if
              end if
    
    
     
           else if(n_KL.eq.1) then ! D_{p barq}, n_K>n_L, barn_K<barn_L
              ne4=1; tmpck=1
              neactA = nael(ndetK) ! <K|e_p^\dagger
              do ne3=1,norb,2 ! odd
                 if(activeorbA(ne3,ndetK).eq.activeorbA(ne3,ndetL)) then
                    if(activeorbA(ne3,ndetK).eq.1)  ne4 = ne4 + 1 
                 else
                    if(tmpck.eq.2) cycle lpK ! <2>
                    denssub1(1) = ne3  !denssub1(1)->p for D_{barp q} norb_index
                    denssub2(1) = ne4  !denssub2(1)->p for D_{barp q} neact_index
                    tmpck = tmpck +1
                 end if
                 if((tmpck.eq.2).and.(ne4.eq.neactA)) exit
              end do
     
              ne4=1; tmpck=1
              neactB = nbel(ndetL)! e_barq|L>
              do ne3=2,norb,2! even for loop kramers partners
                 if(activeorbB(ne3,ndetK).eq.activeorbB(ne3,ndetL)) then
                    if(activeorbB(ne3,ndetK).eq.1)  ne4 = ne4 + 1 
                 else
                    if(tmpck.eq.2) cycle lpK  ! <2>
                    denssub1(2) = ne3  !denssub1(2)->barq for D_{barp q} norb_index
                    denssub2(2) = ne4  !denssub2(2)->barq for D_{barp q} neact_index
                    tmpck = tmpck +1
                 end if
                 if((tmpck.eq.2).and.(ne4.eq.neactB)) exit
              end do

              if((activeorbA(denssub1(1),ndetK).eq.1).and.(activeorbB(denssub1(2),ndetL).eq.1)) then 
                 if(mod(denssub2(1)+denssub2(2),2).eq.0) then
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) = &
                         & tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) + dconjg(occdet(ndetK))*occdet(ndetL)
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) = &
                         & tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) + dconjg(occdet(ndetL))*occdet(ndetK)
                 else
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) = &
                           tmpdnstmtx(omp_get_thread_num()+1,denssub1(1),denssub1(2)) - dconjg(occdet(ndetK))*occdet(ndetL)
                    tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) = &
                         & tmpdnstmtx(omp_get_thread_num()+1,denssub1(2),denssub1(1)) - dconjg(occdet(ndetL))*occdet(ndetK)
                 end if
              end if
    
!           else 
!              exit
           end if
        end do lpK
     end do
     !$omp end parallel do
 
     write(*,*)'ck3'
  
 
     !$omp parallel do private(ndetL,ne1a,ne1b)
     do ndetL=1,detcnt
        do ne1a=1,nael(ndetL)
           tmpdnstmtx(omp_get_thread_num()+1,gasA(ne1a,ndetL),gasA(ne1a,ndetL)) = &
                  tmpdnstmtx(omp_get_thread_num()+1,gasA(ne1a,ndetL),gasA(ne1a,ndetL)) + dconjg(occdet(ndetL))*occdet(ndetL)
        end do
        do ne1b=1,nbel(ndetL)
           tmpdnstmtx(omp_get_thread_num()+1,gasB(ne1b,ndetL),gasB(ne1b,ndetL)) = &
                  tmpdnstmtx(omp_get_thread_num()+1,gasB(ne1b,ndetL),gasB(ne1b,ndetL)) + dconjg(occdet(ndetL))*occdet(ndetL)
        end do
     end do
     !$omp end parallel do

     write(*,*) 'Make Density Matrix 4'
      dnstmtx(:,:)  = (0.0d0,0.0d0)
      !$omp parallel do private(k,i,j)
      do j=1, norb
        do i=1, norb
          do k = 1, omp_get_max_threads()
               dnstmtx(i,j) = dnstmtx(i,j) + tmpdnstmtx(k,i,j)
          end do
!          write(*,*)j,i,dnstmtx(i,j)
        end do
      end do
      !$omp end parallel do

       call cpu_time(timeb)

      write(*,*) ' '
      write(*,*) 'Finish : Density Matrix '
      write(*,*) 'CPU time is ', timeb - timea
!      write(*,*) 'Products of determinant :', itmpk
      write(*,*) 'Make Natural Orbital '
      write(*,*) ' '
       
         cdnsty = (0._dp,0._dp)
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

          call FDATE(today)
          write(*,*)'Start diagonalizing density matrix :', today

!          call dsyevd('V', 'U', norb, dnstmtx, norb, eigvals, &
!                      eigwork, itemp, eigiwork, itemp2, itemp3) 
          call zheevd('V', 'U', norb, dnstmtx, norb, eigvals, &
                      eigwork, itemp, rwork, itemp3, eigiwork, itemp2, itemp4) 

         call cpu_time(timea)

         write(*,*) ' '
         write(*,*) 'Finish : Natural Orbital '
         write(*,*) 'CPU time is ', timea - timeb
         write(*,*) ' '

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
                  end do 
               end do 
               if(i.eq.j) write(100,*) tmp2, eigvals(i) - tmp2, (eigvals(i) - tmp2)/tmp2
               if(abs(tmp2).gt.10.0d-15) write(14,*) i,j,tmp2
             end do 
          end do 
!        stop

     ! generate Kramers partners' coef
     if(SYMM.eq.2) then 
        do i=1,norb
          call index_from_Qmat(i+coreorb,p,a)
          call pm12(a,aa)
          call copy_DiracOutput_cp(p,aa,tmpc_natu(:,:,i))
        end do
     else ! D2h symmetry
        do i=1,norbg
          call index_from_Qmat(i+coreorbg,p,a)
          call pm12(a,aa)
          call copy_DiracOutput_cp(p,aa,tmpc_natu(:,:,i))
        end do
        itmp = NBS_Eg*2 + coreorbu
!        write(*,*)'NBS_Eg, coreorbg, coreorbu, norbg',NBS_Eg, coreorbg, coreorbu, norbg
!        stop
        do i=1,norbu
          call index_from_Qmat(i+itmp,p,a)
          call pm12(a,aa)
          call copy_DiracOutput_cp(p,aa,tmpc_natu(:,:,i+norbg))
        end do
     end if

!     open(unit=53,file='iorbconv.dat')
!        do l=1,norb/2
!           read(53,*) iorbnew, iorbold
!           write(*,*) iorbnew, iorbold
!           tmpc_natu2(:,:,2*iorbnew-1) = tmpc_natu(:,:,2*iorbold-1)
!           tmpc_natu2(:,:,2*iorbnew) = tmpc_natu(:,:,2*iorbold)
!        end do
!     close(unit=53)
   
     do l=1,norb
        do k=1,norb
          do j=1,NBS_L
             c_natu(1,j,k) = c_natu(1,j,k) + dconjg(dnstmtx(l,k)) * tmpc_natu(1,j,l) ! Is the "dconjg" needed?
             c_natu(2,j,k) = c_natu(2,j,k) + dconjg(dnstmtx(l,k)) * tmpc_natu(2,j,l)
!             c_natu(1,j,k) = c_natu(1,j,k) + dconjg(dnstmtx(l,k)) * tmpc_natu2(1,j,l) ! Is the "dconjg" needed?
!             c_natu(2,j,k) = c_natu(2,j,k) + dconjg(dnstmtx(l,k)) * tmpc_natu2(2,j,l)
          end do
          do j=1,NBS_S
             c_natu(3,j,k) = c_natu(3,j,k) + dconjg(dnstmtx(l,k)) * tmpc_natu(3,j,l)
             c_natu(4,j,k) = c_natu(4,j,k) + dconjg(dnstmtx(l,k)) * tmpc_natu(4,j,l)
!             c_natu(3,j,k) = c_natu(3,j,k) + dconjg(dnstmtx(l,k)) * tmpc_natu2(3,j,l)
!             c_natu(4,j,k) = c_natu(4,j,k) + dconjg(dnstmtx(l,k)) * tmpc_natu2(4,j,l)
          end do
        end do
     end do

  else if(DMrw=='r') then
     write(*,*)'# read density matrix'
  else
     write(*,*)'DMrw should be "r" or "w".'
  end if

  write(*,*)'# NBS_Eg, coreorbg, coreorbu, norbg',NBS_Eg, coreorbg, coreorbu, norbg

  open(unit=50,file='c_natu.dat')
  open(unit=52,file='occup_natu.dat')
  if(DMrw=='w') then
     do k=1,norb
        do j=1,NBS_L
           do i=1,2
             write(50,'(2e22.14)') c_natu(i,j,k)
           end do 
        end do 
        do j=1,NBS_S
           do i=3,4,1
             write(50,'(2e22.14)') c_natu(i,j,k)
           end do 
        end do 
     end do  
     do k=1,norb
        occup_natu(k) = eigvals(k)
        write(52,'(e22.14)') occup_natu(k)
     end do
  else if(DMrw=='r') then
     do k=1,norb
        do j=1,NBS_L
           do i=1,2
             read(50,'(2e22.14)') c_natu(i,j,k)
           end do 
        end do 
        do j=1,NBS_S
           do i=3,4,1
             read(50,'(2e22.14)') c_natu(i,j,k)
           end do 
        end do 
     end do  
     do k=1,norb
        read(52,'(e22.14)') occup_natu(k)
     end do
  else 
     write(*,*)'error DMrw'
     stop
  end if
  close(unit=50)
  close(unit=52)

!  stop

  deallocate(eigwork)
  deallocate(eigiwork)
  deallocate(eigvals)
  deallocate(occdet)
  deallocate(gasA,gasB)
  deallocate(nael,nbel)
  deallocate(dnstmtx)
  deallocate(dnstmtx2)
  deallocate(tmpc_natu)
!  deallocate(tmpc_natu2)
  deallocate(tmpdnstmtx)
  deallocate(activeorbA)
  deallocate(activeorbB)

  return

end subroutine set_KRCI_coef

