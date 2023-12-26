! Last Change:07-Nov-2014.
!==============================================================================
subroutine readVectors_KRCI(bnum_L,bnum_S,bnum_LT,bnum_ST,bnum_LTA,bnum_STA,cpg_L,cpg_S,n,cont)
!==============================================================================
  use DiracOutput

  implicit none

  character(len=3) MOJI,L,S
  character(len=4) Angu
  character(len=5) DUMMY
  character(len=6) VECTOR
  integer,intent(in):: bnum_L(100),bnum_S(100)
  integer,intent(in):: bnum_LT,bnum_ST
  integer,intent(in):: bnum_LTA(NAT),bnum_STA(NAT)
  integer i,j,ip,ie,kk,k0,bnum,num,inat,m(NAT),m_total
  integer num1, num2
  integer nt,Pns(NBS0,4),Ens(NBS0,4)  !---------------checker of n(NBS0,4) 
  integer tmpPns(NAT,NBS0,4),tmpEns(NAT,NBS0,4)
  integer tmpcont
  integer mcont(NAT)
  integer,intent(in):: cont(NBS00,NAT)
  integer,intent(in):: n(NBS0,4)
  double precision eig
  double precision,intent(in):: cpg_L(NBS_L),cpg_S(NBS_S)
  double precision renorm
  complex(kind(0d0)) tmpA, tmpB
  complex(kind(0d0)) tmpd_La(NAT,NBS_L), tmpd_Lb(NAT,NBS_L),tmpd_Sa(NAT,NBS_S), tmpd_Sb(NAT,NBS_S)
  complex(kind(0d0)) tmpc_La(NAT,NBS_L), tmpc_Lb(NAT,NBS_L),tmpc_Sa(NAT,NBS_S), tmpc_Sb(NAT,NBS_S)
  integer tmpNBS


  read(42,'(a6)') VECTOR
  if (VECTOR.eq.'VECTOR') then
   write(*,*) 'start reading Vectors.txt'
  end if

 !----------------single atom--------------------------------------------------
  if (SYMM.eq.0) then
    write(*,*)"SYMM = 0 is not available."
    write(*,*)"SYMM = 0 is symmetry calculation for atom"
    stop
  !-------------------------molecule-------------------------------------------
  else    !--------SYMM=1 => molecule using D2h symmmetry calculation, SYMM=2 => the other calculation, SYMM=3 => professional option
   ip=0 ! counter of positronic orbitals
   ie=0 ! counter of electronic orbitals

   if(KBAL.eq.0) tmpNBS = NBS_TOTAL
   if(KBAL.eq.1) tmpNBS = NBS*2
   write(*,*)'tmpNBS=',tmpNBS

   do kk=1,tmpNBS
!    read (42,'(A3,5x,E18.10)') MOJI,eig
    read (42,'(A3,5x,E22.14)') MOJI,eig
!    write(*, '(A3,5x,E22.14)') MOJI,eig
 !-------read PIV---------------------------------------------------------------------------
     if (MOJI.eq.'PIV') then
      ip=ip+1
      p_eig(ip) = eig

      do inat=1,NAT
        m(inat)=0  ! counter of (inat)th atom's primitive gaussian number
        mcont(inat)=0 ! counter of (inat)th atom's contracted primitive gaussian number
      end do
      read(42,'(A5)')DUMMY
 !-------------read large part of PIV--------------------------------------------------------
      do i=1,bnum_LT
      
        if(SYMM.eq.1) then
           read (42,'(i4,a3,4x,i3,a4,2x,4e22.14)') bnum, L, num, Angu, tmpA, tmpB
        else if(SYMM.eq.2) then
           read (42,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, L, num, Angu, tmpA, tmpB
!           write (*,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, L, num, Angu, tmpA, tmpB
        else if(SYMM.eq.0) then
           read (42,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, L, num, Angu, tmpA, tmpB
        else if(SYMM.eq.3) then
           read (42,'(i4,a3,2x,i2,i3,a4,2x,4e22.14)') bnum, L, num1, num2, Angu, tmpA, tmpB
           num = num1 + num2 -1
        end if

         if(bnum.ne.i) then
            write(*,*)'error read vectors PIV Large'
            write(*,*)'i=',i,'bnum=',bnum
            stop
         end if

        call angular(nt, Angu)

        do inat=1,NAT
          if(num.eq.inat) then
            mcont(inat) = mcont(inat) + 1
            tmpcont = cont(mcont(inat),inat)
            do j=1,tmpcont
              m(inat) = m(inat)+1 
              tmpPns(inat,m(inat),1) = nt
              tmpPns(inat,m(inat),2) = nt
              tmpd_La(inat,m(inat)) = tmpA
              tmpd_Lb(inat,m(inat)) = tmpB
	    end do
          end if
        end do

      end do

         m_total=0
       do inat=1,NAT
          do i=1,m(inat)
           Pns(i+m_total,1) = tmpPns(inat,i,1)
           Pns(i+m_total,2) = tmpPns(inat,i,2)
           d_La(i+m_total,ip) = tmpd_La(inat,i)*cpg_L(i+m_total)*renorm(nx_L(i+m_total),ny_L(i+m_total),nz_L(i+m_total))   !---------coefficient of La
!           d_La(i+m_total,ip) = tmpd_La(inat,i)*cpg_L(i+m_total)   !---------coefficient of La
           d_Lb(i+m_total,ip) = tmpd_Lb(inat,i)*cpg_L(i+m_total)*renorm(nx_L(i+m_total),ny_L(i+m_total),nz_L(i+m_total))   !---------coefficient of La
!           d_Lb(i+m_total,ip) = tmpd_Lb(inat,i)*cpg_L(i+m_total)   !---------coefficient of Lb
          end do
         m_total=m_total+m(inat)
       end do
 !--------------read small part of PIV---------------------------------------------------------
      do inat=1,NAT
        m(inat)=0  ! counter of (inat)th atom's primitive gaussian number
        mcont(inat)=0 ! counter of (inat)th atom's contracted primitive gaussian number
      end do

      do i=1,bnum_ST
        if(SYMM.eq.1) then
          read (42,'(i4,a3,4x,i3,a4,2x,4e22.14)') bnum, S, num, Angu, tmpA, tmpB
!          write(*,'(i4,a3,4x,i3,a4,2x,4e22.14)') bnum, S, num, Angu, tmpA, tmpB
        else if(SYMM.eq.2) then
          read (42,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, S, num, Angu, tmpA, tmpB
!          write (*,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, S, num, Angu, tmpA, tmpB
        else if(SYMM.eq.0) then
          read (42,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, S, num, Angu, tmpA, tmpB
        else if(SYMM.eq.3) then
           read (42,'(i4,a3,2x,i2,i3,a4,2x,4e22.14)') bnum, S, num1, num2, Angu, tmpA, tmpB
           num = num1 + num2 -1
        end if
         if(bnum.ne.i+bnum_LT) then
          write(*,*)'error read vectors PIV Small'
          write(*,*)'i+bnum_LT=',i+bnum_LT,'bnum=',bnum
          stop
         end if
        call angular(nt, Angu)
        do inat=1,NAT
          if(num.eq.inat) then
            mcont(inat) = mcont(inat) + 1
            tmpcont = cont(mcont(inat)+bnum_LTA(inat),inat)
            do j=1,tmpcont
              m(inat) = m(inat)+1 
              tmpPns(inat,m(inat),3) = nt
              tmpPns(inat,m(inat),4) = nt
              tmpd_Sa(inat,m(inat)) = tmpA
              tmpd_Sb(inat,m(inat)) = tmpB
	    end do
          end if
        end do
      end do

         m_total=0
       do inat=1,NAT
          do i=1,m(inat)
           Pns(i+m_total,3) = tmpPns(inat,i,3)
           Pns(i+m_total,4) = tmpPns(inat,i,4)
           d_Sa(i+m_total,ip) = tmpd_Sa(inat,i)*cpg_S(i+m_total)*renorm(nx_S(i+m_total),ny_S(i+m_total),nz_S(i+m_total))   !---------coefficient of La
!           d_Sa(i+m_total,ip) = tmpd_Sa(inat,i)*cpg_S(i+m_total)   !---------coefficient of La
           d_Sb(i+m_total,ip) = tmpd_Sb(inat,i)*cpg_S(i+m_total)*renorm(nx_S(i+m_total),ny_S(i+m_total),nz_S(i+m_total))   !---------coefficient of Lb
!           d_Sb(i+m_total,ip) = tmpd_Sb(inat,i)*cpg_S(i+m_total)   !---------coefficient of Lb
          end do
         m_total=m_total+m(inat)
       end do
 !-------read EIV---------------------------------------------------------------------------
     else if (MOJI.eq.'EIV') then
      ie=ie+1
      e_eig(ie) = eig

      do inat=1,NAT
        m(inat)=0  ! counter of (inat)th atom's primitive gaussian number
        mcont(inat)=0 ! counter of (inat)th atom's contracted primitive gaussian number
      end do

      read(42,'(A5)')DUMMY
 !---------------read large part of EIV--------------------------------------------

      do i=1,bnum_LT
        if(SYMM.eq.1) then
          read (42,'(i4,a3,4x,i3,a4,2x,4e22.14)') bnum, L, num, Angu, tmpA, tmpB
        else if(SYMM.eq.2) then
          read (42,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, L, num, Angu, tmpA, tmpB
!        write (*,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, L, num, Angu, tmpA, tmpB
        else if(SYMM.eq.0) then
          read (42,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, L, num, Angu, tmpA, tmpB
        else if(SYMM.eq.3) then
           read (42,'(i4,a3,2x,i2,i3,a4,2x,4e22.14)') bnum, L, num1, num2, Angu, tmpA, tmpB
           num = num1 + num2 -1
        end if
         if(bnum.ne.i) then
          write(*,*)'error read vectors EIV Large'
          write(*,*)'i=',i,'bnum=',bnum
          stop
         end if
        call angular(nt, Angu)

        do inat=1,NAT
          if(num.eq.inat) then
            mcont(inat) = mcont(inat) + 1
            tmpcont = cont(mcont(inat),inat)
            do j=1,tmpcont
              m(inat) = m(inat)+1 
              tmpEns(inat,m(inat),1) = nt
              tmpEns(inat,m(inat),2) = nt
              tmpc_La(inat,m(inat)) = tmpA
              tmpc_Lb(inat,m(inat)) = tmpB
	    end do
          end if
        end do
      end do

         m_total=0
       do inat=1,NAT
          do i=1,m(inat)
           Ens(i+m_total,1) = tmpEns(inat,i,1)
           Ens(i+m_total,2) = tmpEns(inat,i,2)
           c_La(i+m_total,ie) = tmpc_La(inat,i)*cpg_L(i+m_total)*renorm(nx_L(i+m_total),ny_L(i+m_total),nz_L(i+m_total))   !---------coefficient of La
!           c_La(i+m_total,ie) = tmpc_La(inat,i)*cpg_L(i+m_total)   !---------coefficient of La
           c_Lb(i+m_total,ie) = tmpc_Lb(inat,i)*cpg_L(i+m_total)*renorm(nx_L(i+m_total),ny_L(i+m_total),nz_L(i+m_total))   !---------coefficient of Lb
!           c_Lb(i+m_total,ie) = tmpc_Lb(inat,i)*cpg_L(i+m_total)   !---------coefficient of Lb
          end do
         m_total=m_total+m(inat)
       end do

 !--------------read small part of EIV---------------------------------------------------------
      do inat=1,NAT
        m(inat)=0  ! counter of (inat)th atom's primitive gaussian number
        mcont(inat)=0 ! counter of (inat)th atom's contracted primitive gaussian number
      end do

      do i=1,bnum_ST
        if(SYMM.eq.1) then
        read (42,'(i4,a3,4x,i3,a4,2x,4e22.14)') bnum, S, num, Angu, tmpA, tmpB
        else if(SYMM.eq.2) then
        read (42,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, S, num, Angu, tmpA, tmpB
        else if(SYMM.eq.0) then
        read (42,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, S, num, Angu, tmpA, tmpB
        else if(SYMM.eq.3) then
           read (42,'(i4,a3,2x,i2,i3,a4,2x,4e22.14)') bnum, S, num1, num2, Angu, tmpA, tmpB
           num = num1 + num2 -1
        end if
         if(bnum.ne.i+bnum_LT) then
          write(*,*)'error read vectors EIV Small'
          write(*,*)'i+bnum_LT=',i+bnum_LT,'bnum=',bnum
          stop
         end if
        call angular(nt, Angu)

        do inat=1,NAT
          if(num.eq.inat) then
            mcont(inat) = mcont(inat) + 1
            tmpcont = cont(mcont(inat)+bnum_LTA(inat),inat)
            do j=1,tmpcont
              m(inat) = m(inat)+1 
              tmpEns(inat,m(inat),3) = nt
              tmpEns(inat,m(inat),4) = nt
              tmpc_Sa(inat,m(inat)) = tmpA
              tmpc_Sb(inat,m(inat)) = tmpB
	    end do
          end if
        end do !inat
      end do

         m_total=0
       do inat=1,NAT
          do i=1,m(inat)
           Ens(i+m_total,3) = tmpEns(inat,i,3)
           Ens(i+m_total,4) = tmpEns(inat,i,4)
           c_Sa(i+m_total,ie) = tmpc_Sa(inat,i)*cpg_S(i+m_total)*renorm(nx_S(i+m_total),ny_S(i+m_total),nz_S(i+m_total))   !---------coefficient of La
!           c_Sa(i+m_total,ie) = tmpc_Sa(inat,i)*cpg_S(i+m_total)   !---------coefficient of La
           c_Sb(i+m_total,ie) = tmpc_Sb(inat,i)*cpg_S(i+m_total)*renorm(nx_S(i+m_total),ny_S(i+m_total),nz_S(i+m_total))   !---------coefficient of Lb
!           c_Sb(i+m_total,ie) = tmpc_Sb(inat,i)*cpg_S(i+m_total)   !---------coefficient of Lb
          end do
         m_total=m_total+m(inat)
       end do

     end if
   end do
  end if  !SYMM

!--------------error check-------------------------
  do i=1,NBS_L
     do k0=1,2
        if(Pns(i,k0).ne.n(i,k0)) then
           write(*,*)'error Pns(',i,',',k0,')=',Pns(i,k0),'n=',n(i,k0)
           stop
        end if
        if(Ens(i,k0).ne.n(i,k0)) then
           write(*,*)'error Ens(',i,',',k0,')=',Ens(i,k0),'n=',n(i,k0)
           stop
        end if
     end do
  end do

  do i=1,NBS_S
     do k0=3,4
        if(Pns(i,k0).ne.n(i,k0)) then
           write(*,*)'error Pns(',i,',',k0,')=',Pns(i,k0),'n=',n(i,k0)
           stop
        end if
        if(Ens(i,k0).ne.n(i,k0)) then
           write(*,*)'error Ens(',i,',',k0,')=',Ens(i,k0),'n=',n(i,k0)
           stop
        end if
     end do
  end do
!-------------arrange orbitals to eigenvalue--------
!  call arrangeEigenvalueE(e_eig)
!  call arrangeEigenvalueP(p_eig)
!  call arrangeEigenvalue(e_eig,p_eig)
  write(*,*)'finish reading vectors.txt'

end subroutine readVectors_KRCI

