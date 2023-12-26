! Last Change:08-Jun-2014.
!!$!============================================================
!!$function intEeff_EDM_ele_pt_Qmat(nn,mm,pp,qq)
!!$! nn,mm,pp,qq : labels for molecular orbitals (KP) . 1~2*NBS
!!$!============================================================
!!$  use Precision
!!$  use DefineTypes
!!$  use DiracOutput
!!$  use Constants
!!$  use PTcoef, only : c_psi
!!$  implicit none
!!$
!!$  integer,intent(in) :: nn,mm,pp,qq
!!$  integer :: i,j,k,l ! p.g. indice
!!$  integer :: NLS !function
!!$  integer :: NL,NS
!!$  integer :: a,b ! spinor indices 1,2->L, 3,4->S
!!$
!!$  real(kind=dp) :: inttwoelegrad_pg !function
!!$  type(primitive_gaussian) :: pg
!!$  complex(kind=dp) :: tmp, tmptot
!!$  complex(kind=dp) :: intEeff_EDM_ele_pt_Qmat
!!$
!!$  call set_GammaMatrix
!!$  call copy_DiracOutput_pg(pg)
!!$  NL=NBS_L; NS=NBS_S
!!$
!!$    tmptot=(0._dp,0._dp)
!!$    tmp=(0._dp,0._dp)
!!$    !m=1
!!$    !LLLL
!!$    a=1; b=1
!!$    do i=1,NL
!!$       do j=1,NL
!!$          do k=1,NL
!!$             do l=1,NL
!!$                tmpij11 = conjg(c_psi(1,i,nn))*c_psi(1,j,mm)
!!$                tmpij12 = conjg(c_psi(1,i,nn))*c_psi(2,j,mm)
!!$                tmpij21 = conjg(c_psi(2,i,nn))*c_psi(1,j,mm)
!!$                tmpij22 = conjg(c_psi(2,i,nn))*c_psi(2,j,mm)
!!$                tmpkl11 = conjg(c_psi(1,k,pp))*c_psi(1,l,qq)
!!$                tmpkl22 = conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
!!$                tmpkl1122 = tmpkl11 + tmpkl22
!!$                !m=1
!!$!                tmp =  conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)& 
!!$!                     &+conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
!!$!                     &+conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!!$!                     &+conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
!!$                tmp = (tmpij21+tmpij12)*tmpkl1122
!!$                tmptot = tmptot + tmp*inttwoelegrad_pg(1,i,j,k,l,1,1,1,1,NL,NS,pg)
!!$                !m=2
!!$!                tmp =  IU*conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!!$!                     &+IU*conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
!!$!                     &-IU*conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!!$!                     &-IU*conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
!!$                tmp = IU*(tmpij21-tmpij12)*tmpkl1122
!!$                tmptot = tmptot + tmp*inttwoelegrad_pg(2,i,j,k,l,1,1,1,1,NL,NS,pg)
!!$                !m=3
!!$!                tmp =  conjg(c_psi(1,i,nn))*c_psi(1,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!!$!                     &+conjg(c_psi(1,i,nn))*c_psi(1,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
!!$!                     &-conjg(c_psi(2,i,nn))*c_psi(2,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!!$!                     &-conjg(c_psi(2,i,nn))*c_psi(2,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
!!$                tmp = (tmpij11-tmpij22)*tmpkl1122
!!$                tmptot = tmptot + tmp*inttwoelegrad_pg(3,i,j,k,l,1,1,1,1,NL,NS,pg)
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    !LLSS
!!$    a=1; b=3
!!$    do i=1,NL
!!$       do j=1,NL
!!$          do k=1,NS
!!$             do l=1,NS
!!$                tmpij11 = conjg(c_psi(1,i,nn))*c_psi(1,j,mm)
!!$                tmpij12 = conjg(c_psi(1,i,nn))*c_psi(2,j,mm)
!!$                tmpij21 = conjg(c_psi(2,i,nn))*c_psi(1,j,mm)
!!$                tmpij22 = conjg(c_psi(2,i,nn))*c_psi(2,j,mm)
!!$                tmpkl33 = conjg(c_psi(3,k,pp))*c_psi(3,l,qq)
!!$                tmpkl44 = conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
!!$                tmpkl3344 = tmpkl33 + tmpkl44
!!$                !m=1
!!$!                tmp =  conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!!$!                     &+conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
!!$!                     &+conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!!$!                     &+conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
!!$                tmp = (tmpij21 + tmpij12)*tmpkl3344
!!$                tmptot = tmptot + tmp*inttwoelegrad_pg(1,i,j,k,l,1,1,3,3,NL,NS,pg)
!!$                !m=2
!!$!                tmp =  IU*conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!!$!                     &+IU*conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
!!$!                     &-IU*conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!!$!                     &-IU*conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
!!$                tmp = IU*(tmpij21 - tmpij12)*tmpkl3344
!!$                tmptot = tmptot + tmp*inttwoelegrad_pg(2,i,j,k,l,1,1,3,3,NL,NS,pg)
!!$                !m=3
!!$!                tmp =  conjg(c_psi(1,i,nn))*c_psi(1,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!!$!                     &+conjg(c_psi(1,i,nn))*c_psi(1,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
!!$!                     &-conjg(c_psi(2,i,nn))*c_psi(2,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!!$!                     &-conjg(c_psi(2,i,nn))*c_psi(2,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
!!$                tmp = (tmpij11 - tmpij22)*tmpkl3344
!!$                tmptot = tmptot + tmp*inttwoelegrad_pg(3,i,j,k,l,1,1,3,3,NL,NS,pg)
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    !SSLL
!!$    a=3; b=1
!!$    do i=1,NS
!!$       do j=1,NS
!!$          do k=1,NL
!!$             do l=1,NL
!!$                tmpij33 = conjg(c_psi(3,i,nn))*c_psi(3,j,mm)
!!$                tmpij34 = conjg(c_psi(3,i,nn))*c_psi(4,j,mm)
!!$                tmpij43 = conjg(c_psi(4,i,nn))*c_psi(3,j,mm)
!!$                tmpij44 = conjg(c_psi(4,i,nn))*c_psi(4,j,mm)
!!$                tmpkl11 = conjg(c_psi(1,k,pp))*c_psi(1,l,qq)
!!$                tmpkl22 = conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
!!$                tmpkl1122 = tmpkl11 + tmpkl22
!!$                !m=1
!!$!                tmp = -conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!!$!                     &-conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
!!$!                     &-conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!!$!                     &-conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
!!$                tmp = -(tmpij43 + tmpij34)*tmpkl1122
!!$                tmptot = tmptot +  tmp*inttwoelegrad_pg(1,i,j,k,l,3,3,1,1,NL,NS,pg)
!!$                !m=2
!!$!                tmp = -IU*conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!!$!                     &-IU*conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
!!$!                     &+IU*conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!!$!                     &+IU*conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
!!$                tmp = -IU*(tmpij43 - tmpij34)*tmpkl1122
!!$                tmptot = tmptot + tmp*inttwoelegrad_pg(2,i,j,k,l,3,3,1,1,NL,NS,pg)
!!$                !m=3
!!$!                tmp = -conjg(c_psi(3,i,nn))*c_psi(3,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!!$!                     &-conjg(c_psi(3,i,nn))*c_psi(3,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
!!$!                     &+conjg(c_psi(4,i,nn))*c_psi(4,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!!$!                     &+conjg(c_psi(4,i,nn))*c_psi(4,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
!!$                tmp = -(tmpij33 - tmpij44)*tmpkl1122
!!$                tmptot = tmptot + tmp*inttwoelegrad_pg(3,i,j,k,l,3,3,1,1,NL,NS,pg)
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    !SSSS
!!$    a=3; b=3
!!$    do i=1,NS
!!$       do j=1,NS
!!$          do k=1,NS
!!$             do l=1,NS
!!$                tmpij33 = conjg(c_psi(3,i,nn))*c_psi(3,j,mm)
!!$                tmpij34 = conjg(c_psi(3,i,nn))*c_psi(4,j,mm)
!!$                tmpij43 = conjg(c_psi(4,i,nn))*c_psi(3,j,mm)
!!$                tmpij44 = conjg(c_psi(4,i,nn))*c_psi(4,j,mm)
!!$                tmpkl33 = conjg(c_psi(3,k,pp))*c_psi(3,l,qq)
!!$                tmpkl44 = conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
!!$                tmpkl3344 = tmpkl33 + tmpkl44
!!$                !m=1
!!$!                tmp = -conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!!$!                     &-conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
!!$!                     &-conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!!$!                     &-conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
!!$                tmp = -(tmpij43+tmpij34)*tmpkl3344
!!$                tmptot = tmptot + tmp*inttwoelegrad_pg(1,i,j,k,l,3,3,3,3,NL,NS,pg)
!!$                !m=2
!!$!                tmp = -IU*conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!!$!                     &-IU*conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
!!$!                     &+IU*conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!!$!                     &+IU*conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
!!$                tmp = -IU*(tmpij43-tmpij43)*tmpkl3344
!!$                tmptot = tmptot + tmp*inttwoelegrad_pg(2,i,j,k,l,3,3,3,3,NL,NS,pg)
!!$                !m=3
!!$!                tmp = -conjg(c_psi(3,i,nn))*c_psi(3,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!!$!                     &-conjg(c_psi(3,i,nn))*c_psi(3,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
!!$!                     &+conjg(c_psi(4,i,nn))*c_psi(4,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!!$!                     &+conjg(c_psi(4,i,nn))*c_psi(4,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
!!$                tmp = -(tmpij33-tmpij44)*tmpkl3344
!!$                tmptot = tmptot + tmp*inttwoelegrad_pg(3,i,j,k,l,3,3,3,3,NL,NS,pg)
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    intEeff_EDM_ele_pt_Qmat = Ze*tmptot
!!$
!!$    return
!!$end function intEeff_EDM_ele_pt_Qmat
!!$
!============================================================
subroutine coef_intEeff_EDM_ele(nn,mm,pp,qq,i,j,k,l,a,b,coef)
! nn,mm,pp,qq : labels for molecular orbitals (KP) . 1~2*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use PTcoef, only : c_psi
  implicit none

  integer,intent(in) :: nn,mm,pp,qq
  integer,intent(in) :: i,j,k,l ! p.g. indice
  integer,intent(in) :: a,b ! spinor indices 1,2->L, 3,4->S
  complex(kind=dp),intent(out) :: coef(3)

  complex(kind=dp) :: tmpij11,tmpij12,tmpij21,tmpij22
  complex(kind=dp) :: tmpij33,tmpij34,tmpij43,tmpij44
  complex(kind=dp) :: tmpkl11,tmpkl22,tmpkl33,tmpkl44
  complex(kind=dp) :: tmpkl1122,tmpkl3344

  call set_GammaMatrix

    !m=1
    !LLLL
    
      if(a==1.and.b==1) then
         tmpij11 = conjg(c_psi(1,i,nn))*c_psi(1,j,mm)
         tmpij12 = conjg(c_psi(1,i,nn))*c_psi(2,j,mm)
         tmpij21 = conjg(c_psi(2,i,nn))*c_psi(1,j,mm)
         tmpij22 = conjg(c_psi(2,i,nn))*c_psi(2,j,mm)
         tmpkl11 = conjg(c_psi(1,k,pp))*c_psi(1,l,qq)
         tmpkl22 = conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
         tmpkl1122 = tmpkl11 + tmpkl22
         !m=1
!         tmp =  conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)& 
!              &+conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
!              &+conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!              &+conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
         coef(1) = (tmpij21+tmpij12)*tmpkl1122
         !m=2
!         tmp =  IU*conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!              &+IU*conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
!              &-IU*conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!              &-IU*conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
         coef(2) = IU*(tmpij21-tmpij12)*tmpkl1122
         !m=3
!         tmp =  conjg(c_psi(1,i,nn))*c_psi(1,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!              &+conjg(c_psi(1,i,nn))*c_psi(1,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
!              &-conjg(c_psi(2,i,nn))*c_psi(2,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!              &-conjg(c_psi(2,i,nn))*c_psi(2,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
         coef(3) = (tmpij11-tmpij22)*tmpkl1122

    !LLSS
      else if(a==1.and.b==3) then
         tmpij11 = conjg(c_psi(1,i,nn))*c_psi(1,j,mm)
         tmpij12 = conjg(c_psi(1,i,nn))*c_psi(2,j,mm)
         tmpij21 = conjg(c_psi(2,i,nn))*c_psi(1,j,mm)
         tmpij22 = conjg(c_psi(2,i,nn))*c_psi(2,j,mm)
         tmpkl33 = conjg(c_psi(3,k,pp))*c_psi(3,l,qq)
         tmpkl44 = conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
         tmpkl3344 = tmpkl33 + tmpkl44
         !m=1
!         tmp =  conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!              &+conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
!              &+conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!              &+conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
         coef(1) = (tmpij21 + tmpij12)*tmpkl3344
         !m=2
!         tmp =  IU*conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!              &+IU*conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
!              &-IU*conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!              &-IU*conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
         coef(2) = IU*(tmpij21 - tmpij12)*tmpkl3344
         !m=3
!         tmp =  conjg(c_psi(1,i,nn))*c_psi(1,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!              &+conjg(c_psi(1,i,nn))*c_psi(1,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
!              &-conjg(c_psi(2,i,nn))*c_psi(2,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!              &-conjg(c_psi(2,i,nn))*c_psi(2,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
         coef(3) = (tmpij11 - tmpij22)*tmpkl3344

    !SSLL
      else if(a==3.and.b==1) then
         tmpij33 = conjg(c_psi(3,i,nn))*c_psi(3,j,mm)
         tmpij34 = conjg(c_psi(3,i,nn))*c_psi(4,j,mm)
         tmpij43 = conjg(c_psi(4,i,nn))*c_psi(3,j,mm)
         tmpij44 = conjg(c_psi(4,i,nn))*c_psi(4,j,mm)
         tmpkl11 = conjg(c_psi(1,k,pp))*c_psi(1,l,qq)
         tmpkl22 = conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
         tmpkl1122 = tmpkl11 + tmpkl22
         !m=1
!         tmp = -conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!              &-conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
!              &-conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!              &-conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
         coef(1) = -(tmpij43 + tmpij34)*tmpkl1122
         !m=2
!         tmp = -IU*conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!              &-IU*conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
!              &+IU*conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!              &+IU*conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
         coef(2) = -IU*(tmpij43 - tmpij34)*tmpkl1122
         !m=3
!         tmp = -conjg(c_psi(3,i,nn))*c_psi(3,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!              &-conjg(c_psi(3,i,nn))*c_psi(3,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
!              &+conjg(c_psi(4,i,nn))*c_psi(4,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!              &+conjg(c_psi(4,i,nn))*c_psi(4,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
         coef(3) = -(tmpij33 - tmpij44)*tmpkl1122

    !SSSS
      else if(a==3.and.b==3) then
         tmpij33 = conjg(c_psi(3,i,nn))*c_psi(3,j,mm)
         tmpij34 = conjg(c_psi(3,i,nn))*c_psi(4,j,mm)
         tmpij43 = conjg(c_psi(4,i,nn))*c_psi(3,j,mm)
         tmpij44 = conjg(c_psi(4,i,nn))*c_psi(4,j,mm)
         tmpkl33 = conjg(c_psi(3,k,pp))*c_psi(3,l,qq)
         tmpkl44 = conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
         tmpkl3344 = tmpkl33 + tmpkl44
         !m=1
!         tmp = -conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!              &-conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
!              &-conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!              &-conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
         coef(1) = -(tmpij43+tmpij34)*tmpkl3344
         !m=2
!         tmp = -IU*conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!              &-IU*conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
!              &+IU*conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!              &+IU*conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
         coef(2) = -IU*(tmpij43-tmpij43)*tmpkl3344
         !m=3
!         tmp = -conjg(c_psi(3,i,nn))*c_psi(3,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!              &-conjg(c_psi(3,i,nn))*c_psi(3,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
!              &+conjg(c_psi(4,i,nn))*c_psi(4,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!              &+conjg(c_psi(4,i,nn))*c_psi(4,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
         coef(3) = -(tmpij33-tmpij44)*tmpkl3344
      end if

    return
end subroutine coef_intEeff_EDM_ele

!============================================================
subroutine calc_intEeff_EDM_ele_anz2_pt(actorb0,actorb,occ,intEeff_EDM_ele_pt)
! nn,mm,pp,qq : labels for molecular orbitals (KP) . 1~2*NBS
! anzatz 2
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use PTcoef, only : c_psi
  implicit none

  integer,intent(in) :: actorb0,actorb
  real(kind=dp),intent(in) :: occ(actorb)
  complex(kind=dp),intent(out) :: intEeff_EDM_ele_pt

  integer :: nn,mm,pp,qq
  integer :: i,j,k,l,m ! p.g. indice
  integer :: NLS !function
  integer :: NL,NS
  integer :: a,b ! spinor indices 1,2->L, 3,4->S

  real(kind=dp) :: inttwoelegrad_pg !function
  type(primitive_gaussian) :: pg
  complex(kind=dp) :: tmp, tmptot, sum, sumall
  complex(kind=dp),allocatable :: sumomp(:)
  complex(kind=dp) :: tmpcoef(3),coefNNMM(3),coefNMMN(3)
  integer :: itmp

  call set_GammaMatrix
  call copy_DiracOutput_pg(pg)
  NL=NBS_L; NS=NBS_S

    sumall=(0._dp,0._dp)
    do a=1,3,2
    do b=1,3,2
    do k=1,NLS(b,NL,NS)
       itmp=NLS(b,NL,NS)
       allocate(sumomp(itmp))
       sumomp(:)=(0._dp,0._dp)
       !$omp parallel do private(l,i,j,m,nn,mm,tmpcoef,coefNNMM,coefNMMN,sum)
       do l=1,NLS(b,NL,NS)
          sum=(0._dp,0._dp)
          do i=1,NLS(a,NL,NS)
             do j=1,i

                tmpcoef(:)=(0._dp,0._dp)
                if(i==j) then
                  do nn=actorb0,actorb
                    do mm=actorb0,nn-1
                      call coef_intEeff_EDM_ele(nn,nn,mm,mm,i,j,k,l,a,b,coefNNMM)
                      call coef_intEeff_EDM_ele(nn,mm,mm,nn,i,j,k,l,a,b,coefNMMN)
                      do m=1,3
                        tmpcoef(m) = tmpcoef(m) + occ(nn)*occ(mm)*(coefNNMM(m)-2.d0*dble(coefNMMN(m)))
                      end do
                    end do
                    if(nn.ne.actorb) then
                       do mm=nn+1,actorb
                         do m=1,3
                           tmpcoef(m) = tmpcoef(m) + occ(nn)*occ(mm)*coefNNMM(m)
                         end do
                       end do
                    end if
                  end do!nn
                  do m=1,3
                    sum = sum + tmpcoef(m)*inttwoelegrad_pg(m,i,j,k,l,a,a,b,b,NL,NS,pg)
                  end do

                else !i>j
                  do nn=actorb0,actorb
                    do mm=actorb0,nn-1
                      call coef_intEeff_EDM_ele(nn,nn,mm,mm,i,j,k,l,a,b,coefNNMM)
                      call coef_intEeff_EDM_ele(nn,mm,mm,nn,i,j,k,l,a,b,coefNMMN)
                      do m=1,3
                        tmpcoef(m) = tmpcoef(m) + 2.d0*occ(nn)*occ(mm)*(dble(coefNNMM(m))-dble(coefNMMN(m)))
                      end do
                    end do
                    if(nn.ne.actorb) then
                       do mm=nn+1,actorb
                         call coef_intEeff_EDM_ele(nn,nn,mm,mm,i,j,k,l,a,b,coefNNMM)
                         call coef_intEeff_EDM_ele(nn,mm,mm,nn,i,j,k,l,a,b,coefNMMN)
                         do m=1,3
                           tmpcoef(m) = tmpcoef(m) + 2.d0*occ(nn)*occ(mm)*(dble(coefNNMM(m))-dble(coefNMMN(m)))
                         end do
                       end do
                    end if
                  end do!nn
                  do m=1,3
                    sum = sum + tmpcoef(m)*inttwoelegrad_pg(m,i,j,k,l,a,a,b,b,NL,NS,pg)
                  end do
                end if

             end do!j
          end do!i
          sumomp(l) = sum
       end do!l
       !$omp end parallel do

       do l=1,NLS(b,NL,NS)
         sumall = sumall + sumomp(l)
       end do!l
       deallocate(sumomp)
       write(*,*) a,b,k,'/', NLS(b,NL,NS)
    end do!k
    end do!b
    end do!a

    intEeff_EDM_ele_pt = Ze*sumall

    return
end subroutine calc_intEeff_EDM_ele_anz2_pt

!============================================================
subroutine calc_intEeff_EDM_ele_anz1_pt(actorb0,actorb,occ,intEeff_EDM_ele_pt)
! nn,mm,pp,qq : labels for molecular orbitals (KP) . 1~2*NBS
! anzatz 1
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use PTcoef, only : c_psi
  implicit none

  integer,intent(in) :: actorb0,actorb
  real(kind=dp),intent(in) :: occ(actorb)
  complex(kind=dp),intent(out) :: intEeff_EDM_ele_pt

  integer :: nn,mm,pp,qq
  integer :: i,j,k,l,m ! p.g. indice
  integer :: NLS !function
  integer :: NL,NS
  integer :: a,b ! spinor indices 1,2->L, 3,4->S

  real(kind=dp) :: inttwoelegrad_pg !function
  type(primitive_gaussian) :: pg
  complex(kind=dp) :: tmp, tmptot, sum, sumall
  complex(kind=dp),allocatable :: sumomp(:)
  complex(kind=dp) :: tmpcoef(3),coefNNMM(3),coefNMMN(3)
  integer :: itmp

  call set_GammaMatrix
  call copy_DiracOutput_pg(pg)
  NL=NBS_L; NS=NBS_S

    sumall=(0._dp,0._dp)
    do a=1,3,2
    do b=1,3,2
    do k=1,NLS(b,NL,NS)
       itmp=NLS(b,NL,NS)
       allocate(sumomp(itmp))
       sumomp(:)=(0._dp,0._dp)
       !$omp parallel do private(l,i,j,m,nn,mm,tmpcoef,coefNNMM,coefNMMN,sum)
       do l=1,NLS(b,NL,NS)
          sum=(0._dp,0._dp)
          do i=1,NLS(a,NL,NS)
             do j=1,i

                tmpcoef(:)=(0._dp,0._dp)
                if(i==j) then
                  do nn=actorb0,actorb
                    do mm=actorb0, actorb
                      call coef_intEeff_EDM_ele(nn,nn,mm,mm,i,j,k,l,a,b,coefNNMM)
                      do m=1,3
                        tmpcoef(m) = tmpcoef(m) + occ(nn)*occ(mm)*(coefNNMM(m))
                      end do
                    end do
                  end do!nn
                  do m=1,3
                    sum = sum + tmpcoef(m)*inttwoelegrad_pg(m,i,j,k,l,a,a,b,b,NL,NS,pg)
                  end do

                else !i>j
                  do nn=actorb0,actorb
                    do mm=actorb0,actorb
                      call coef_intEeff_EDM_ele(nn,nn,mm,mm,i,j,k,l,a,b,coefNNMM)
                      do m=1,3
                        tmpcoef(m) = tmpcoef(m) + 2.d0*occ(nn)*occ(mm)*(dble(coefNNMM(m)))
                      end do
                    end do
                  end do!nn
                  do m=1,3
                    sum = sum + tmpcoef(m)*inttwoelegrad_pg(m,i,j,k,l,a,a,b,b,NL,NS,pg)
                  end do
                end if

             end do!j
          end do!i
          sumomp(l) = sum
       end do!l
       !$omp end parallel do

       do l=1,NLS(b,NL,NS)
         sumall = sumall + sumomp(l)
       end do!l
       deallocate(sumomp)
       write(*,*) a,b,k,'/', NLS(b,NL,NS)
    end do!k
    end do!b
    end do!a

    intEeff_EDM_ele_pt = Ze*sumall

    return
end subroutine calc_intEeff_EDM_ele_anz1_pt
