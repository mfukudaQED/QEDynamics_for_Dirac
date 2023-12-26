!============================================================
function intEeff_EDM_ele_pt_Qmat(nn,mm,pp,qq)
! nn,mm,pp,qq : labels for molecular orbitals (KP) . 1~2*NBS
! use (ij ^|^ kl)^k = - (kl ^|^ ij)^k
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use GammaMatrix
  use PTcoef, only : c_psi
  implicit none

  integer,intent(in) :: nn,mm,pp,qq
  integer :: i,j,k,l ! p.g. indice
  integer :: NLS !function
  integer :: NL,NS
  integer :: a,b ! spinor indices 1,2->L, 3,4->S

  real(kind=dp) :: inttwoelegrad_pg !function
  type(primitive_gaussian) :: pg
  complex(kind=dp) :: tmp, tmptot
  complex(kind=dp) :: intEeff_EDM_ele_pt_Qmat

  call set_GammaMatrix
  call copy_DiracOutput_pg(pg)
  NL=NBS_L; NS=NBS_S

    tmptot=(0._dp,0._dp)
    tmp=(0._dp,0._dp)
    !m=1
    !LLLL
    a=1; b=1
    do i=1,NL
       do j=1,NL
          do k=1,NL
             do l=1,NL
                tmpij11 = conjg(c_psi(1,i,nn))*c_psi(1,j,mm)
                tmpij12 = conjg(c_psi(1,i,nn))*c_psi(2,j,mm)
                tmpij21 = conjg(c_psi(2,i,nn))*c_psi(1,j,mm)
                tmpij22 = conjg(c_psi(2,i,nn))*c_psi(2,j,mm)
                tmpkl11 = conjg(c_psi(1,k,pp))*c_psi(1,l,qq)
                tmpkl22 = conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
                tmpkl1122 = tmpkl11 + tmpkl22
                !m=1
!                tmp =  conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)& 
!                     &+conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
!                     &+conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!                     &+conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
                tmp = (tmpij21+tmpij12)*tmpkl1122
                tmptot = tmptot + tmp*inttwoelegrad_pg(1,i,j,k,l,1,1,1,1,NL,NS,pg)
                !m=2
!                tmp =  IU*conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!                     &+IU*conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
!                     &-IU*conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!                     &-IU*conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
                tmp = IU*(tmpij21-tmpij12)*tmpkl1122
                tmptot = tmptot + tmp*inttwoelegrad_pg(2,i,j,k,l,1,1,1,1,NL,NS,pg)
                !m=3
!                tmp =  conjg(c_psi(1,i,nn))*c_psi(1,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!                     &+conjg(c_psi(1,i,nn))*c_psi(1,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
!                     &-conjg(c_psi(2,i,nn))*c_psi(2,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!                     &-conjg(c_psi(2,i,nn))*c_psi(2,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
                tmp = (tmpij11-tmpij22)*tmpkl1122
                tmptot = tmptot + tmp*inttwoelegrad_pg(3,i,j,k,l,1,1,1,1,NL,NS,pg)
             end do
          end do
       end do
    end do

    !LLSS
    a=1; b=3
    do i=1,NL
       do j=1,NL
          do k=1,NS
             do l=1,NS
                tmpij11 = conjg(c_psi(1,i,nn))*c_psi(1,j,mm)
                tmpij12 = conjg(c_psi(1,i,nn))*c_psi(2,j,mm)
                tmpij21 = conjg(c_psi(2,i,nn))*c_psi(1,j,mm)
                tmpij22 = conjg(c_psi(2,i,nn))*c_psi(2,j,mm)
                tmpkl33 = conjg(c_psi(3,k,pp))*c_psi(3,l,qq)
                tmpkl44 = conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
                tmpkl3344 = tmpkl33 + tmpkl44
                !m=1
!                tmp =  conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!                     &+conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
!                     &+conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!                     &+conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
                tmp = (tmpij21 + tmpij12)*tmpkl3344
                tmptot = tmptot + tmp*inttwoelegrad_pg(1,i,j,k,l,1,1,3,3,NL,NS,pg)
                !m=2
!                tmp =  IU*conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!                     &+IU*conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
!                     &-IU*conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!                     &-IU*conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
                tmp = IU*(tmpij21 - tmpij12)*tmpkl3344
                tmptot = tmptot + tmp*inttwoelegrad_pg(2,i,j,k,l,1,1,3,3,NL,NS,pg)
                !m=3
!                tmp =  conjg(c_psi(1,i,nn))*c_psi(1,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!                     &+conjg(c_psi(1,i,nn))*c_psi(1,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
!                     &-conjg(c_psi(2,i,nn))*c_psi(2,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!                     &-conjg(c_psi(2,i,nn))*c_psi(2,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
                tmp = (tmpij11 - tmpij22)*tmpkl3344
                tmptot = tmptot + tmp*inttwoelegrad_pg(3,i,j,k,l,1,1,3,3,NL,NS,pg)
             end do
          end do
       end do
    end do

    !SSLL
    a=3; b=1
    do k=1,NS
       do l=1,NS
          do i=1,NL
             do j=1,NL
                tmpkl33 = conjg(c_psi(3,k,nn))*c_psi(3,l,mm)
                tmpkl34 = conjg(c_psi(3,k,nn))*c_psi(4,l,mm)
                tmpkl43 = conjg(c_psi(4,k,nn))*c_psi(3,l,mm)
                tmpkl44 = conjg(c_psi(4,k,nn))*c_psi(4,l,mm)
                tmpij11 = conjg(c_psi(1,i,pp))*c_psi(1,j,qq)
                tmpij22 = conjg(c_psi(2,i,pp))*c_psi(2,j,qq)
                tmpij1122 = tmpij11 + tmpij22
                !m=1
!                tmp = -conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!                     &-conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
!                     &-conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!                     &-conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
                tmp = -(tmpkl43 + tmpkl34)*tmpij1122
                tmptot = tmptot - tmp*inttwoelegrad_pg(1,i,j,k,l,1,1,3,3,NL,NS,pg)
                !m=2
!                tmp = -IU*conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!                     &-IU*conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
!                     &+IU*conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!                     &+IU*conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
                tmp = -IU*(tmpkl43 - tmpkl34)*tmpij1122
                tmptot = tmptot - tmp*inttwoelegrad_pg(2,i,j,k,l,1,1,3,3,NL,NS,pg)
                !m=3
!                tmp = -conjg(c_psi(3,i,nn))*c_psi(3,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!                     &-conjg(c_psi(3,i,nn))*c_psi(3,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
!                     &+conjg(c_psi(4,i,nn))*c_psi(4,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
!                     &+conjg(c_psi(4,i,nn))*c_psi(4,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
                tmp = -(tmpkl33 - tmpkl44)*tmpij1122
                tmptot = tmptot - tmp*inttwoelegrad_pg(3,i,j,k,l,1,1,3,3,NL,NS,pg)
             end do
          end do
       end do
    end do

    !SSSS
    a=3; b=3
    do i=1,NS
       do j=1,NS
          do k=1,NS
             do l=1,NS
                tmpij33 = conjg(c_psi(3,i,nn))*c_psi(3,j,mm)
                tmpij34 = conjg(c_psi(3,i,nn))*c_psi(4,j,mm)
                tmpij43 = conjg(c_psi(4,i,nn))*c_psi(3,j,mm)
                tmpij44 = conjg(c_psi(4,i,nn))*c_psi(4,j,mm)
                tmpkl33 = conjg(c_psi(3,k,pp))*c_psi(3,l,qq)
                tmpkl44 = conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
                tmpkl3344 = tmpkl33 + tmpkl44
                !m=1
!                tmp = -conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!                     &-conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
!                     &-conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!                     &-conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
                tmp = -(tmpij43+tmpij34)*tmpkl3344
                tmptot = tmptot + tmp*inttwoelegrad_pg(1,i,j,k,l,3,3,3,3,NL,NS,pg)
                !m=2
!                tmp = -IU*conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!                     &-IU*conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
!                     &+IU*conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!                     &+IU*conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
                tmp = -IU*(tmpij43-tmpij43)*tmpkl3344
                tmptot = tmptot + tmp*inttwoelegrad_pg(2,i,j,k,l,3,3,3,3,NL,NS,pg)
                !m=3
!                tmp = -conjg(c_psi(3,i,nn))*c_psi(3,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!                     &-conjg(c_psi(3,i,nn))*c_psi(3,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
!                     &+conjg(c_psi(4,i,nn))*c_psi(4,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
!                     &+conjg(c_psi(4,i,nn))*c_psi(4,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
                tmp = -(tmpij33-tmpij44)*tmpkl3344
                tmptot = tmptot + tmp*inttwoelegrad_pg(3,i,j,k,l,3,3,3,3,NL,NS,pg)
             end do
          end do
       end do
    end do

    intEeff_EDM_ele_pt_Qmat = Ze*tmptot

    return
end function intEeff_EDM_ele_pt_Qmat

