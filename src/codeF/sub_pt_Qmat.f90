! Last Change:28-Sep-2015.
!====================================================================
! integrals
!   function intN_pt_Qmat(nn,mm)
!   function intSpin_pt_Qmat(k,nn,mm)
!   function intjmoment_pt_Qmat(i,j,nn,mm)
!   function intspol_pt_Qmat(k,nn,mm)
!   function intgamk_pt_Qmat(k,nn,mm)
!   function intEDMtorqE_pt_Qmat(i,nn,mm)
!   function intEDMtorqB_pt_Qmat(i,nn,mm)
!   function intEeff_EDM_ele_pt_Qmat(nn,mm,pp,qq) written 13/11/08
!   function intEeff_EDM_nuc_pt_Qmat(nn,mm)
!   function intEeff_EDM_ob_pt_Qmat(nn,mm)
!   subroutine calc_intOrbang_pt_Qmat(nn,mm,intOrbang_Qmat)
!   subroutine calc_intOrbang_AM_pt_Qmat(nn,mm,intOrbang_AM_Qmat)
!   subroutine calc_intdensmom_pt_Qmat(nn,mm,intdensmom_mat)
!   function intMagHyp_pt_Qmat(posR,nn,mm)

! densities
!   function Hso_nuc_pt_Qmat(x,y,z,nn,mm)
!   function js_pt_Qmat(i,j,x,y,z,nn,mm)
!   function zeta_pt_Qmat(i,x,y,z,nn,mm)   :zeta force density
!   function t_pt_Qmat(i,x,y,z,nn,mm)      :spintorque density
!   function tau_pt_Qmat(k,l,x,y,z,nn,mm)  :stress tensor density
!   function divj_pt_Qmat(x,y,z,nn,mm)     :div of current density
!   function j_pt_Qmat(i,x,y,z,nn,mm)      :current density
!   function s_pt_Qmat(i,x,y,z,nn,mm)      :spin density
!   function rho_pt_Qmat(x,y,z,nn,mm)      :charge density
!   function density_pt_Qmat(x,y,z,nn,mm)  :position probability density
!   function t_AM_pt_Qmat(t,i,x,y,z,nn,mm) :AM term of spin torque 
!   function pi_pt_Qmat(i,x,y,z,nn,mm)     :kinetic momentum density
!   function orbang_pt_Qmat(i,x,y,z,nn,mm) :orbital angular momentum density
!   function AA_pt_Qmat(i,x,y,z,nn,mm)     :vector potential A_A
!   function tau_AA_pt_Qmat(k,l,x,y,z,nn,mm,pp,qq)
!   function t_AA_pt_Qmat(i,x,y,z,nn,mm,pp,qq)
!====================================================================


!============================================================
function intN_pt_Qmat(nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef
  implicit none
 
  complex(kind=dp) :: intN_pt_Qmat
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: intN_pq
  integer :: i
  
  call copy_DiracOutput_pg(pg)
  
  intN_pt_Qmat = (0._dp,0._dp)
  do i=1,4
     call calc_intN_pq(i,i,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intN_pq) 
!     write(*,*) i,intN_pq
     intN_pt_Qmat = intN_pt_Qmat +intN_pq
  end do
  return
end function intN_pt_Qmat

!============================================================
function intSpina_pt_Qmat(nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use PTcoef
  implicit none
 
  complex(kind=dp) :: intSpina_pt_Qmat
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: intSpin_pq_11,intSpin_pq_22,intSpin_pq_33,intSpin_pq_44
  complex(kind=dp) :: intSpin_pq_12,intSpin_pq_21,intSpin_pq_34,intSpin_pq_43
  
  call set_GammaMatrix
  call copy_DiracOutput_pg(pg)
  
  intSpina_pt_Qmat = (0._dp,0._dp)
     call calc_intN_pq(1,1,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intSpin_pq_11) 
     call calc_intN_pq(3,3,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intSpin_pq_33) 
    
     intSpina_pt_Qmat = (intSpin_pq_11 +intSpin_pq_33 )*0.5d0

  return

end function intSpina_pt_Qmat

!============================================================
function intSpinb_pt_Qmat(nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use PTcoef
  implicit none
 
  complex(kind=dp) :: intSpinb_pt_Qmat
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: intSpin_pq_11,intSpin_pq_22,intSpin_pq_33,intSpin_pq_44
  complex(kind=dp) :: intSpin_pq_12,intSpin_pq_21,intSpin_pq_34,intSpin_pq_43
  
  call set_GammaMatrix
  call copy_DiracOutput_pg(pg)
  
  intSpinb_pt_Qmat = (0._dp,0._dp)
     call calc_intN_pq(2,2,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intSpin_pq_22) 
     call calc_intN_pq(4,4,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intSpin_pq_44) 
    
     intSpinb_pt_Qmat = (intSpin_pq_22 + intSpin_pq_44 )*0.5d0

  return

end function intSpinb_pt_Qmat

!============================================================
function intSpin_pt_Qmat(k,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use PTcoef
  implicit none
 
  complex(kind=dp) :: intSpin_pt_Qmat
  integer,intent(in) :: k
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: intSpin_pq_11,intSpin_pq_22,intSpin_pq_33,intSpin_pq_44
  complex(kind=dp) :: intSpin_pq_12,intSpin_pq_21,intSpin_pq_34,intSpin_pq_43
  
  call set_GammaMatrix
  call copy_DiracOutput_pg(pg)
  
  intSpin_pt_Qmat = (0._dp,0._dp)
  if(k==1) then
     call calc_intN_pq(1,2,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intSpin_pq_12) 
     call calc_intN_pq(2,1,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intSpin_pq_21) 
     call calc_intN_pq(3,4,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intSpin_pq_34) 
     call calc_intN_pq(4,3,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intSpin_pq_43) 
    
     intSpin_pt_Qmat = (intSpin_pq_12 +intSpin_pq_21 +intSpin_pq_34 +intSpin_pq_43 )*0.5d0

  else if(k==2) then
     call calc_intN_pq(1,2,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intSpin_pq_12) 
     call calc_intN_pq(2,1,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intSpin_pq_21) 
     call calc_intN_pq(3,4,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intSpin_pq_34) 
     call calc_intN_pq(4,3,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intSpin_pq_43) 
    
     intSpin_pt_Qmat = (-IU*intSpin_pq_12 +IU*intSpin_pq_21 -IU*intSpin_pq_34 +IU*intSpin_pq_43 )*0.5d0

  else if(k==3) then
     call calc_intN_pq(1,1,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intSpin_pq_11) 
     call calc_intN_pq(2,2,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intSpin_pq_22) 
     call calc_intN_pq(3,3,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intSpin_pq_33) 
     call calc_intN_pq(4,4,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intSpin_pq_44) 
    
     intSpin_pt_Qmat = (intSpin_pq_11 -intSpin_pq_22 +intSpin_pq_33 -intSpin_pq_44 )*0.5d0
   end if

  return
end function intSpin_pt_Qmat

!============================================================
function int_torqAM(i,nn,mm)
!============================================================
  use Precision
  use Constants
  use param_AM
  implicit none

  integer,intent(in) :: i,nn,mm
  complex(kind=dp) :: int_torqAM

  complex(kind=dp) :: intjmoment_pt_Qmat
  complex(kind=dp) :: int_tauAM(2)

  if (i.eq.1) then

    int_tauAM(1) = Bmag(3)*intjmoment_pt_Qmat(1,3,nn,mm) - Bmag(1)*intjmoment_pt_Qmat(3,3,nn,mm)
    int_tauAM(2) = Bmag(1)*intjmoment_pt_Qmat(2,2,nn,mm) - Bmag(2)*intjmoment_pt_Qmat(1,2,nn,mm)
    int_tauAM = 0.5_dp*Ze*int_tauAM
    int_torqAM = -int_tauAM(1) +int_tauAM(2)

  else if (i.eq.2) then

    int_tauAM(1) = Bmag(1)*intjmoment_pt_Qmat(2,1,nn,mm) - Bmag(2)*intjmoment_pt_Qmat(1,1,nn,mm)
    int_tauAM(2) = Bmag(2)*intjmoment_pt_Qmat(3,3,nn,mm) - Bmag(3)*intjmoment_pt_Qmat(2,3,nn,mm)
    int_tauAM = 0.5_dp*Ze*int_tauAM
    int_torqAM = -int_tauAM(1) +int_tauAM(2)

  else if (i.eq.3) then

    int_tauAM(1) = Bmag(2)*intjmoment_pt_Qmat(3,2,nn,mm) - Bmag(3)*intjmoment_pt_Qmat(2,2,nn,mm)
    int_tauAM(2) = Bmag(3)*intjmoment_pt_Qmat(1,1,nn,mm) - Bmag(1)*intjmoment_pt_Qmat(3,1,nn,mm)
    int_tauAM = 0.5_dp*Ze*int_tauAM
    int_torqAM = -int_tauAM(1) +int_tauAM(2)

  end if

  return
end function int_torqAM

!============================================================
function intjmoment_pt_Qmat(i,j,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
! \psi^\dagger x^i \gamma^0 \gamma^j \psi
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef, only : c_psi
  use Constants  ! use IU
  implicit none

  complex(kind=8) :: intjmoment_pt_Qmat
  integer,intent(in) :: i,j
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: intmom_pq_14(3),intmom_pq_23(3),intmom_pq_32(3),intmom_pq_41(3)
  complex(kind=dp) :: intmom_pq_13(3),intmom_pq_24(3),intmom_pq_31(3),intmom_pq_42(3)

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in intjmoment_pt_Qmat."
     stop
  end if

  if(j.lt.1 .or. j.gt.3) then
     write(*,*) "j should be 1-3 in intjmoment_pt_Qmat."
     stop
  end if
  
  call set_GammaMatrix
  call copy_DiracOutput_pg(pg)
  
  intjmoment_pt_Qmat = (0._dp,0._dp)
  if(j==1) then
     call calc_intmoment_pq(1,4,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intmom_pq_14) 
     call calc_intmoment_pq(2,3,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intmom_pq_23) 
     call calc_intmoment_pq(3,2,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intmom_pq_32) 
     call calc_intmoment_pq(4,1,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intmom_pq_41) 
    
     intjmoment_pt_Qmat = intmom_pq_14(i) +intmom_pq_23(i) +intmom_pq_32(i) +intmom_pq_41(i) 

  else if(j==2) then
     call calc_intmoment_pq(1,4,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intmom_pq_14) 
     call calc_intmoment_pq(2,3,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intmom_pq_23) 
     call calc_intmoment_pq(3,2,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intmom_pq_32) 
     call calc_intmoment_pq(4,1,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intmom_pq_41) 
    
     intjmoment_pt_Qmat = -IU*intmom_pq_14(i) +IU*intmom_pq_23(i) -IU*intmom_pq_32(i) +IU*intmom_pq_41(i)

  else if(j==3) then
     call calc_intmoment_pq(1,3,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intmom_pq_13) 
     call calc_intmoment_pq(2,4,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intmom_pq_24) 
     call calc_intmoment_pq(3,1,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intmom_pq_31) 
     call calc_intmoment_pq(4,2,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intmom_pq_42) 
    
     intjmoment_pt_Qmat = intmom_pq_13(i) -intmom_pq_24(i) +intmom_pq_31(i) -intmom_pq_42(i)
   end if

  return
end function intjmoment_pt_Qmat

!============================================================
function intspol_pt_Qmat(k,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
! int d^3x (\psi^\dagger \gamma^0 \Sigma^k \psi)
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use PTcoef
  implicit none
 
  complex(kind=dp) :: intspol_pt_Qmat
  integer,intent(in) :: k
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: intspol_pq_11,intspol_pq_22,intspol_pq_33,intspol_pq_44
  complex(kind=dp) :: intspol_pq_12,intspol_pq_21,intspol_pq_34,intspol_pq_43
  
  call set_GammaMatrix
  call copy_DiracOutput_pg(pg)
  
  intspol_pt_Qmat = (0._dp,0._dp)
  if(k==1) then
     call calc_intN_pq(1,2,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intspol_pq_12) 
     call calc_intN_pq(2,1,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intspol_pq_21) 
     call calc_intN_pq(3,4,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intspol_pq_34) 
     call calc_intN_pq(4,3,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intspol_pq_43) 
    
     intspol_pt_Qmat = (intspol_pq_12 +intspol_pq_21 -intspol_pq_34 -intspol_pq_43 )*0.5d0

  else if(k==2) then
     call calc_intN_pq(1,2,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intspol_pq_12) 
     call calc_intN_pq(2,1,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intspol_pq_21) 
     call calc_intN_pq(3,4,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intspol_pq_34) 
     call calc_intN_pq(4,3,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intspol_pq_43) 
    
     intspol_pt_Qmat = (-IU*intspol_pq_12 +IU*intspol_pq_21 +IU*intspol_pq_34 -IU*intspol_pq_43 )*0.5d0

  else if(k==3) then
     call calc_intN_pq(1,1,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intspol_pq_11) 
     call calc_intN_pq(2,2,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intspol_pq_22) 
     call calc_intN_pq(3,3,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intspol_pq_33) 
     call calc_intN_pq(4,4,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intspol_pq_44) 
    
     intspol_pt_Qmat = (intspol_pq_11 -intspol_pq_22 -intspol_pq_33 +intspol_pq_44 )*0.5d0
   end if

  return
end function intspol_pt_Qmat

!============================================================
function intgamk_pt_Qmat(k,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
! int d^3x (\psi^\dagger \gamma^k \psi)
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use PTcoef
  implicit none
 
  complex(kind=dp) :: intgamk_pt_Qmat
  integer,intent(in) :: k
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: intgamk_pq_14,intgamk_pq_23,intgamk_pq_32,intgamk_pq_41
  complex(kind=dp) :: intgamk_pq_13,intgamk_pq_24,intgamk_pq_31,intgamk_pq_42
  
  call set_GammaMatrix
  call copy_DiracOutput_pg(pg)
  
  intgamk_pt_Qmat = (0._dp,0._dp)
  if(k==1) then
     call calc_intN_pq(1,4,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intgamk_pq_14) 
     call calc_intN_pq(2,3,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intgamk_pq_23) 
     call calc_intN_pq(3,2,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intgamk_pq_32) 
     call calc_intN_pq(4,1,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intgamk_pq_41) 
    
     intgamk_pt_Qmat = intgamk_pq_14 +intgamk_pq_23 -intgamk_pq_32 -intgamk_pq_41

  else if(k==2) then
     call calc_intN_pq(1,4,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intgamk_pq_14) 
     call calc_intN_pq(2,3,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intgamk_pq_23) 
     call calc_intN_pq(3,2,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intgamk_pq_32) 
     call calc_intN_pq(4,1,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intgamk_pq_41) 
    
     intgamk_pt_Qmat = -IU*intgamk_pq_14 +IU*intgamk_pq_23 +IU*intgamk_pq_32 -IU*intgamk_pq_41

  else if(k==3) then
     call calc_intN_pq(1,3,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intgamk_pq_13) 
     call calc_intN_pq(2,4,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intgamk_pq_24) 
     call calc_intN_pq(3,1,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intgamk_pq_31) 
     call calc_intN_pq(4,2,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intgamk_pq_42) 
    
     intgamk_pt_Qmat = intgamk_pq_13 -intgamk_pq_24 -intgamk_pq_31 +intgamk_pq_42 
   end if

  return
end function intgamk_pt_Qmat

!============================================================
function intEDMtorqE_pt_Qmat(i,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef, only : c_psi
  use param_AM, only : Eelec
  implicit none

  complex(kind=8) :: intEDMtorqE_pt_Qmat
  complex(kind=8) :: intspol_pt_Qmat !function
  integer,intent(in) :: i
  integer,intent(in) :: nn,mm

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in intEDMtorqE_pt_Qmat."
     stop
  end if
  
  if(i.eq.1) then
    intEDMtorqE_pt_Qmat = intspol_pt_Qmat(2,nn,mm)*Eelec(3) - intspol_pt_Qmat(3,nn,mm)*Eelec(2)
  else if(i.eq.2) then
    intEDMtorqE_pt_Qmat = intspol_pt_Qmat(3,nn,mm)*Eelec(1) - intspol_pt_Qmat(1,nn,mm)*Eelec(3)
  else if(i.eq.3) then
    intEDMtorqE_pt_Qmat = intspol_pt_Qmat(1,nn,mm)*Eelec(2) - intspol_pt_Qmat(2,nn,mm)*Eelec(1)
  end if

  intEDMtorqE_pt_Qmat = intEDMtorqE_pt_Qmat*2._dp

  return
end function intEDMtorqE_pt_Qmat

!============================================================
function intEDMtorqB_pt_Qmat(i,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef, only : c_psi
  use param_AM, only : Bmag
  use Constants  ! use IU
  implicit none

  complex(kind=8) :: intEDMtorqB_pt_Qmat
  complex(kind=8) :: intgamk_pt_Qmat !function
  integer,intent(in) :: i
  integer,intent(in) :: nn,mm

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in intEDMtorqB_pt_Qmat."
     stop
  end if
  
  if(i.eq.1) then
    intEDMtorqB_pt_Qmat = intgamk_pt_Qmat(2,nn,mm)*Bmag(3) - intgamk_pt_Qmat(3,nn,mm)*Bmag(2)
  else if(i.eq.2) then
    intEDMtorqB_pt_Qmat = intgamk_pt_Qmat(3,nn,mm)*Bmag(1) - intgamk_pt_Qmat(1,nn,mm)*Bmag(3)
  else if(i.eq.3) then
    intEDMtorqB_pt_Qmat = intgamk_pt_Qmat(1,nn,mm)*Bmag(2) - intgamk_pt_Qmat(2,nn,mm)*Bmag(1)
  end if

  intEDMtorqB_pt_Qmat = IU*intEDMtorqB_pt_Qmat

  return
end function intEDMtorqB_pt_Qmat

!============================================================
function intEeff_EDM_ele_pt_Qmat(nn,mm,pp,qq)
! nn,mm,pp,qq : labels for molecular orbitals (KP) . 1~2*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
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
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
          do k=1,NLS(b,NL,NS)
             do l=1,NLS(b,NL,NS)
                tmp =  conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)& 
                     &+conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
                     &+conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
                     &+conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
                tmptot = tmptot + tmp*inttwoelegrad_pg(1,i,j,k,l,a,a,b,b,NL,NS,pg)
             end do
          end do
       end do
    end do

    !LLSS
    a=1; b=3
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
          do k=1,NLS(b,NL,NS)
             do l=1,NLS(b,NL,NS)
                tmp =  conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
                     &+conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
                     &+conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
                     &+conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
                tmptot = tmptot + tmp*inttwoelegrad_pg(1,i,j,k,l,a,a,b,b,NL,NS,pg)
             end do
          end do
       end do
    end do

    !SSLL
    a=3; b=1
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
          do k=1,NLS(b,NL,NS)
             do l=1,NLS(b,NL,NS)
                tmp = -conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
                     &-conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
                     &-conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
                     &-conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
                tmptot = tmptot +  tmp*inttwoelegrad_pg(1,i,j,k,l,a,a,b,b,NL,NS,pg)
             end do
          end do
       end do
    end do

    !SSSS
    a=3; b=3
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
          do k=1,NLS(b,NL,NS)
             do l=1,NLS(b,NL,NS)
                tmp = -conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
                     &-conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
                     &-conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
                     &-conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
                tmptot = tmptot + tmp*inttwoelegrad_pg(1,i,j,k,l,a,a,b,b,NL,NS,pg)
             end do
          end do
       end do
    end do

    !m=2
    !LLLL
    a=1; b=1
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
          do k=1,NLS(b,NL,NS)
             do l=1,NLS(b,NL,NS)
                tmp =  IU*conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
                     &+IU*conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
                     &-IU*conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
                     &-IU*conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
                tmptot = tmptot + tmp*inttwoelegrad_pg(2,i,j,k,l,a,a,b,b,NL,NS,pg)
             end do
          end do
       end do
    end do

    !LLSS
    a=1; b=3
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
          do k=1,NLS(b,NL,NS)
             do l=1,NLS(b,NL,NS)
                tmp =  IU*conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
                     &+IU*conjg(c_psi(2,i,nn))*c_psi(1,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
                     &-IU*conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
                     &-IU*conjg(c_psi(1,i,nn))*c_psi(2,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
                tmptot = tmptot + tmp*inttwoelegrad_pg(2,i,j,k,l,a,a,b,b,NL,NS,pg)
             end do
          end do
       end do
    end do

    !SSLL
    a=3; b=1
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
          do k=1,NLS(b,NL,NS)
             do l=1,NLS(b,NL,NS)
                tmp = -IU*conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
                     &-IU*conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
                     &+IU*conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
                     &+IU*conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
                tmptot = tmptot + tmp*inttwoelegrad_pg(2,i,j,k,l,a,a,b,b,NL,NS,pg)
             end do
          end do
       end do
    end do

    !SSSS
    a=3; b=3
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
          do k=1,NLS(b,NL,NS)
             do l=1,NLS(b,NL,NS)
                tmp = -IU*conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
                     &-IU*conjg(c_psi(4,i,nn))*c_psi(3,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
                     &+IU*conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
                     &+IU*conjg(c_psi(3,i,nn))*c_psi(4,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
                tmptot = tmptot + tmp*inttwoelegrad_pg(2,i,j,k,l,a,a,b,b,NL,NS,pg)
             end do
          end do
       end do
    end do

    !m=3
    !LLLL
    a=1; b=1
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
          do k=1,NLS(b,NL,NS)
             do l=1,NLS(b,NL,NS)
                tmp =  conjg(c_psi(1,i,nn))*c_psi(1,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
                     &+conjg(c_psi(1,i,nn))*c_psi(1,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
                     &-conjg(c_psi(2,i,nn))*c_psi(2,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
                     &-conjg(c_psi(2,i,nn))*c_psi(2,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
                tmptot = tmptot + tmp*inttwoelegrad_pg(3,i,j,k,l,a,a,b,b,NL,NS,pg)
             end do
          end do
       end do
    end do

    !LLSS
    a=1; b=3
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
          do k=1,NLS(b,NL,NS)
             do l=1,NLS(b,NL,NS)
                tmp =  conjg(c_psi(1,i,nn))*c_psi(1,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
                     &+conjg(c_psi(1,i,nn))*c_psi(1,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
                     &-conjg(c_psi(2,i,nn))*c_psi(2,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
                     &-conjg(c_psi(2,i,nn))*c_psi(2,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
                tmptot = tmptot + tmp*inttwoelegrad_pg(3,i,j,k,l,a,a,b,b,NL,NS,pg)
             end do
          end do
       end do
    end do

    !SSLL
    a=3; b=1
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
          do k=1,NLS(b,NL,NS)
             do l=1,NLS(b,NL,NS)
                tmp = -conjg(c_psi(3,i,nn))*c_psi(3,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
                     &-conjg(c_psi(3,i,nn))*c_psi(3,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)&
                     &+conjg(c_psi(4,i,nn))*c_psi(4,j,mm)*conjg(c_psi(1,k,pp))*c_psi(1,l,qq)&
                     &+conjg(c_psi(4,i,nn))*c_psi(4,j,mm)*conjg(c_psi(2,k,pp))*c_psi(2,l,qq)
                tmptot = tmptot + tmp*inttwoelegrad_pg(3,i,j,k,l,a,a,b,b,NL,NS,pg)
             end do
          end do
       end do
    end do

    !SSSS
    a=3; b=3
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
          do k=1,NLS(b,NL,NS)
             do l=1,NLS(b,NL,NS)
                tmp = -conjg(c_psi(3,i,nn))*c_psi(3,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
                     &-conjg(c_psi(3,i,nn))*c_psi(3,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)&
                     &+conjg(c_psi(4,i,nn))*c_psi(4,j,mm)*conjg(c_psi(3,k,pp))*c_psi(3,l,qq)&
                     &+conjg(c_psi(4,i,nn))*c_psi(4,j,mm)*conjg(c_psi(4,k,pp))*c_psi(4,l,qq)
                tmptot = tmptot + tmp*inttwoelegrad_pg(3,i,j,k,l,a,a,b,b,NL,NS,pg)
             end do
          end do
       end do
    end do

    intEeff_EDM_ele_pt_Qmat = Ze*tmptot

    return
end function intEeff_EDM_ele_pt_Qmat

!============================================================
function intEeff_EDM_nuc_pt_Qmat(nn,mm)
! nn,mm : labels for molecular orbitals (KP) . 1~2*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use PTcoef, only : c_psi
  implicit none

  integer,intent(in) :: nn,mm
  integer :: i,j,k,l ! p.g. indice
  integer :: NLS !function
  integer :: NL,NS
  integer :: a,b ! spinor indices 1,2->L, 3,4->S

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: tmp, tmptot
  complex(kind=dp) :: intEeff_EDM_nuc_pt_Qmat
  real(kind=dp) :: posR(3)  
  complex(kind=dp) :: intE_pg !function

  call set_GammaMatrix
  call copy_DiracOutput_pg(pg)
  NL=NBS_L; NS=NBS_S

    intEeff_EDM_nuc_pt_Qmat=(0._dp,0._dp)
    tmp=(0._dp,0._dp)
  do k=1,NAT !number of atoms
    posR(1)=xc(k) ! position of nucleus
    posR(2)=yc(k)
    posR(3)=zc(k)
!    write(*,*)'posR',k, posR(1), posR(2), posR(3)
    tmptot=(0._dp,0._dp)
    !m=1
    !LL
    a=1
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
         tmp =  conjg(c_psi(2,i,nn))*c_psi(1,j,mm)& 
              &+conjg(c_psi(1,i,nn))*c_psi(2,j,mm)
         tmptot = tmptot + tmp*intE_pg(1,posR,i,j,a,a,NL,NS,pg)
       end do
    end do

    !SS
    a=3
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
         tmp = -conjg(c_psi(4,i,nn))*c_psi(3,j,mm)& 
              &-conjg(c_psi(3,i,nn))*c_psi(4,j,mm)
         tmptot = tmptot + tmp*intE_pg(1,posR,i,j,a,a,NL,NS,pg)
       end do
    end do

    !m=2
    !LL
    a=1
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
         tmp =  IU*conjg(c_psi(2,i,nn))*c_psi(1,j,mm)& 
              &-IU*conjg(c_psi(1,i,nn))*c_psi(2,j,mm)
         tmptot = tmptot + tmp*intE_pg(2,posR,i,j,a,a,NL,NS,pg)
       end do
    end do

    !SS
    a=3
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
         tmp = -IU*conjg(c_psi(4,i,nn))*c_psi(3,j,mm)& 
              &+IU*conjg(c_psi(3,i,nn))*c_psi(4,j,mm)
         tmptot = tmptot + tmp*intE_pg(2,posR,i,j,a,a,NL,NS,pg)
       end do
    end do

    !m=3
    !LL
    a=1
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
         tmp =  conjg(c_psi(1,i,nn))*c_psi(1,j,mm)& 
              &-conjg(c_psi(2,i,nn))*c_psi(2,j,mm)
         tmptot = tmptot + tmp*intE_pg(3,posR,i,j,a,a,NL,NS,pg)
       end do
    end do

    !SS
    a=3
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(a,NL,NS)
         tmp = -conjg(c_psi(3,i,nn))*c_psi(3,j,mm)& 
              &+conjg(c_psi(4,i,nn))*c_psi(4,j,mm)
         tmptot = tmptot + tmp*intE_pg(3,posR,i,j,a,a,NL,NS,pg)
       end do
    end do

    intEeff_EDM_nuc_pt_Qmat = intEeff_EDM_nuc_pt_Qmat + cn(k)*tmptot ! cn -> nuclear charge
!    write(*,*)'NAT=',k,cn(k), nn,mm, tmptot
  end do !k=1,NAT

  intEeff_EDM_nuc_pt_Qmat = - intEeff_EDM_nuc_pt_Qmat ! minus for intE_pg , added in 20140204
    
    return
end function intEeff_EDM_nuc_pt_Qmat

!============================================================
function intEeff_EDM_ob_pt_Qmat(nn,mm)
! nn,mm : labels for molecular orbitals (KP) . 1~2*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use PTcoef, only : c_psi
  implicit none

  integer,intent(in) :: nn,mm
  integer :: i,j,k,l ! p.g. indice
  integer :: NLS !function
  integer :: NL,NS
  integer :: a,b ! spinor indices 1,2->L, 3,4->S

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: tmp, tmptot
  complex(kind=dp) :: intEeff_EDM_ob_pt_Qmat
  real(kind=dp) :: posR(3)  
  real(kind=dp) :: intEeff_EDM_ob_pg !function


  call set_GammaMatrix
  call copy_DiracOutput_pg(pg)
  NL=NBS_L; NS=NBS_S

    intEeff_EDM_ob_pt_Qmat=(0._dp,0._dp)
    tmptot=(0._dp,0._dp)
    tmp=(0._dp,0._dp)
    !LS
    a=1; b=3
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(b,NL,NS)
         tmp =  conjg(c_psi(1,i,nn))*c_psi(3,j,mm)& 
              &+conjg(c_psi(2,i,nn))*c_psi(4,j,mm)
         tmptot = tmptot + tmp*intEeff_EDM_ob_pg(i,j,a,b,NL,NS,pg)
       end do
    end do

    !SL
    a=3; b=1
    do i=1,NLS(a,NL,NS)
       do j=1,NLS(b,NL,NS)
         tmp = -conjg(c_psi(3,i,nn))*c_psi(1,j,mm)& 
              &-conjg(c_psi(4,i,nn))*c_psi(2,j,mm)
         tmptot = tmptot + tmp*intEeff_EDM_ob_pg(i,j,a,b,NL,NS,pg)
       end do
    end do

    intEeff_EDM_ob_pt_Qmat = 2._dp*IU*CCC*tmptot/Ze

    return
end function intEeff_EDM_ob_pt_Qmat

!============================================================
subroutine calc_intOrbang_pt_Qmat(nn,mm,intOrbang_Qmat)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
! intOrbang_mat = -i \hbar \int d^3r \psi_p^a x \times \nabla \psi_q^b
! The vector potential term is NOT included.
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use PTcoef
  implicit none
 
  complex(kind=dp),intent(out) :: intOrbang_Qmat(3)
  integer,intent(in) :: nn,mm

  integer :: k
  type(primitive_gaussian) :: pg
  complex(kind=dp) :: intOrbang_pq_11(3),intOrbang_pq_22(3),intOrbang_pq_33(3),intOrbang_pq_44(3)
  
  call set_GammaMatrix
  call copy_DiracOutput_pg(pg)
  
     call calc_intorbang_pq(1,1,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intOrbang_pq_11) 
     call calc_intorbang_pq(2,2,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intOrbang_pq_22) 
     call calc_intorbang_pq(3,3,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intOrbang_pq_33) 
     call calc_intorbang_pq(4,4,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intOrbang_pq_44) 
    
  intOrbang_Qmat(:) = (0._dp,0._dp)
     do k=1,3
        intOrbang_Qmat(k) = intOrbang_pq_11(k) +intOrbang_pq_22(k) +intOrbang_pq_33(k) +intOrbang_pq_44(k)
     end do

  return
end subroutine calc_intOrbang_pt_Qmat

!============================================================
subroutine calc_intOrbang_AM_pt_Qmat(nn,mm,intOrbang_AM_Qmat)
! intOrbang_AM_mat = -(Zee/2c)* \int d^3r \psi_p^a [x^j x^j Bmag(i) - x^i x^j Bmag(j)] \psi_q^b
! AM = 1/2 * B \times r
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use PTcoef
  use param_AM
  implicit none
 
  complex(kind=dp),intent(out) :: intOrbang_AM_Qmat(3)
  integer,intent(in) :: nn,mm

  integer :: k,l
  type(primitive_gaussian) :: pg
  complex(kind=dp) :: int2moment_pq_11(6),int2moment_pq_22(6),int2moment_pq_33(6),int2moment_pq_44(6)
  complex(kind=dp) :: tmp(6), tmpmat(3,3)
  
  call set_GammaMatrix
  call copy_DiracOutput_pg(pg)
  
     call calc_int_2nd_moment_pq(1,1,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),int2moment_pq_11) 
     call calc_int_2nd_moment_pq(2,2,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),int2moment_pq_22) 
     call calc_int_2nd_moment_pq(3,3,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),int2moment_pq_33) 
     call calc_int_2nd_moment_pq(4,4,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),int2moment_pq_44) 

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
    
  intOrbang_AM_Qmat(:) = (0._dp,0._dp)
     do k=1,3
        do l=1,3
          intOrbang_AM_Qmat(k) = tmpmat(l,l)*Bmag(k) - tmpmat(k,l)*Bmag(l)
        end do
        intOrbang_AM_Qmat(k) = - Ze*intOrbang_AM_Qmat(k)/CCC/2.d0
     end do

!     do k=1,3
!        write(*,*) k, intOrbang_AM_Qmat(k)
!     end do
!     stop

  return
end subroutine calc_intOrbang_AM_pt_Qmat

!============================================================
subroutine calc_intdensmom_pt_Qmat(nn,mm,intdensmom_mat)
! intdensmom_mat(i) = \int d^3r \psi_nn x^i  \psi_mm
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use PTcoef
  implicit none
 
  complex(kind=dp),intent(out) :: intdensmom_mat(3)
  integer,intent(in) :: nn,mm

  integer :: k
  type(primitive_gaussian) :: pg
  complex(kind=dp) :: intmom_pq_11(3),intmom_pq_22(3),intmom_pq_33(3),intmom_pq_44(3)
  
  call set_GammaMatrix
  call copy_DiracOutput_pg(pg)
  
     call calc_intmoment_pq(1,1,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intmom_pq_11) 
     call calc_intmoment_pq(2,2,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intmom_pq_22) 
     call calc_intmoment_pq(3,3,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intmom_pq_33) 
     call calc_intmoment_pq(4,4,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intmom_pq_44) 
    
  intdensmom_mat(:) = (0._dp,0._dp)
     do k=1,3
        intdensmom_mat(k) = intmom_pq_11(k) +intmom_pq_22(k) +intmom_pq_33(k) +intmom_pq_44(k)
     end do

  return
end subroutine calc_intdensmom_pt_Qmat

!============================================================
subroutine calc_intdensmom_wpos_pt_Qmat(posR,nn,mm,intdensmom_mat)
! intdensmom_mat(i) = \int d^3r \psi_nn x^i  \psi_mm
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use PTcoef
  implicit none
 
  complex(kind=dp),intent(out) :: intdensmom_mat(3)
  integer,intent(in) :: nn,mm
  real(kind=dp),intent(in) :: posR(3)

  integer :: k
  type(primitive_gaussian) :: pg
  complex(kind=dp) :: intmom_pq_11(3),intmom_pq_22(3),intmom_pq_33(3),intmom_pq_44(3)
  
  call set_GammaMatrix
  call copy_DiracOutput_pg(pg)
  
     call calc_intmoment_wpos_pq(posR,1,1,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intmom_pq_11) 
     call calc_intmoment_wpos_pq(posR,2,2,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intmom_pq_22) 
     call calc_intmoment_wpos_pq(posR,3,3,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intmom_pq_33) 
     call calc_intmoment_wpos_pq(posR,4,4,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intmom_pq_44) 
    
  intdensmom_mat(:) = (0._dp,0._dp)
     do k=1,3
        intdensmom_mat(k) = intmom_pq_11(k) +intmom_pq_22(k) +intmom_pq_33(k) +intmom_pq_44(k)
     end do

  return
end subroutine calc_intdensmom_wpos_pt_Qmat

!============================================================
function intMagHyp_pt_Qmat(posR,nn,mm)
! Z_ee \mu^i epsilon_{ijk} \int d^3r \psi^\dagger \gamma^0 \gamma^j r^k/|\vec{r}|^3 \psi
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants  ! use IU
  use PTcoef
  use param_AM !magdip
  implicit none

  complex(kind=8) :: intMagHyp_pt_Qmat
  integer,intent(in) :: nn,mm
  real(kind=dp),intent(in) :: posR(3)

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
  call copy_DiracOutput_pg(pg)
  
  tmp_mat(:,:) = (0._dp,0._dp)
     call calc_intE_pq(1,4,posR,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intE_pq_14) 
     call calc_intE_pq(2,3,posR,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intE_pq_23) 
     call calc_intE_pq(3,2,posR,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intE_pq_32) 
     call calc_intE_pq(4,1,posR,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intE_pq_41) 
     call calc_intE_pq(1,3,posR,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intE_pq_13) 
     call calc_intE_pq(2,4,posR,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intE_pq_24) 
     call calc_intE_pq(3,1,posR,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intE_pq_31) 
     call calc_intE_pq(4,2,posR,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),intE_pq_42) 
    
     do k=1,3
        tmp_mat(1,k) = intE_pq_14(k) +intE_pq_23(k) +intE_pq_32(k) +intE_pq_41(k) 
        tmp_mat(2,k) = -IU*intE_pq_14(k) +IU*intE_pq_23(k) -IU*intE_pq_32(k) +IU*intE_pq_41(k)
        tmp_mat(3,k) = intE_pq_13(k) -intE_pq_24(k) +intE_pq_31(k) -intE_pq_42(k)
     end do

     intMagHyp_pt_Qmat = magdip(1)*( tmp_mat(2,3)-tmp_mat(3,2) )&
                       &+magdip(2)*( tmp_mat(3,1)-tmp_mat(1,3) )&
                       &+magdip(3)*( tmp_mat(1,2)-tmp_mat(2,1) )

     intMagHyp_pt_Qmat = Ze*intMagHyp_pt_Qmat

  return
end function intMagHyp_pt_Qmat

!============================================================
function Hso_ele_pt_Qmat(x,y,z,pp,qq,rr,ss)
! Hso_ele = \frac{i Z_ee \hbar^2}{4m^2c^2} \epsilon_{ijk} \psi^\dagger {E}^i \Sigma^k \partial_j \psi 
! contribution from eleleus
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use DefineTypes
  use DiracOutput
  implicit none
 
  complex(kind=dp) :: Hso_ele_pt_Qmat
  real(kind=dp),intent(in) :: x,y,z
  integer,intent(in) :: pp,qq,rr,ss

  integer :: i,k
  complex(kind=dp) :: Eele_pt_Qmat(3)
  complex(kind=8) :: js_pt_Qmat !function

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_r(4,NMAX_PG),c_s(4,NMAX_PG)

  call calc_Eele_pt_Qmat(x,y,z,rr,ss,Eele_pt_Qmat)

  Hso_ele_pt_Qmat = Eele_pt_Qmat(1)*(js_pt_Qmat(2,3,x,y,z,pp,qq)-js_pt_Qmat(3,2,x,y,z,pp,qq)) &
                  &+Eele_pt_Qmat(2)*(js_pt_Qmat(3,1,x,y,z,pp,qq)-js_pt_Qmat(1,3,x,y,z,pp,qq)) &
                  &+Eele_pt_Qmat(3)*(js_pt_Qmat(1,2,x,y,z,pp,qq)-js_pt_Qmat(2,1,x,y,z,pp,qq))
  
  return
end function Hso_ele_pt_Qmat

!============================================================
function Hso_nuc_pt_Qmat(x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef
  implicit none

  complex(kind=8) :: Hso_nuc_pt_Qmat
  complex(kind=8) :: js_pt_Qmat !function
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm
  real(kind=dp) :: vecR(3)
  real(kind=8) :: posRA(3)
  real(kind=dp) :: vecD(3) ! vecR-vecRA
  real(kind=dp) :: norm !  |vecR-vecRA|
  integer k,iN
  complex(kind=dp) :: Enuc(3)

  vecR(1) = x; vecR(2) = y; vecR(3) = z

  Enuc(:) = (0._dp,0._dp)
  do k=1,3
     do iN=1,NAT
        posRA(1) = xc(iN); posRA(2) = yc(iN); posRA(3) = zc(iN)
        vecD(:) = vecR(:)-posRA(:)
        norm = sqrt(vecD(1)**2+vecD(2)**2+vecD(3)**2)
        Enuc(k) = Enuc(k) +cn(iN)*vecD(k)/norm**3
!        write(*,*) iN,posRA(:),vecD(:),norm
     end do
  end do
!  write(*,*) x,y,z,Enuc(1),Enuc(2),Enuc(3)

  Hso_nuc_pt_Qmat = Enuc(1)*(js_pt_Qmat(2,3,x,y,z,nn,mm)-js_pt_Qmat(3,2,x,y,z,nn,mm)) &
                  &+Enuc(2)*(js_pt_Qmat(3,1,x,y,z,nn,mm)-js_pt_Qmat(1,3,x,y,z,nn,mm)) &
                  &+Enuc(3)*(js_pt_Qmat(1,2,x,y,z,nn,mm)-js_pt_Qmat(2,1,x,y,z,nn,mm))

  return
end function Hso_nuc_pt_Qmat

!============================================================
function js_pt_Qmat(i,j,x,y,z,nn,mm)
! ij-th component of js(ab)_pq , i,j=1-3
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
! spin current = -(i/2) \psi^\dagger \Sigma^i \partial_j \psi
! non Hermite
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef
  implicit none

  complex(kind=8) :: js_pt_Qmat
  integer,intent(in) :: i,j
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg
  

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in js_pt_mat."
     stop
  end if
  if(j.lt.1 .or. j.gt.3) then
     write(*,*) "j should be 1-3 in js_pt_mat."
     stop
  end if

  call copy_DiracOutput_pg(pg)
  call calc_js_pq(i,j,x,y,z,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),js_pt_Qmat)
 
  return
end function js_pt_Qmat

!============================================================
function zeta_pt_Qmat(i,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef
  implicit none

  complex(kind=8) :: zeta_pt_Qmat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg
  
  call copy_DiracOutput_pg(pg)
  call calc_zeta_pq(i,x,y,z,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),zeta_pt_Qmat)

  return
end function zeta_pt_Qmat

!============================================================
function t_pt_Qmat(i,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef
  implicit none

  complex(kind=8) :: t_pt_Qmat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  integer :: n,m ! index for 1~2*NBS (+ or - specified)

  type(primitive_gaussian) :: pg

  complex(kind=8) :: tau(3,3)
  integer :: k,l
  
  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in t_pt_Qmat."
     stop
  end if

  call copy_DiracOutput_pg(pg)
  do k=1,3
     do l=1,3
        call calc_tau_pq(k,l,x,y,z,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),tau(k,l))
     end do
  end do
  if(i.eq.1) then
     t_pt_Qmat = -tau(2,3)+tau(3,2)
  elseif(i.eq.2) then
     t_pt_Qmat =  tau(1,3)-tau(3,1)
  elseif(i.eq.3) then
     t_pt_Qmat = -tau(1,2)+tau(2,1)
  else
     write(*,*) "i should be 1-3 in t_pt_Qmat."
     stop
  end if

  return
end function t_pt_Qmat

!============================================================
function tau_pt_Qmat(k,l,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef
  implicit none

  complex(kind=8) :: tau_pt_Qmat
  integer,intent(in) :: k,l
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg

  if(k.lt.1 .or. k.gt.3 .or. l.lt.1 .or. l.gt.3) then
     write(*,*) "k and l should be 1-3 in tau_pt_Qmat."
     stop
  end if
  
  call copy_DiracOutput_pg(pg)
  call calc_tau_pq(k,l,x,y,z,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),tau_pt_Qmat)
  

  return
end function tau_pt_Qmat

!============================================================
function divj_pt_Qmat(x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef
  implicit none

  complex(kind=8) :: divj_pt_Qmat
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg
  
  call copy_DiracOutput_pg(pg)
  call calc_divj_pq(x,y,z,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),divj_pt_Qmat)

  return
end function divj_pt_Qmat

!============================================================
function j_pt_Qmat(i,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef
  implicit none

  complex(kind=8) :: j_pt_Qmat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in j_pt_Qmat."
     stop
  end if
  
  call copy_DiracOutput_pg(pg)
  call calc_j_pq(i,x,y,z,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),j_pt_Qmat)
  
  return
end function j_pt_Qmat

!============================================================
function s_pt_Qmat(i,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef
  implicit none

  complex(kind=8) :: s_pt_Qmat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in s_pt_Qmat."
     stop
  end if
  
  call copy_DiracOutput_pg(pg)
  call calc_s_pq(i,x,y,z,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),s_pt_Qmat)

  return
end function s_pt_Qmat

!============================================================
function ssmall_pt_Qmat(i,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef
  implicit none

  complex(kind=8) :: ssmall_pt_Qmat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in ssmall_pt_Qmat."
     stop
  end if
  
  call copy_DiracOutput_pg(pg)
  call calc_ssmall_pq(i,x,y,z,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),ssmall_pt_Qmat)

  return
end function ssmall_pt_Qmat

!============================================================
subroutine calc_ssmall_pq(i,x,y,z,NL,NS,pg,c_p,c_q,ssmall_pq_i)
! calculate i-th component of ssmall_pq
!============================================================
  use DefineTypes
  use Constants
  implicit none
  
  integer,intent(in) :: i   !s^i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=8),intent(out) :: ssmall_pq_i

  complex(kind=8) :: psi_p(4),psi_q(4)

  call calc_psi(x,y,z,NL,NS,pg,c_p,psi_p)
  call calc_psi(x,y,z,NL,NS,pg,c_q,psi_q)

  ssmall_pq_i = (0.d0,0.d0)
  if(i==1) then
!    ssmall_pq_i = conjg(psi_p(4))*psi_q(3) + conjg(psi_p(4))*psi_q(3)
    ssmall_pq_i = conjg(psi_p(4))*psi_q(3) + conjg(psi_p(3))*psi_q(4) !2014/11/1 corrected
  else if(i==2) then 
!    ssmall_pq_i = -IU*conjg(psi_p(4))*psi_q(3) + IU*conjg(psi_p(4))*psi_q(3)
    ssmall_pq_i =  IU*conjg(psi_p(4))*psi_q(3) - IU*conjg(psi_p(3))*psi_q(4) !2014/11/1 corrected
  else if(i==3) then 
    ssmall_pq_i = conjg(psi_p(3))*psi_q(3) - conjg(psi_p(4))*psi_q(4)
  else
     write(*,*) "i should be 1-3 in calc_ssmall_pq."
  end if

  ssmall_pq_i = ssmall_pq_i/2.d0

  return
end subroutine calc_ssmall_pq

!============================================================
function curl_spin_pt_Qmat(i,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef
  implicit none

  complex(kind=8) :: curl_spin_pt_Qmat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  integer :: n,m ! index for 1~2*NBS (+ or - specified)

  type(primitive_gaussian) :: pg

  complex(kind=8) :: grad_spin(3,3)
  integer :: k,l
  
  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in curl_spin_pt_Qmat."
     stop
  end if

  call copy_DiracOutput_pg(pg)
  do k=1,3
     do l=1,3
        call calc_grad_spin_pq(k,l,x,y,z,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),grad_spin(k,l))
     end do
  end do
  if(i.eq.1) then
     curl_spin_pt_Qmat = grad_spin(2,3) - grad_spin(3,2)
  else if(i.eq.2) then
     curl_spin_pt_Qmat = grad_spin(3,1) - grad_spin(1,3)
  else if(i.eq.3) then
     curl_spin_pt_Qmat = grad_spin(1,2) - grad_spin(2,1)
  else
     write(*,*) "i should be 1-3 in curl_spin_pt_Qmat."
     stop
  end if

  return
end function curl_spin_pt_Qmat

!============================================================
subroutine calc_grad_spin_pq(k,l,x,y,z,NL,NS,pg,c_p,c_q,grad_spin_pq_kl)
! calculate (k,l) component of \partial_k spin^l
!============================================================
  use DefineTypes
  use Constants
  implicit none
  
  integer,intent(in) :: k,l  ! partial_k spin^l
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=8),intent(out) :: grad_spin_pq_kl

  complex(kind=8) :: psi_p(4),psi_q(4)
  complex(kind=8) :: dpsi_p(3,4),dpsi_q(3,4)  ! 3:x,y,z
  integer :: m1,m2 ! spinor indice

  call calc_psi(x,y,z,NL,NS,pg,c_p,psi_p)
  call calc_psi(x,y,z,NL,NS,pg,c_q,psi_q)
  call calc_dpsi(x,y,z,NL,NS,pg,c_p,dpsi_p)
  call calc_dpsi(x,y,z,NL,NS,pg,c_q,dpsi_q)

  grad_spin_pq_kl = (0.d0,0.d0)
  do m1=1,4
     do m2=1,4
        grad_spin_pq_kl = grad_spin_pq_kl &
             +conjg(psi_p(m1))*Sigma(l,m1,m2)*dpsi_q(k,m2) &
             +conjg(dpsi_p(k,m1))*Sigma(l,m1,m2)*psi_q(m2) 
     end do
  end do

  grad_spin_pq_kl = grad_spin_pq_kl/2.d0

  return
end subroutine calc_grad_spin_pq

!============================================================
function grad_rho_pt_Qmat(i,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef
  implicit none

  complex(kind=8) :: grad_rho_pt_Qmat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in grad_rho_pt_Qmat."
     stop
  end if
  
  call copy_DiracOutput_pg(pg)
  call calc_grad_rho_pq(i,x,y,z,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),grad_rho_pt_Qmat)

  return
end function grad_rho_pt_Qmat

!============================================================
subroutine calc_grad_rho_pq(i,x,y,z,NL,NS,pg,c_p,c_q,grad_rho_pq_i)
! calculate i-th component of grad rho
!============================================================
  use DefineTypes
  use Constants
  implicit none
  
  integer,intent(in) :: i  ! \partial_i rho
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=8),intent(out) :: grad_rho_pq_i

  complex(kind=8) :: psi_p(4),psi_q(4)
  complex(kind=8) :: dpsi_p(3,4),dpsi_q(3,4)  ! 3:x,y,z
  integer :: m1 ! spinor indice

  call calc_psi(x,y,z,NL,NS,pg,c_p,psi_p)
  call calc_psi(x,y,z,NL,NS,pg,c_q,psi_q)
  call calc_dpsi(x,y,z,NL,NS,pg,c_p,dpsi_p)
  call calc_dpsi(x,y,z,NL,NS,pg,c_q,dpsi_q)

  grad_rho_pq_i = (0.d0,0.d0)
  do m1=1,4
    grad_rho_pq_i = grad_rho_pq_i +conjg(psi_p(m1))*dpsi_q(i,m1) +conjg(dpsi_p(i,m1))*psi_q(m1)
  end do

  return
end subroutine calc_grad_rho_pq

!============================================================
function t5_pt_Qmat(i,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
! i-th component of \epsilon_{ijk} \tau^{jk}_5
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef
  implicit none

  complex(kind=8) :: t5_pt_Qmat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  integer :: n,m ! index for 1~2*NBS (+ or - specified)

  type(primitive_gaussian) :: pg

  complex(kind=8) :: tau5(3,3)
  integer :: k,l
  
  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in t5_pt_Qmat."
     stop
  end if

  call copy_DiracOutput_pg(pg)
  do k=1,3
     do l=1,3
        call calc_tau5_pq(k,l,x,y,z,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),tau5(k,l))
     end do
  end do
  if(i.eq.1) then
     t5_pt_Qmat = tau5(2,3) - tau5(3,2)
  else if(i.eq.2) then
     t5_pt_Qmat = tau5(3,1) - tau5(1,3)
  else if(i.eq.3) then
     t5_pt_Qmat = tau5(1,2) - tau5(2,1)
  else
     write(*,*) "i should be 1-3 in t5_pt_Qmat."
     stop
  end if

  return
end function t5_pt_Qmat

!============================================================
subroutine calc_tau5_pq(k,l,x,y,z,NL,NS,pg,c_p,c_q,tau5_pq_kl)
! calculate (k,l) component of (c/2)*(\psi^\dagger \Sigma^l IU*\hbar \partial_k psi + h.c.)
!============================================================
  use DefineTypes
  use Constants
  implicit none
  
  integer,intent(in) :: k,l  ! (c/2)*(\psi^\dagger \Sigma^k IU*\hbar \partial_l psi + h.c.)
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=8),intent(out) :: tau5_pq_kl

  complex(kind=8) :: psi_p(4),psi_q(4)
  complex(kind=8) :: dpsi_p(3,4),dpsi_q(3,4)  ! 3:x,y,z
  integer :: m1,m2 ! spinor indice

  call calc_psi(x,y,z,NL,NS,pg,c_p,psi_p)
  call calc_psi(x,y,z,NL,NS,pg,c_q,psi_q)
  call calc_dpsi(x,y,z,NL,NS,pg,c_p,dpsi_p)
  call calc_dpsi(x,y,z,NL,NS,pg,c_q,dpsi_q)

  tau5_pq_kl = (0.d0,0.d0)
  do m1=1,4
     do m2=1,4
        tau5_pq_kl = tau5_pq_kl&
            & +conjg(psi_p(m1))*Sigma(l,m1,m2)*dpsi_q(k,m2) &
            & -conjg(dpsi_p(k,m1))*Sigma(l,m1,m2)*psi_q(m2) 
     end do
  end do

  tau5_pq_kl = IU*CCC*tau5_pq_kl/2.d0

  return
end subroutine calc_tau5_pq

!============================================================
function ipdgamip_pt_Qmat(i,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef
  implicit none

  complex(kind=8) :: ipdgamip_pt_Qmat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in ipdgamip_pt_Qmat."
     stop
  end if
  
  call copy_DiracOutput_pg(pg)
  call calc_ipdgamip_pq(i,x,y,z,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),ipdgamip_pt_Qmat)

  return
end function ipdgamip_pt_Qmat

!============================================================
subroutine calc_ipdgamip_pq(i,x,y,z,NL,NS,pg,c_p,c_q,ipdgamip_pq_i)
! calculate i-th component of IU *\psi^\dagger \gamma^i \psi
!============================================================
  use DefineTypes
  use Constants
  implicit none
  
  integer,intent(in) :: i  ! IU *\psi^\dagger \gamma^i \psi
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=8),intent(out) :: ipdgamip_pq_i

  complex(kind=8) :: psi_p(4),psi_q(4)
  integer :: m1,m2 ! spinor indice

  call calc_psi(x,y,z,NL,NS,pg,c_p,psi_p)
  call calc_psi(x,y,z,NL,NS,pg,c_q,psi_q)

  ipdgamip_pq_i = (0.d0,0.d0)
  do m1=1,4
    do m2=1,4
      ipdgamip_pq_i = ipdgamip_pq_i +conjg(psi_p(m1))*Gam(i,m1,m2)*psi_q(m2)
    end do
  end do

  ipdgamip_pq_i = IU *ipdgamip_pq_i

  return
end subroutine calc_ipdgamip_pq

!============================================================
function rho_pt_Qmat(x,y,z,nn,mm)
! Nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef
  implicit none

  complex(kind=8) :: rho_pt_Qmat
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg
  
  call copy_DiracOutput_pg(pg)
  call calc_rho_pq(x,y,z,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),rho_pt_Qmat) 
  

  return
end function rho_pt_Qmat

!============================================================
function density_pt_Qmat(x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef
  implicit none

  complex(kind=8) :: density_pt_Qmat
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg
  
  call copy_DiracOutput_pg(pg)
  call calc_density_pq(x,y,z,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),density_pt_Qmat) 

  return
end function density_pt_Qmat

!============================================================
function chiral_pt_Qmat(x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef
  implicit none

  complex(kind=8) :: chiral_pt_Qmat
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg
  
  call copy_DiracOutput_pg(pg)
  call calc_chiral_pq(x,y,z,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),chiral_pt_Qmat)

  return
end function chiral_pt_Qmat

!============================================================
function tau_AM_pt_Qmat(t,k,l,x,y,z,nn,mm)
! t : time = timestep*DeltaT
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================

  use Constants  ! use CCC, Nph
  use DefineTypes
  use DiracOutput
  use PTcoef
  implicit none

  complex(kind=8) :: tau_AM_pt_Qmat
  real(kind=8),intent(in) :: t  ! time
  integer,intent(in) :: k,l
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=8) :: j_pt_Qmat,tau_pt_Qmat,func_AMvec

  tau_AM_pt_Qmat = j_pt_Qmat(l,x,y,z,nn,mm)*func_AMvec(t,k,x,y,z)/CCC

  return
end function tau_AM_pt_Qmat

!============================================================
function t_AM_pt_Qmat(t,i,x,y,z,nn,mm)
! t : time = timestep*DeltaT
! i-th component of t(ab)_pq, i=1,3
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use DefineTypes
  use DiracOutput
  use Constants  ! use Nph
  implicit none
 
  complex(kind=8) :: t_AM_pt_Qmat
  real(kind=8),intent(in) :: t  ! time
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=8) :: tau_AM_pt_Qmat

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in t_AM_pt_Qmat."
     stop
  end if

  if(i.eq.1) then
     t_AM_pt_Qmat = -tau_AM_pt_Qmat(t,2,3,x,y,z,nn,mm) +tau_AM_pt_Qmat(t,3,2,x,y,z,nn,mm)
  elseif(i.eq.2) then
     t_AM_pt_Qmat =  tau_AM_pt_Qmat(t,1,3,x,y,z,nn,mm) -tau_AM_pt_Qmat(t,3,1,x,y,z,nn,mm)
  elseif(i.eq.3) then
     t_AM_pt_Qmat = -tau_AM_pt_Qmat(t,1,2,x,y,z,nn,mm) +tau_AM_pt_Qmat(t,2,1,x,y,z,nn,mm)
  else
     write(*,*) "i should be 1-3 in t_AM_pt_Qmat."
     stop
  end if
  
  return
end function t_AM_pt_Qmat

!!$!============================================================
!!$function pi_pt_Qmat(i,x,y,z,nn,mm)
!!$! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!!$!============================================================
!!$  use Precision
!!$  use DefineTypes
!!$  use DiracOutput
!!$  use PTcoef
!!$  implicit none
!!$
!!$  complex(kind=8) :: pi_pt_Qmat
!!$  integer,intent(in) :: i
!!$  real(kind=8),intent(in) :: x,y,z
!!$  integer,intent(in) :: nn,mm
!!$
!!$  type(primitive_gaussian) :: pg
!!$
!!$  if(i.lt.1 .or. i.gt.3) then
!!$     write(*,*) "i should be 1-3 in pi_pt_Qmat."
!!$     stop
!!$  end if
!!$  
!!$  call copy_DiracOutput_pg(pg)
!!$  call calc_pi_pq(i,x,y,z,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),pi_pt_Qmat)
!!$
!!$  return
!!$end function pi_pt_Qmat
!!$
!!$!============================================================
!!$function orbang_pt_Qmat(i,x,y,z,nn,mm)
!!$! i-th component of le_pq, i=1,3
!!$! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!!$!============================================================
!!$  use DefineTypes
!!$  use DiracOutput
!!$  implicit none
!!$ 
!!$  complex(kind=8) :: orbang_pt_Qmat
!!$  complex(kind=8) :: pi_pt_Qmat
!!$  integer,intent(in) :: i
!!$  real(kind=8),intent(in) :: x,y,z
!!$  integer,intent(in) :: nn,mm
!!$
!!$  integer :: k,l
!!$  
!!$  if(i.lt.1 .or. i.gt.3) then
!!$     write(*,*) "i should be 1-3 in orbang_pt_Qmat."
!!$     stop
!!$  end if
!!$
!!$  if(i.eq.1) then
!!$     orbang_pt_Qmat = y *pi_pt_Qmat(3,x,y,z,nn,mm) - z *pi_pt_Qmat(2,x,y,z,nn,mm)
!!$  elseif(i.eq.2) then
!!$     orbang_pt_Qmat = z *pi_pt_Qmat(1,x,y,z,nn,mm) - x *pi_pt_Qmat(3,x,y,z,nn,mm)
!!$  elseif(i.eq.3) then
!!$     orbang_pt_Qmat = x *pi_pt_Qmat(2,x,y,z,nn,mm) - y *pi_pt_Qmat(1,x,y,z,nn,mm)
!!$  else
!!$     write(*,*) "i should be 1-3 in orbang_pt_Qmat."
!!$     stop
!!$  end if
!!$  
!!$  return
!!$end function orbang_pt_Qmat

!============================================================
function AA_pt_Qmat(i,x,y,z,nn,mm)
! A_A^i = Z_ee Int ds psi_nn^dagger gam^0 gam^k psi_mm /|r-s|
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef, only : c_psi
  use Constants
  implicit none

  complex(kind=8) :: AA_pt_Qmat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg

  real(kind=dp) :: posRA(3)  ! position of nucleus
!!$  complex(kind=dp) :: AA_pq
!!$  integer ip,iq,ir
  complex(kind=dp) :: AA_pq_14, AA_pq_23, AA_pq_32, AA_pq_41
  complex(kind=dp) :: AA_pq_13, AA_pq_24, AA_pq_31, AA_pq_42

  call set_GammaMatrix ! for IU

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in AA_pt_Qmat."
     stop
  end if

  posRA(1)=x; posRA(2)=y; posRA(3)=z
  
  call copy_DiracOutput_pg(pg)
  AA_pt_Qmat = (0._dp,0._dp)
!!$  ! psi^a gam^ab gam^bc psi^c
!!$  do ip=1,4
!!$    do iq=1,4
!!$      do ir=1,4
!!$        call calc_AA_pq(ip,ir,posRA,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),AA_pq)
!!$        AA_pt_Qmat = AA_pt_Qmat + AA_pq*Gam0(ip,iq)*Gam(i,iq,ir)
!!$      end do
!!$    end do
!!$  end do
  if(i==1) then
    call calc_AA_pq(1,4,posRA,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),AA_pq_14)
    call calc_AA_pq(2,3,posRA,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),AA_pq_23)
    call calc_AA_pq(3,2,posRA,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),AA_pq_32)
    call calc_AA_pq(4,1,posRA,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),AA_pq_41)
      AA_pt_Qmat = AA_pq_14 + AA_pq_23 + AA_pq_32 + AA_pq_41
  else if(i==2) then
    call calc_AA_pq(1,4,posRA,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),AA_pq_14)
    call calc_AA_pq(2,3,posRA,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),AA_pq_23)
    call calc_AA_pq(3,2,posRA,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),AA_pq_32)
    call calc_AA_pq(4,1,posRA,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),AA_pq_41)
      AA_pt_Qmat = IU*(-AA_pq_14 + AA_pq_23 - AA_pq_32 + AA_pq_41)
  else if(i==3) then
    call calc_AA_pq(1,3,posRA,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),AA_pq_13)
    call calc_AA_pq(2,4,posRA,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),AA_pq_24)
    call calc_AA_pq(3,1,posRA,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),AA_pq_31)
    call calc_AA_pq(4,2,posRA,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),AA_pq_42)
      AA_pt_Qmat = AA_pq_13 - AA_pq_24 + AA_pq_31 - AA_pq_42
  end if

  AA_pt_Qmat = AA_pt_Qmat*Ze

  return
end function AA_pt_Qmat

!============================================================
function tau_AA_pt_Qmat(k,l,x,y,z,nn,mm,pp,qq)
! nn,mm,pp,qq : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================

  use Constants  ! use CCC, Nph
  use DefineTypes
  use DiracOutput
  use PTcoef, only : c_psi
  implicit none

  complex(kind=8) :: tau_AA_pt_Qmat
  integer,intent(in) :: k,l
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm,pp,qq

  complex(kind=8) :: j_pt_Qmat,AA_pt_Qmat !function

  tau_AA_pt_Qmat = j_pt_Qmat(l,x,y,z,nn,mm)*AA_pt_Qmat(k,x,y,z,pp,qq)/CCC

  return
end function tau_AA_pt_Qmat

!============================================================
function t_AA_pt_Qmat(i,x,y,z,nn,mm,pp,qq)
! i-th component of t(ab)_pq, i=1,3
! nn,mm,pp,qq : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use DefineTypes
  use DiracOutput
  use Constants  ! use Nph
  implicit none
 
  complex(kind=8) :: t_AA_pt_Qmat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm,pp,qq

  complex(kind=8) :: tau_AA_pt_Qmat

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in t_AA_pt_Qmat."
     stop
  end if

  if(i.eq.1) then
     t_AA_pt_Qmat = -tau_AA_pt_Qmat(2,3,x,y,z,nn,mm,pp,qq) +tau_AA_pt_Qmat(3,2,x,y,z,nn,mm,pp,qq)
  elseif(i.eq.2) then
     t_AA_pt_Qmat =  tau_AA_pt_Qmat(1,3,x,y,z,nn,mm,pp,qq) -tau_AA_pt_Qmat(3,1,x,y,z,nn,mm,pp,qq)
  elseif(i.eq.3) then
     t_AA_pt_Qmat = -tau_AA_pt_Qmat(1,2,x,y,z,nn,mm,pp,qq) +tau_AA_pt_Qmat(2,1,x,y,z,nn,mm,pp,qq)
  else
     write(*,*) "i should be 1-3 in t_AA_pt_Qmat."
     stop
  end if
  
  return
end function t_AA_pt_Qmat

!============================================================
subroutine calc_Eele_pt_Qmat(x,y,z,nn,mm,Eele_pt_Qmat)
! E^i(r) = Z_ee Int ds psi_nn^dagger  psi_mm (r-s)^i/|r-s|^3
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef, only : c_psi
  use Constants
  implicit none

  complex(kind=8),intent(out) :: Eele_pt_Qmat(3)
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg

  real(kind=dp) :: posRA(3)  ! position of nucleus
  complex(kind=dp) :: Eele_pq_11(3), Eele_pq_22(3), Eele_pq_33(3), Eele_pq_44(3)
  integer :: k

  call set_GammaMatrix ! for IU

  posRA(1)=x; posRA(2)=y; posRA(3)=z
  
  call copy_DiracOutput_pg(pg)
  Eele_pt_Qmat(:) = (0._dp,0._dp)
    call calc_intE_pq(1,1,posRA,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),Eele_pq_11)
    call calc_intE_pq(2,2,posRA,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),Eele_pq_22)
    call calc_intE_pq(3,3,posRA,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),Eele_pq_33)
    call calc_intE_pq(4,4,posRA,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),Eele_pq_44)

  do k=1,3
    Eele_pt_Qmat(k) = Eele_pq_11(k) + Eele_pq_22(k) + Eele_pq_33(k) + Eele_pq_44(k)
    Eele_pt_Qmat(k) = Eele_pt_Qmat(k)*Ze
  end do

  return
end subroutine calc_Eele_pt_Qmat

!============================================================
function spol_pt_Qmat(i,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef, only : c_psi
  implicit none

  complex(kind=8) :: spol_pt_Qmat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in spol_pt_Qmat."
     stop
  end if
  
  call copy_DiracOutput_pg(pg)
  call calc_spol_pq(i,x,y,z,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),spol_pt_Qmat)

  return
end function spol_pt_Qmat

!============================================================
function gamk_pt_Qmat(i,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef, only : c_psi
  implicit none

  complex(kind=8) :: gamk_pt_Qmat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in gamk_pt_Qmat."
     stop
  end if
  
  call copy_DiracOutput_pg(pg)
  call calc_gamk_pq(i,x,y,z,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),gamk_pt_Qmat)

  return
end function gamk_pt_Qmat

!============================================================
function EDMtorqE_pt_Qmat(i,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef, only : c_psi
  use param_AM, only : Eelec
  implicit none

  complex(kind=8) :: EDMtorqE_pt_Qmat !function
  complex(kind=8) :: spol_pt_Qmat !function
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in EDMtorqE_pt_Qmat."
     stop
  end if
  
  if(i.eq.1) then
    EDMtorqE_pt_Qmat = spol_pt_Qmat(2,x,y,z,nn,mm)*Eelec(3) - spol_pt_Qmat(3,x,y,z,nn,mm)*Eelec(2)
  else if(i.eq.2) then
    EDMtorqE_pt_Qmat = spol_pt_Qmat(3,x,y,z,nn,mm)*Eelec(1) - spol_pt_Qmat(1,x,y,z,nn,mm)*Eelec(3)
  else if(i.eq.3) then
    EDMtorqE_pt_Qmat = spol_pt_Qmat(1,x,y,z,nn,mm)*Eelec(2) - spol_pt_Qmat(2,x,y,z,nn,mm)*Eelec(1)
  end if

  EDMtorqE_pt_Qmat = EDMtorqE_pt_Qmat*2._dp

  return
end function EDMtorqE_pt_Qmat

!============================================================
function EDMtorqB_pt_Qmat(i,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef, only : c_psi
  use param_AM, only : Bmag
  use Constants  ! use IU
  implicit none

  complex(kind=8) :: EDMtorqB_pt_Qmat !function
  complex(kind=8) :: gamk_pt_Qmat !function
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  type(primitive_gaussian) :: pg

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in EDMtorqB_pt_Qmat."
     stop
  end if
  
  if(i.eq.1) then
    EDMtorqB_pt_Qmat = gamk_pt_Qmat(2,x,y,z,nn,mm)*Bmag(3) - gamk_pt_Qmat(3,x,y,z,nn,mm)*Bmag(2)
  else if(i.eq.2) then
    EDMtorqB_pt_Qmat = gamk_pt_Qmat(3,x,y,z,nn,mm)*Bmag(1) - gamk_pt_Qmat(1,x,y,z,nn,mm)*Bmag(3)
  else if(i.eq.3) then
    EDMtorqB_pt_Qmat = gamk_pt_Qmat(1,x,y,z,nn,mm)*Bmag(2) - gamk_pt_Qmat(2,x,y,z,nn,mm)*Bmag(1)
  end if

  EDMtorqB_pt_Qmat = IU*EDMtorqB_pt_Qmat

  return
end function EDMtorqB_pt_Qmat

!============================================================
function Eeff_EDM_ele_pt_Qmat(x,y,z,nn,mm,pp,qq)
! nn,mm,pp,qq : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================

  use Precision
  use Constants  ! use CCC, Nph
  use DefineTypes
  use DiracOutput
  use PTcoef, only : c_psi
  implicit none

  complex(kind=8) :: Eeff_EDM_ele_pt_Qmat
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm,pp,qq

  integer :: i
  complex(kind=8) :: Eele_pt_Qmat(3)
  complex(kind=8) :: spol_pt_Qmat !function

  Eeff_EDM_ele_pt_Qmat = (0._dp,0._dp)
  call calc_Eele_pt_Qmat(x,y,z,pp,qq,Eele_pt_Qmat)
  do i=1,3
    Eeff_EDM_ele_pt_Qmat = Eeff_EDM_ele_pt_Qmat + spol_pt_Qmat(i,x,y,z,nn,mm)*Eele_pt_Qmat(i) *2.d0 ! for spol is (1/2)psi^+ gam0 Sigma^i psi
  end do

  return
end function Eeff_EDM_ele_pt_Qmat

!============================================================
function Eeff_EDM_nuc_pt_Qmat(x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================

  use Precision
  use Constants  ! use CCC, Nph
  use DefineTypes
  use DiracOutput
  use PTcoef, only : c_psi
  implicit none

  complex(kind=8) :: Eeff_EDM_nuc_pt_Qmat
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  integer :: i,iN
  complex(kind=8) :: Enuc(3)
  real(kind=dp) :: vecR(3)
  real(kind=8) :: posRA(3)
  real(kind=dp) :: vecD(3) ! vecR-vecRA
  real(kind=dp) :: norm !  |vecR-vecRA|
  complex(kind=8) :: spol_pt_Qmat !function

  Eeff_EDM_nuc_pt_Qmat = (0._dp,0._dp)
     Enuc(:) = (0._dp,0._dp)
     do iN=1,NAT
        posRA(1) = xc(iN); posRA(2) = yc(iN); posRA(3) = zc(iN)
        vecR(1) = x; vecR(2) = y; vecR(3) = z
        vecD(:) = posRA(:)-vecR(:)
        norm = sqrt(vecD(1)**2+vecD(2)**2+vecD(3)**2)
        do i=1,3
           Enuc(i) = Enuc(i) +cn(iN)*vecD(i)/norm**3
        end do
     end do
!        write(*,*) iN,posRA(:),vecD(:),norm

  do i=1,3
    Eeff_EDM_nuc_pt_Qmat = Eeff_EDM_nuc_pt_Qmat + spol_pt_Qmat(i,x,y,z,nn,mm)*Enuc(i) *2._dp ! for spol is (1/2)psi^+ gam0 Sigma^i psi
  end do

  return
end function Eeff_EDM_nuc_pt_Qmat

!============================================================
function Eeff_EDM_nuc2_pt_Qmat(x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================

  use Precision
  use Constants  ! use CCC, Nph
  use DefineTypes
  use DiracOutput
  use PTcoef, only : c_psi
  implicit none

  complex(kind=dp) :: Eeff_EDM_nuc2_pt_Qmat
  real(kind=dp),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  integer :: i,iN
  complex(kind=dp) :: Enuc(3)
  real(kind=dp) :: vecR(3)
  real(kind=dp) :: posRA(3)
  real(kind=dp) :: vecD(3) ! vecR-vecRA
  real(kind=dp) :: norm !  |vecR-vecRA|
  complex(kind=dp) :: spol_pt_Qmat !function
  complex(kind=dp) :: ssmall_pt_Qmat !function

  Eeff_EDM_nuc2_pt_Qmat = (0._dp,0._dp)
     Enuc(:) = (0._dp,0._dp)
     do iN=1,NAT
        posRA(1) = xc(iN); posRA(2) = yc(iN); posRA(3) = zc(iN)
        vecR(1) = x; vecR(2) = y; vecR(3) = z
        vecD(:) = posRA(:)-vecR(:)
        norm = sqrt(vecD(1)**2+vecD(2)**2+vecD(3)**2)
        do i=1,3
           Enuc(i) = Enuc(i) +cn(iN)*vecD(i)/norm**3
        end do
     end do
!        write(*,*) iN,posRA(:),vecD(:),norm

  do i=1,3
    Eeff_EDM_nuc2_pt_Qmat = Eeff_EDM_nuc2_pt_Qmat - ssmall_pt_Qmat(i,x,y,z,nn,mm)*Enuc(i) *4._dp
  end do

  return
end function Eeff_EDM_nuc2_pt_Qmat

!============================================================
subroutine calc_tau_AA_ansatz1(x,y,z,actorb0,actorb,occ,tau_AA_pt)
! nn,mm,pp,qq : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================

  use Precision
  use Constants  ! use CCC, Nph
  use DefineTypes
  use DiracOutput
  use PTcoef, only : c_psi
  implicit none

  complex(kind=dp), intent(out) :: tau_AA_pt(3,3)
  real(kind=dp),intent(in) :: x,y,z
  integer,intent(in) :: actorb0,actorb
  real(kind=dp),intent(in) :: occ(actorb)

  type(primitive_gaussian) :: pg

  integer :: nn,mm,pp,qq
  integer :: i,j,k,l,ip,iq

  complex(kind=dp) :: j_pt_Qmat !function
  real(kind=dp) :: posR(3)  
  integer :: NLS !function
  
  real(kind=dp) :: overlap_normpg,ef_normpg(3) ! --> dummy variables here 
  real(kind=dp) :: nucatt_normpg
  real(kind=dp) :: nucatt_pg(4,4,NBS0,NBS0)
  complex(kind=dp) :: AAptQmat(3,actorb,actorb)

  call set_GammaMatrix ! for IU, CCC
  posR(1)=x; posR(2)=y; posR(3)=z

  call copy_DiracOutput_pg(pg)

  nucatt_pg(:,:,:,:) = 0._dp

!  write(*,*)"ck1"
  do ip=1,3,2
    do iq=1,3,2
      do i=1,NLS(ip,NBS_L,NBS_S)
        do j=1,NLS(iq,NBS_L,NBS_S)
          call save_gauss_int(ip,iq,i,j,posR,NBS_L,NBS_S,pg,overlap_normpg,nucatt_normpg,ef_normpg)
!          write(*,*)'ip iq i j',ip,iq,i,j
          do k=1,2
            nucatt_pg(ip+k-1,iq+k-1,i,j)=nucatt_normpg
          end do
        end do
      end do
    end do
  end do
!  write(*,*)'nucatt_pg2'

!!$  !LL
!!$  ip=1; iq=1
!!$  do i=1,NLS(ip,NL,NS)
!!$    do j=1,NLS(iq,NL,NS)
!!$      call save_gauss_int(ip,iq,i,j,posR,NL,NS,pg,overlap_normpg,nucatt_normpg,ef_normpg)
!!$      do k=0,1
!!$        nucatt_pg(ip+k,iq+k,i,j)=nucatt_normpg
!!$      end do
!!$    end do
!!$  end do
!!$  !LS
!!$  ip=1; iq=3
!!$  do i=1,NLS(ip,NL,NS)
!!$    do j=1,NLS(iq,NL,NS)
!!$      call save_gauss_int(ip,iq,i,j,posR,NL,NS,pg,overlap_normpg,nucatt_normpg,ef_normpg)
!!$      do k=0,1
!!$        nucatt_pg(ip+k,iq+k,i,j)=nucatt_normpg
!!$      end do
!!$    end do
!!$  end do
!!$  !SL
!!$  ip=3; iq=1
!!$  do i=1,NLS(ip,NL,NS)
!!$    do j=1,NLS(iq,NL,NS)
!!$      call save_gauss_int(ip,iq,i,j,posR,NL,NS,pg,overlap_normpg,nucatt_normpg,ef_normpg)
!!$      do k=0,1
!!$        nucatt_pg(ip+k,iq+k,i,j)=nucatt_normpg
!!$      end do
!!$    end do
!!$  end do
!!$  !SS
!!$  ip=3; iq=3
!!$  do i=1,NLS(ip,NL,NS)
!!$    do j=1,NLS(iq,NL,NS)
!!$      call save_gauss_int(ip,iq,i,j,posR,NL,NS,pg,overlap_normpg,nucatt_normpg,ef_normpg)
!!$      do k=0,1
!!$        nucatt_pg(ip+k,iq+k,i,j)=nucatt_normpg
!!$      end do
!!$    end do
!!$  end do

!  write(*,*)'start calc_AAptQmat'
  call calc_AAptQmat2(actorb0,actorb,pg,nucatt_pg,AAptQmat)
!  call calc_AAptQmat(actorb0,actorb,posR,AAptQmat)
!  write(*,*)'end calc_AAptQmat'

  tau_AA_pt(:,:) = (0._dp,0._dp)
  do nn=actorb0,actorb
    do mm=actorb0,actorb
      do k=1,3
        do l=1,3
          tau_AA_pt(k,l) = tau_AA_pt(k,l) &
                         &  + j_pt_Qmat(l,x,y,z,nn,nn)*AAptQmat(k,mm,mm)*occ(nn)*occ(mm)/CCC
!          tau_AA_pt(k,l) = tau_AA_pt(k,l) &
!                         &  - j_pt_Qmat(l,x,y,z,nn,mm)*AAptQmat(k,mm,nn)*occ(nn)*occ(mm)/CCC
        end do
      end do
    end do
  end do

end subroutine calc_tau_AA_ansatz1

!============================================================
subroutine calc_tau_AA_ansatz2(x,y,z,actorb0,actorb,occ,tau_AA_pt)
! nn,mm,pp,qq : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================

  use Precision
  use Constants  ! use CCC, Nph
  use DefineTypes
  use DiracOutput
  use PTcoef, only : c_psi
  implicit none

  complex(kind=dp), intent(out) :: tau_AA_pt(3,3)
  real(kind=dp),intent(in) :: x,y,z
  integer,intent(in) :: actorb0,actorb
  real(kind=dp),intent(in) :: occ(actorb)

  type(primitive_gaussian) :: pg

  integer :: nn,mm,pp,qq
  integer :: i,j,k,l,ip,iq

  complex(kind=dp) :: j_pt_Qmat !function
  real(kind=dp) :: posR(3)  
  integer :: NLS !function
  
  real(kind=dp) :: overlap_normpg,ef_normpg(3) ! --> dummy variables here 
  real(kind=dp) :: nucatt_normpg
  real(kind=dp) :: nucatt_pg(4,4,NBS0,NBS0)
  complex(kind=dp) :: AAptQmat(3,actorb,actorb)

  call set_GammaMatrix ! for IU, CCC
  posR(1)=x; posR(2)=y; posR(3)=z

  call copy_DiracOutput_pg(pg)

  nucatt_pg(:,:,:,:) = 0._dp

!  write(*,*)"ck1"
  do ip=1,3,2
    do iq=1,3,2
      do i=1,NLS(ip,NBS_L,NBS_S)
        do j=1,NLS(iq,NBS_L,NBS_S)
          call save_gauss_int(ip,iq,i,j,posR,NBS_L,NBS_S,pg,overlap_normpg,nucatt_normpg,ef_normpg)
!          write(*,*)'ip iq i j',ip,iq,i,j
          do k=1,2
            nucatt_pg(ip+k-1,iq+k-1,i,j)=nucatt_normpg
          end do
        end do
      end do
    end do
  end do
!  write(*,*)'nucatt_pg2'

!!$  !LL
!!$  ip=1; iq=1
!!$  do i=1,NLS(ip,NL,NS)
!!$    do j=1,NLS(iq,NL,NS)
!!$      call save_gauss_int(ip,iq,i,j,posR,NL,NS,pg,overlap_normpg,nucatt_normpg,ef_normpg)
!!$      do k=0,1
!!$        nucatt_pg(ip+k,iq+k,i,j)=nucatt_normpg
!!$      end do
!!$    end do
!!$  end do
!!$  !LS
!!$  ip=1; iq=3
!!$  do i=1,NLS(ip,NL,NS)
!!$    do j=1,NLS(iq,NL,NS)
!!$      call save_gauss_int(ip,iq,i,j,posR,NL,NS,pg,overlap_normpg,nucatt_normpg,ef_normpg)
!!$      do k=0,1
!!$        nucatt_pg(ip+k,iq+k,i,j)=nucatt_normpg
!!$      end do
!!$    end do
!!$  end do
!!$  !SL
!!$  ip=3; iq=1
!!$  do i=1,NLS(ip,NL,NS)
!!$    do j=1,NLS(iq,NL,NS)
!!$      call save_gauss_int(ip,iq,i,j,posR,NL,NS,pg,overlap_normpg,nucatt_normpg,ef_normpg)
!!$      do k=0,1
!!$        nucatt_pg(ip+k,iq+k,i,j)=nucatt_normpg
!!$      end do
!!$    end do
!!$  end do
!!$  !SS
!!$  ip=3; iq=3
!!$  do i=1,NLS(ip,NL,NS)
!!$    do j=1,NLS(iq,NL,NS)
!!$      call save_gauss_int(ip,iq,i,j,posR,NL,NS,pg,overlap_normpg,nucatt_normpg,ef_normpg)
!!$      do k=0,1
!!$        nucatt_pg(ip+k,iq+k,i,j)=nucatt_normpg
!!$      end do
!!$    end do
!!$  end do

!  write(*,*)'start calc_AAptQmat'
  call calc_AAptQmat2(actorb0,actorb,pg,nucatt_pg,AAptQmat)
!  call calc_AAptQmat(actorb0,actorb,posR,AAptQmat)
!  write(*,*)'end calc_AAptQmat'

  tau_AA_pt(:,:) = (0._dp,0._dp)
  do nn=actorb0,actorb
    do mm=actorb0,actorb
      do k=1,3
        do l=1,3
          tau_AA_pt(k,l) = tau_AA_pt(k,l) &
                         &  + j_pt_Qmat(l,x,y,z,nn,nn)*AAptQmat(k,mm,mm)*occ(nn)*occ(mm)/CCC
          tau_AA_pt(k,l) = tau_AA_pt(k,l) &
                         &  - j_pt_Qmat(l,x,y,z,nn,mm)*AAptQmat(k,mm,nn)*occ(nn)*occ(mm)/CCC
        end do
      end do
    end do
  end do

end subroutine calc_tau_AA_ansatz2

!============================================================
subroutine calc_tau_AA_noaprox(x,y,z,actorb0,actorb,occ,tau_AA_pt)
! This subroutine DO NOT written completely !!
! nn,mm,pp,qq : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================

  use Precision
  use Constants  ! use CCC, Nph
  use DefineTypes
  use DiracOutput
  use PTcoef, only : c_psi
  use AA_nucatt_pg, only : calEnpqm
  implicit none

  complex(kind=dp), intent(out) :: tau_AA_pt(3,3)
  real(kind=dp),intent(in) :: x,y,z
  integer,intent(in) :: actorb0,actorb
  real(kind=dp),intent(in) :: occ(actorb)

  type(primitive_gaussian) :: pg

  integer :: nn,mm,pp,qq
  integer :: i,j,k,l,ip,iq

  complex(kind=dp) :: j_pt_Qmat !function
  real(kind=dp) :: posR(3)  
  integer :: NLS !function
  
  real(kind=dp) :: overlap_normpg,ef_normpg(3) ! --> dummy variables here 
  real(kind=dp) :: nucatt_normpg
  real(kind=dp) :: nucatt_pg(4,4,NBS0,NBS0)
  complex(kind=dp) :: AAptQmat(3,actorb,actorb)

  call set_GammaMatrix ! for IU, CCC
  posR(1)=x; posR(2)=y; posR(3)=z

  call copy_DiracOutput_pg(pg)

  nucatt_pg(:,:,:,:) = 0._dp

!  write(*,*)"ck1"
  do ip=1,3,2
    do iq=1,3,2
      do i=1,NLS(ip,NBS_L,NBS_S)
        do j=1,NLS(iq,NBS_L,NBS_S)
          call save_gauss_int(ip,iq,i,j,posR,NBS_L,NBS_S,pg,overlap_normpg,nucatt_normpg,ef_normpg)
!          write(*,*)'ip iq i j',ip,iq,i,j
          do k=1,2
            nucatt_pg(ip+k-1,iq+k-1,i,j)=nucatt_normpg
          end do
        end do
      end do
    end do
  end do
!  write(*,*)'nucatt_pg2'

!  write(*,*)'start calc_AAptQmat'
  call calc_AAptQmat2(actorb0,actorb,pg,nucatt_pg,AAptQmat)
!  call calc_AAptQmat(actorb0,actorb,posR,AAptQmat)
!  write(*,*)'end calc_AAptQmat'

  tau_AA_pt(:,:) = (0._dp,0._dp)
  do nn=actorb0,actorb
    do mm=actorb0,actorb
      do pp=actorb0,actorb
        do qq=actorb0,actorb
          do k=1,3
            do l=1,3
              tau_AA_pt(k,l) = tau_AA_pt(k,l) &
                             &  + j_pt_Qmat(l,x,y,z,nn,mm)*AAptQmat(k,pp,qq)*calEnpqm(nn,pp,qq,mm)/CCC
            end do
          end do
        end do
      end do
    end do
  end do

end subroutine calc_tau_AA_noaprox

!============================================================
subroutine calc_AAptQmat3(x,y,z,actorb0,actorb,occ,AA_pt)
! This subroutine DO NOT written completely !!
! nn,mm,pp,qq : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================

  use Precision
  use Constants  ! use CCC, Nph
  use DefineTypes
  use DiracOutput
  use PTcoef, only : c_psi
  use AA_nucatt_pg, only : calEnpqm
  implicit none

  complex(kind=dp), intent(out) :: AA_pt(3)
  real(kind=dp),intent(in) :: x,y,z
  integer,intent(in) :: actorb0,actorb
  real(kind=dp),intent(in) :: occ(actorb)

  type(primitive_gaussian) :: pg

  integer :: nn,mm,pp,qq
  integer :: i,j,k,l,ip,iq

  complex(kind=dp) :: j_pt_Qmat !function
  real(kind=dp) :: posR(3)  
  integer :: NLS !function
  
  real(kind=dp) :: overlap_normpg,ef_normpg(3) ! --> dummy variables here 
  real(kind=dp) :: nucatt_normpg
  real(kind=dp) :: nucatt_pg(4,4,NBS0,NBS0)
  complex(kind=dp) :: AAptQmat(3,actorb,actorb)

  call set_GammaMatrix ! for IU, CCC
  posR(1)=x; posR(2)=y; posR(3)=z

  call copy_DiracOutput_pg(pg)

  nucatt_pg(:,:,:,:) = 0._dp

!  write(*,*)"ck1"
  do ip=1,3,2
    do iq=1,3,2
      do i=1,NLS(ip,NBS_L,NBS_S)
        do j=1,NLS(iq,NBS_L,NBS_S)
          call save_gauss_int(ip,iq,i,j,posR,NBS_L,NBS_S,pg,overlap_normpg,nucatt_normpg,ef_normpg)
!          write(*,*)'ip iq i j',ip,iq,i,j
          do k=1,2
            nucatt_pg(ip+k-1,iq+k-1,i,j)=nucatt_normpg
          end do
        end do
      end do
    end do
  end do
!  write(*,*)'nucatt_pg2'

!  write(*,*)'start calc_AAptQmat'
  call calc_AAptQmat2(actorb0,actorb,pg,nucatt_pg,AAptQmat)
!  call calc_AAptQmat(actorb0,actorb,posR,AAptQmat)
!  write(*,*)'end calc_AAptQmat'

  AA_pt(:) = (0._dp,0._dp)
  do k=1,3
    do nn=actorb0,actorb
      AA_pt(k) = AA_pt(k) + AAptQmat(k,nn,nn)*occ(nn)
    end do
  end do

end subroutine calc_AAptQmat3

!============================================================
subroutine calc_AAptQmat2(actorb0,actorb,pg,nucatt_pg,AAptQmat)
! A_A^i = Z_ee Int ds psi_nn^dagger gam^0 gam^k psi_mm /|r-s|
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef, only : c_psi
  use Constants
  implicit none

  real(kind=dp),intent(in) :: nucatt_pg(4,4,NBS0,NBS0)
  integer,intent(in) :: actorb0,actorb
  complex(kind=dp),intent(out) :: AAptQmat(3,actorb,actorb)

  type(primitive_gaussian),intent(in) :: pg

  integer :: nn,mm
  complex(kind=dp) :: AA_pq_14, AA_pq_23, AA_pq_32, AA_pq_41
  complex(kind=dp) :: AA_pq_13, AA_pq_24, AA_pq_31, AA_pq_42
  integer :: k

  call set_GammaMatrix ! for IU
  
!  call copy_DiracOutput_pg(pg)
  AAptQmat(:,:,:) = (0._dp,0._dp)

  do nn=actorb0,actorb
    do mm=actorb0,actorb
      call calc_AA_pq2(1,4,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),NBS0,nucatt_pg,AA_pq_14)
      call calc_AA_pq2(2,3,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),NBS0,nucatt_pg,AA_pq_23)
      call calc_AA_pq2(3,2,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),NBS0,nucatt_pg,AA_pq_32)
      call calc_AA_pq2(4,1,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),NBS0,nucatt_pg,AA_pq_41)
      call calc_AA_pq2(1,3,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),NBS0,nucatt_pg,AA_pq_13)
      call calc_AA_pq2(2,4,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),NBS0,nucatt_pg,AA_pq_24)
      call calc_AA_pq2(3,1,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),NBS0,nucatt_pg,AA_pq_31)
      call calc_AA_pq2(4,2,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),NBS0,nucatt_pg,AA_pq_42)
      AAptQmat(1,nn,mm) = AA_pq_14 + AA_pq_23 + AA_pq_32 + AA_pq_41
      AAptQmat(2,nn,mm) = IU*(-AA_pq_14 + AA_pq_23 - AA_pq_32 + AA_pq_41)
      AAptQmat(3,nn,mm) = AA_pq_13 - AA_pq_24 + AA_pq_31 - AA_pq_42

      do k=1,3
        AAptQmat(k,nn,mm) = AAptQmat(k,nn,mm)*Ze
      end do
!      write(*,*)nn,AAptQmat(1,nn,1),AAptQmat(2,nn,1),AAptQmat(3,nn,1)
    end do
  end do

  return
end subroutine calc_AAptQmat2

!============================================================
subroutine calc_AAptQmat(actorb0,actorb,posR,nucatt_pg,AAptQmat)
! A_A^i = Z_ee Int ds psi_nn^dagger gam^0 gam^k psi_mm /|r-s|
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use PTcoef, only : c_psi
  use Constants
  implicit none

  real(kind=dp),intent(in) :: nucatt_pg(4,4,NBS0,NBS0)
  integer,intent(in) :: actorb0,actorb
  real(kind=dp),intent(in) :: posR(3)
  complex(kind=dp),intent(out) :: AAptQmat(3,actorb,actorb)

  type(primitive_gaussian) :: pg

  integer :: nn,mm
  complex(kind=dp) :: AA_pq_14, AA_pq_23, AA_pq_32, AA_pq_41
  complex(kind=dp) :: AA_pq_13, AA_pq_24, AA_pq_31, AA_pq_42
  integer :: k

  call set_GammaMatrix ! for IU
  
  call copy_DiracOutput_pg(pg)
  AAptQmat(:,:,:) = (0._dp,0._dp)

  do nn=actorb0,actorb
    do mm=actorb0,actorb
      call calc_AA_pq(1,4,posR,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),AA_pq_14)
      call calc_AA_pq(2,3,posR,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),AA_pq_23)
      call calc_AA_pq(3,2,posR,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),AA_pq_32)
      call calc_AA_pq(4,1,posR,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),AA_pq_41)
      call calc_AA_pq(1,3,posR,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),AA_pq_13)
      call calc_AA_pq(2,4,posR,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),AA_pq_24)
      call calc_AA_pq(3,1,posR,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),AA_pq_31)
      call calc_AA_pq(4,2,posR,NBS_L,NBS_S,pg,c_psi(:,:,nn),c_psi(:,:,mm),AA_pq_42)
      AAptQmat(1,nn,mm) = AA_pq_14 + AA_pq_23 + AA_pq_32 + AA_pq_41
      AAptQmat(2,nn,mm) = IU*(-AA_pq_14 + AA_pq_23 - AA_pq_32 + AA_pq_41)
      AAptQmat(3,nn,mm) = AA_pq_13 - AA_pq_24 + AA_pq_31 - AA_pq_42

      do k=1,3
        AAptQmat(k,nn,mm) = AAptQmat(k,nn,mm)*Ze
      end do
!      write(*,*)nn,AAptQmat(1,nn,1),AAptQmat(2,nn,1),AAptQmat(3,nn,1)
    end do
  end do

!  write(*,*)'end end calc_AAptQmat'

end subroutine calc_AAptQmat

!============================================================
subroutine calc_t_AA_ci(x,y,z,actorb0,actorb,occ,t_AA_pt)
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  implicit none
 
  real(kind=dp),intent(in) :: x,y,z
  integer,intent(in) :: actorb0,actorb
  real(kind=dp),intent(in) :: occ(actorb)
  complex(kind=dp), intent(out) :: t_AA_pt(3)

  complex(kind=dp) :: tau_AA_pt(3,3)

!  write(*,*)'call calc_tau_AA_ci'
!  call calc_tau_AA_ansatz1(x,y,z,actorb0,actorb,occ,tau_AA_pt)
!  call calc_tau_AA_ansatz2(x,y,z,actorb0,actorb,occ,tau_AA_pt)
  call calc_tau_AA_noaprox(x,y,z,actorb0,actorb,occ,tau_AA_pt)

  t_AA_pt(1) = -tau_AA_pt(2,3) +tau_AA_pt(3,2)
  t_AA_pt(2) =  tau_AA_pt(1,3) -tau_AA_pt(3,1)
  t_AA_pt(3) = -tau_AA_pt(1,2) +tau_AA_pt(2,1)
  
end subroutine calc_t_AA_ci
