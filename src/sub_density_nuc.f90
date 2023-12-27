!============================================================
function rho_nuc_mat(x,y,z,ii,jj)
! atomic nuclei charge density 
!============================================================
  use Constants
  use NucBasis
  implicit none
 
  complex(kind=8) :: rho_nuc_mat
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: ii,jj

  real(kind=8) :: pos(3) ! position of center of gaussian 
  real(kind=8) :: func_pg  ! primitive gaussian function (normalized)

  complex(kind=8) :: chi_i,chi_j
  integer :: i,j,k,s_i,s_j

  chi_i = (0.d0,0.d0)
  chi_j = (0.d0,0.d0)

  call convert_index_nuc(NUCTYPE,ii,i,s_i)
  call convert_index_nuc(NUCTYPE,jj,j,s_j)

  ! use analytic expressions
  if(abs(z).le.(L_well/2.d0)) then
     chi_i = cmplx( sqrt(2.d0*alpha_N/PI)*exp(-alpha_N*y**2)*exp(-alpha_N*x**2)*sqrt(2/L_well)*sin(i*PI/L_well*(z+L_well/2.d0)) )
     chi_j = cmplx( sqrt(2.d0*alpha_N/PI)*exp(-alpha_N*y**2)*exp(-alpha_N*x**2)*sqrt(2/L_well)*sin(j*PI/L_well*(z+L_well/2.d0)) )
  else
     chi_i = (0.d0,0.d0)
     chi_j = (0.d0,0.d0)
  end if
!!$
!!$  ! use expression expanded by gaussian functions (may not be good outside the well)
!!$  do k=1,N_PGN
!!$     pos(:) = vecR_N(:,k)
!!$     chi_i = chi_i +c_nuc(i,k)*func_pg(x,y,z,pos,alpha_N,0,0,0)  ! 1s type
!!$     chi_j = chi_j +c_nuc(j,k)*func_pg(x,y,z,pos,alpha_N,0,0,0)  ! 1s type
!!$  end do

  rho_nuc_mat = Z_N*conjg(chi_i)*chi_j

  return
end function rho_nuc_mat

