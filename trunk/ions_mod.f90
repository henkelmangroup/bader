!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module with ion data structure
!-----------------------------------------------------------------------------------!

MODULE ions_mod
  USE kind_mod
  IMPLICIT NONE

TYPE :: ions_obj
  REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: r_car,r_dir,r_lat
  REAL(q2),ALLOCATABLE,DIMENSION(:) :: ion_chg
  REAL(q2),DIMENSION(3,3) :: lattice,dir2car,car2dir
  INTEGER,ALLOCATABLE,DIMENSION(:) :: num_ion
  INTEGER :: niontypes,nions
END TYPE

  PRIVATE
  PUBLIC :: ions_obj

!-----------------------------------------------------------------------------------!

END MODULE ions_mod
