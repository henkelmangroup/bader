!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!    Module with ionic information
!
! By Andri Arnaldson and Graeme Henkelman
! Last modified
!-----------------------------------------------------------------------------------!

MODULE ions_mod
  USE kind_mod , ONLY : q2 
  IMPLICIT NONE

TYPE :: ions_obj
  REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: r_car,r_dir
  INTEGER,ALLOCATABLE,DIMENSION(:) :: ion_chg
  REAL(q2),DIMENSION(3,3) :: lattice
  REAL(q2),DIMENSION(3) :: corner
  INTEGER,ALLOCATABLE,DIMENSION(:) :: num_ion
  INTEGER :: niontypes,nions
END TYPE

  PRIVATE
  PUBLIC :: ions_obj

!------------------------------------------------------------------------------------!

END MODULE ions_mod
