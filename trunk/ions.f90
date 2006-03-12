!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!    Module with ionic information
!
! By Andri Arnaldson and Graeme Henkelman
! Last modified
!-----------------------------------------------------------------------------------!

MODULE ions
  USE vars , ONLY : q2 
  USE matrix , ONLY : transpose_matrix,matrix_vector
  IMPLICIT NONE

  REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: r_car,r_dir
  REAL(q2),DIMENSION(3,3) :: lattice
  INTEGER :: natypes,na

  PRIVATE
  PUBLIC :: r_car,r_dir,lattice,natypes,na

!------------------------------------------------------------------------------------!

END MODULE ions
