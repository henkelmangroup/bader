!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!    Module with ionic information
!
! By Andri Arnaldson and Graeme Henkelman
! Last modified
!-----------------------------------------------------------------------------------!

MODULE ions_mod
  USE vars_mod , ONLY : q2 
!  USE matrix , ONLY : transpose_matrix,matrix_vector
  IMPLICIT NONE

type, public :: ions_obj
  REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: r_car,r_dir
  REAL(q2),DIMENSION(3) :: corner
  REAL(q2),DIMENSION(3,3) :: lattice
  INTEGER :: niontypes,na,nions
  INTEGER,ALLOCATABLE,DIMENSION(:) :: num_ions,nel
end type

  PRIVATE
!  PUBLIC :: r_car,r_dir,lattice,natypes,na,num_atom,nel

!------------------------------------------------------------------------------------!

END MODULE ions_mod
