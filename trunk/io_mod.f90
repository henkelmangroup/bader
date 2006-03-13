!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!    Module for reading and writing charge density data
!
! By Andri Arnaldsson and Graeme Henkelman
! Last modified by AA on Aug. 02 2005
!-----------------------------------------------------------------------------------!

MODULE io_mod
  USE kind_mod , ONLY : q1,q2
  USE matrix_mod , ONLY : transpose_matrix,matrix_vector
  USE options_mod 
  USE ions_mod
  USE charge_mod
  USE chgcar_mod
  USE cube_mod
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_charge,write_charge

  CONTAINS

!-----------------------------------------------------------------------------------!
! read_charge: Reads the charge density from a file in vasp or Gaussian cube format,
!    by first reading the header, and then charges
!-----------------------------------------------------------------------------------!

  SUBROUTINE read_charge(ions,chg,opts,chargefile)

    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg
    TYPE(options_obj) :: opts
    CHARACTER(LEN=120) :: chargefile

!    INTEGER :: cr,count_max,t1,t2,nx,ny,nz
         
!    CALL system_clock(t1,cr,count_max)
        
!    vasp=.false.
        
!    OPEN(100,FILE=filename,STATUS='old',ACTION='read',BLANK='null',PAD='yes')
        
!    WRITE(*,'(/,1A11,1A20)') 'OPEN ... ',chargefile
!    READ(100,'(6/,1A7)') text 
!    REWIND(100)
       
    IF (opts%llist%li_chgcar) CALL read_charge_chgcar(ions,chg,chargefile)
    IF (opts%llist%li_cube) CALL read_charge_cube(ions,chg,chargefile)
    
  RETURN
  END SUBROUTINE read_charge

!-----------------------------------------------------------------------------------!
! read_charge: Reads the charge density from a file in vasp or Gaussian cube format,
!    by first reading the header, and then charges
!-----------------------------------------------------------------------------------!

  SUBROUTINE write_charge(ions,chg,opts,chargefile)

    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg
    TYPE(options_obj) :: opts
    CHARACTER(LEN=120) :: chargefile

!    INTEGER :: i
!    INTEGER :: cr,count_max,t1,t2,nx,ny,nz
         
!    CALL system_clock(t1,cr,count_max)
        
!    vasp=.false.
        
!    OPEN(100,FILE=filename,STATUS='old',ACTION='read',BLANK='null',PAD='yes')
        
!    WRITE(*,'(/,1A11,1A20)') 'OPEN ... ',chargefile
!    READ(100,'(6/,1A7)') text 
!    REWIND(100)
       
    IF (opts%llist%li_chgcar) CALL write_charge_chgcar(ions,chg,chargefile)
    IF (opts%llist%li_cube) CALL write_charge_cube(ions,chg,chargefile)
    
  RETURN
  END SUBROUTINE write_charge

!------------------------------------------------------------------------------------!

END MODULE io_mod
