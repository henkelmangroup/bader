!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!    Module for reading and writing charge density data
!
! By Andri Arnaldson and Graeme Henkelman
! Last modified by AA on Aug. 02 2005
!-----------------------------------------------------------------------------------!

MODULE io_mod
  USE vars_mod , ONLY : q1,q2
  USE options_mod 
  USE matrix_mod , ONLY : transpose_matrix,matrix_vector
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

  SUBROUTINE read_charge(ions,chg,chargefile)

    type(ions_obj) :: ions
    type(charge_obj) :: chg
    CHARACTER(15) :: chargefile

    REAL(q2),DIMENSION(3,3) :: B
    REAL(q2),DIMENSION(3) :: box,v
    REAL(q2) :: side,vol,t
    INTEGER :: i
    INTEGER,DIMENSION(110) :: elements
    CHARACTER(LEN=7) :: text
    INTEGER :: cr,count_max,t1,t2,nx,ny,nz
         
    CALL system_clock(t1,cr,count_max)
        
!    vasp=.false.
        
!    OPEN(100,FILE=filename,STATUS='old',ACTION='read',BLANK='null',PAD='yes')
        
!    WRITE(*,'(/,1A11,1A20)') 'OPEN ... ',chargefile
!    READ(100,'(6/,1A7)') text 
!    REWIND(100)
       
    IF (opts%li_chgcar) CALL read_charge_chgcar(ions,chg,chargefile)
    IF (opts%li_cube) CALL read_charge_cube(ions,chg,chargefile)
    
  RETURN
  END SUBROUTINE read_charge

!-----------------------------------------------------------------------------------!
! read_charge: Reads the charge density from a file in vasp or Gaussian cube format,
!    by first reading the header, and then charges
!-----------------------------------------------------------------------------------!

  SUBROUTINE write_charge(ions,chg,chargefile)

    type(ions_obj) :: ions
    type(charge_obj) :: chg
    CHARACTER(15) :: chargefile

    REAL(q2),DIMENSION(3,3) :: B
    REAL(q2),DIMENSION(3) :: box,v
    REAL(q2) :: side,vol,t
    INTEGER :: i
    INTEGER,DIMENSION(110) :: elements
    CHARACTER(LEN=7) :: text
    INTEGER :: cr,count_max,t1,t2,nx,ny,nz
         
    CALL system_clock(t1,cr,count_max)
        
!    vasp=.false.
        
!    OPEN(100,FILE=filename,STATUS='old',ACTION='read',BLANK='null',PAD='yes')
        
!    WRITE(*,'(/,1A11,1A20)') 'OPEN ... ',chargefile
!    READ(100,'(6/,1A7)') text 
!    REWIND(100)
       
    IF (opts%li_chgcar) CALL read_charge_chgcar(ions,chg,chargefile)
    IF (opts%li_cube) CALL read_charge_cube(ions,chg,chargefile)
    
  RETURN
  END SUBROUTINE write_charge

!------------------------------------------------------------------------------------!

END MODULE io_mod
