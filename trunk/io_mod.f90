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

    CHARACTER(LEN=7) :: text
    INTEGER :: cr,count_max,t1,t2

    CALL system_clock(t1,cr,count_max)

    IF (.NOT. (opts%llist%li_chgcar .OR. opts%llist%li_cube)) THEN
      ! Try to guess the file type
      OPEN(100,FILE=chargefile,STATUS='old',ACTION='read',BLANK='null',PAD='yes')
      READ(100,'(6/,1A7)') text
      CLOSE(100)
      IF (text == 'Direct') THEN
        opts%llist%li_chgcar=.TRUE.
      ELSE 
        opts%llist%li_chgcar=.TRUE.
      ENDIF
    ENDIF

    write(*,*) ' chgcar ??? ',opts%llist%li_chgcar
    write(*,*) ' cube   ??? ',opts%llist%li_cube

    IF (opts%llist%li_chgcar) THEN
      CALL read_charge_chgcar(ions,chg,chargefile)
    ELSEIF(opts%llist%li_cube) THEN
      CALL read_charge_cube(ions,chg,chargefile)
    ENDIF

    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(1A12,1F6.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'
   
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
