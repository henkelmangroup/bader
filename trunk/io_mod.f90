!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!    Module for reading and writing charge density data
!
! By Andri Arnaldsson and Graeme Henkelman
! Last modified by GH on Apr. 23 2006
!-----------------------------------------------------------------------------------!

MODULE io_mod
  USE kind_mod , ONLY : q1,q2
  USE matrix_mod
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

  SUBROUTINE read_charge(ions,chg,opts)

    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg
    TYPE(options_obj) :: opts

    CHARACTER(LEN=120) :: chargefile
    CHARACTER(LEN=7) :: text
    INTEGER :: cr,count_max,t1,t2

    CALL system_clock(t1,cr,count_max)

    chargefile=opts%chargefile
    IF ( opts%in_opt == opts%in_auto ) THEN
      ! Try to guess the file type
      OPEN(100,FILE=chargefile,STATUS='old',ACTION='read',BLANK='null',PAD='yes')
      READ(100,'(6/,1A7)') text
      CLOSE(100)
      IF (text == 'Direct') THEN
        opts%in_opt=opts%in_chgcar
      ELSE 
        opts%in_opt=opts%in_cube
      ENDIF
    ENDIF

    IF (opts%in_opt == opts%in_chgcar) THEN
      CALL read_charge_chgcar(ions,chg,chargefile)
    ELSEIF (opts%in_opt == opts%in_cube) THEN
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

    IF (opts%out_opt == opts%out_chgcar) THEN
      CALL write_charge_chgcar(ions,chg,chargefile)
    ENDIF
    IF (opts%out_opt == opts%out_cube) THEN 
      CALL write_charge_cube(ions,chg,chargefile)
    ENDIF
    
  RETURN
  END SUBROUTINE write_charge

!------------------------------------------------------------------------------------!

END MODULE io_mod
