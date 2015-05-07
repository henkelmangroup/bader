MODULE brute_mod

  USE kind_mod
  USE matrix_mod
  USE bader_mod
  USE charge_mod 
  USE options_mod
  USE ions_mod
  USE io_mod
  USE chgcar_mod
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: brute_force

  CONTAINS

  SUBROUTINE brute_force(bdr,ions,chgval,opts)
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chgval,chgtemp
    TYPE(ions_obj) :: ions,ionstemp
    TYPE(options_obj) :: opts
    ! first read in the entire density grid
    IF (opts%ref_flag) THEN
      CALL read_charge_ref(ionstemp,chgtemp,opts)
    ELSE
      chgtemp = chgval
    END IF
    ALLOCATE(bdr%volnum(chgtemp%npts(1),chgtemp%npts(2),chgtemp%npts(3)))
    ALLOCATE(bdr%nnion(bdr%nvols), bdr%iondist(bdr%nvols),bdr%ionchg(ions%nions))
  END SUBROUTINE
END MODULE brute_mod
