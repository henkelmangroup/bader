!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
! Version 0.22a, 03/07 neargrid and reference charge
!
! Authors:
!   Andri Arnaldsson, Wenjie Tang, and Graeme Henkelman
!
! Based on algorithms described in the following publications:
!
!   A fast and robust algorithm for Bader decomposition of charge density,
!   G. Henkelman, A. Arnaldsson, and H. Jonsson,
!   Comput. Mater. Sci. 36 254-360 (2006).
!
!   An improved grid-based algorithm for Bader charge allocation
!   E. Sanville, S. Kenny, R. Smith, and G. Henkelman
!   J. Comput. Chem. 28 899-908 (2007).
!
!   (no title yet, but this will describe the default near-grid algorithm used here)
!   W. Tang, E. Sanville, and G. Henkelman
!
!-----------------------------------------------------------------------------------!

  PROGRAM Charge

     USE options_mod
     USE ions_mod
     USE charge_mod
     USE io_mod
     USE bader_mod 
     USE voronoi_mod
     USE chgcar_mod
     
     IMPLICIT NONE

     ! Variables
     TYPE(options_obj) :: opts
     TYPE(ions_obj) :: ions
     TYPE(charge_obj) :: chgval
     TYPE(bader_obj) :: bdr
     TYPE(voronoi_obj) :: vor
 
    ! Write the version number
     WRITE(*,'(/,2X,A)') 'GRID BASED BADER ANALYSIS  (v0.22a 03/07)'

     ! Get the control variables
     CALL get_options(opts)

     ! Call the read routines from io_mod
     CALL read_charge(ions,chgval,opts)

     IF (opts%bader_flag) THEN
       CALL bader_calc(bdr,ions,chgval,opts)
       CALL bader_mindist(bdr,ions,chgval)
       CALL bader_output(bdr,ions,chgval)
       IF (opts%print_opt==opts%print_all_bader) THEN
         CALL write_all_bader(bdr,opts,ions,chgval)
       ELSEIF (opts%print_opt==opts%print_all_atom) THEN
         CALL write_all_atom(bdr,opts,ions,chgval)
       ENDIF
     ENDIF
!     IF (opts%dipole_flag) CALL multipole()
     IF (opts%voronoi_flag) CALL voronoi(vor,ions,chgval)

    WRITE(*,*)
  END PROGRAM Charge

