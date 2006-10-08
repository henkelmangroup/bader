!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
! Version 0.2, 09/28/26 near-grid algorithm 
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
!   (submitted)
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

     IMPLICIT NONE

     ! Variables
     TYPE(options_obj) :: opts
     TYPE(ions_obj) :: ions
     TYPE(charge_obj) :: chg
     TYPE(bader_obj) :: bdr
     TYPE(voronoi_obj) :: vor

     ! Write the version number
     WRITE(*,'(/,2X,A)') 'GRID BASED BADER ANALYSIS  (v0.20 10/08/06)'

     ! Get the control variables
     CALL get_options(opts)

     ! Call the read routines from io_mod
     CALL read_charge(ions,chg,opts)

     ! Calculations
     IF (opts%bader_flag) THEN
       CALL bader_calc(bdr,ions,chg,opts)
       CALL bader_mindist(bdr,ions,chg)
       CALL bader_output(bdr,ions,chg)
       IF (opts%print_opt==opts%print_all_bader) THEN
         CALL write_all_bader(bdr,opts,ions,chg)
       ELSEIF (opts%print_opt==opts%print_all_atom) THEN
         CALL write_all_atom(bdr,opts,ions,chg)
       ENDIF
     ENDIF
!     IF (opts%dipole_flag) CALL multipole()
     IF (opts%voronoi_flag) CALL voronoi(vor,ions,chg)

      ! Write density files
!      IF (opts.print_opt==opts.print_all_bader) THEN
!        CALL write_all_bader(bdr,ions,chg)
!      ELSEIF (opts.print_opt==opts.print_all_bader) THEN
!        CALL write_all_atom(bdr,ions,chg)
!      ENDIF

!GH: need to add write selected bader and atomic regions


!    INTEGER :: i,choose
!    CHARACTER(LEN=20) :: cc
!
!    IF (n == 1 .OR. n == 2) THEN
!      WRITE(*,'(/,2X,A)') 'CHOOSE OUTPUT:'
!      WRITE(*,'(4X,A)') '1 ... DO NOT WRITE OUT ANY BADER VOLUMES'
!      WRITE(*,'(4X,A)') '2 ... WRITE ALL BADER VOLUMES TO FILES'
!      WRITE(*,'(4X,A)') '3 ... WRITE BADER VOLUME(S) FOR EACH ATOM TO FILE'
!      WRITE(*,'(4X,A)') '4 ... WRITE SPECIFIED BADER VOLUMES TO FILE'
!      WRITE(*,'(2X,A,$)') 'CHOOSE 1,2,3 or 4 : '
!      READ(*,*) choose
!    ELSE IF (n > 2) THEN
!      CALL GETARG(2,cc)
!      READ(cc,*) choose
!    END IF
!    IF (choose == 3) THEN
!      WRITE(*,'(2X,A,$)') 'NUMBER OF VOLUMES TO ADD UP: '
!      READ(*,*) na
!      ALLOCATE(addup(na))
!      WRITE(*,'(2X,A,$)') 'INPUT THE VOLUMES TO ADD UP: '   
!      READ(*,*) addup
!    END IF
!
!! ECHO INPUT
!    WRITE(*,'(/,2X,A)') 'ECHO INPUT'
!    WRITE(*,'(2X,A,A)') 'CHARGE DENSITY FILE: ',chargefile
!    IF (choose == 1) WRITE(*,'(2X,A)') 'OUTPUT: WRITING ALL BADER VOLUMES TO FILES'
!    IF (choose == 2) WRITE(*,'(2X,A)') 'OUTPUT: WRITING BADER VOLUME FOR EACH ATOM TO FILE'
!    IF (choose == 3) THEN 
!      WRITE(*,'(2X,A)') 'OUTPUT: WRITING SPECIFIED BADER VOLUMES TO FILE'
!      WRITE(*,*) ' WITH VOLUMES: ',addup
!    END IF
!    IF (choose == 4) WRITE(*,'(2X,A)') 'OUTPUT: NOT WRITING OUT ANY BADER VOLUMES'
!
!    CALL read_charge()
!    CALL bader()
!    CALL multipole()
!    CALL mindist()
!    CALL voronoi()
!
!    SELECT CASE (choose)
!      CASE (1)
!        IF(vasp) THEN
!          CALL Vallvolume()
!        ELSE
!          CALL Gallvolume()
!        END IF
!        CALL write_bader_num()
!      CASE (2)
!        IF(vasp) THEN
!          CALL Vatomvolume()
!        ELSE
!          CALL Gatomvolume()
!        END IF
!        CALL write_bader_num()
!      CASE (3)
!        IF(vasp) THEN
!          CALL Vspecvolume()
!        ELSE
!          CALL Gspecvolume()
!        END IF
!      CASE (4)
!!       DO NOTHING
!      CASE DEFAULT
!        WRITE(*,'(2X,A)') 'PLEASE CHOOSE 1,2,3 or 4'
!        STOP
!    END SELECT
!
!    CALL output()

    WRITE(*,*)
    STOP
  END PROGRAM Charge

