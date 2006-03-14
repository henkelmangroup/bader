!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
! Version 0.09c mindist
!
! By Andri Arnaldsson and Graeme Henkelman
! Last modified by GH on June 14 2004
!-----------------------------------------------------------------------------------!
  PROGRAM Charge
     USE options_mod , ONLY : get_options,options_obj
     USE ions_mod
     USE charge_mod
     USE io_mod
     USE bader_mod 
     USE voronoi_mod
!     USE multipole_mod
     USE results_mod , ONLY : output     

     IMPLICIT NONE

! Variables
     TYPE(options_obj) :: opts
     TYPE(ions_obj) :: ions
     TYPE(charge_obj) :: chg
     TYPE(bader_obj) :: bdr
     TYPE(voronoi_obj) :: vor
     INTEGER :: i
     LOGICAL :: ertil

! Get the control variables
     CALL get_options(opts)

! Make sure that that charge density file exists here
     i=LEN_TRIM(ADJUSTL(opts%chargefile))
     INQUIRE(FILE=opts%chargefile(1:i),EXIST=ertil)
     IF (.NOT. ertil) THEN
       WRITE(*,'(2X,A,A)') opts%chargefile(1:i),' DOES NOT EXIST IN THIS DIRECTORY'
       STOP
     ENDIF

     write(*,*) opts

! Call the read routines  .... from io.f90
     CALL read_charge(ions,chg,opts,opts%chargefile)

     write(*,*) ' in main again '
     write(*,*) ' if bader   ',opts%llist%lc_bader
     write(*,*) ' if voronoi ',opts%llist%lc_voronoi

! Choose bader if no other calculation method specified
     IF (.NOT.(opts%llist%lc_bader .OR. opts%llist%lc_voronoi)) THEN
       opts%llist%lc_bader=.TRUE.
     ENDIF
! Calculate
     IF (opts%llist%lc_bader) CALL bader(bdr,ions,chg)
!     IF (opts%llist%lc_dipole) CALL multipole()
! call mindist()
     IF (opts%llist%lc_voronoi) CALL voronoi(vor,ions,chg)

! Write out the volumes !!!!

     CALL output(bdr,vor,ions,chg)

!    USE varsM , ONLY : q2,addup,chargefile,vasp,na
!    USE ChargeM, ONLY : bader,multipole,voronoi,mindist
!    USE ChargeIOM 
!    IMPLICIT NONE
!
!    INTEGER :: i,n,iargc,choose
!    LOGICAL :: ertil
!    CHARACTER(LEN=20) :: cc
!
!    n=IARGC()+1
!    IF (n > 1) THEN
!      CALL GETARG(1,chargefile)
!    ELSE
!      WRITE(*,'(2x,A,$)') 'ENTER CHARGEFILE NAME: '
!      READ(*,*) chargefile
!    END IF
!    i=LEN_TRIM(chargefile)
!    INQUIRE(FILE=chargefile,EXIST=ertil)
!    IF (.NOT. ertil) THEN 
!      WRITE(*,'(2X,A,A)') chargefile(1:i),' DOES NOT EXIST IN THIS DIRECTORY' 
!      STOP
!    END IF
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

  STOP
  END PROGRAM Charge

