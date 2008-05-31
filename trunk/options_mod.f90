!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for specifying input options
!-----------------------------------------------------------------------------------!
  MODULE options_mod
    USE kind_mod , ONLY : q2
    IMPLICIT NONE

    TYPE :: options_obj
      CHARACTER(LEN=128) :: chargefile,refchgfile
      REAL(q2) :: badertol, stepsize
      INTEGER :: print_opt, print_none = 0, print_all_bader = 1, print_all_atom = 2, &
                 &  print_sel_bader = 3, print_sel_atom = 4
      INTEGER :: out_opt, out_auto = 0, out_cube = 1, out_chgcar = 2
      INTEGER :: in_opt, in_auto = 0, in_cube = 1, in_chgcar = 2
      INTEGER :: bader_opt, bader_offgrid = 0, bader_ongrid = 1, bader_neargrid = 2
      INTEGER :: quit_opt, quit_max = 0, quit_known = 1
      INTEGER :: refine_edge_itrs
      INTEGER :: selnum
      INTEGER,ALLOCATABLE,DIMENSION(:) :: selvol
      LOGICAL :: bader_flag, voronoi_flag, dipole_flag, ldos_flag
      LOGICAL :: verbose_flag,ref_flag
    END TYPE options_obj

    PRIVATE
    PUBLIC :: get_options,options_obj

    CONTAINS

!-----------------------------------------------------------------------------------!
! get_options: Read any input flags and the charge density file name
!-----------------------------------------------------------------------------------!

    SUBROUTINE get_options(opts)

      TYPE(options_obj) :: opts
      LOGICAL :: existflag
      LOGICAL :: readchgflag
      INTEGER :: n,iargc,i,ip,m,it,ini
      INTEGER :: j, sel
      REAL(q2) :: temp
      CHARACTER(LEN=128) :: p
      CHARACTER*128 :: inc
      INTEGER :: COMMAND_ARGUMENT_COUNT

! Default values
      opts%out_opt = opts%out_chgcar
      opts%in_opt = opts%in_auto
      opts%print_opt = opts%print_none
      opts%bader_opt = opts%bader_neargrid
      opts%quit_opt = opts%quit_known
      opts%refine_edge_itrs = -1
      opts%bader_flag = .TRUE.
      opts%voronoi_flag = .FALSE.
      opts%dipole_flag = .FALSE.
      opts%ldos_flag = .FALSE.
      opts%verbose_flag = .FALSE.
      opts%badertol = 1.0e-4_q2
      opts%stepsize = 0.0_q2
      opts%ref_flag=.FALSE.

!      n=IARGC()
      n=COMMAND_ARGUMENT_COUNT()
      IF (n == 0) THEN
        call write_options()
        STOP
      END IF

      ! Loop over all arguments
      m=0
      readchgflag = .FALSE.
      readopts: DO WHILE(m<n)
        m=m+1
!        CALL GETARG(m,p)
        CALL GET_COMMAND_ARGUMENT(m,p)
        p=ADJUSTL(p)
        ip=LEN_TRIM(p)
        i=INDEX(p,'-')
        IF (i /= 1) THEN
          ! Not a flag, so read the charge density file name
          IF (readchgflag) THEN
            WRITE(*,'(A,A,A)') ' Option "',p(1:ip),'" is not valid'
            CALL write_options()
            STOP
          END IF
          opts%chargefile=p
          INQUIRE(FILE=opts%chargefile,EXIST=existflag)
          IF (.NOT. existflag) THEN
            WRITE(*,'(2X,A,A)') opts%chargefile(1:ip),' does not exist'
            STOP
          END IF
          readchgflag=.TRUE.
        ! Help
        ELSEIF (p(1:ip) == '-h') THEN
          CALL write_help()
          STOP
        ! Verbose
        ELSEIF (p(1:ip) == '-v') THEN
          opts%verbose_flag = .TRUE.
        ! Bader options
        ELSEIF (p(1:ip) == '-b') THEN
          m=m+1
!          CALL GETARG(m,inc)
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc=ADJUSTL(inc)
          it=LEN_TRIM(inc)
          IF (inc(1:it) == 'OFFGRID' .OR. inc(1:it) == 'offgrid') THEN
            opts%bader_opt = opts%bader_offgrid
          ELSEIF (inc(1:it) == 'ONGRID' .OR. inc(1:it) == 'ongrid') THEN
            opts%bader_opt = opts%bader_ongrid
          ELSEIF (inc(1:it) == 'NEARGRID' .OR. inc(1:it) == 'neargrid') THEN
            opts%bader_opt = opts%bader_neargrid
          ELSE
            WRITE(*,'(A,A,A)') ' Unknown option "',inc(1:it),'"'
            STOP
          END IF
        ! Quit options
        ELSEIF (p(1:ip) == '-m') THEN
          m=m+1
!          CALL GETARG(m,inc)
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc=ADJUSTL(inc)
          it=LEN_TRIM(inc)
          IF (inc(1:it) == 'MAX' .OR. inc(1:it) == 'max') THEN
            opts%quit_opt = opts%quit_max
          ELSEIF (inc(1:it) == 'KNOWN' .OR. inc(1:it) == 'known') THEN
            opts%quit_opt = opts%quit_known
          ELSE
            WRITE(*,'(A,A,A)') ' Unknown option "',inc(1:it),'"'
            STOP
          END IF
        ! Print options
        ELSEIF (p(1:ip) == '-p') THEN
          m=m+1
!          CALL GETARG(m,inc)
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc=ADJUSTL(inc)
          it=LEN_TRIM(inc)
          IF (inc(1:it) == 'NONE' .OR. inc(1:it) == 'none') THEN
            opts%print_opt = opts%print_none
          ELSEIF (inc(1:it) == 'ALL_BADER' .OR. inc(1:it) == 'all_bader') THEN
            opts%print_opt = opts%print_all_bader
          ELSEIF (inc(1:it) == 'ALL_ATOM' .OR. inc(1:it) == 'all_atom') THEN
            opts%print_opt = opts%print_all_atom
          ELSEIF (inc(1:it) == 'SEL_BADER' .OR. inc(1:it) == 'sel_bader') THEN
            opts%print_opt = opts%print_sel_bader
          ELSEIF (inc(1:it) == 'SEL_ATOM' .OR. inc(1:it) == 'sel_atom') THEN
            opts%print_opt = opts%print_sel_atom
          ELSE
            WRITE(*,'(A,A,A)') ' Unknown option "',inc(1:it),'"'
            STOP
          END IF
          IF (opts%print_opt == opts%print_sel_bader .OR. opts.print_opt == opts%print_sel_atom) THEN
            ALLOCATE(opts%selvol(n))
            opts%selnum=0
            DO
              m=m+1
              CALL GET_COMMAND_ARGUMENT(m,inc)
              inc=ADJUSTL(inc)
              it=LEN_TRIM(inc)
              READ (inc(1:it),'(I10)',ERR=110) sel
              opts%selnum=opts%selnum+1
              opts%selvol(opts%selnum)=sel
              IF(m==n) EXIT readopts
              CYCLE
   110        m=m-1
              EXIT
            END DO
          END IF
        ! Output file type
        ELSEIF (p(1:ip) == '-o') THEN
          m=m+1
!          CALL GETARG(m,inc)
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc=ADJUSTL(inc)
          it=LEN_TRIM(inc)
          IF (inc(1:it) == 'CUBE' .OR. inc(1:it) == 'cube') THEN
            opts%out_opt = opts%out_cube
          ELSEIF (inc(1:it) == 'CHGCAR' .OR. inc(1:it) == 'chgcar') THEN
            opts%out_opt = opts%out_chgcar
          ELSE
            WRITE(*,'(A,A,A)') ' Unknown option "',inc(1:it),'"'
            STOP
          END IF  
        ! Calculation options
        ELSEIF (p(1:ip) == '-c') THEN
          m=m+1
!          CALL GETARG(m,inc)
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc=ADJUSTL(inc)
          it=LEN_TRIM(inc)
          IF (inc(1:it) == 'ALL' .OR. inc(1:it) == 'all') THEN
            opts%bader_flag = .TRUE.
            opts%voronoi_flag = .TRUE.
            opts%dipole_flag = .TRUE.
            opts%ldos_flag = .TRUE.
          ELSEIF (inc(1:it) == 'BADER' .OR. inc(1:it) == 'bader') THEN
            opts%bader_flag = .TRUE.
          ELSEIF (inc(1:it) == 'VORONOI' .OR. inc(1:it) == 'voronoi') THEN
            opts%voronoi_flag = .TRUE.
          ELSEIF (inc(1:it) == 'DIPOLE' .OR. inc(1:it) == 'dipole') THEN
            opts%dipole_flag = .TRUE.
          ELSEIF (inc(1:it) == 'LDOS' .OR. inc(1:it) == 'ldos') THEN
            opts%ldos_flag = .TRUE.
          ELSE
            WRITE(*,'(A,A,A)') ' Unknown option "',inc(1:it),'"'
            STOP
          END IF
        ELSEIF (p(1:ip) == '-n') THEN
          m=m+1
!          CALL GETARG(m,inc)
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc=ADJUSTL(inc)
          it=LEN_TRIM(inc)
          IF (inc(1:it) == 'ALL' .OR. inc(1:it) == 'all') THEN
            opts%bader_flag = .FALSE.
            opts%voronoi_flag = .FALSE.
            opts%dipole_flag = .FALSE.
            opts%ldos_flag = .FALSE.
          ELSEIF (inc(1:it) == 'BADER' .OR. inc(1:it) == 'bader') THEN
            opts%bader_flag = .FALSE.
          ELSEIF (inc(1:it) == 'VORONOI' .OR. inc(1:it) == 'voronoi') THEN
            opts%voronoi_flag = .FALSE.
          ELSEIF (inc(1:it) == 'DIPOLE' .OR. inc(1:it) == 'dipole') THEN
            opts%dipole_flag = .FALSE.
          ELSEIF (inc(1:it) == 'LDOS' .OR. inc(1:it) == 'ldos') THEN
            opts%ldos_flag = .FALSE.
          ELSE
            WRITE(*,'(A,A,A)') ' Unknown option "',inc(1:it),'"'
            STOP
          ENDIF
        ! Input file type
        ELSEIF (p(1:ip) == '-i') THEN
          m=m+1
!          CALL GETARG(m,inc)
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc=ADJUSTL(inc)
          it=LEN_TRIM(inc)
          IF (inc(1:it) == 'CUBE' .OR. inc(1:it) == 'cube') THEN
            opts%out_opt=opts%out_cube
          ELSEIF (inc(1:it) == 'CHGCAR' .OR. inc(1:it) == 'chgcar') THEN
            opts%out_opt=opts%out_chgcar
          ELSE
            WRITE(*,'(A,A,A)') ' Unknown option "',inc(1:it),'"'
            STOP
          END IF
        ! Bader tolerance
        ELSEIF (p(1:ip) == '-t') THEN
          m=m+1
!          CALL GETARG(m,inc)
          CALL GET_COMMAND_ARGUMENT(m,inc)
          read(inc,*) opts%badertol
        ! Refine edge iterations  -- change this to a flag once working
        ELSEIF (p(1:ip) == '-r') THEN
          m=m+1
!          CALL GETARG(m,inc)
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc=ADJUSTL(inc)
          it=LEN_TRIM(inc) 
          IF (inc(1:it) == 'AUTO' .OR. inc(1:it) == 'auto') THEN
            opts%refine_edge_itrs=-1
          ELSE
            read(inc,*) opts%refine_edge_itrs
          END IF
        ! Step size
        ELSEIF (p(1:ip) == '-s') THEN
          m=m+1
!          CALL GETARG(m,inc)
          CALL GET_COMMAND_ARGUMENT(m,inc)
          read(inc,*) opts%stepsize
        !doing analysis with ref charge
        ELSEIF (p(1:ip) == '-ref') THEN
          m=m+1
!          CALL GETARG(m,inc)
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc=ADJUSTL(inc)
          it=LEN_TRIM(inc)
          IF (inc(1:it) == 'NONE' .OR. inc(1:it) == 'none') THEN
            opts%ref_flag = .FALSE.
          ELSE
            opts%ref_flag = .TRUE.
            opts%refchgfile = inc(1:it)
          END IF
        ! Unknown flag
        ELSE
          WRITE(*,'(A,A,A)') ' Unknown option flag "',p(1:ip),'"'
          STOP
        END IF

      END DO readopts

    ! If no file name, we die
    IF (.NOT. readchgflag) THEN
      WRITE(*,*) ' ERROR: Did not read a charge file name in the arguments'
      CALL write_options()
      STOP
    ENDIF

    ! Default to no edge refinement for the ongrid algorithm
    IF (opts%bader_opt==opts%bader_ongrid) THEN
      opts%refine_edge_itrs=0
    END IF
   
    IF (opts%print_opt==opts%print_sel_atom .AND. opts%selnum==0) THEN
      WRITE(*,'(/,A)') 'NO ATOMIC VOLUMES SELECTED'
      STOP
    END IF

    IF (opts%print_opt==opts%print_sel_bader .AND. opts%selnum==0) THEN
      WRITE(*,'(/,A)') 'NO BADER VOLUMES SELECTED'
      STOP
    END IF

    RETURN
    END SUBROUTINE get_options

!-----------------------------------------------------------------------------------!
! write_opts: write flag options
!-----------------------------------------------------------------------------------!

    SUBROUTINE write_options()

      WRITE(*,*) ''
      WRITE(*,*) 'Usage:'
      WRITE(*,*) '   bader [ -c bader | voronoi ]'
      WRITE(*,*) '         [ -n bader | voronoi ]'
      WRITE(*,*) '         [ -b neargrid | ongrid ]'
      WRITE(*,*) '         [ -r refine_edge_iterations ]'
      WRITE(*,*) '         [ -ref reference_charge ]'
      WRITE(*,*) '         [ -m known | max ]'
      WRITE(*,*) '         [ -p none | all_atom | all_bader | sel_atom | sel_bader ] [ volume list]'
      WRITE(*,*) '         [ -i cube | chgcar ]'
      WRITE(*,*) '         [ -o cube | chgcar ]'
      WRITE(*,*) '         [ -h ] [ -v ]'
      WRITE(*,*) '         chargefile'
      WRITE(*,*) ''

    END SUBROUTINE write_options

!-----------------------------------------------------------------------------------!
! write_help: write help
!-----------------------------------------------------------------------------------!

    SUBROUTINE write_help()

      WRITE(*,*) ''
      WRITE(*,*) 'Description of flags'
      WRITE(*,*) ''
!      WRITE(*,*) '   -c | -n  < bader | voronoi | dipole | ldos >'
      WRITE(*,*) '   -c | -n  < bader | voronoi >'
      WRITE(*,*) '        Turn on [-c] or off [-n] the following calculations'
      WRITE(*,*) '           bader: Bader atoms in molecules (default)'
      WRITE(*,*) '           voronoi: population analysis based on distance'
!      WRITE(*,*) '           dipole: multiple moments in Bader volumes'
!      WRITE(*,*) '           ldos: local density of states in Bader volumes'
      WRITE(*,*) ''
      WRITE(*,*) '   -b < neargrid | ongrid >'
      WRITE(*,*) '        Use the default near-grid bader partitioning or the'
      WRITE(*,*) '        original on-grid based algorithm.'
      WRITE(*,*) ''
!      WRITE(*,*) '   -s < stepsize >'
!      WRITE(*,*) '        Steepest asent trajectory step size.  This parameter is'
!      WRITE(*,*) '        (only) used for the default offgrid Bader analysis.  If'
!      WRITE(*,*) '        not specified, the stepsize is set to the minimum distance'
!      WRITE(*,*) '        between charge density grid points.'
!      WRITE(*,*) ''
      WRITE(*,*) '   -r < iterations >'
      WRITE(*,*) '        Number of times the bader edges are refined.'
      WRITE(*,*) ''
      WRITE(*,*) '   -ref < reference_charge >'
      WRITE(*,*) '        Use a reference charge file to do the Bader partitioning.'
      WRITE(*,*) '        This is the recommended way to analyze vasp output files:'
      WRITE(*,*) '           bader CHGCAR -ref CHGCAR_total'
      WRITE(*,*) '        where CHGCAR_total is the sum of AECCAR0 and AECCAR2.'
      WRITE(*,*) ''
      WRITE(*,*) '   -m < known | max >'
      WRITE(*,*) '        Determines how trajectories terminate'
      WRITE(*,*) '           known: stop when a point is surrounded by known points' 
      WRITE(*,*) '           max: stop only when a charge density maximum is reached '
      WRITE(*,*) ''
      WRITE(*,*) '   -p < none | all_atom | all_bader | sel_atom | sel_bader > <volume list>'
      WRITE(*,*) '        Print options for calculated Bader volumes'
      WRITE(*,*) '           none: do not output any Bader volumes'
      WRITE(*,*) '           all_atom: all atomic volumes'
      WRITE(*,*) '           all_bader: all Bader volumes'
      WRITE(*,*) '           sel_atom: atomic volume(s) around the selected atom(s)'
      WRITE(*,*) '           sel_bader: selected Bader volumes'
      WRITE(*,*) ''
      WRITE(*,*) '   -i < cube | chgcar >'
      WRITE(*,*) '        Input charge density file type.  If not specified, the'
      WRITE(*,*) '        program will try to determine the charge density file type'
      WRITE(*,*) '        automatically.'
      WRITE(*,*) ''
      WRITE(*,*) '   -o < cube | chgcar >'
      WRITE(*,*) '        Output charge density file type.'
      WRITE(*,*) ''
      WRITE(*,*) '   -h'
      WRITE(*,*) '        Help.'
      WRITE(*,*) ''
      WRITE(*,*) '   -v'
      WRITE(*,*) '        Verbose output.'
      WRITE(*,*) ''
    
    END SUBROUTINE write_help

!-----------------------------------------------------------------------------------!

  END module options_mod
