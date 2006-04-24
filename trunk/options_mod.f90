  MODULE options_mod
    USE kind_mod , ONLY : q2
    IMPLICIT NONE

    TYPE :: opts_obj
      CHARACTER(LEN=120) :: chargefile
      REAL(q2) :: badertol
      INTEGER :: print_opt, print_none = 0, print_all = 1, print_atom = 2
      INTEGER :: out_opt, out_cube = 1, out_chgcar = 2
      INTEGER :: in_opt, in_cube = 1, in_chgcar = 2
      LOGICAL :: bader_flag, voronoi_flag, dipole_flag, ldos_flag
      LOGICAL :: verbose_flag
    END TYPE opts_obj

    PRIVATE
    PUBLIC :: get_opts,opts_obj

    CONTAINS
!-----------------------------------------------------------------------------------!
! get_opts: Read any input flags and the charge density file name
!-----------------------------------------------------------------------------------!
    SUBROUTINE get_opts(opts)

      TYPE(opts_obj) :: opts
      LOGICAL :: existflag, inl
      LOGICAL :: readchgflag
      INTEGER :: n,iargc,i,ip,m,it,ini
      REAL(q2) :: inr,temp
      CHARACTER(LEN=120) :: p
      CHARACTER(LEN=30) :: inc

! Default values
      opts%out_opt = 0
      opts%in_opt = 0
      opts%print_opt = 0
      opts%bader_flag = .TRUE.
      opts%voronoi_flag = .FALSE.
      opts%dipole_flag = .FALSE.
      opts%ldos_flag = .FALSE.
      opts%verbose_flag = .FALSE.
      opts%badertol=1.0e-4_q2

      n=IARGC()
      IF (n == 0) THEN
        call write_opts()
        STOP
      END IF

! Loop over all arguments
      m=0
      readchgflag=.FALSE.
      DO WHILE(m<n)
        m=m+1
        CALL GETARG(m,p)
        p=ADJUSTL(p)
        ip=LEN_TRIM(p)
        i=INDEX(p,'-')
        IF (i /= 1) THEN
! Not a flag, so read the charge density file name
          IF (readchgflag) THEN
            WRITE(*,'(A,A,A)') ' Option "',p(1:ip),'" is not valid'
            CALL write_opts()
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
! Print options
        ELSEIF (p(1:ip) == '-p') THEN
          m=m+1
          CALL GETARG(m,inc)
          inc=ADJUSTL(inc)
          it=LEN_TRIM(inc)
          IF (inc(1:it) == 'ALL' .OR. inc(1:it) == 'all') THEN
            opts%print_opt = opts%print_all
          ELSEIF (inc(1:it) == 'ATOM' .OR. inc(1:it) == 'atom') THEN
            opts%print_opt = opts%print_atom
          ELSEIF (inc(1:it) == 'NONE' .OR. inc(1:it) == 'none') THEN
            opts%print_opt = opts%print_none
          ELSE
            WRITE(*,'(A,A,A)') ' Unknown option "',inc(1:it),'"'
            STOP
          END IF
! Output file type
        ELSEIF (p(1:ip) == '-o') THEN
          m=m+1
          CALL GETARG(m,inc)
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
          CALL GETARG(m,inc)
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
          CALL GETARG(m,inc)
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
! Output file types
        ELSEIF (p(1:ip) == '-i') THEN
          m=m+1
          CALL GETARG(m,inc)
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
          CALL GETARG(m,inr)
          opts%badertol=inr
        ELSE
          WRITE(*,'(A,A,A)') ' Unknown option flag "',p(1:ip),'"'
          STOP
        END IF 

      END DO

! If no file name, we die
    IF (.NOT. readchgflag) THEN
      WRITE(*,*) 'Did not read a charge file name in the arguments'
      CALL write_opts()
      STOP
    ENDIF
 
    RETURN
    END SUBROUTINE get_opts

!-----------------------------------------------------------------------------------!
! write_opts: write flag options
!-----------------------------------------------------------------------------------!
    SUBROUTINE write_opts()

      WRITE(*,*) ''
      WRITE(*,*) 'Usage:'
      WRITE(*,*) '   bader [ -c bader|voronoi|dipole|ldos ]'
      WRITE(*,*) '         [ -n bader|voronoi|dipole|ldos ]'
      WRITE(*,*) '         [ -i cube|chgcar ] [ -o cube|chgcar ]'
      WRITE(*,*) '         [ -p none|atom|all ] [ -h ] [ -v ]'
      WRITE(*,*) '         chargefile'
      WRITE(*,*) ''

    END SUBROUTINE write_opts


!-----------------------------------------------------------------------------------!
! write_help: write help
!-----------------------------------------------------------------------------------!
    SUBROUTINE write_help()

      WRITE(*,*) ''
      WRITE(*,*) 'Description of flags'
      WRITE(*,*) ''
    
    END SUBROUTINE write_help
!-----------------------------------------------------------------------------------!

  END module options_mod
