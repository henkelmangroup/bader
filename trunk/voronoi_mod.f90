!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for analyzing the charge with a voronoi analysis
!
! By Andri Arnaldsson and Graeme Henkelman
! Last modified by 
!
!-----------------------------------------------------------------------------------!
MODULE voronoi_mod
  USE kind_mod , ONLY : q2
  USE matrix_mod
  USE ions_mod
  USE charge_mod
  IMPLICIT NONE

! Public, allocatable variables
  TYPE voronoi_obj
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: vorchg
  END TYPE

  PRIVATE
  PUBLIC :: voronoi_obj
  PUBLIC :: voronoi
  CONTAINS

!-----------------------------------------------------------------------------------!
!  voronoi:  Calculate the charge density populations based upon Voronoi polyhedra.
!    In this scheme each element of charge density is associated with the atom that
!    it is closest to.
!-----------------------------------------------------------------------------------!
  SUBROUTINE voronoi(vor,ions,chg)

    TYPE(voronoi_obj) :: vor
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

!    REAL(q2),DIMENSION(ndim,3) :: r_atm
    REAL(q2),DIMENSION(3) :: r_cur,dr,ngf,ngf_2
    REAL(q2) :: dist,min_dist,shift
    INTEGER :: i,nx,ny,nz,closest,tenths_done,cr,count_max,t1,t2
    INTEGER :: nxf,nyf,nzf

    CALL system_clock(t1,cr,count_max)

    WRITE(*,'(/,2x,A)') 'CALCULATING VORONOI CHARGE DISTRIBUTION'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'

    ALLOCATE(vor%vorchg(ions%nions))

    nxf=chg%npts(1)
    nyf=chg%npts(2)
    nzf=chg%npts(3)

    ngf(1)=REAL(nxf,q2)
    ngf(2)=REAL(nyf,q2)
    ngf(3)=REAL(nzf,q2)
    ngf_2=REAL(ngf,q2)/2.0_q2

    shift=1.0_q2
    IF (chg%halfstep) shift=0.5_q2

!    r_atm(:,1)=r_dir(:,1)*ngf(1)+shift
!    r_atm(:,2)=r_dir(:,2)*ngf(2)+shift
!    r_atm(:,3)=r_dir(:,3)*ngf(3)+shift

    vor%vorchg=0.0_q2
    tenths_done=0
    DO nx=1,nxf
      r_cur(1)=REAL(nx,q2)
      IF ((nx*10/nxf) > tenths_done) THEN
        tenths_done=(nx*10/nxf)
        WRITE(*,'(A,$)') '**'
      END IF
      DO ny=1,nyf
        r_cur(2)=REAL(ny,q2)
        DO nz=1,nzf
          r_cur(3)=REAL(nz,q2)
          closest=1
          dr=r_cur-ions%r_car(1,:)
          CALL dpbc(dr,ngf,ngf_2)
          min_dist=DOT_PRODUCT(dr,dr)
          DO i=2,ions%nions
            dr=r_cur-ions%r_car(i,:)
            CALL dpbc(dr,ngf,ngf_2)
            dist=DOT_PRODUCT(dr,dr)
            IF (dist < min_dist) THEN
              min_dist=dist
              closest=i
            END IF
          END DO
          vor%vorchg(closest)=vor%vorchg(closest)+rho_val(chg,nx,ny,nz)
        END DO
      END DO
    END DO
    WRITE(*,*)
! Don't have this normalization for MONDO
    vor%vorchg(:)=vor%vorchg(:)/REAL(chg%nrho,q2)

    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(1A12,1F6.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE voronoi

!-----------------------------------------------------------------------------------!

END MODULE voronoi_mod
