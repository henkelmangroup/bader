!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for analyzing the charge with a voronoi analysis
!
! By Andri Arnaldsson and Graeme Henkelman
! Last modified by 
!
!-----------------------------------------------------------------------------------!
MODULE voronoi_mod
  USE vars_mod , ONLY : q2
  USE matrix_mod
  USE ions_mod
  USE charge_mod
  IMPLICIT NONE

! Public, allocatable variables
  REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: voronoi_charge

!  INTEGER :: ndim,bdim,nrho,wdim
  INTEGER :: bdim,nrho,wdim

  PRIVATE
  PUBLIC :: voronoi
  CONTAINS

!-----------------------------------------------------------------------------------!
!  vornoi:  Calculate the charge density populations based upon Vornoi polyhedia.
!    In this scheme each element of charge density is associated with the atom that
!    it is closest to.
!-----------------------------------------------------------------------------------!
  SUBROUTINE voronoi(ions,chg)

    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    REAL(q2),DIMENSION(ndim,3) :: r_atm
    REAL(q2),DIMENSION(3) :: r_cur,dr,ngf,ngf_2
    REAL(q2) :: dist,min_dist,shift
    INTEGER :: i,nx,ny,nz,closest,tenths_done,cr,count_max,t1,t2

    CALL system_clock(t1,cr,count_max)

    WRITE(*,'(/,2x,A)') 'CALCULATING VORONOI CHARGE DISTRIBUTION'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'

    wdim=ndim
    ALLOCATE(voronoi_charge(wdim,4))

    ngf(1)=REAL(nxf,q2)
    ngf(2)=REAL(nyf,q2)
    ngf(3)=REAL(nzf,q2)
    ngf_2=REAL(ngf,q2)/2.0_q2

    shift=0.5_q2
    IF (vasp) shift=1.0_q2

    r_atm(:,1)=r_dir(:,1)*ngf(1)+shift
    r_atm(:,2)=r_dir(:,2)*ngf(2)+shift
    r_atm(:,3)=r_dir(:,3)*ngf(3)+shift

    voronoi_charge=0.0_q2
    tenths_done=0
    DO nx=1,nxf
      r_cur(1)=REAL(nx,q2)
      IF ((nx*10/nxf) > tenths_done) THEN
        tenths_done=(nx*10/nxf)
!        WRITE(*,'(1X,1I4,1A6)') (tenths_done*10),'% done'
        WRITE(*,'(A,$)') '**'
      END IF
      DO ny=1,nyf
        r_cur(2)=REAL(ny,q2)
        DO nz=1,nzf
          r_cur(3)=REAL(nz,q2)
          closest=1
          dr=r_cur-r_atm(1,:)
          CALL dpbc(dr,ngf,ngf_2)
          min_dist=DOT_PRODUCT(dr,dr)
          DO i=2,wdim
            dr=r_cur-r_atm(i,:)
            CALL dpbc(dr,ngf,ngf_2)
            dist=DOT_PRODUCT(dr,dr)
            IF (dist < min_dist) THEN
              min_dist=dist
              closest=i
            END IF
          END DO
          voronoi_charge(closest,4)=voronoi_charge(closest,4)+                    &
                               rho_val(chg,nx,ny,nz,nxf,nyf,nzf)
        END DO
      END DO
    END DO
    WRITE(*,*)
! Don't have this normalization for MONDO
    voronoi_charge(:,4)=voronoi_charge(:,4)/REAL(nrho,q2)
    voronoi_charge(:,1:3)=r_car

    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(1A12,1F6.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE voronoi

!-----------------------------------------------------------------------------------!

END MODULE voronoi_mod
