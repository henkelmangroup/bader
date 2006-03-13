!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for analyzing the charge
!
! By Andri Arnaldsson and Graeme Henkelman
! Last modified by
!-----------------------------------------------------------------------------------!
MODULE charge_mod
  USE kind_mod , ONLY : q2
  USE matrix_mod
  USE ions_mod , ONLY : ions_obj
  IMPLICIT NONE

! Public, allocatable variables

  TYPE :: charge_obj
    REAL(q2),ALLOCATABLE,DIMENSION(:,:,:) :: rho
    REAL(q2),DIMENSION(3,3) :: lattice
    REAL(q2),DIMENSION(3) :: corner,steps
    INTEGER :: nxf,nyf,nzf,nrho
    LOGICAL :: halfstep
  END TYPE

  PRIVATE
  PUBLIC :: charge_obj
  PUBLIC :: rho_val,pbc,dpbc_dir,dpbc,pbcq2,verticies
  CONTAINS

!-----------------------------------------------------------------------------------!
!  rho_val:  Return the density at the point (px,py,pz) taking into account the
!    boundary conditions.  This function is used to address points outside the
!    charge density array without a bunch of if statements at the place the value
!    is needed.
!-----------------------------------------------------------------------------------!
  FUNCTION rho_val(chg,px,py,pz)
    TYPE(charge_obj) :: chg
    INTEGER,INTENT(IN) :: px,py,pz
    REAL(q2) :: rho_val

    INTEGER :: pxt,pyt,pzt

    pxt=px
    pyt=py
    pzt=pz
    DO
      IF(pxt >= 1) EXIT
      pxt=pxt+chg%nxf
    END DO
    DO
      IF(pxt <= chg%nxf) EXIT
      pxt=pxt-chg%nxf
    END DO
    DO
      IF(pyt >= 1) EXIT
      pyt=pyt+chg%nyf
    END DO
    DO
      IF(pyt <= chg%nyf) EXIT
      pyt=pyt-chg%nyf
    END DO
    DO
      IF(pzt >= 1) EXIT
      pzt=pzt+chg%nzf
    END DO
    DO
      IF(pzt <= chg%nzf) EXIT
      pzt=pzt-chg%nzf
    END DO

    rho_val=chg%rho(pxt,pyt,pzt)

  RETURN
  END FUNCTION rho_val

!-----------------------------------------------------------------------------------!
! pbc: Wrap the point (px,py,pz) to the boundary conditions [0,nf].
!-----------------------------------------------------------------------------------!
  SUBROUTINE pbc(px,py,pz,nxf,nyf,nzf)

    INTEGER,INTENT(INOUT) :: px,py,pz
    INTEGER,INTENT(IN) :: nxf,nyf,nzf

    DO
      IF(px > 0) EXIT
      px=px+nxf
    END DO
    DO
      IF(px <= nxf) EXIT
      px=px-nxf
    END DO
    DO
      IF(py > 0) EXIT
      py=py+nyf
    END DO
    DO
      IF(py <= nyf) EXIT
      py=py-nyf
    END DO
    DO
      IF(pz > 0) EXIT
      pz=pz+nzf
    END DO
    DO
      IF(pz <= nzf) EXIT
      pz=pz-nzf
    END DO

  RETURN
  END SUBROUTINE pbc

!-----------------------------------------------------------------------------------!
! dpbc_dir:  Wrap the vector dr to the boundary conditions [-1/2,1/2].
!-----------------------------------------------------------------------------------!
  SUBROUTINE dpbc_dir(dr)
    REAL(q2),INTENT(INOUT),DIMENSION(3) :: dr

    INTEGER :: i

    DO i=1,3
      DO
        IF(dr(i) > -0.5_q2) EXIT
        dr(i)=dr(i)+1.0_q2
      END DO
      DO
        IF(dr(i) < 0.5_q2) EXIT
        dr(i)=dr(i)-1.0_q2
      END DO
    END DO
  RETURN
  END SUBROUTINE dpbc_dir

!-----------------------------------------------------------------------------------!
! dpbc:  Wrap the vector dr to the boundary conditions [-ngf/2,ngf/2].
!-----------------------------------------------------------------------------------!
  SUBROUTINE dpbc(dr,nf,nf_2)
    REAL(q2),INTENT(IN),DIMENSION(3) :: nf,nf_2
    REAL(q2),INTENT(INOUT),DIMENSION(3) :: dr

    INTEGER :: i

    DO i=1,3
      DO
        IF(dr(i) > -nf_2(i)) EXIT
        dr(i)=dr(i)+nf(i)
      END DO
      DO
        IF(dr(i) < nf_2(i)) EXIT
        dr(i)=dr(i)-nf(i)
      END DO
    END DO

  RETURN
  END SUBROUTINE dpbc

!-----------------------------------------------------------------------------------!
! pbc: Wrap the REAL(q2) point (px,py,pz) to the boundary conditions [0,nf].
!-----------------------------------------------------------------------------------!
  SUBROUTINE pbcq2(px,py,pz,rnxf,rnyf,rnzf)
    REAL(q2),INTENT(INOUT) :: px,py,pz
    REAL(q2),INTENT(IN) :: rnxf,rnyf,rnzf

    DO
      IF(px >= 0.0_q2) EXIT
      px=px+rnxf
    END DO
    DO
      IF(px < rnxf) EXIT
      px=px-rnxf
    END DO
    DO
      IF(py >= 0.0_q2) EXIT
      py=py+rnyf
    END DO
    DO
      IF(py < rnyf) EXIT
      py=py-rnyf
    END DO
    DO
      IF(pz >= 0.0_q2) EXIT
      pz=pz+rnzf
    END DO
    DO
      IF(pz < rnzf) EXIT
      pz=pz-rnzf
    END DO

  RETURN
  END SUBROUTINE pbcq2

!-----------------------------------------------------------------------------------!
! verticies: Return the points (ix,iy,iz) and (ixm,iym,izm) which are the verticies
!   closest to the REAL(q2) point (px,py,pz)
!-----------------------------------------------------------------------------------!
  SUBROUTINE verticies(px,py,pz,nxf,nyf,nzf,ix,iy,iz,ixm,iym,izm)

    REAL(q2),INTENT(INOUT) :: px,py,pz
    INTEGER,INTENT(IN) :: nxf,nyf,nzf
    INTEGER,INTENT(OUT) :: ix,iy,iz,ixm,iym,izm

    ix=FLOOR(px)
    IF (ix == 0) THEN
      ix=nxf
      ixm=1
    ELSE
      ixm=ix+1
    END IF
    iy=FLOOR(py)
    IF (iy == 0) THEN
      iy=nyf
      iym=1
    ELSE
      iym=iy+1
    END IF
    iz=FLOOR(pz)
    IF (iz == 0) THEN
      iz=nzf
      izm=1
    ELSE
      izm=iz+1
    END IF

  RETURN
  END SUBROUTINE verticies

!-----------------------------------------------------------------------------------!

END MODULE charge_mod
