!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for analyzing the charge
!
! By Andri Arnaldsson and Graeme Henkelman
! Last modified by GH on Apr 24 2005
!
!-----------------------------------------------------------------------------------!
MODULE charge
  USE vars , ONLY : q2
  USE matrix
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: rho_value,pbc,dpbc_dir,dpbc,pbcq2,verticies
  CONTAINS

!-----------------------------------------------------------------------------------!
!  rho_value:  Return the density at the point (px,py,pz) taking into account the
!    boundary conditions.  This function is used to address points outside the
!    charge density array without a bunch of if statements at the place the value
!    is needed.
!-----------------------------------------------------------------------------------!
  FUNCTION rho_value(px,py,pz,ngxf,ngyf,ngzf)
    INTEGER,INTENT(IN) :: px,py,pz
    INTEGER,INTENT(IN) :: ngxf,ngyf,ngzf
    REAL(q2) :: rho_value

    INTEGER :: pxt,pyt,pzt

    pxt=px
    pyt=py
    pzt=pz
    DO
      IF(pxt >= 1) EXIT
      pxt=pxt+ngxf
    END DO
    DO
      IF(pxt <= ngxf) EXIT
      pxt=pxt-ngxf
    END DO
    DO
      IF(pyt >= 1) EXIT
      pyt=pyt+ngyf
    END DO
    DO
      IF(pyt <= ngyf) EXIT
      pyt=pyt-ngyf
    END DO
    DO
      IF(pzt >= 1) EXIT
      pzt=pzt+ngzf
    END DO
    DO
      IF(pzt <= ngzf) EXIT
      pzt=pzt-ngzf
    END DO

    rho_value=rho(pxt,pyt,pzt)

  RETURN
  END FUNCTION rho_value

!-----------------------------------------------------------------------------------!
! pbc: Wrap the point (px,py,pz) to the boundary conditions [0,ngf].
!-----------------------------------------------------------------------------------!
  SUBROUTINE pbc(px,py,pz)

    INTEGER,INTENT(INOUT) :: px,py,pz

    DO
      IF(px > 0) EXIT
      px=px+ngxf
    END DO
    DO
      IF(px <= ngxf) EXIT
      px=px-ngxf
    END DO
    DO
      IF(py > 0) EXIT
      py=py+ngyf
    END DO
    DO
      IF(py <= ngyf) EXIT
      py=py-ngyf
    END DO
    DO
      IF(pz > 0) EXIT
      pz=pz+ngzf
    END DO
    DO
      IF(pz <= ngzf) EXIT
      pz=pz-ngzf
    END DO

  RETURN
  END SUBROUTINE pbc

!-----------------------------------------------------------------------------------!
! dpbc_dir:  Wrap the vector dR to the boundary conditions [-1/2,1/2].
!-----------------------------------------------------------------------------------!
  SUBROUTINE dpbc_dir(dR)
    REAL(q2),INTENT(INOUT),DIMENSION(3) :: dR

    INTEGER :: i

    DO i=1,3
      DO
        IF(dR(i) > -0.5_q2) EXIT
        dR(i)=dR(i)+1.0_q2
      END DO
      DO
        IF(dR(i) < 0.5_q2) EXIT
        dR(i)=dR(i)-1.0_q2
      END DO
    END DO
  RETURN
  END SUBROUTINE dpbc_dir

!-----------------------------------------------------------------------------------!
! dpbc:  Wrap the vector dR to the boundary conditions [-ngf/2,ngf/2].
!-----------------------------------------------------------------------------------!
  SUBROUTINE dpbc(dR,ngf,ngf_2)
    REAL(q2),INTENT(IN),DIMENSION(3) :: ngf,ngf_2
    REAL(q2),INTENT(INOUT),DIMENSION(3) :: dR

    INTEGER :: i

    DO i=1,3
      DO
        IF(dR(i) > -ngf_2(i)) EXIT
        dR(i)=dR(i)+ngf(i)
      END DO
      DO
        IF(dR(i) < ngf_2(i)) EXIT
        dR(i)=dR(i)-ngf(i)
      END DO
    END DO

  RETURN
  END SUBROUTINE dpbc

!-----------------------------------------------------------------------------------!
! pbc: Wrap the REAL(q2) point (px,py,pz) to the boundary conditions [0,ngf].
!-----------------------------------------------------------------------------------!
  SUBROUTINE pbcq2(px,py,pz,rngxf,rngyf,rngzf)
    REAL(q2),INTENT(INOUT) :: px,py,pz
    REAL(q2),INTENT(IN) :: rngxf,rngyf,rngzf

    DO
      IF(px >= 0.0_q2) EXIT
      px=px+rngxf
    END DO
    DO
      IF(px < rngxf) EXIT
      px=px-rngxf
    END DO
    DO
      IF(py >= 0.0_q2) EXIT
      py=py+rngyf
    END DO
    DO
      IF(py < rngyf) EXIT
      py=py-rngyf
    END DO
    DO
      IF(pz >= 0.0_q2) EXIT
      pz=pz+rngzf
    END DO
    DO
      IF(pz < rngzf) EXIT
      pz=pz-rngzf
    END DO

  RETURN
  END SUBROUTINE pbcq2

!-----------------------------------------------------------------------------------!
! verticies: Return the points (ix,iy,iz) and (ixm,iym,izm) which are the verticies
!   closest to the REAL(q2) point (px,py,pz)
!-----------------------------------------------------------------------------------!
  SUBROUTINE verticies(px,py,pz,ix,iy,iz,ixm,iym,izm)

    REAL(q2),INTENT(INOUT) :: px,py,pz
    INTEGER,INTENT(OUT) :: ix,iy,iz,ixm,iym,izm

    ix=FLOOR(px)
    IF (ix == 0) THEN
      ix=ngxf
      ixm=1
    ELSE
      ixm=ix+1
    END IF
    iy=FLOOR(py)
    IF (iy == 0) THEN
      iy=ngyf
      iym=1
    ELSE
      iym=iy+1
    END IF
    iz=FLOOR(pz)
    IF (iz == 0) THEN
      iz=ngzf
      izm=1
    ELSE
      izm=iz+1
    END IF

  RETURN
  END SUBROUTINE verticies

!-----------------------------------------------------------------------------------!

END MODULE charge
