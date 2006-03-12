!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for analyzing the charge
!
! By Andri Arnaldsson and Graeme Henkelman
! Last modified by
!-----------------------------------------------------------------------------------!
MODULE charge
  USE vars , ONLY : q2
  USE matrix
  USE ions , ONLY : r_car,r_dir,lattice,natypes,na
  USE chgcar
  USE cube
  IMPLICIT NONE

! Public, allocatable variables
  REAL(q2),ALLOCATABLE,DIMENSION(:,:,:) :: rho
  INTEGER :: ngxf,ngyf,ngzf,nrho
  LOGICAL :: halfstep

! Public, static variables
  INTEGER :: rho,nxf,nyf,nzf,nrho,halfstep

  PRIVATE
!  PUBLIC :: read_charge,write_all_bader,write_sel_bader,write_all_atom,write_sel_atom
  PUBLIC :: rho_value,pbc,dpbc_dir,dpbc,pbcq2,verticies
  CONTAINS

!-----------------------------------------------------------------------------------!
!  rho_value:  Return the density at the point (px,py,pz) taking into account the
!    boundary conditions.  This function is used to address points outside the
!    charge density array without a bunch of if statements at the place the value
!    is needed.
!-----------------------------------------------------------------------------------!
  FUNCTION rho_value(px,py,pz)
    INTEGER,INTENT(IN) :: px,py,pz
    REAL(q2) :: rho_value

    INTEGER :: pxt,pyt,pzt

    pxt=px
    pyt=py
    pzt=pz
    DO
      IF(pxt >= 1) EXIT
      pxt=pxt+nxf
    END DO
    DO
      IF(pxt <= nxf) EXIT
      pxt=pxt-nxf
    END DO
    DO
      IF(pyt >= 1) EXIT
      pyt=pyt+nyf
    END DO
    DO
      IF(pyt <= nyf) EXIT
      pyt=pyt-nyf
    END DO
    DO
      IF(pzt >= 1) EXIT
      pzt=pzt+nzf
    END DO
    DO
      IF(pzt <= nzf) EXIT
      pzt=pzt-nzf
    END DO

    rho_value=rho(pxt,pyt,pzt)

  RETURN
  END FUNCTION rho_value

!-----------------------------------------------------------------------------------!
! pbc: Wrap the point (px,py,pz) to the boundary conditions [0,nf].
!-----------------------------------------------------------------------------------!
  SUBROUTINE pbc(px,py,pz)

    INTEGER,INTENT(INOUT) :: px,py,pz

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
      py=py+ngyf
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
  SUBROUTINE dpbc(dR,nf,nf_2)
    REAL(q2),INTENT(IN),DIMENSION(3) :: nf,nf_2
    REAL(q2),INTENT(INOUT),DIMENSION(3) :: dR

    INTEGER :: i

    DO i=1,3
      DO
        IF(dR(i) > -nf_2(i)) EXIT
        dR(i)=dR(i)+nf(i)
      END DO
      DO
        IF(dR(i) < nf_2(i)) EXIT
        dR(i)=dR(i)-nf(i)
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
