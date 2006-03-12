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
  PUBLIC :: read_charge,write_all_bader,write_sel_bader,write_all_atom,write_sel_atom
  PUBLIC :: rho_value,pbc,dpbc_dir,dpbc,pbcq2,verticies
  CONTAINS

!-----------------------------------------------------------------------------------!
! read_charge: Reads the charge density from a file in vasp or Gaussian cube format,
!    by first reading the header, and then charges
!-----------------------------------------------------------------------------------!

  SUBROUTINE read_charge()
      
    REAL(q2),DIMENSION(3,3) :: B
    REAL(q2),DIMENSION(3) :: box,v
    REAL(q2) :: side,vol,t
    INTEGER :: i
    INTEGER,DIMENSION(110) :: elements
    CHARACTER(LEN=7) :: text
    INTEGER :: cr,count_max,t1,t2,nx,ny,nz
      
    CALL system_clock(t1,cr,count_max)
      
    vasp=.false. 

    OPEN(100,FILE=chargefile,STATUS='old',ACTION='read',BLANK='null',PAD='yes')
    WRITE(*,'(/,1A11,1A20)') 'OPEN ... ',chargefile
    READ(100,'(6/,1A7)') text
    REWIND(100)
    IF (text == 'Direct') THEN
      elements=0
      vasp=.true.
      WRITE(*,'(2x,A)') 'VASP-STYLE INPUT FILE'
      READ(100,'(/,1F20.16)') side
      READ(100,'(3F13.6)') (lattice(i,1:3) , i=1,3)
      READ(100,'(110I4)') elements
      READ(100,*)
      DO i=1,110
        if(elements(i).eq.0) exit
      ENDDO
      natypes=i-1
      ALLOCATE(num_atom(natypes))
      DO i=1,natypes
        num_atom(i)=elements(i)
      END DO
      ndim=SUM(elements)
      lattice=side*lattice
      CALL transpose_matrix(lattice,B,3,3)
      wdim=ndim
!      ALLOCATE(r_car(ndim,3),r_dir(ndim,3),voronoi_charge(wdim,4))
! Do not allocate voronoi_charge here
      ALLOCATE(r_car(ndim,3),r_dir(ndim,3))
      DO i=1,ndim
!   Shouldn't r_dir be multiplied by side?
        READ(100,'(3(2X,1F8.6))') r_dir(i,:)
        CALL matrix_vector(B,r_dir(i,:),v,3,3)
        r_car(i,:)=v
      END DO
      READ(100,*)
      READ(100,*) ngxf,ngyf,ngzf
      ALLOCATE(max_rho(ngxf,ngyf,ngzf))
    ELSE
! temp. readin
      WRITE(*,'(1A27)') 'GAUSSIAN-STYLE INPUT FILE'
! Skip the first two lines
      READ(100,*) text
      READ(100,*) text
      READ(100,*) ndim,box
      corner=box
      wdim=ndim
!      ALLOCATE(r_car(ndim,3),r_dir(ndim,3),voronoi_charge(wdim,4),nel(ndim))
! Do not allocate voronoi_charge here
      ALLOCATE(r_car(ndim,3),r_dir(ndim,3),nel(ndim))
      READ(100,*) ngxf,steps(1),t,t
      READ(100,*) ngyf,t,steps(2),t
      READ(100,*) ngzf,t,t,steps(3)
      ALLOCATE(max_rho(ngxf,ngyf,ngzf))
      IF(ngxf<0) ngxf=(-1)*ngxf  ! This should really indicate the units (Bohr/Ang)
      box(1)=REAL((ngxf),q2)*steps(1)
      box(2)=REAL((ngyf),q2)*steps(2)
      box(3)=REAL((ngzf),q2)*steps(3)
      lattice=0.0_q2
      DO i=1,3
        lattice(i,i)=box(i)
      END DO
      vol=volume(lattice)
      DO i=1,3
        lattice(i,i)=lattice(i,i)-steps(i)
      END DO
      CALL transpose_matrix(lattice,B,3,3)
      DO i=1,ndim
        READ(100,*) nel(i),t,r_dir(i,:)
        r_dir(i,:)=(r_dir(i,:)-corner)/(box-steps)
        CALL matrix_vector(B,r_dir(i,:),v,3,3)
        r_car(i,:)=v
      END DO
    END IF
    nrho=ngxf*ngyf*ngzf
    ALLOCATE(rho(ngxf,ngyf,ngzf))
    IF (vasp) THEN
      READ(100,*) (((rho(nx,ny,nz),nx=1,ngxf),ny=1,ngyf),nz=1,ngzf)
    ELSE
      READ(100,*) (((rho(nx,ny,nz),nz=1,ngzf),ny=1,ngyf),nx=1,ngxf)
      rho=rho*vol
    END IF
    WRITE(*,'(1A12,1I5,1A2,1I4,1A2,1I4)') 'FFT-grid: ',ngxf,'x',ngyf,'x',ngzf
    WRITE(*,'(2x,A,1A20)') 'CLOSE ... ', chargefile
    CLOSE(100)
    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(1A12,1F6.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE read_charge


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
