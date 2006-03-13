!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!    Module for reading and writing charge density data
!
! By Andri Arnaldson and Graeme Henkelman
! Last modified by 
!-----------------------------------------------------------------------------------!

MODULE cube_mod
  USE kind_mod , ONLY : q1,q2
  USE matrix_mod , ONLY : transpose_matrix,matrix_vector,volume
  USE ions_mod
  USE charge_mod
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_charge_cube, write_charge_cube
  CONTAINS

!-----------------------------------------------------------------------------------!
! read_charge_cube: Reads the charge density from a file in vasp or Gaussian cube 
!   format, by first reading the header, and then charges
!-----------------------------------------------------------------------------------!

  SUBROUTINE read_charge_cube(ions,chg,chargefile)

    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg
    CHARACTER(LEN=120) :: chargefile

    REAL(q2),DIMENSION(3,3) :: B
    REAL(q2),DIMENSION(3) :: box,v,steps
    REAL(q2) :: side,vol,tmp
    INTEGER :: i
    INTEGER,DIMENSION(110) :: nionlist=0
    CHARACTER(LEN=7) :: text
    INTEGER :: nx,ny,nz

    OPEN(100,FILE=chargefile(1:LEN_TRIM(ADJUSTL(chargefile))),STATUS='old',ACTION='read')
    WRITE(*,'(/,1A11,1A20)') 'OPEN ... ',chargefile
    WRITE(*,'(1A27)') 'GAUSSIAN-STYLE INPUT FILE'
! Skip the first two lines
    READ(100,'(2/)') 
    READ(100,*) ions%nions,ions%corner
    chg%corner=ions%corner
    ! GH: do we need this?
    ALLOCATE(ions%r_car(ions%nions,3),ions%r_dir(ions%nions,3),ions%ion_chg(ions%nions))
    READ(100,*) chg%nxf,steps(1),tmp,tmp
    READ(100,*) chg%nyf,tmp,steps(2),tmp
    READ(100,*) chg%nzf,tmp,tmp,steps(3)
    IF(chg%nxf<0) chg%nxf=(-1)*chg%nxf  ! This should really indicate the units (Bohr/Ang)
    box(1)=REAL((chg%nxf),q2)*steps(1)
    box(2)=REAL((chg%nyf),q2)*steps(2)
    box(3)=REAL((chg%nzf),q2)*steps(3)
    ions%lattice=0.0_q2
    DO i=1,3
      ions%lattice(i,i)=box(i)
    END DO
    !GH: this -steps(i) is to account for having points at the edge of the cube
    DO i=1,3
      ions%lattice(i,i)=ions%lattice(i,i)-steps(i)
    END DO
    chg%lattice=ions%lattice
    vol=volume(ions%lattice)
    CALL transpose_matrix(ions%lattice,B,3,3)
    DO i=1,ions%nions
      READ(100,*) tmp,ions%ion_chg(i),ions%r_dir(i,:)
      ions%r_dir(i,:)=(ions%r_dir(i,:)-ions%corner)/(box-steps)
      CALL matrix_vector(B,ions%r_dir(i,:),v,3,3)
      ions%r_car(i,:)=v
    END DO
    chg%nrho=chg%nxf*chg%nyf*chg%nzf
    ALLOCATE(chg%rho(chg%nxf,chg%nyf,chg%nzf))
    READ(100,*) (((chg%rho(nx,ny,nz),nx=1,chg%nxf),ny=1,chg%nyf),nz=1,chg%nzf)
    chg%rho=chg%rho*vol 
    chg%halfstep=.TRUE.
    WRITE(*,'(1A12,1I5,1A2,1I4,1A2,1I4)') 'FFT-grid: ',chg%nxf,'x',chg%nyf,'x',chg%nzf
    WRITE(*,'(2x,A,1A20)') 'CLOSE ... ', chargefile
    CLOSE(100)

  RETURN
  END SUBROUTINE read_charge_cube

!------------------------------------------------------------------------------------!
! write_charge_cube: Write out a Gaussian cube type file
!------------------------------------------------------------------------------------!

  SUBROUTINE write_charge_cube(ions,chg,chargefile)

    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg
    CHARACTER(LEN=120) :: chargefile
    
    REAL(q2),DIMENSION(3) :: box
    INTEGER :: i,nx,ny,nz

    INTEGER,DIMENSION(3) :: nxyz
    REAL(q2) :: vol
     
    DO i=1,3
      box(i)=ions%lattice(i,i)
    END DO
    vol=volume(ions%lattice)

    nxyz=(/chg%nxf,chg%nyf,chg%nzf/)

    OPEN(100,FILE=chargefile(1:LEN_TRIM(ADJUSTL(chargefile))),STATUS='replace')

    WRITE(100,*) 'Gaussian cube file'
    WRITE(100,*) 'Bader charge'
    WRITE(100,'(1I5,3(3X,1F9.6))') ions%nions,chg%corner
    WRITE(100,'(1I5,3(3X,1F9.6))') chg%nxf,chg%steps(1),0.000000,0.000000
    WRITE(100,'(1I5,3(3X,1F9.6))') chg%nyf,0.000000,chg%steps(2),0.000000
    WRITE(100,'(1I5,3(3X,1F9.6))') chg%nzf,0.000000,0.000000,chg%steps(3)
    DO i=1,ions%nions
    WRITE(100,'(1I5,3X,1F9.6,3(3X,1F9.6))')                                    & 
    &         INT(ions%ion_chg(i)),ions%ion_chg(i),                            &
    &         chg%steps*REAL((nxyz-1),q2)*ions%r_dir(i,:)+ions%corner
    END DO
    WRITE(100,'(6E13.5)') (((chg%rho(nx,ny,nz),nz=1,chg%nzf),ny=1,chg%nyf),nx=1,chg%nxf)
    CLOSE(100)
                
  RETURN
  END SUBROUTINE write_charge_cube

END MODULE cube_mod
