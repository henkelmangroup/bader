!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for reading and writing VASP CHGCAR files
!-----------------------------------------------------------------------------------!

MODULE chgcar_mod
  USE kind_mod
  USE matrix_mod
  USE ions_mod 
  USE charge_mod 
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_charge_chgcar,write_charge_chgcar

  CONTAINS

!-----------------------------------------------------------------------------------!
! read_charge_chgcar: Reads the charge density from a file in vasp format
!-----------------------------------------------------------------------------------!

  SUBROUTINE read_charge_chgcar(ions,chg,chargefile)

    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg
    CHARACTER(LEN=128) :: chargefile

    REAL(q2) :: scalefactor
    REAL(q2),DIMENSION(3) :: dlat,dcar
    INTEGER :: i,n1,n2,n3,d1,d2,d3
    INTEGER,DIMENSION(110) :: nionlist=0

    OPEN(100,FILE=chargefile,STATUS='old',ACTION='read',BLANK='null',PAD='yes')
    WRITE(*,'(/,1A11,1A20)') 'OPEN ... ',chargefile
    WRITE(*,'(2x,A)') 'VASP-STYLE INPUT FILE'
    READ(100,'(/,1F20.16)') scalefactor
    READ(100,'(3F13.6)') (ions%lattice(i,1:3) , i=1,3)
    READ(100,'(110I4)') nionlist
    READ(100,*)
    DO i=1,110
     if(nionlist(i).eq.0) exit
    ENDDO
    ions%niontypes=i-1
    ALLOCATE(ions%num_ion(ions%niontypes))
    DO i=1,ions%niontypes
      ions%num_ion(i)=nionlist(i)
    END DO
    ions%nions=SUM(nionlist)
    ions%lattice=scalefactor*ions%lattice
    CALL matrix_transpose(ions%lattice,ions%dir2car)
    CALL matrix_3x3_inverse(ions%dir2car,ions%car2dir)
    ALLOCATE(ions%r_dir(ions%nions,3),ions%r_car(ions%nions,3))
    DO i=1,ions%nions
      READ(100,'(3(2X,1F8.6))') ions%r_dir(i,:)
      CALL matrix_vector(ions%dir2car,ions%r_dir(i,:),ions%r_car(i,:))
    END DO
    READ(100,*) 
    READ(100,*) chg%npts(:)
    chg%nrho=PRODUCT(chg%npts(:))
    chg%i_npts=1.0_q2/REAL(chg%npts,q2)
    ALLOCATE(chg%rho(chg%npts(1),chg%npts(2),chg%npts(3)))
    READ(100,*) (((chg%rho(n1,n2,n3), &
    &  n1=1,chg%npts(1)),n2=1,chg%npts(2)),n3=1,chg%npts(3))
    WRITE(*,'(1A12,1I5,1A2,1I4,1A2,1I4)') &
    &  'FFT-grid: ',chg%npts(1),'x',chg%npts(2),'x',chg%npts(3)
    WRITE(*,'(2x,A,1A20)') 'CLOSE ... ', chargefile
    CLOSE(100)
    DO i=1,3
      chg%lat2car(:,i)=ions%dir2car(:,i)/REAL(chg%npts(i),q2)
      chg%car2lat(:,i)=ions%car2dir(:,i)*REAL(chg%npts(i),q2)
    END DO

    ! origin of the lattice is at chg(1,1,1)
    chg%org_lat=(/1._q2,1._q2,1._q2/)
    chg%org_dir=(/0._q2,0._q2,0._q2/)
    chg%org_car=(/0._q2,0._q2,0._q2/)

    ! ion positions in grid points
    ALLOCATE(ions%r_lat(ions%nions,3))
    DO i=1,ions%nions
      ions%r_lat(i,:)=dir2lat(chg,ions%r_dir(i,:))
    END DO

    ! distance between neighboring points
    DO d1=-1,1
      dlat(1)=REAL(d1,q2)
      DO d2=-1,1
        dlat(2)=REAL(d2,q2)
        DO d3=-1,1
          dlat(3)=REAL(d3,q2)
          CALL matrix_vector(chg%lat2car,dlat,dcar)
          chg%lat_dist(d1,d2,d3)=SQRT(SUM(dcar*dcar))
          IF ((d1 == 0).AND.(d2 == 0).AND.(d3 == 0)) THEN
            chg%lat_i_dist(d1,d2,d3)=0._q2
          ELSE
            chg%lat_i_dist(d1,d2,d3)=1._q2/chg%lat_dist(d1,d2,d3)
          END IF
        END DO
      END DO
    END DO

  RETURN
  END SUBROUTINE read_charge_chgcar

!-----------------------------------------------------------------------------------!
! write_charge_chgcar: Write the charge density from a file in vasp format
!-----------------------------------------------------------------------------------!
    
  SUBROUTINE write_charge_chgcar(ions,chg,chargefile)
    
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg
    CHARACTER(128) :: chargefile

    INTEGER :: i,n1,n2,n3

    OPEN(100,FILE=chargefile(1:LEN_TRIM(ADJUSTL(chargefile))),STATUS='replace')
    WRITE(100,*)'Bader charge density file'
    WRITE(100,*)'1.00'
    WRITE(100,'(3F13.6)') (ions%lattice(i,1:3) , i=1,3)
    WRITE(100,'(110I4)') ions%num_ion
    WRITE(100,*)'DIRECT'
    WRITE(100,'(3(2X,1F8.6))') (ions%r_dir(i,:) , i=1,ions%nions)
    WRITE(100,*)
    WRITE(100,*) chg%npts(1),chg%npts(2),chg%npts(3)
    WRITE(100,'(5E18.11)') (((chg%rho(n1,n2,n3), &
    &  n1=1,chg%npts(1)),n2=1,chg%npts(2)),n3=1,chg%npts(3))
    CLOSE(100)

  RETURN
  END SUBROUTINE write_charge_chgcar

!-----------------------------------------------------------------------------------!

END MODULE chgcar_mod
