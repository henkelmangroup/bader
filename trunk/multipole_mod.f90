!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for doing a multipole expansion
!-----------------------------------------------------------------------------------!
MODULE multipole
  USE kind
  USE matrix
  IMPLICIT NONE

! Public, allocatable variables
  REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: dipole

  PRIVATE
  PUBLIC :: mutipole
  CONTAINS

!-----------------------------------------------------------------------------------!
! multipole: Calculate the multipole moments of atomic volumes
!-----------------------------------------------------------------------------------!
  SUBROUTINE multipole()

    REAL(q2),DIMENSION(3,3) :: B
    REAL(q2),DIMENSION(3) :: dv,v,ringf,shift
    INTEGER :: i,atom,nx,ny,nz,tenths_done
    INTEGER :: cr,count_max,t1,t2

    CALL system_clock(t1,cr,count_max)

    WRITE(*,'(/,2x,A)') 'CALCULATING DIPOLE MOMENTS'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'

    ALLOCATE(dipole(ndim,3))
    dipole=0
    shift=0.5_q2
    IF (vasp) shift=1.0_q2
    ringf(1)=1.0_q2/REAL(nxf,q2)
    ringf(2)=1.0_q2/REAL(nyf,q2) 
    ringf(3)=1.0_q2/REAL(nzf,q2)
    tenths_done=0
    DO nx=1,nxf
      IF ((nx*10/nxf) > tenths_done) THEN
        tenths_done=(nx*10/nxf)
        WRITE(*,'(A,$)') '**'
      END IF
      DO ny=1,nyf
        DO nz=1,nzf
          atom=bader_atom(bader_num(nx,ny,nz))
          v(1:3)=(/nx,ny,nz/)
          dv=(v-shift)*ringf-r_dir(atom,:)
          CALL dpbc_dir(dv)
          dipole(atom,:)=dipole(atom,:)-dv*rho(nx,ny,nz)
        END DO
      END DO
    END DO

    CALL matrix_transpose(lattice,B)
    DO i=1,ndim 
      CALL matrix_vector(B,dipole(i,1:3),v)
      dipole(i,1:3)=v/REAL(nrho,q2)
    END DO
    WRITE(*,*)

    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(1A12,1F7.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE multipole

!-----------------------------------------------------------------------------------!

END MODULE multipole
