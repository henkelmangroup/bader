!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for analyzing the charge
!
! By Andri Arnaldsson and Graeme Henkelman
! Last modified by 
!
! we need to add volume information
!
!-----------------------------------------------------------------------------------!
MODULE bader
  USE vars , ONLY : q2
  USE matrix
  USE charge
  IMPLICIT NONE

! Public parameters
  REAL(q2),PARAMETER :: bader_tol=1.0e-4_q2

! Public, allocatable variables
  INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: max_rho
  REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: bader_charge,min_dist
  REAL(q2),ALLOCATABLE,DIMENSION(:) :: bader_dist,bader_achg
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: path
  INTEGER,ALLOCATABLE,DIMENSION(:) :: bader_atom
!  INTEGER,ALLOCATABLE,DIMENSION(:) :: num_atom,nel,addup

  INTEGER :: ndim,bdim,nrho,wdim

  PRIVATE
  PUBLIC :: max_rho,bader_charge,min_dist,bader_dist,bader_achg,path,bader_atom
  PUBLIC :: bader,mindist

  CONTAINS
!-----------------------------------------------------------------------------------!
! bader: Calculate the Bader volumes and integrate to give the total charge in each
!   volume.
!-----------------------------------------------------------------------------------!
  SUBROUTINE bader()

    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: tmp
    REAL(q2),DIMENSION(3,3) :: B
    REAL(q2),DIMENSION(3) :: v,rnf
    INTEGER :: nx,ny,nz,px,py,pz,i,pdim,pnum,known_max,bnum,p,tenths_done
    INTEGER :: cr,count_max,t1,t2

    CALL system_clock(t1,cr,count_max)

    WRITE(*,'(/,2x,A)')   'CALCULATING BADER CHARGE DISTRIBUTION'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'
    bdim=64
    pdim=64
    ALLOCATE(bader_charge(bdim,4))
    ALLOCATE(path(pdim,3))
    bader_charge=0.0_q2
    max_rho=0
    bnum=0
    tenths_done=0
    DO nx=1,nxf
      IF ((nx*10/nxf) > tenths_done) THEN
        tenths_done=(nx*10/nxf)
        WRITE(*,'(A,$)') '**'
      END IF
      DO ny=1,nyf
        DO nz=1,nzf
          px=nx
          py=ny
          pz=nz
          IF(max_rho(px,py,pz) == 0) THEN
            CALL maximize(px,py,pz,pdim,pnum)
            CALL pbc(px,py,pz)
            known_max=max_rho(px,py,pz)
            IF (known_max == 0) THEN
              bnum=bnum+1
              known_max=bnum
              IF (bnum > bdim) THEN
                ALLOCATE(tmp(bdim,4))
                tmp=bader_charge
                DEALLOCATE(bader_charge)
                bdim=2*bdim
                ALLOCATE(bader_charge(bdim,4))
                bader_charge=0.0_q2
                bader_charge(1:bnum-1,1:4)=tmp
                DEALLOCATE(tmp)
              END IF
              bader_charge(bnum,1)=REAL(px,q2)
              bader_charge(bnum,2)=REAL(py,q2)
              bader_charge(bnum,3)=REAL(pz,q2)
            END IF
            DO p=1,pnum
              max_rho(path(p,1),path(p,2),path(p,3))=known_max
            END DO
          END IF
        END DO
      END DO
    END DO
    WRITE(*,*)

    DO nx=1,nxf
      DO ny=1,nyf
        DO nz=1,nzf
          bader_charge(max_rho(nx,ny,nz),4)=bader_charge(max_rho(nx,ny,nz),4)+      &
  &                    rho(nx,ny,nz)
        END DO
      END DO
    END DO

    ALLOCATE(tmp(bdim,4))
    tmp=bader_charge
    DEALLOCATE(bader_charge)
    ALLOCATE(bader_charge(bnum,4))
    bader_charge=tmp(1:bnum,:)
    bdim=bnum
    ALLOCATE(bader_atom(bdim),bader_dist(bdim),bader_achg(ndim))

! Don't have this normalization in MONDO
    bader_charge(:,4)=bader_charge(:,4)/REAL(nrho,q2)

    rnf(1)=REAL(nxf,q2)
    rnf(2)=REAL(nyf,q2)
    rnf(3)=REAL(nzf,q2)
    DO i=1,bdim
      bader_charge(i,1:3)=(bader_charge(i,1:3)-1.0_q2)/rnf
    END DO

    CALL charge2atom()

    CALL transpose_matrix(lattice,B,3,3)
    DO i=1,bdim
      CALL matrix_vector(B,bader_charge(i,1:3),v,3,3)
      bader_charge(i,1:3)=v
    END DO

    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(1A12,1F6.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE bader

!-----------------------------------------------------------------------------------!
! maximize:  From the point (px,py,pz) do a maximization on the charge density grid
!   and assign the maximum found to the max_rho array.
!-----------------------------------------------------------------------------------!
  SUBROUTINE maximize(px,py,pz,pdim,pnum)
    INTEGER,INTENT(INOUT) :: px,py,pz,pdim,pnum

    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: tmp

    pnum=1
    path(pnum,1:3)=(/px,py,pz/)

    DO
      IF(max_neighbour(px,py,pz,nxf,nyf,nzf)) THEN
        pnum=pnum+1
        IF (pnum > pdim) THEN
          ALLOCATE(tmp(pdim,3))
          tmp=path
          DEALLOCATE(path)
          pdim=2*pdim
          ALLOCATE(path(pdim,3))
          path=0.0_q2
          path(1:pnum-1,:)=tmp
          DEALLOCATE(tmp)
        END IF
        CALL pbc(px,py,pz)
        path(pnum,1:3)=(/px,py,pz/)
        IF(max_rho(px,py,pz) /= 0) EXIT
      ELSE
        EXIT
      END IF
    END DO

  RETURN
  END SUBROUTINE maximize

!-----------------------------------------------------------------------------------!
!  max_neighbour:  Do a single iteration of a maximization on the charge density 
!    grid from the point (px,py,pz).  Return a logical indicating if the current
!    point is a charge density maximum.
!-----------------------------------------------------------------------------------!
  FUNCTION max_neighbour(px,py,pz)
    INTEGER,INTENT(INOUT) :: px,py,pz
    LOGICAL :: max_neighbour

    REAL(q2) :: rho_max,rho_tmp,rho_ctr
    INTEGER :: dx,dy,dz,pxt,pyt,pzt,pxm,pym,pzm
    REAL(q2),DIMENSION(-1:1,-1:1,-1:1),SAVE :: w=RESHAPE((/           &
    &    0.5773502691896_q2,0.7071067811865_q2,0.5773502691896_q2,    &
    &    0.7071067811865_q2,1.0000000000000_q2,0.7071067811865_q2,    &
    &    0.5773502691896_q2,0.7071067811865_q2,0.5773502691896_q2,    &
    &    0.7071067811865_q2,1.0000000000000_q2,0.7071067811865_q2,    &
    &    1.0000000000000_q2,0.0000000000000_q2,1.0000000000000_q2,    &
    &    0.7071067811865_q2,1.0000000000000_q2,0.7071067811865_q2,    &
    &    0.5773502691896_q2,0.7071067811865_q2,0.5773502691896_q2,    &
    &    0.7071067811865_q2,1.0000000000000_q2,0.7071067811865_q2,    &
    &    0.5773502691896_q2,0.7071067811865_q2,0.5773502691896_q2     &
    &    /),(/3,3,3/))

    rho_max=0.0_q2
    pxm=px
    pym=py
    pzm=pz
    rho_ctr=rho_value(px,py,pz,nxf,nyf,nzf)
    DO dx=-1,1
      pxt=px+dx
      DO dy=-1,1
        pyt=py+dy
        DO dz=-1,1
          pzt=pz+dz
          rho_tmp=rho_value(pxt,pyt,pzt,nxf,nyf,nzf)
          rho_tmp=rho_ctr+w(dx,dy,dz)*(rho_tmp-rho_ctr)
          IF (rho_tmp > rho_max) THEN
            rho_max=rho_tmp
            pxm=pxt
            pym=pyt
            pzm=pzt
          END IF
        END DO
      END DO
    END DO

    max_neighbour=((pxm /= px) .or. (pym /= py) .or. (pzm /= pz))
    IF (max_neighbour) THEN
      px=pxm
      py=pym
      pz=pzm
    END IF

  RETURN
  END FUNCTION max_neighbour

!-----------------------------------------------------------------------------------!
! charge2atom: Assign an element of charge to a Bader atom.
!-----------------------------------------------------------------------------------!
  SUBROUTINE charge2atom()

    REAL(q2),DIMENSION(3,3) :: B
    REAL(q2),DIMENSION(3) :: dv,v
    REAL(q2) :: dsq,dminsq
    INTEGER :: i,j,dindex

    bader_achg=0.0_q2
    CALL transpose_matrix(lattice,B,3,3)
    DO i=1,bdim
      dv=bader_charge(i,1:3)-r_dir(1,:)
      CALL dpbc_dir(dv)
      CALL matrix_vector(B,dv,v,3,3)
      dminsq=DOT_PRODUCT(v,v)
      dindex=1
      DO j=2,ndim
        dv=bader_charge(i,1:3)-r_dir(j,:)
        CALL dpbc_dir(dv)
        CALL matrix_vector(B,dv,v,3,3)
        dsq=DOT_PRODUCT(v,v)
        IF (dsq < dminsq) THEN
          dminsq=dsq
          dindex=j
        END IF
      END DO
      bader_dist(i)=SQRT(dminsq)
      bader_atom(i)=dindex
      bader_achg(dindex)=bader_achg(dindex)+bader_charge(i,4)
    END DO
    
  RETURN
  END SUBROUTINE charge2atom

!-----------------------------------------------------------------------------------!
! mindist: Find the minimum distance from the surface of each volume to an atom
!-----------------------------------------------------------------------------------!
  SUBROUTINE mindist()

    REAL(q2),DIMENSION(3,3) :: B
    REAL(q2),DIMENSION(3) :: dv,v,ringf,shift
    REAL :: dist
    INTEGER :: i,atom,atom_tmp,nx,ny,nz,tenths_done
    INTEGER :: cr,count_max,t1,t2
    INTEGER :: dx,dy,dz,nxt,nyt,nzt
    LOGICAL :: bader_edge
 
    CALL system_clock(t1,cr,count_max)

    WRITE(*,'(/,2x,A)') 'CALCULATING MINIMUM DISTANCES TO ATOMS'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'

!   Store the minimum distance and the vector
    ALLOCATE(min_dist(ndim,4))
    min_dist=0.0_q2
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

!         Check to see if this is at the edge of an atomic volume
          atom=bader_atom(max_rho(nx,ny,nz))
          bader_edge=.false.
          neighbourloop: DO dx=-1,1
            nxt=nx+dx
            DO dy=-1,1
              nyt=ny+dy
              DO dz=-1,1
                nzt=nz+dz
                CALL pbc(nxt,nyt,nzt)
                atom_tmp=bader_atom(max_rho(nxt,nyt,nzt))
                IF (atom_tmp /= atom ) THEN
                  bader_edge=.true.
                  EXIT neighbourloop
                END IF
              END DO
            END DO
          END DO neighbourloop

!         If this is an edge cell, check if it is the closest to the atom so far
          IF (bader_edge) THEN
            v(1:3)=(/nx,ny,nz/)
            dv=(v-shift)*ringf-r_dir(atom,:)
            CALL dpbc_dir(dv)
            dist=DOT_PRODUCT(dv,dv)
            IF ((min_dist(atom,4)==0.0_q2) .OR. (dist < min_dist(atom,4))) THEN
              min_dist(atom,4)=dist
              min_dist(atom,1:3)=dv
            END IF
          END IF

        END DO
      END DO
    END DO

    CALL transpose_matrix(lattice,B,3,3)
    DO i=1,ndim 
      CALL matrix_vector(B,min_dist(i,1:3),v,3,3)
      min_dist(i,1:3)=v
      min_dist(i,4)=sqrt(DOT_PRODUCT(v,v))
!      write(*,*) min_dist(i,4)
    END DO
    WRITE(*,*)
    
    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(1A12,1F6.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'
  
  RETURN
  END SUBROUTINE mindist

!-----------------------------------------------------------------------------------!

END MODULE bader
