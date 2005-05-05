!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for analyzing the charge
!
! By Andri Arnaldsson and Graeme Henkelman
! Last modified by GH on Apr 24 2005
!-----------------------------------------------------------------------------------!
MODULE ChargeM
  USE varsM , ONLY : q2,max_rho,rho,Rdir,bader_charge,dipole,Rcar,voronoi_charge,   &
  &                  bader_dist,bader_achg,bader_atom,path,lattice,ngxf,ngyf,ngzf,  &
  &                  nrho,bdim,ndim,vasp,wdim,min_dist
  USE MatrixM
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: bader,voronoi,multipole,mindist
  CONTAINS
!-----------------------------------------------------------------------------------!
! bader: Calculate the Bader volumes and integrate to give the total charge in each
!   volume.
!-----------------------------------------------------------------------------------!
  SUBROUTINE bader()

    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: tmp
    REAL(q2),DIMENSION(3,3) :: B
    REAL(q2),DIMENSION(3) :: v,rngf
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
    DO nx=1,ngxf
      IF ((nx*10/ngxf) > tenths_done) THEN
        tenths_done=(nx*10/ngxf)
        WRITE(*,'(A,$)') '**'
      END IF
      DO ny=1,ngyf
        DO nz=1,ngzf
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

    DO nx=1,ngxf
      DO ny=1,ngyf
        DO nz=1,ngzf
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

    rngf(1)=REAL(ngxf,q2)
    rngf(2)=REAL(ngyf,q2)
    rngf(3)=REAL(ngzf,q2)
    DO i=1,bdim
      bader_charge(i,1:3)=(bader_charge(i,1:3)-1.0_q2)/rngf
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
      IF(max_neighbour(px,py,pz,ngxf,ngyf,ngzf)) THEN
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
  FUNCTION max_neighbour(px,py,pz,ngxf,ngyf,ngzf)
    INTEGER,INTENT(IN) :: ngxf,ngyf,ngzf
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
    rho_ctr=rho_value(px,py,pz,ngxf,ngyf,ngzf)
    DO dx=-1,1
      pxt=px+dx
      DO dy=-1,1
        pyt=py+dy
        DO dz=-1,1
          pzt=pz+dz
          rho_tmp=rho_value(pxt,pyt,pzt,ngxf,ngyf,ngzf)
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
!  vornoi:  Calculate the charge density populations based upon Vornoi polyhedia.
!    In this scheme each element of charge density is associated with the atom that
!    it is closest to.
!-----------------------------------------------------------------------------------!
  SUBROUTINE voronoi()

    REAL(q2),DIMENSION(ndim,3) :: Ratm
    REAL(q2),DIMENSION(3) :: Rcur,dR,ngf,ngf_2
    REAL(q2) :: dist,min_dist,shift
    INTEGER :: i,nx,ny,nz,closest,tenths_done,cr,count_max,t1,t2

    CALL system_clock(t1,cr,count_max)

    WRITE(*,'(/,2x,A)') 'CALCULATING VORONOI CHARGE DISTRIBUTION'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'

    wdim=ndim
    ALLOCATE(voronoi_charge(wdim,4))

    ngf(1)=REAL(ngxf,q2)
    ngf(2)=REAL(ngyf,q2)
    ngf(3)=REAL(ngzf,q2)
    ngf_2=REAL(ngf,q2)/2.0_q2

    shift=0.5_q2
    IF (vasp) shift=1.0_q2

    Ratm(:,1)=Rdir(:,1)*ngf(1)+shift
    Ratm(:,2)=Rdir(:,2)*ngf(2)+shift
    Ratm(:,3)=Rdir(:,3)*ngf(3)+shift

    voronoi_charge=0.0_q2
    tenths_done=0
    DO nx=1,ngxf
      Rcur(1)=REAL(nx,q2)
      IF ((nx*10/ngxf) > tenths_done) THEN
        tenths_done=(nx*10/ngxf)
!        WRITE(*,'(1X,1I4,1A6)') (tenths_done*10),'% done'
        WRITE(*,'(A,$)') '**'
      END IF
      DO ny=1,ngyf
        Rcur(2)=REAL(ny,q2)
        DO nz=1,ngzf
          Rcur(3)=REAL(nz,q2)
          closest=1
          dR=Rcur-Ratm(1,:)
          CALL dpbc(dR,ngf,ngf_2)
          min_dist=DOT_PRODUCT(dR,dR)
          DO i=2,wdim
            dR=Rcur-Ratm(i,:)
            CALL dpbc(dR,ngf,ngf_2)
            dist=DOT_PRODUCT(dR,dR)
            IF (dist < min_dist) THEN
              min_dist=dist
              closest=i
            END IF
          END DO
          voronoi_charge(closest,4)=voronoi_charge(closest,4)+                    &
                               rho_value(nx,ny,nz,ngxf,ngyf,ngzf)
        END DO
      END DO
    END DO
    WRITE(*,*)
! Don't have this normalization for MONDO
    voronoi_charge(:,4)=voronoi_charge(:,4)/REAL(nrho,q2)
    voronoi_charge(:,1:3)=Rcar

    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(1A12,1F6.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE voronoi

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
      dv=bader_charge(i,1:3)-Rdir(1,:)
      CALL dpbc_dir(dv)
      CALL matrix_vector(B,dv,v,3,3)
      dminsq=DOT_PRODUCT(v,v)
      dindex=1
      DO j=2,ndim
        dv=bader_charge(i,1:3)-Rdir(j,:)
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
    ringf(1)=1.0_q2/REAL(ngxf,q2)
    ringf(2)=1.0_q2/REAL(ngyf,q2) 
    ringf(3)=1.0_q2/REAL(ngzf,q2)
    tenths_done=0
    DO nx=1,ngxf
      IF ((nx*10/ngxf) > tenths_done) THEN
        tenths_done=(nx*10/ngxf)
        WRITE(*,'(A,$)') '**'
      END IF
      DO ny=1,ngyf
        DO nz=1,ngzf
          atom=bader_atom(max_rho(nx,ny,nz))
          v(1:3)=(/nx,ny,nz/)
          dv=(v-shift)*ringf-Rdir(atom,:)
          CALL dpbc_dir(dv)
          dipole(atom,:)=dipole(atom,:)-dv*rho(nx,ny,nz)
        END DO
      END DO
    END DO

    CALL transpose_matrix(lattice,B,3,3)
    DO i=1,ndim 
      CALL matrix_vector(B,dipole(i,1:3),v,3,3)
      dipole(i,1:3)=v/REAL(nrho,q2)
    END DO
    WRITE(*,*)

    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(1A12,1F6.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE multipole

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
    ringf(1)=1.0_q2/REAL(ngxf,q2)
    ringf(2)=1.0_q2/REAL(ngyf,q2)
    ringf(3)=1.0_q2/REAL(ngzf,q2)
    tenths_done=0
    DO nx=1,ngxf
      IF ((nx*10/ngxf) > tenths_done) THEN
        tenths_done=(nx*10/ngxf)
        WRITE(*,'(A,$)') '**'
      END IF
      DO ny=1,ngyf
        DO nz=1,ngzf

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
            dv=(v-shift)*ringf-Rdir(atom,:)
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

END MODULE ChargeM
