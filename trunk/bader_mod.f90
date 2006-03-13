!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for analyzing the charge
!
! By Andri Arnaldsson and Graeme Henkelman
! Last modified by 
!-----------------------------------------------------------------------------------!
MODULE bader_mod
  USE kind_mod , ONLY : q2
  USE matrix_mod
  USE charge_mod
  USE io_mod
  IMPLICIT NONE

! Public parameters
  REAL(q2),PARAMETER :: bader_tol=1.0e-4_q2

! Public, allocatable variables
  INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: bader_num
  REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: bader_charge,min_dist
  REAL(q2),ALLOCATABLE,DIMENSION(:) :: bader_dist,bader_achg
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: path
  INTEGER,ALLOCATABLE,DIMENSION(:) :: bader_atom
!  INTEGER,ALLOCATABLE,DIMENSION(:) :: num_atom,nel,addup

  INTEGER :: ndim,bdim,nrho,wdim

  PRIVATE
  PUBLIC :: bader_num,bader_charge,min_dist,bader_dist,bader_achg,path,bader_atom
  PUBLIC :: bader,mindist

  CONTAINS
!-----------------------------------------------------------------------------------!
! bader: Calculate the Bader volumes and integrate to give the total charge in each
!   volume.
!-----------------------------------------------------------------------------------!
  SUBROUTINE bader(chg)
    TYPE(charge_obj) :: chg
    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: tmp
    REAL(q2),DIMENSION(3,3) :: B
    REAL(q2),DIMENSION(3) :: v,rnf
    INTEGER :: nx,ny,nz,px,py,pz,i,pdim,pnum,known_max,bnum,p,tenths_done
    INTEGER :: cr,count_max,t1,t2

    CALL system_clock(t1,cr,count_max)

    WRITE(*,'(/,2x,A)')   'CALCULATING BADER CHARGE DISTRIBUTION'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'
    nxf=chg%nxf
    nyf=chg%nyf
    nzf=chg%nzf
    bdim=64
    pdim=64
    ALLOCATE(bader_charge(bdim,4))
    ALLOCATE(path(pdim,3))
    bader_charge=0.0_q2
    bader_num=0
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
          IF(bader_num(px,py,pz) == 0) THEN
            CALL maximize(px,py,pz,pdim,pnum)
            CALL pbc(px,py,pz,chg%nxf,chg%nyf,chg%nzf)
            known_max=bader_num(px,py,pz)
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
              bader_num(path(p,1),path(p,2),path(p,3))=known_max
            END DO
          END IF
        END DO
      END DO
    END DO
    WRITE(*,*)

    DO nx=1,nxf
      DO ny=1,nyf
        DO nz=1,nzf
          bader_charge(bader_num(nx,ny,nz),4)=bader_charge(bader_num(nx,ny,nz),4)+      &
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
!   and assign the maximum found to the bader_num array.
!-----------------------------------------------------------------------------------!
  SUBROUTINE maximize(chg,px,py,pz,pdim,pnum)
    TYPE(charge_obj) :: chg
    INTEGER,INTENT(INOUT) :: px,py,pz,pdim,pnum

    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: tmp

    pnum=1
    path(pnum,1:3)=(/px,py,pz/)

    DO
      IF(max_neighbour(px,py,pz,chg%nxf,chg%nyf,chg%nzf)) THEN
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
        CALL pbc(px,py,pz,chg%nxf,chg%nyf,chg%nzf)
        path(pnum,1:3)=(/px,py,pz/)
        IF(bader_num(px,py,pz) /= 0) EXIT
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

  FUNCTION max_neighbour(chg,px,py,pz)

    TYPE(charge_obj) :: chg
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
    rho_ctr=rho_val(chg,px,py,pz)
    DO dx=-1,1
      pxt=px+dx
      DO dy=-1,1
        pyt=py+dy
        DO dz=-1,1
          pzt=pz+dz
          rho_tmp=rho_val(chg,pxt,pyt,pzt)
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

  SUBROUTINE mindist(chg)

    TYPE(charge_obj) :: chg

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
          atom=bader_atom(bader_num(nx,ny,nz))
          bader_edge=.false.
          neighbourloop: DO dx=-1,1
            nxt=nx+dx
            DO dy=-1,1
              nyt=ny+dy
              DO dz=-1,1
                nzt=nz+dz
                CALL pbc(nxt,nyt,nzt,chg%nxf,chg%nyf,chg%nzf)
                atom_tmp=bader_atom(bader_num(nxt,nyt,nzt))
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

!------------------------------------------------------------------------------------!
! output: Write out a summary of the bader analysis.
!         AtomVolumes.dat: Stores the 'significant' Bader volumes associated with
!                          each atom.
!         ACF.dat        : Stores the main output to the screen.
!         BCF.dat        : Stores 'significant' Bader volumes, their coordinates and
!                          charge, atom associated and distance to it. 
!------------------------------------------------------------------------------------!

  SUBROUTINE output()

    REAL(q2) :: sum_achg
    INTEGER :: i,bdimsig,mib,mab,cc,j,nmax
    INTEGER,DIMENSION(bdim) :: rck

    mab=MAXVAL(bader_atom)
    mib=MINVAL(bader_atom)
    OPEN(100,FILE='AtomVolumes.dat',STATUS='replace',ACTION='write')
    WRITE(100,'(A)') '   Atom                     Volume(s)   '
    WRITE(100,'(A,A)') '-----------------------------------------------------------',&
  &                    '-------------'

    DO i=mib,mab
      cc=0
      rck=0
      nmax=0
      DO j=1,bdim
        IF (bader_charge(j,4) > bader_tol) THEN
          nmax=nmax+1
          IF(bader_atom(j) == i) THEN
            cc=cc+1
            rck(cc)=nmax
          END IF
        END IF
      END DO
      IF (cc == 0) CYCLE
      WRITE(100,'(2X,1I4,2X,A,2X,10000I5)') i,' ... ',rck(1:cc)
    END DO
    CLOSE(100)

    OPEN(100,FILE='ACF.dat',STATUS='replace',ACTION='write')
    WRITE(*,555) '#','X','Y','Z','VORONOI','BADER','%','MIN DIST'
    WRITE(100,555) '#','X','Y','Z','VORONOI','BADER','%','MIN DIST'
    555 FORMAT(/,4X,1A,9X,1A1,2(11X,1A1),8X,1A7,5X,1A5,9X,1A1,6X,1A10)
    WRITE(*,666) '----------------------------------------------------------------', &
  &              '------------------------------'
    WRITE(100,667) '---------------------------------------------------------------',&
  &              '------------------------------'
    666 FORMAT(1A66,1A26)
    667 FORMAT(1A65,1A27)

    
    sum_achg=SUM(bader_achg)
    DO i=1,wdim
      WRITE(*,'(1I5,7F12.4)') i,voronoi_charge(i,1:4),bader_achg(i),                 &
  &                           100.*bader_achg(i)/sum_achg,min_dist(i,4)              
      WRITE(100,'(1I5,7F12.4)') i,voronoi_charge(i,1:4),bader_achg(i),               &
  &                           100.*bader_achg(i)/sum_achg,min_dist(i,4)
    END DO
    CLOSE(100)

    bdimsig=0
    OPEN(200,FILE='BCF.dat',STATUS='replace',ACTION='write')

    WRITE(200,556) '#','X','Y','Z','CHARGE','NEAREST ATOM','DISTANCE'
    556 FORMAT(/,4X,1A1,13X,1A1,2(16X,1A1),13X,1A7,2X,1A14,2X,1A10)

    WRITE(200,668) '---------------------------------------------------------------',&
  &                '----------------------------------------'
    668 FORMAT(1A65,1A37)

    DO i=1,bdim
        IF(bader_charge(i,4) > bader_tol) THEN
           bdimsig=bdimsig+1
           WRITE(200,777) bdimsig,bader_charge(i,:),bader_atom(i),bader_dist(i)
           777 FORMAT(1I5,4(5X,F12.4),5X,1I5,5X,1F12.4)
        END IF
    END DO
    CLOSE(200)

    OPEN(300,FILE='dipole.dat',STATUS='replace',ACTION='write')
    WRITE(300,557) '#','X','Y','Z','MAGNITUDE'
    557 FORMAT(/,4X,1A1,10X,1A1,2(15X,1A1),10X,1A10)
    WRITE(300,*) '--------------------------------------------------------------------'
    DO i=1,ndim
      WRITE(300,888) i,dipole(i,:)*4.803_q2,                                         &
  &                  sqrt(DOT_PRODUCT(dipole(i,:),dipole(i,:)))*4.803_q2
!      888 FORMAT(1I5,4ES16.5)
      888 FORMAT(1I5,4F16.6)
    END DO
    CLOSE(300)

    WRITE(*,'(/,2x,A,6X,1I8)')     'NUMBER OF BADER MAXIMA FOUND: ',bdim
    WRITE(*,'(2x,A,6X,1I8)')       '    SIGNIFICANT MAXIMA FOUND: ',bdimsig
    WRITE(*,'(2x,A,2X,1F12.5,/)')  '         NUMBER OF ELECTRONS: ',                 &
  &                                          SUM(bader_charge(1:bdim,4))

  RETURN
  END SUBROUTINE output

!------------------------------------------------------------------------------------!
! write_badernum: Write out a CHGCAR type file with each entry containing an integer
!    indicating the associated Bader maximum.
! GH: change this to write the appropriate type of file
!------------------------------------------------------------------------------------!

  SUBROUTINE write_badernum()

!    INTEGER nx,ny,nz
!
!    OPEN(100,FILE='bader_rho.dat',STATUS='replace',ACTION='write')        
!    WRITE(*,'(2x,A)') 'WRITING BADER VOLUMES TO BADER_RHO.DAT'
!    WRITE(100,'(5E18.11)')
!  &              (((REAL(bader_num(nx,ny,nz),q2),nx=1,nxf),ny=1,nyf),nz=1,nzf)
!    CLOSE(100)

  RETURN
  END SUBROUTINE write_badernum

!------------------------------------------------------------------------------------!
! write_all_bader: Write out a CHGCAR type file for each of the Bader volumes found.
!------------------------------------------------------------------------------------!

  SUBROUTINE write_all_bader()

    INTEGER :: nx,ny,nz,i,AtomNum,BaderCur,tenths_done,t1,t2,cr,count_max
    CHARACTER(15) :: AtomFileName,AtomNumText
    REAL(q1),dimension(nxf,nyf,nzf) :: rho_tmp
    
    WRITE(*,'(/,2x,A)') 'WRITING BADER VOLUMES'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'
    CALL system_clock(t1,cr,count_max)
    AtomNum=0
    tenths_done=0
    DO BaderCur=1,bdim  
      DO WHILE ((BaderCur*10/bdim) > tenths_done)
        tenths_done=tenths_done+1
        WRITE(*,'(A,$)') '**'
      ENDDO
      IF(Bader_charge(BaderCur,4)>bader_tol) THEN
        AtomNum=AtomNum+1
        IF(AtomNum<10) THEN
          WRITE(AtomNumText,'(1A3,I1)') '000',AtomNum
        ELSE IF(AtomNum<100) THEN
          WRITE(AtomNumText,'(1A2,I2)') '00',AtomNum
        ELSE IF(AtomNum<1000) THEN
          WRITE(AtomNumText,'(1A1,I3)') '0',AtomNum
        ELSE
          WRITE(AtomNumText,'(I4)') AtomNum
        END IF
        AtomFileName = "Bvol"//Trim(AtomNumText(1:))//".dat"
        OPEN(100,FILE=AtomFileName) 
        WRITE(100,*)'Bader volume of Atom',AtomNum
        WRITE(100,*)'1.00'
        WRITE(100,'(3F13.6)') (lattice(i,1:3) , i=1,3)
        WRITE(100,'(110I4)') num_atom
        WRITE(100,*)'DIRECT'
        WRITE(100,'(3(2X,1F8.6))') (r_dir(i,:) , i=1,ndim)
        WRITE(100,*)
        WRITE(100,*) nxf,nyf,nzf
        rho_tmp=0.0_q1
        WHERE(bader_num == BaderCur) rho_tmp=rho
        WRITE(100,'(5E18.11)') (((rho_tmp(nx,ny,nz),nx=1,nxf),ny=1,nyf),nz=1,nzf)
        CLOSE(100)
      END IF
    END DO
    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(/,1A12,1F6.2,1A8,/)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE write_all_bader

!------------------------------------------------------------------------------------!
! write_all_atom: Write out a CHGCAR type file for each atom where all Bader volumes
!              assigned to that atom are added together. This is only done if the 
!              atoms has any 'significant' bader volumes associated with it.
!------------------------------------------------------------------------------------!

  SUBROUTINE write_all_bader()

    INTEGER :: nx,ny,nz,i,j,b,mab,mib,ik,sc,cc,tenths_done,t1,t2,cr,count_max
    INTEGER,DIMENSION(bdim) :: rck
    CHARACTER(15) :: AtomFileName,AtomNumText
    REAL(q1),dimension(nxf,nyf,nzf) :: rho_tmp

    CALL system_clock(t1,cr,count_max)
    WRITE(*,'(/,2x,A)') 'WRITING BADER VOLUMES '
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'
    tenths_done=0
    mab=MAXVAL(bader_atom)
    mib=MINVAL(bader_atom)
    sc=0
    DO ik=mib,mab
      cc=0
      rck=0
      DO j=1,bdim
        IF (bader_charge(j,4) < bader_tol) CYCLE
        IF (bader_atom(j) == ik) THEN
          cc=cc+1
          rck(cc)=j
        END IF
      END DO
      sc=sc+cc
      DO WHILE ((ik*10/(mab-mib+1)) > tenths_done)
        tenths_done=tenths_done+1
        WRITE(*,'(A,$)') '**'
      ENDDO
      IF(cc == 0) CYCLE
      IF(ik < 10) THEN
        WRITE(AtomNumText,'(1A3,I1)') '000',ik
      ELSE IF(ik<100) THEN
        WRITE(AtomNumText,'(1A2,I2)') '00',ik
      ELSE IF(ik<1000) THEN
        WRITE(AtomNumText,'(1A1,I3)') '0',ik
      ELSE
        WRITE(AtomNumText,'(I4)') ik
      END IF
      AtomFileName = "BvAt"//Trim(AtomNumText(1:))//".dat"
      OPEN(100,FILE=AtomFileName)
      WRITE(100,*)'Sum of Bader Volumes'
      WRITE(100,*)'1.00'
      WRITE(100,'(3F13.6)') (lattice(i,1:3) , i=1,3)
      WRITE(100,'(110I4)') num_atom
      WRITE(100,*)'DIRECT'
      WRITE(100,'(3(2X,1F8.6))') (r_dir(i,:) , i=1,ndim)
      WRITE(100,*)
      WRITE(100,*) nxf,nyf,nzf
      rho_tmp=0.0_q1
      DO b=1,cc
        WHERE(bader_num == rck(b)) rho_tmp=rho
      END DO
      WRITE(100,'(5E18.11)') (((rho_tmp(nx,ny,nz),nx=1,nxf),ny=1,nyf),nz=1,nzf)
      CLOSE(100)
    END DO
    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(/,1A12,1F6.2,1A8,/)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE write_all_bader

!------------------------------------------------------------------------------------!
! write_sel_bader: Write out a CHGCAR type file for specified Bader volumes by the user.
!              Volumes associated with a atom can be read from AtomVolumes.dat
!------------------------------------------------------------------------------------!

  SUBROUTINE write_sel_bader()

    CHARACTER(15) :: AtomFileName
    INTEGER,DIMENSION(bdim,2) :: volsig
    INTEGER,DIMENSION(na) :: vols
    REAL(q1),dimension(nxf,nyf,nzf) :: rho_tmp

    CALL system_clock(t1,cr,count_max)
! Correlate the number for each 'significant' bader volume to its real number
    bdimsig=0
    volsig=0
    DO i=1,bdim
      IF (bader_charge(i,4) > bader_tol) THEN
        bdimsig=bdimsig+1
        volsig(bdimsig,1)=bdimsig
        volsig(bdimsig,2)=i
      END IF
    END DO
    vols=volsig(addup,2)
    WRITE(*,'(/,2x,A)') 'WRITING SPECIFIED BADER VOLUMES '
    AtomFileName = "Bvsm.dat"
    OPEN(100,FILE=AtomFileName)
    WRITE(100,*)'Sum of Bader Volumes'
    WRITE(100,*)'1.00'
    WRITE(100,'(3F13.6)') (lattice(i,1:3) , i=1,3)
    WRITE(100,'(110I4)') num_atom
    WRITE(100,*)'DIRECT'
    WRITE(100,'(3(2X,1F8.6))') (r_dir(i,:) , i=1,ndim)
    WRITE(100,*)
    WRITE(100,*) nxf,nyf,nzf
    rho_tmp=0.0_q2
    DO b=1,na
      WHERE(bader_num == vols(b)) rho_tmp=rho
    END DO
    WRITE(100,'(5E18.11)') (((rho_tmp(nx,ny,nz),nx=1,nxf),ny=1,nyf),nz=1,nzf)
    CLOSE(100)
    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(1A12,1F6.2,1A8,/)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE write_sel_bader


!-----------------------------------------------------------------------------------!

END MODULE bader
