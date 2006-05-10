!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for analyzing the charge with the Bader atom in molecules approach
!
! By Andri Arnaldsson and Graeme Henkelman
! Last modified by 
!-----------------------------------------------------------------------------------!
MODULE bader_mod
  USE kind_mod
  USE matrix_mod
  USE options_mod
  USE ions_mod
  USE charge_mod
  USE io_mod
  IMPLICIT NONE

! Public parameters

! volnum: Bader volume number for each grid point
! volpos: position of maximum in each Bader volume
! colchg: integrated charge of each Bader volume
! ionchg: integrated charge over all Bader volumes associated with each ion
! iondist: distance from each Bader maximum to the nearest ion
! nnion: index of the nearst ion used to calculate iondist
! path: array of points along the current charge density maximization
! minsurfdist: minimum distance from the Bader surface to the included ion

! Public, allocatable variables
  TYPE bader_obj
    REAL(q2) :: tol
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: volnum
    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: volpos_lat,volpos_car,volpos_dir
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: volchg,iondist,ionchg,minsurfdist
    INTEGER,ALLOCATABLE,DIMENSION(:) :: nnion
    REAL(q2) :: stepsize
    INTEGER nvols
  END TYPE

  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: path

  PRIVATE
  PUBLIC :: bader_obj
  PUBLIC :: bader_calc,bader_mindist,bader_output

  CONTAINS
!-----------------------------------------------------------------------------------!
! bader_calc: Calculate the Bader volumes and integrate to give the total
!   charge in each volume.
!-----------------------------------------------------------------------------------!
  SUBROUTINE bader_calc(bdr,ions,chg,opts)

    TYPE(bader_obj) :: bdr
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg
    TYPE(options_obj) :: opts

    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: tmpvolpos
    REAL(q2),DIMENSION(3) :: v
    INTEGER,DIMENSION(3) :: p
    INTEGER :: n1,n2,n3,i,known_max,pn,tenths_done
    INTEGER :: pdim,pnum,bdim,bnum
    INTEGER :: cr,count_max,t1,t2

    REAL(q2),DIMENSION(3) :: grad,voxlen
    REAL(q2) :: rho

    CALL system_clock(t1,cr,count_max)
    
    WRITE(*,'(/,2x,A)')   'CALCULATING BADER CHARGE DISTRIBUTION'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'

! copy bader variable from opts
    bdr%tol=opts%badertol
    bdr%stepsize=opts%stepsize
    IF (opts%stepsize == 0) THEN  ! check for transpose error
      DO i=1,3
       voxlen(i)=SQRT(SUM(chg%lat2car(:,i)*chg%lat2car(:,i)))
      END DO
      bdr%stepsize=MINVAL(voxlen(:))
    END IF

    bdim=64  ! temporary number of bader volumes, will be expanded as needed
    pdim=64  ! temporary path length, also expanded as needed
    ALLOCATE(bdr%volpos_lat(bdim,3))
    ALLOCATE(path(pdim,3))
    ALLOCATE(bdr%volnum(chg%npts(1),chg%npts(2),chg%npts(3)))
    bdr%volchg=0.0_q2
    bdr%volnum=0
    bnum=0
    bdr%nvols=0  ! True number of Bader volumes

    tenths_done=0
    DO n1=1,chg%npts(1)
      IF ((n1*10/chg%npts(1)) > tenths_done) THEN
        tenths_done=(n1*10/chg%npts(1))
        WRITE(*,'(A,$)') '**'
      END IF
      DO n2=1,chg%npts(2)
        DO n3=1,chg%npts(3)
          p=(/n1,n2,n3/)
          IF(bdr%volnum(p(1),p(2),p(3)) == 0) THEN
            IF(opts%bader_opt == opts%bader_offgrid) THEN
              CALL max_offgrid(bdr,chg,p,pdim,pnum)
            ELSE
              CALL max_ongrid(bdr,chg,p,pdim,pnum)
            END IF
            known_max=bdr%volnum(p(1),p(2),p(3))
            IF (known_max == 0) THEN
              bnum=bnum+1
              known_max=bnum
              IF (bnum > bdim) THEN
                ALLOCATE(tmpvolpos(bdim,3))
                tmpvolpos=bdr%volpos_lat
                DEALLOCATE(bdr%volpos_lat)
                bdim=2*bdim
                ALLOCATE(bdr%volpos_lat(bdim,3))
                bdr%volpos_lat(1:bnum-1,:) = tmpvolpos
                DEALLOCATE(tmpvolpos)
              END IF
              bdr%volpos_lat(bnum,:) = REAL(p,q2)
            END IF
            DO i=1,pnum
              bdr%volnum(path(i,1),path(i,2),path(i,3)) = known_max
            END DO
          END IF
        END DO
      END DO
    END DO
    
! Sum up the charge included in each volume
    bdr%nvols=bnum
    ALLOCATE(bdr%volchg(bdr%nvols))
    bdr%volchg=0.0_q2
    DO n1=1,chg%npts(1)
      DO n2=1,chg%npts(2)
        DO n3=1,chg%npts(3)
          bdr%volchg(bdr%volnum(n1,n2,n3)) = &
          &  bdr%volchg(bdr%volnum(n1,n2,n3)) + chg%rho(n1,n2,n3)
        END DO
      END DO
    END DO

    ALLOCATE(tmpvolpos(bdim,3))
    tmpvolpos=bdr%volpos_lat
    DEALLOCATE(bdr%volpos_lat)
    ALLOCATE(bdr%volpos_lat(bdr%nvols,3))
    bdr%volpos_lat=tmpvolpos(1:bdr%nvols,:)
    DEALLOCATE(tmpvolpos) 
 
    ALLOCATE(bdr%nnion(bdr%nvols),bdr%iondist(bdr%nvols),bdr%ionchg(ions%nions))

    bdr%volchg=bdr%volchg/REAL(chg%nrho,q2)

    ! volpos in direct coordinates
    ALLOCATE(bdr%volpos_dir(bdr%nvols,3))
    ALLOCATE(bdr%volpos_car(bdr%nvols,3))
    DO i=1,bdr%nvols
      bdr%volpos_dir(i,:)=lat2dir(chg,bdr%volpos_lat(i,:))
      bdr%volpos_car(i,:)=lat2car(chg,bdr%volpos_lat(i,:))
    END DO

    CALL assign_chg2atom(bdr,ions,chg)

    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(2/,1A12,1F6.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE bader_calc

!-----------------------------------------------------------------------------------!
! max_offgrid:  From the point (px,py,pz) do a maximization in the charge density
!-----------------------------------------------------------------------------------!
  SUBROUTINE max_offgrid(bdr,chg,p,pdim,pnum)
    
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3),INTENT(INOUT) :: p
    INTEGER,INTENT(INOUT) :: pdim,pnum
    
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: tmp
    REAL :: rx,ry,rz
    
    pnum=1
    path(pnum,1:3)=p
    DO
      ! GH: replace this with a steepest ascent trajectory
      IF(step_offgrid(chg,rx,ry,rz)) THEN

        ! Only if this is a new point
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
!        CALL pbc(px,py,pz,nxf,nyf,nzf)
!        path(pnum,1:3)=(/px,py,pz/)
!        IF(bdr%volnum(px,py,pz) /= 0) EXIT

      ELSE
        EXIT
      END IF
    END DO
    
  RETURN
  END SUBROUTINE max_offgrid

!-----------------------------------------------------------------------------------!
!  step_offgrid:
!-----------------------------------------------------------------------------------!

  FUNCTION step_offgrid(chg,rx,ry,rz)

    TYPE(charge_obj) :: chg
    REAL,INTENT(INOUT) :: rx,ry,rz
    LOGICAL :: step_offgrid

! get grad_rho
! move rx -> rx+step*grad_rho

!GH: change this for non-orthogonal cells
    REAL(q2) :: rho_max,rho_tmp,rho_ctr

!   step_offgrid=.TRUE.
!   step_offgrid=is_max(p)
    step_offgrid=.TRUE.
  RETURN
  END FUNCTION step_offgrid

!-----------------------------------------------------------------------------------!
! max_ongrid:  From the point (p(1),p(2),p(3)) do a maximization on the charge
!   density grid and assign the maximum found to the volnum array.
!-----------------------------------------------------------------------------------!
  SUBROUTINE max_ongrid(bdr,chg,p,pdim,pnum)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3),INTENT(INOUT) :: p
    INTEGER,INTENT(INOUT) :: pdim,pnum

    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: tmp

    pnum=1
    path(pnum,1:3)=p
    DO
      IF(step_ongrid(chg,p(1),p(2),p(3))) THEN
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
        CALL pbc(p,chg%npts)
        path(pnum,1:3)=p
        IF(bdr%volnum(p(1),p(2),p(3)) /= 0) EXIT
      ELSE
        EXIT
      END IF
    END DO

  RETURN
  END SUBROUTINE max_ongrid

!-----------------------------------------------------------------------------------!
!  step_ongrid:  Do a single iteration of a maximization on the charge density 
!    grid from the point (px,py,pz).  Return a logical indicating if the current
!    point is a charge density maximum.
!-----------------------------------------------------------------------------------!

  FUNCTION step_ongrid(chg,px,py,pz)

    TYPE(charge_obj) :: chg
    INTEGER,INTENT(INOUT) :: px,py,pz
    LOGICAL :: step_ongrid

    REAL(q2) :: rho_max,rho_tmp,rho_ctr
    INTEGER :: dx,dy,dz,pxt,pyt,pzt,pxm,pym,pzm

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
          rho_tmp=rho_ctr+(rho_tmp-rho_ctr)*chg%lat_i_dist(dx,dy,dz)
          IF (rho_tmp > rho_max) THEN
            rho_max=rho_tmp
            pxm=pxt
            pym=pyt
            pzm=pzt
          END IF
        END DO
      END DO
    END DO

    step_ongrid=((pxm /= px) .or. (pym /= py) .or. (pzm /= pz))
    IF (step_ongrid) THEN
      px=pxm
      py=pym
      pz=pzm
    END IF

  RETURN
  END FUNCTION step_ongrid

!-----------------------------------------------------------------------------------!
! assign_chg2atom: Assign an element of charge to a Bader atom.
!-----------------------------------------------------------------------------------!

  SUBROUTINE assign_chg2atom(bdr,ions,chg)

    TYPE(bader_obj) :: bdr
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    REAL(q2),DIMENSION(3) :: dv,v
    REAL(q2) :: dsq,dminsq
    INTEGER :: i,j,dindex

    bdr%ionchg=0.0_q2
    DO i=1,bdr%nvols
      dv=bdr%volpos_dir(i,:)-ions%r_dir(1,:)
      CALL dpbc_dir(dv)
      CALL matrix_vector(ions%dir2car,dv,v)
      dminsq=DOT_PRODUCT(v,v)
      dindex=1
      DO j=2,ions%nions
        dv=bdr%volpos_dir(i,:)-ions%r_dir(j,:)
        CALL dpbc_dir(dv)
        CALL matrix_vector(ions%dir2car,dv,v)
        dsq=DOT_PRODUCT(v,v)
        IF (dsq < dminsq) THEN
          dminsq=dsq
          dindex=j
        END IF
      END DO
      bdr%iondist(i)=SQRT(dminsq)
      bdr%nnion(i)=dindex
      bdr%ionchg(dindex)=bdr%ionchg(dindex)+bdr%volchg(i)
    END DO
 
  RETURN
  END SUBROUTINE assign_chg2atom

!-----------------------------------------------------------------------------------!
! bader_mindist: Find the minimum distance from the surface of each volume to 
!-----------------------------------------------------------------------------------!

  SUBROUTINE bader_mindist(bdr,ions,chg)

    TYPE(bader_obj) :: bdr
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    REAL(q2),DIMENSION(3) :: shift,v,dv_dir,dv_car
    INTEGER,DIMENSION(3) :: n
    REAL :: dist
    INTEGER :: i,atom,atom_tmp,n1,n2,n3,d1,d2,d3
    INTEGER :: cr,count_max,t1,t2,tenths_done
    LOGICAL :: surfflag

    CALL system_clock(t1,cr,count_max)

    WRITE(*,'(/,2x,A)') 'CALCULATING MINIMUM DISTANCES TO ATOMS'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'

!   Store the minimum distance and the vector
    ALLOCATE(bdr%minsurfdist(ions%nions))
    bdr%minsurfdist=0.0_q2

    tenths_done=0
    DO n1=1,chg%npts(1)
      IF ((n1*10/chg%npts(1)) > tenths_done) THEN
        tenths_done=(n1*10/chg%npts(1))
        WRITE(*,'(A,$)') '**'
      END IF
      DO n2=1,chg%npts(2)
        DO n3=1,chg%npts(3)

!         Check to see if this is at the edge of an atomic volume
          atom=bdr%nnion(bdr%volnum(n1,n2,n3))
          surfflag=.FALSE.
          neighbourloop: DO d1=-1,1
            n(1)=n1+d1
            DO d2=-1,1
              n(2)=n2+d2
              DO d3=-1,1
                n(3)=n3+d3
                CALL pbc(n,chg%npts)
                atom_tmp=bdr%nnion(bdr%volnum(n(1),n(2),n(3)))
                IF (atom_tmp /= atom ) THEN
                  surfflag=.TRUE.
                  EXIT neighbourloop
                END IF
              END DO
            END DO
          END DO neighbourloop

!         If this is an edge cell, check if it is the closest to the atom so far
          IF (surfflag) THEN
            v(1:3)=(/n1,n2,n3/)
            dv_dir=(v-chg%org_lat)/REAL(chg%npts,q2)-ions%r_dir(atom,:)
            CALL dpbc_dir(dv_dir)
            CALL matrix_vector(ions%dir2car,dv_dir,dv_car)
            dist=DOT_PRODUCT(dv_car,dv_car)
            IF ((bdr%minsurfdist(atom) == 0.0_q2) .OR.  &
            &   (bdr%minsurfdist(atom) > dist)) THEN
              bdr%minsurfdist(atom)=dist
            END IF
          END IF

        END DO
      END DO
    END DO

    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(2/,1A12,1F6.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'
  
  RETURN
  END SUBROUTINE bader_mindist

!------------------------------------------------------------------------------------!
! write_volnum: Write out a CHGCAR type file with each entry containing an integer
!    indicating the associated Bader maximum.
!------------------------------------------------------------------------------------!

  SUBROUTINE write_volnum(bdr,opts,ions,chg)

     TYPE(bader_obj) :: bdr
     TYPE(options_obj) :: opts
     TYPE(ions_obj) :: ions
     TYPE(charge_obj) :: chg

     TYPE(charge_obj) :: tmp
     INTEGER :: nx,ny,nz
     CHARACTER(LEN=120) :: filename

     tmp=chg
     tmp%rho=bdr%volnum
     
     filename='VOLUME_INDEX'
     CALL write_charge(ions,chg,opts,filename)

  RETURN
  END SUBROUTINE write_volnum

!------------------------------------------------------------------------------------!
! write_all_bader: Write out a CHGCAR type file for each of the Bader volumes found.
!------------------------------------------------------------------------------------!

  SUBROUTINE write_all_bader(bdr,opts,ions,chg)

    TYPE(bader_obj) :: bdr
    TYPE(options_obj) :: opts
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    TYPE(charge_obj) :: tmp
    INTEGER :: nx,ny,nz,i,atomnum,badercur,tenths_done,t1,t2,cr,count_max
!    CHARACTER(LEN=120) :: atomfilename,atomnumtext
    CHARACTER(LEN=120) :: atomfilename
    
    tmp=chg

    WRITE(*,'(/,2x,A)') 'WRITING BADER VOLUMES'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'
    CALL system_clock(t1,cr,count_max)
    atomnum=0
    tenths_done=0

    DO badercur=1,bdr%nvols
      DO WHILE ((badercur*10/bdr%nvols) > tenths_done)
        tenths_done=tenths_done+1
        WRITE(*,'(A,$)') '**'
      ENDDO
      IF(bdr%volchg(badercur) > bdr%tol) THEN
        atomnum=atomnum+1
        WRITE(atomfilename,'(A4,I4.4,A4)') "Bvol",atomnum,".dat"
!        IF(atomnum<10) THEN
!          WRITE(atomnumtext,'(1A3,I1)') '000',atomnum
!        ELSE IF(atomnum<100) THEN
!          WRITE(atomnumtext,'(1A2,I2)') '00',atomnum
!        ELSE IF(atomnum<1000) THEN
!          WRITE(atomnumtext,'(1A1,I3)') '0',atomnum
!        ELSE
!          WRITE(atomnumtext,'(I4)') atomnum
!        END IF
!        atomfilename = "Bvol"//Trim(atomnumtext(1:))//".dat"

        tmp%rho=0.0_q2
        WHERE(bdr%volnum == badercur) tmp%rho=chg%rho
        CALL write_charge(ions,chg,opts,atomfilename)

      END IF
    END DO

    DEALLOCATE(tmp%rho)

    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(2/,1A12,1F6.2,1A8,/)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE write_all_bader

!------------------------------------------------------------------------------------!
! write_all_atom: Write out a CHGCAR type file for each atom where all Bader volumes
!              assigned to that atom are added together. This is only done if the 
!              atoms has any 'significant' bader volumes associated with it.
!------------------------------------------------------------------------------------!

  SUBROUTINE write_all_atom(bdr,opts,ions,chg)

    TYPE(bader_obj) :: bdr
    TYPE(options_obj) :: opts
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    TYPE(charge_obj) :: tmp

    INTEGER :: nx,ny,nz,i,j,b,mab,mib,ik,sc,cc,tenths_done,t1,t2,cr,count_max
    INTEGER,DIMENSION(bdr%nvols) :: rck
!    CHARACTER(LEN=120) :: atomfilename,atomnumtext
    CHARACTER(LEN=120) :: atomfilename

    CALL system_clock(t1,cr,count_max)

    tmp=chg

    WRITE(*,'(/,2x,A)') 'WRITING BADER VOLUMES '
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'
    tenths_done=0
    mab=MAXVAL(bdr%nnion)
    mib=MINVAL(bdr%nnion)
    sc=0

    DO ik=mib,mab
      cc=0
      rck=0
      DO j=1,bdr%nvols
        IF (bdr%volchg(j) < bdr%tol) CYCLE
        IF (bdr%nnion(j) == ik) THEN
          cc=cc+1
          rck(cc)=j
        END IF
      END DO
      sc=sc+cc
      DO WHILE ((ik*10/(mab-mib+1)) > tenths_done)
        tenths_done=tenths_done+1
        WRITE(*,'(A,$)') '**'
      END DO
      IF(cc == 0) CYCLE
!      IF(ik < 10) THEN
!        WRITE(atomnumtext,'(1A3,I1)') '000',ik
!      ELSE IF(ik<100) THEN
!        WRITE(atomnumtext,'(1A2,I2)') '00',ik
!      ELSE IF(ik<1000) THEN
!        WRITE(atomnumtext,'(1A1,I3)') '0',ik
!      ELSE
!        WRITE(atomnumtext,'(I4)') ik
!      END IF
!      atomfilename = "BvAt"//Trim(atomnumtext(1:))//".dat"
      WRITE(atomfilename,'(A4,I4.4,A4)') "BvAt",ik,".dat"

      tmp%rho=0.0_q2
      DO b=1,cc
        WHERE(bdr%volnum == rck(b)) tmp%rho=chg%rho
      END DO 
      CALL write_charge(ions,chg,opts,atomfilename)

    END DO
    DEALLOCATE(tmp%rho)

    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(2/,1A12,1F6.2,1A8,/)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE write_all_atom

!------------------------------------------------------------------------------------!
! write_sel_bader: Write out a CHGCAR type file for the user specified Bader volumes.
!              Volumes associated with a atom can be read from AtomVolumes.dat
!------------------------------------------------------------------------------------!

  SUBROUTINE write_sel_bader(bdr,opts,ions,chg)

    TYPE(bader_obj) :: bdr
    TYPE(options_obj) :: opts
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    TYPE(charge_obj) :: tmp
    CHARACTER(LEN=120) :: atomfilename
    INTEGER,DIMENSION(bdr%nvols,2) :: volsig
!    INTEGER,DIMENSION(na) :: vols
    INTEGER :: cr,count_max,t1,t2,i,bdimsig

    CALL system_clock(t1,cr,count_max)

    tmp=chg

! Correlate the number for each 'significant' bader volume to its real number
    bdimsig=0
    volsig=0

    DO i=1,bdr%nvols
      IF (bdr%volchg(i) > bdr%tol) THEN
        bdimsig=bdimsig+1
        volsig(bdimsig,1)=bdimsig
        volsig(bdimsig,2)=i
      END IF
    END DO
!    vols=volsig(addup,2)
    WRITE(*,'(/,2x,A)') 'WRITING SPECIFIED BADER VOLUMES '
    atomfilename = "Bvsm.dat"

    tmp%rho=0.0_q2
! fix this when we get na input through options
!    DO b=1,na
!      WHERE(bdr%volnum == vols(b)) tmp%rho=chg%rho
!    END DO
    CALL write_charge(ions,chg,opts,atomfilename)

    DEALLOCATE(tmp%rho)

    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(2/,1A12,1F6.2,1A8,/)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE write_sel_bader

!------------------------------------------------------------------------------------!
! bader_output: Write out a summary of the bader analysis.
!         AtomVolumes.dat: Stores the 'significant' Bader volumes associated with
!                          each atom.
!         ACF.dat        : Stores the main output to the screen.
!         BCF.dat        : Stores 'significant' Bader volumes, their coordinates and
!                          charge, atom associated and distance to it. 
!------------------------------------------------------------------------------------!

  SUBROUTINE bader_output(bdr,ions,chg)

    TYPE(bader_obj) :: bdr
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

  
    REAL(q2) :: sum_ionchg
    INTEGER :: i,bdimsig,mib,mab,cc,j,nmax
    INTEGER,DIMENSION(bdr%nvols) :: rck
  
    mab=MAXVAL(bdr%nnion)
    mib=MINVAL(bdr%nnion)
    OPEN(100,FILE='AVF.dat',STATUS='replace',ACTION='write')
    WRITE(100,'(A)') '   Atom                     Volume(s)   '
    WRITE(100,'(A,A)') '-----------------------------------------------------------',&
    &                  '-------------'

    DO i=mib,mab
      cc=0
      rck=0
      nmax=0
      DO j=1,bdr%nvols
        IF (bdr%volchg(j) > bdr%tol) THEN
          nmax=nmax+1
          IF(bdr%nnion(j) == i) THEN
            cc=cc+1
            rck(cc)=nmax
          END IF
        END IF
      END DO 
      IF (cc == 0) CYCLE
      WRITE(100,'(2X,1I4,2X,A,2X,10000I5)') i,' ... ',rck(1:cc)
    END DO
    CLOSE(100)
    
    WRITE(*,'(/,A41)') 'WRITING BADER ATOMIC CHARGES TO ACF.dat'
    WRITE(*,'(A41,/)') 'WRITING BADER VOLUME CHARGES TO BCF.dat'

    OPEN(100,FILE='ACF.dat',STATUS='replace',ACTION='write')
    WRITE(100,555) '#','X','Y','Z','CHARGE','MIN DIST'
    555 FORMAT(4X,1A,9X,1A1,2(11X,1A1),8X,1A6,5X,1A8)
    WRITE(100,666) '----------------------------------------------------------------'
    666 FORMAT(1A66)
    
    sum_ionchg=SUM(bdr%ionchg)
    DO i=1,ions%nions
      WRITE(100,'(1I5,6F12.4)') i,ions%r_car(i,:),bdr%ionchg(i),bdr%minsurfdist(i)
    END DO
    CLOSE(100)

    bdimsig=0
    OPEN(200,FILE='BCF.dat',STATUS='replace',ACTION='write')

    WRITE(200,556) '#','X','Y','Z','CHARGE','ATOM','DISTANCE'
    556 FORMAT(4X,1A1,9X,1A1,2(11X,1A1),8X,1A6,5X,1A4,4X,1A8)
    
    WRITE(200,668) '---------------------------------------------------------------',&
    &              '----------'
    668 FORMAT(1A65,1A10)
   
    DO i=1,bdr%nvols
        IF(bdr%volchg(i) > bdr%tol) THEN
           bdimsig=bdimsig+1
           WRITE(200,777) bdimsig,bdr%volpos_car(i,:),bdr%volchg(i), &
           &              bdr%nnion(i),bdr%iondist(i)
           777 FORMAT(1I5,4F12.4,3X,1I5,1F12.4)
        END IF
    END DO
    CLOSE(200)

    WRITE(*,'(2x,A,6X,1I8)')       'NUMBER OF BADER MAXIMA FOUND: ',bdr%nvols
    WRITE(*,'(2x,A,6X,1I8)')       '    SIGNIFICANT MAXIMA FOUND: ',bdimsig
    WRITE(*,'(2x,A,2X,1F12.5,/)')  '         NUMBER OF ELECTRONS: ', &
    &                                        SUM(bdr%volchg(1:bdr%nvols))

  RETURN
  END SUBROUTINE bader_output

!-----------------------------------------------------------------------------------!

END MODULE bader_mod
