!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for analyzing the charge with the Bader atom in molecules approach
!
! By Andri Arnaldsson and Graeme Henkelman
! Last modified by GH on May 11, 2009
!-----------------------------------------------------------------------------------!

MODULE bader_mod

  USE kind_mod
  USE matrix_mod
  USE options_mod
  USE ions_mod
  USE charge_mod
  USE io_mod
  IMPLICIT NONE

! Variable descriptions:

! volnum: Bader volume number for each grid point
! volpos: position of maximum in each Bader volume
! colchg: integrated charge of each Bader volume
! ionchg: integrated charge over all Bader volumes associated with each ion
! iondist: distance from each Bader maximum to the nearest ion
! nnion: index of the nearst ion used to calculate iondist
! path: array of points along the current charge density maximization
! minsurfdist: minimum distance from the Bader surface to the included ion

  TYPE bader_obj
    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: volpos_lat,volpos_car,volpos_dir
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: volchg,iondist,ionchg,minsurfdist
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: volnum
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: path
    INTEGER,ALLOCATABLE,DIMENSION(:) :: nnion
    REAL(q2) :: stepsize, tol
    INTEGER nvols,pnum,bnum,pdim,bdim
  END TYPE

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
    INTEGER :: n1,n2,n3,i,path_volnum,pn,tenths_done
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

    bdr%bdim=64
    bdr%pdim=64
    ALLOCATE(bdr%volpos_lat(bdr%bdim,3)) ! will be expanded as needed
    ALLOCATE(bdr%path(bdr%pdim,3))
    ALLOCATE(bdr%volnum(chg%npts(1),chg%npts(2),chg%npts(3)))
    bdr%volchg=0.0_q2
    bdr%volnum=0
    bdr%bnum=0
    bdr%nvols=0  ! True number of Bader volumes
    
    DO i=1,3
      print*,'ions%car:',ions%dir2car(:,i)
    END DO  
    tenths_done=0
    DO n1=1,chg%npts(1)
      IF ((n1*10/chg%npts(1)) > tenths_done) THEN
        tenths_done=(n1*10/chg%npts(1))
        WRITE(*,'(A,$)') '**'
      END IF
      DO n2=1,chg%npts(2)
        DO n3=1,chg%npts(3)
          p=(/n1,n2,n3/)
          !write(*,*) ''
          !write(*,'(A,3I4)') ' grd: ',p
          IF(bdr%volnum(p(1),p(2),p(3)) == 0) THEN
            IF(opts%bader_opt == opts%bader_offgrid) THEN
              CALL max_offgrid(bdr,chg,p)
            ELSE IF(opts%bader_opt == opts%bader_ongrid) THEN
              CALL max_ongrid(bdr,chg,p)
            ELSE 
              CALL max_neargrid(bdr,chg,p)
            END IF
            path_volnum=bdr%volnum(p(1),p(2),p(3))
            IF (path_volnum == 0) THEN
              IF (bdr%bnum >= bdr%bdim) THEN
                CALL reallocate_volpos(bdr,bdr%bdim*2)
              END IF
              bdr%bnum=bdr%bnum+1
              !write(*,*) '   new max',bdr%bnum
              path_volnum=bdr%bnum
              bdr%volpos_lat(bdr%bnum,:) = REAL(p,q2)
            END IF
            DO i=1,bdr%pnum
     IF ( (bdr%volnum(bdr%path(i,1),bdr%path(i,2),bdr%path(i,3)) /= 0) .AND. &
     &    (bdr%volnum(bdr%path(i,1),bdr%path(i,2),bdr%path(i,3)) /= path_volnum) ) &
     &    write(*,*) '   overwrite along path'
              bdr%volnum(bdr%path(i,1),bdr%path(i,2),bdr%path(i,3)) = path_volnum
            END DO
          END IF
        END DO
      END DO
    END DO
    
!    print*,'bdr%bnum at the end of normal loop',bdr%bnum
    print*,'opts%checkep_times',opts%checkep_times
    DO i=1,opts%checkep_times
      CALL checkep(bdr,chg,opts)
      print*,'times',i
    END DO      
!    print*,'finish call checkedp'
 
    ! The total number of bader volumes is now known
    bdr%nvols=bdr%bnum
    CALL reallocate_volpos(bdr,bdr%nvols)
    ALLOCATE(bdr%volpos_dir(bdr%nvols,3))
    ALLOCATE(bdr%volpos_car(bdr%nvols,3))
    DO i=1,bdr%nvols
      bdr%volpos_dir(i,:)=lat2dir(chg,bdr%volpos_lat(i,:))
      bdr%volpos_car(i,:)=lat2car(chg,bdr%volpos_lat(i,:))
    END DO
!    print*,'end of cal total bader vol'


    ! Sum up the charge included in each volume
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
    bdr%volchg=bdr%volchg/REAL(chg%nrho,q2)
    print*,'end of sum up charge'

    ALLOCATE(bdr%nnion(bdr%nvols),bdr%iondist(bdr%nvols),bdr%ionchg(ions%nions))
    CALL assign_chg2atom(bdr,ions,chg)

    DEALLOCATE(bdr%path)

    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(2/,1A12,1F6.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE bader_calc

!-----------------------------------------------------------------------------------!
! max_offgrid:  From the point p, do a maximization in the charge density
!-----------------------------------------------------------------------------------!

  SUBROUTINE max_offgrid(bdr,chg,p)
    
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3),INTENT(INOUT) :: p
    INTEGER,DIMENSION(3) :: pm
    REAL(q2) :: cp,cm 
    REAL(q2),DIMENSION(3) :: r
 
    bdr%pnum=1
    bdr%path(bdr%pnum,:)=p
    IF (is_max(chg,p)) THEN
      write(*,*) '   max found (init)'
      RETURN
    END IF
    r=REAL(p,q2)

!    write(*,*) 'max: ',p
    DO
      cp=rho_val(chg,p(1),p(2),p(3))
      CALL step_offgrid(bdr,chg,r)
      pm=to_lat(chg,r)
      cm=rho_val(chg,pm(1),pm(2),pm(3))
     
      IF(cm < cp) EXIT
      p=pm

      ! if the point is new, add it to the path
      IF (.NOT.ALL(p(:)==bdr%path(bdr%pnum,:))) THEN
        IF (bdr%pnum >= bdr%pdim) THEN
          CALL reallocate_path(bdr,bdr%pdim*2)
        ENDIF
        bdr%pnum=bdr%pnum+1
        CALL pbc(p,chg%npts)
        bdr%path(bdr%pnum,:)=p
      END IF

      ! quit if this is a maximum or a knoown volume number
      IF (is_max(chg,p).OR.(known_volnum(bdr,chg,r)/=0)) EXIT

    END DO
    
  RETURN
  END SUBROUTINE max_offgrid

!-----------------------------------------------------------------------------------!
!  step_offgrid: step a distance of StepSize along rho_grad
!-----------------------------------------------------------------------------------!

  SUBROUTINE step_offgrid(bdr,chg,r)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    REAL(q2),DIMENSION(3),INTENT(INOUT) :: r

    REAL(q2),DIMENSION(3) :: grad,dr_car,dr_lat
    REAL(q2) :: rho

!    write(*,'(A,3F12.4)') '   r_before: ',r
    grad = rho_grad(chg,r,rho)
!    write(*,'(A,3F12.4)') '   rho_grad: ',grad
!db    write(*,'(A,F12.4)') '   stepsize: ',bdr%stepsize
    dr_car = grad*bdr%stepsize/SQRT(SUM(grad*grad))
!db    write(*,'(A,3F12.4)') '     dr_car: ',dr_car
    CALL matrix_vector(chg%car2lat,dr_car,dr_lat)
!db    write(*,'(A,3F12.4)') '     dr_lat: ',dr_lat
    r = r + dr_lat
!db    write(*,'(A,3F12.4)') '    r_after: ',r

  RETURN
  END SUBROUTINE step_offgrid

!-----------------------------------------------------------------------------------!
! max_ongrid:  From the point p do a maximization on the charge density grid and
!   assign the maximum found to the volnum array.
!-----------------------------------------------------------------------------------!

  SUBROUTINE max_ongrid(bdr,chg,p)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3),INTENT(INOUT) :: p

    bdr%pnum=1
    bdr%path(bdr%pnum,:)=p
    DO
      CALL step_ongrid(chg,p)
      ! if we didn't move, we're at a maximum
      IF (ALL(p == bdr%path(bdr%pnum,:))) EXIT
      ! otherwise, add point to path
      IF (bdr%pnum >= bdr%pdim) THEN
        CALL reallocate_path(bdr,bdr%pdim*2)
      ENDIF
      bdr%pnum=bdr%pnum+1
      CALL pbc(p,chg%npts)
      bdr%path(bdr%pnum,:)=p
      IF(bdr%volnum(p(1),p(2),p(3)) /= 0) EXIT
    END DO

  RETURN
  END SUBROUTINE max_ongrid

!-----------------------------------------------------------------------------------!
!  step_ongrid:  Do a single iteration of a maximization on the charge density 
!    grid from the point (px,py,pz).  Return a logical indicating if the current
!    point is a charge density maximum.
!-----------------------------------------------------------------------------------!

  SUBROUTINE step_ongrid(chg,p)

    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3),INTENT(INOUT) :: p

    REAL(q2) :: rho_max,rho_tmp,rho_ctr
    INTEGER,DIMENSION(3) :: pt,pm
    INTEGER :: d1,d2,d3

    rho_max=0.0_q2
    pm=p
    rho_ctr=rho_val(chg,p(1),p(2),p(3))
    DO d1=-1,1
!      pt(1)=p(1)+d1
      DO d2=-1,1
!        pt(2)=p(2)+d2
        DO d3=-1,1
!          pt(3)=p(3)+d3
          ! compare the time for this operation to individual assignments
          pt=p+(/d1,d2,d3/)
          rho_tmp=rho_val(chg,pt(1),pt(2),pt(3))
          rho_tmp=rho_ctr+(rho_tmp-rho_ctr)*chg%lat_i_dist(d1,d2,d3)
          IF (rho_tmp > rho_max) THEN
            rho_max=rho_tmp
            pm=pt
          END IF
        END DO
      END DO
    END DO
    p=pm

  RETURN
  END SUBROUTINE step_ongrid

!-----------------------------------------------------------------------------------!
! max_neargrid:  From the point p do a maximization on the charge density grid and
!   assign the maximum found to the volnum array.
!-----------------------------------------------------------------------------------!

  SUBROUTINE max_neargrid(bdr,chg,p)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3),INTENT(INOUT) :: p

    bdr%pnum=1
    bdr%path(bdr%pnum,:)=p
    DO
      CALL step_neargrid(chg,p)
      ! if we didn't move, we're at a maximum
      IF (ALL(p == bdr%path(bdr%pnum,:))) EXIT
      ! otherwise, add point to path
      IF (bdr%pnum >= bdr%pdim) THEN
        CALL reallocate_path(bdr,bdr%pdim*2)
      ENDIF
      bdr%pnum=bdr%pnum+1
      bdr%path(bdr%pnum,:)=p
      IF(bdr%volnum(p(1),p(2),p(3)) /= 0) EXIT
    END DO

  RETURN
  END SUBROUTINE max_neargrid

!-----------------------------------------------------------------------------------!
!  step_neargrid:  Do a single iteration of a maximization on the charge density 
!    grid from the point (px,py,pz).  Return a logical indicating if the current
!-----------------------------------------------------------------------------------!

  SUBROUTINE step_neargrid(chg,p)

    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3),INTENT(INOUT) :: p
    INTEGER,DIMENSION(3) :: pm
    INTEGER :: d1,d2,d3,i

    REAL(q2),DIMENSION(3) :: gradrl,dr=(/0.0_q2,0.0_q2,0.0_q2/)  
    REAL(q2) :: cx,cy,cz,coeff,cp,cm
    SAVE dr

    !print*,'p initial:',p
 
    gradrl=rho_grad_dir(chg,p)
    cx=gradrl(1)*chg%lat_i_dist(1,0,0)
    cy=gradrl(2)*chg%lat_i_dist(0,1,0)
    cz=gradrl(3)*chg%lat_i_dist(0,0,1)
    coeff=1.0_q2/MAX(ABS(cx),ABS(cy),ABS(cz))
    gradrl=coeff*(/cx,cy,cz/)
    !print*,'gradrl:',gradrl
    pm=p+ANINT(gradrl)
    dr=dr+gradrl-ANINT(gradrl)
    !print*,'dr=',dr
        
!    CALL pbc(pm,chg%npts)
!    IF(is_max(chg,pm)) THEN
!       print*, '    max:', p
!       dr=(/0.0_q2,0.0_q2,0.0_q2/)
!       RETURN
!    END IF

    DO i=1,3
      IF (ABS(dr(i)) > 0.5_q2) THEN
        pm(i)=pm(i)+ANINT(dr(i))
        dr(i)=dr(i)-ANINT(dr(i))
      END IF
    END DO

    CALL pbc(pm,chg%npts)
    cp=rho_val(chg,p(1),p(2),p(3))
    cm=rho_val(chg,pm(1),pm(2),pm(3))
    IF (cm < cp) THEN
       print*, '    max:', p
       dr=(/0.0_q2,0.0_q2,0.0_q2/)
       RETURN
    END IF

    p=pm
    !print*,'p end',p
    !print*,' '

  RETURN
  END SUBROUTINE step_neargrid

!-----------------------------------------------------------------------------------!
! is_ep: return .true. if the grid point is on the edge points of a Bader volume.
!-----------------------------------------------------------------------------------!

  SUBROUTINE is_ep(bdr,chg,p,n)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg

    INTEGER,DIMENSION(3),INTENT(IN) :: p
    INTEGER,DIMENSION(3) ::pt
    INTEGER :: d1,d2,d3,path_volnum
    INTEGER,INTENT(INOUT) :: n

    path_volnum=bdr%volnum(p(1),p(2),p(3))
    DO d1=-1,1
      DO d2=-1,1
        DO d3=-1,1
          pt=p+(/d1,d2,d3/)
          CALL pbc(pt,chg%npts)
          IF(bdr%volnum(pt(1),pt(2),pt(3))/=0.AND.bdr%volnum(pt(1),pt(2),pt(3))/=path_volnum) THEN
             bdr%volnum(p(1),p(2),p(3))=0
          END IF
        END DO
      END DO
    END DO

    IF(bdr%volnum(p(1),p(2),p(3))==0) THEN
             n=n+1
    END IF
   
    RETURN
    END SUBROUTINE is_ep

!-----------------------------------------------------------------------------------!
! checkep: re-check the grid point on the edge points of a Bader volume.
!-----------------------------------------------------------------------------------!

  SUBROUTINE checkep(bdr,chg,opts)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    TYPE(options_obj) :: opts

    INTEGER,DIMENSION(3) :: p
    INTEGER :: n1,n2,n3,path_volnum,i,n=0

!    print*,'begin checkep subroutne'
!    print*,'points in each direction',chg%npts(1),chg%npts(2),chg%npts(3)
    DO n1=1,chg%npts(1)
      DO n2=1,chg%npts(2)
        DO n3=1,chg%npts(3)
          p=(/n1,n2,n3/)
          CALL is_ep(bdr,chg,p,n)
!          print*,'check one point',p
        END DO
      END DO
    END DO
    print*,'have resigned the ep', n

!    print*,'bdr%bnum in the beginning of checkep',bdr%bnum
    DO n1=1,chg%npts(1)
      DO n2=1,chg%npts(2)
        DO n3=1,chg%npts(3)
          p=(/n1,n2,n3/)
          IF(bdr%volnum(p(1),p(2),p(3)) == 0) THEN
            IF(opts%bader_opt == opts%bader_offgrid) THEN
              CALL max_offgrid(bdr,chg,p)
            ELSE IF(opts%bader_opt == opts%bader_ongrid) THEN
              CALL max_ongrid(bdr,chg,p)
            ELSE
              CALL max_neargrid(bdr,chg,p)
            END IF
            path_volnum=bdr%volnum(p(1),p(2),p(3))
            IF (path_volnum == 0) THEN
              IF (bdr%bnum >= bdr%bdim) THEN
                CALL reallocate_volpos(bdr,bdr%bdim*2)
              END IF
              bdr%bnum=bdr%bnum+1
              write(*,*) '   new max',bdr%bnum
              path_volnum=bdr%bnum
              bdr%volpos_lat(bdr%bnum,:) = REAL(p,q2)
            END IF
            DO i=1,bdr%pnum
              bdr%volnum(bdr%path(i,1),bdr%path(i,2),bdr%path(i,3))=path_volnum
            END DO
          END IF
        END DO
      END DO
    END DO
  print*,'finish resign'

  RETURN
  END SUBROUTINE checkep

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
        ! GH: this is not correct for non-orthogonal cells, find
        ! proper min length vector by subtracting lattice vectors
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

          ! Check to see if this is at the edge of an atomic volume
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
            v=REAL((/n1,n2,n3/),q2)
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

!-----------------------------------------------------------------------------------!
! write_volnum: Write out a CHGCAR type file with each entry containing an integer
!    indicating the associated Bader maximum.
!-----------------------------------------------------------------------------------!

  SUBROUTINE write_volnum(bdr,opts,ions,chg)

     TYPE(bader_obj) :: bdr
     TYPE(options_obj) :: opts
     TYPE(ions_obj) :: ions
     TYPE(charge_obj) :: chg

     TYPE(charge_obj) :: tmp
     INTEGER :: nx,ny,nz
     CHARACTER(LEN=128) :: filename

     tmp=chg
     tmp%rho=bdr%volnum
     
     filename='VOLUME_INDEX'
     CALL write_charge(ions,chg,opts,filename)

  RETURN
  END SUBROUTINE write_volnum

!-----------------------------------------------------------------------------------!
! write_all_bader: Write out a CHGCAR type file for each of the Bader volumes found.
!-----------------------------------------------------------------------------------!

  SUBROUTINE write_all_bader(bdr,opts,ions,chg)

    TYPE(bader_obj) :: bdr
    TYPE(options_obj) :: opts
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    TYPE(charge_obj) :: tmp
    INTEGER :: nx,ny,nz,i,atomnum,badercur,tenths_done,t1,t2,cr,count_max
    CHARACTER(LEN=128) :: atomfilename
    
    CALL system_clock(t1,cr,count_max)

    WRITE(*,'(/,2x,A)') 'WRITING BADER VOLUMES'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'

    tmp=chg

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

!-----------------------------------------------------------------------------------!
! write_all_atom: Write out a CHGCAR type file for each atom where all Bader volumes
!              assigned to that atom are added together. This is only done if the 
!              atoms has any 'significant' bader volumes associated with it.
!-----------------------------------------------------------------------------------!

  SUBROUTINE write_all_atom(bdr,opts,ions,chg)

    TYPE(bader_obj) :: bdr
    TYPE(options_obj) :: opts
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    TYPE(charge_obj) :: tmp

    INTEGER :: nx,ny,nz,i,j,b,mab,mib,ik,sc,cc,tenths_done,t1,t2,cr,count_max
    INTEGER,DIMENSION(bdr%nvols) :: rck
    CHARACTER(LEN=128) :: atomfilename

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

!-----------------------------------------------------------------------------------!
! write_sel_bader: Write out a CHGCAR type file for the user specified Bader volumes.
!              Volumes associated with a atom can be read from AtomVolumes.dat
!-----------------------------------------------------------------------------------!

  SUBROUTINE write_sel_bader(bdr,opts,ions,chg)

    TYPE(bader_obj) :: bdr
    TYPE(options_obj) :: opts
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    TYPE(charge_obj) :: tmp
    CHARACTER(LEN=128) :: atomfilename
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

!-----------------------------------------------------------------------------------!
! bader_output: Write out a summary of the bader analysis.
!         AtomVolumes.dat: Stores the 'significant' Bader volumes associated with
!                          each atom.
!         ACF.dat        : Stores the main output to the screen.
!         BCF.dat        : Stores 'significant' Bader volumes, their coordinates and
!                          charge, atom associated and distance to it. 
!-----------------------------------------------------------------------------------!

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
    WRITE(100,'(A,A)') '----------------------------------------------------------',&
    &                  '--------------'

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
    
    WRITE(200,668) '--------------------------------------------------------------',&
    &              '-----------'
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
!  rho_val:  Return the density at the point (p1,p2,p3) taking into account the
!    boundary conditions.  This function is used to address points outside the
!    charge density array without a bunch of if statements at the place the value
!    is needed.
!-----------------------------------------------------------------------------------!

  FUNCTION volnum_val(bdr,chg,p1,p2,p3)

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER,INTENT(IN) :: p1,p2,p3
    INTEGER :: volnum_val
    
    INTEGER,DIMENSION(3) :: p
    INTEGER :: i
    
    p=(/p1,p2,p3/)
    DO i=1,3
      DO
        IF(p(i) >= 1) EXIT
        p(i)=p(i)+chg%npts(i)
      END DO
      DO
        IF(p(i) <= chg%npts(i)) EXIT
        p(i)=p(i)-chg%npts(i)
      END DO
    END DO
    
    volnum_val=bdr%volnum(p(1),p(2),p(3))
    
  RETURN
  END FUNCTION volnum_val

!-----------------------------------------------------------------------------------!
! known_volnum: return number of the associated bader volnum if all surrounding
!    grid points are known to be associated with the same bader volnum
!-----------------------------------------------------------------------------------!

  FUNCTION known_volnum(bdr,chg,p)
   
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    REAL(q2),DIMENSION(3),INTENT(IN) :: p
    INTEGER :: known_volnum
   
    INTEGER :: volnum,d1,d2,d3,p1,p2,p3
    LOGICAL :: first_flag

    known_volnum=0
    first_flag=.TRUE.

    p1=FLOOR(p(1))
    p2=FLOOR(p(2))
    p3=FLOOR(p(3))

    DO d1=0,1
      p1=p(1)+d1
      DO d2=0,1
        p2=p(2)+d2
        DO d3=0,1
          p3=p(3)+d3
          IF(first_flag) THEN
            volnum=volnum_val(bdr,chg,p1,p2,p3)
            IF(volnum <= 0) RETURN
            first_flag=.FALSE.
          ELSE
            IF(volnum /= volnum_val(bdr,chg,p1,p2,p3)) RETURN
          END IF
        END DO
      END DO
    END DO
    known_volnum=volnum

  RETURN
  END FUNCTION known_volnum

!-----------------------------------------------------------------------------------!
! reallocate_volpos: 
!-----------------------------------------------------------------------------------!

  SUBROUTINE reallocate_volpos(bdr,newsize)

    TYPE(bader_obj) :: bdr
    INTEGER :: newsize

    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: tmpvolpos

    IF(newsize<bdr%pnum) write(*,*) 'Error: new volpos length too small'

    ALLOCATE(tmpvolpos(bdr%bdim,3))
    tmpvolpos=bdr%volpos_lat
    DEALLOCATE(bdr%volpos_lat)
    bdr%bdim=newsize
    ALLOCATE(bdr%volpos_lat(bdr%bdim,3))
    bdr%volpos_lat(1:bdr%bnum,:) = tmpvolpos(1:bdr%bnum,:)
    DEALLOCATE(tmpvolpos)

  END SUBROUTINE reallocate_volpos

!-----------------------------------------------------------------------------------!
! reallocate_path: 
!-----------------------------------------------------------------------------------!

  SUBROUTINE reallocate_path(bdr,newsize)

    TYPE(bader_obj) :: bdr
    INTEGER :: newsize

    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: tmppath

    IF(newsize<bdr%pnum) write(*,*) 'Error: new path length too small'

    ALLOCATE(tmppath(bdr%pdim,3))
    tmppath=bdr%path
    DEALLOCATE(bdr%path)
    bdr%pdim=newsize
    ALLOCATE(bdr%path(bdr%pdim,3))
    bdr%path(1:bdr%pnum,:) = tmppath(1:bdr%pnum,:)
    DEALLOCATE(tmppath)

  END SUBROUTINE reallocate_path

!-----------------------------------------------------------------------------------!

END MODULE bader_mod
