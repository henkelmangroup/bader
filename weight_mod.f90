  !-----------------------------------------------------------------------------------!
  ! Bader charge density analysis program
  !  Module implementing the weight method by Yu and Trinkle [JCP 134, 064111 (2011)]
  !-----------------------------------------------------------------------------------!

  MODULE weight_mod

    USE kind_mod
    USE matrix_mod
    USE bader_mod
    USE charge_mod 
    USE options_mod
    USE ions_mod
    USE io_mod
    USE chgcar_mod
    IMPLICIT NONE
    PRIVATE

    TYPE weight_obj
      REAL(q2) :: rho
      INTEGER,DIMENSION(3) :: pos
      INTEGER :: volnum !  similar to bader volnum
      REAL(q2),DIMENSION(:),ALLOCATABLE :: volwgt ! weight to the bader volume
      LOGICAL :: isInterior
    END TYPE

    PUBLIC :: weight_obj
    PUBLIC :: bader_weight_calc,quick_sort

  CONTAINS 

  SUBROUTINE bader_weight_calc(bdr,ions,chgval,opts)

    TYPE(weight_obj),ALLOCATABlE,DIMENSION(:) :: chgList, sortedList
    TYPE(weight_obj) :: tempwobj
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chgval,chgtemp
    TYPE(weight_obj) :: wt
    TYPE(ions_obj) :: ions,ionstemp
    TYPE(options_obj) :: opts
    INTEGER(KIND=8) :: totalLength,temp
    INTEGER :: i, j, k, l, n1, n2, n3, walker
    REAL(q2) :: tempflx,totalE,tempw
    LOGICAL :: nbb ! never boundry before
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: indList
    nbb=.TRUE.
    tempw=0
    totalE=0
    bdr%nvols=0
    IF (opts%ref_flag) THEN
      CALL read_charge_ref(ionstemp,chgtemp,opts)
    ELSE
      chgtemp = chgval
    END IF
    totalLength=chgtemp%npts(1)*chgtemp%npts(2)*chgtemp%npts(3)
    ALLOCATE (sortedList(totalLength)) ! sortedList should have the correct size
    PRINT *,'assigned totalLength', totalLength
!    temp=totalLength
    ALLOCATE (chgList(totalLength))
    ALLOCATE (bdr%volnum(chgtemp%npts(1),chgtemp%npts(2),chgtemp%npts(3)))
    ALLOCATE (indlist(chgtemp%npts(1),chgtemp%npts(2),chgtemp%npts(3)))
    walker=1  
    DO n1=1, chgtemp%npts(1)
      DO n2=1, chgtemp%npts(2)
        DO n3=1,chgtemp%npts(3)
          chgList(walker)%rho=chgtemp%rho(n1,n2,n3)
          chgList(walker)%pos=(/n1,n2,n3/)
          chgList(walker)%isInterior=.TRUE.
          chgList(walker)%volnum=0
          walker=walker+1
        END DO
      END DO
    END DO
    PRINT *, 'sorting...'
    CALL quick_sort(chgList)
    ! this is a ascending loop. so future loop starts from the end
    ! or flipping it should be faster
    PRINT *,'DONE, recording indicies'
    DO walker=1,totalLength/2
      tempwobj=chgList(totalLength+1-walker)
      chgList(totalLength+1-walker)=chgList(walker)
      chgList(walker)=tempwobj
    END DO
    PRINT *,'DONE.'
    ! firsst loop,deal with all interior points.
    PRINT *,'looking through for interior points'
    DO walker=1,totalLength
      indList(chgList(walker)%pos(1), &
              chgList(walker)%pos(2), &
              chgList(walker)%pos(3) &
             )=walker
      CALL interiors(bdr,chgtemp,chgList(walker))
    END DO  
    ALLOCATE (bdr%volchg(bdr%nvols))
    DO walker=1,bdr%nvols
      bdr%volchg(walker)=0
    END DO
    !second loop
    PRINT *,'DONE'
    PRINT *,'Running through all boundry points'
    DO walker=1,totalLength
      CALL boundries(bdr,chgtemp,chgList,walker,ions,indList)
    END DO
    PRINT *,'DONE'
    bdr%volchg = bdr%volchg/REAL(chgval%nrho,q2)
    PRINT *,bdr%volchg
    totalE=0
    DO n1=1,bdr%nvols
      totalE=totalE+bdr%volchg(n1)
    END DO
    PRINT *, 'number of electrons:', totalE
  END SUBROUTINE bader_weight_calc

  !-----------------------------------------------------------------------------------!
  ! interiors: this subroutine loops through all grid points, finding out which
  ! are interiors and which are boundries, and number of basins.
  !-----------------------------------------------------------------------------------!
  SUBROUTINE interiors(bdr,chgtemp,cLW) ! chgListWalker
    TYPE(weight_obj) :: cLW
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chgtemp,chgval
    INTEGER :: length,i,n1,n2,n3,temp,tempnvol,tempvolnum,temp2,j
    REAL(q2) :: temprho,denom,rho,areasum,tempw
    INTEGER :: neivolnum,nbd,nin,hdn,temppos,ots ! one to six
    INTEGER,DIMENSION(3) :: xyz,p
    LOGICAL :: ismax
    ismax=.TRUE.
    neivolnum=0
    DO n1=-1,1
      DO n2=-1,1
        DO n3=-1,1
          p=cLW%pos+(/n1,n2,n3/)
          CALL pbc(p,chgtemp%npts)
          IF (cLW%rho < chgtemp%rho(p(1),p(2),p(3))) THEN
            ismax=.FALSE.
            IF (bdr%volnum(p(1),p(2),p(3))==0) THEN 
              ! that neighbor must be a boundary point
              cLW%isInterior=.FALSE.
              cLW%volnum=0
              EXIT
              ! this subroutine can now be stopped
            ELSEIF (neivolnum==0) THEN 
              neivolnum=bdr%volnum(p(1),p(2),p(3))
              cLW%volnum=bdr%volnum(p(1),p(2),p(3))
              bdr%volnum(cLW%pos(1),cLW%pos(2),cLW%pos(3))=bdr%volnum(p(1),p(2),p(3))
              ! if this point turnout to be a boundary point, this will be
              ! erased. 
              CYCLE
            ELSEIF (neivolnum/=bdr%volnum(p(1),p(2),p(3))) THEN
              cLW%isInterior=.FALSE.
              cLW%volnum=0
              EXIT
            END IF
          END IF        
        END DO
      END DO
    END DO
    ! by the end if it gets no basin assignment, it is boundary.
    IF (ismax .AND. cLW%isInterior) THEN
        bdr%nvols=bdr%nvols+1
        cLW%volnum=bdr%nvols
        bdr%volnum(cLW%pos(1),cLW%pos(2),cLW%pos(3))=bdr%nvols
    ELSEIF (neivolnum==0) THEN
      cLW%isInterior=.FALSE.
    END IF
    IF (.NOT. cLW%isInterior) THEN
      cLW%volnum=bdr%volnum(p(1),p(2),p(3))
      bdr%volnum(cLW%pos(1),cLW%pos(2),cLW%pos(3))=bdr%volnum(p(1),p(2),p(3))
    END IF      
  END SUBROUTINE interiors

  ! This subroutine goes through all grid points, adding up rhos of interior
  ! points, calculating flx and weight for boundry points and then adding their
  ! rhos up
  SUBROUTINE boundries(bdr,chgtemp,chgList,walker,ions,indList)
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj):: chgtemp
    TYPE(ions_obj):: ions
    INTEGER         :: walker,i
    INTEGER,DIMENSION(:,:,:) :: indList
!    LOGICAL         :: ismax,isbon
!    REAL(q2)        :: rho,temp
!    INTEGER,DIMENSION(3) :: xyz,p,add
    TYPE(weight_obj),DIMENSION(:) :: chgList
!    REAL(q2),DIMENSION(3,3) :: lat
!      lat(1)=ions%lattice(1,:)
!      lat(2)=ions%lattice(2,:)
!      lat(3)=ions%lattice(3,:)
    IF (chgList(walker)%isInterior) THEN
      bdr%volchg(chgList(walker)%volnum)=bdr%volchg(chgList(walker)%volnum)+chgList(walker)%rho
    ELSE
      CALL flux_weight(bdr,chgtemp,chgList,walker,ions,indList)
    END IF
!    nvols=1
!    DO i=1, SIZE(sortedList)
!      xyz(1)=sortedList(i)%pos(1)
!      xyz(2)=sortedList(i)%pos(2)
!      xyz(3)=sortedList(i)%pos(3)
!      ! First assign weights to maximas. This has to be done to both sorted list
!      ! and chglist
!      rho=sortedList(i)%rho
!      ismax=.TRUE.
!!      ! diagonal neighbors needs to be checked as well
!!      p=xyz+(/1,0,0/)
!!      CALL pbc(p,chgtemp%npts)
!!      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
!!        ismax=.FALSE.
!!        ots=1
!!        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
!!        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
!!      END IF
!!
!!      p=xyz+(/-1,0,0/)
!!      CALL pbc(p,chgtemp%npts)
!!      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
!!        ismax=.FALSE.
!!        ots=2
!!        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
!!        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
!!      END IF
!!
!!      p=xyz+(/0,1,0/)
!!      CALL pbc(p,chgtemp%npts)
!!      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
!!        ismax=.FALSE.
!!        ots=3
!!        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
!!        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
!!      END IF
!!
!!      p=xyz+(/0,-1,0/)
!!      CALL pbc(p,chgtemp%npts)
!!      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
!!        ismax=.FALSE.
!!        ots=4
!!        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
!!        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
!!      END IF
!!
!!      p=xyz+(/0,0,1/)
!!      CALL pbc(p,chgtemp%npts)
!!      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
!!        ismax=.FALSE.
!!        ots=5
!!        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
!!        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
!!      END IF
!!
!!      p=xyz+(/0,0,-1/)
!!      CALL pbc(p,chgtemp%npts)
!!      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
!!        ismax=.FALSE.
!!        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
!!        ots=6
!!        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
!!      END IF
!      ! diagonla neighbors needs to be treated separately. They do not get flux
!      ! as there is no area between them. Dallas treat them as interior points.
!      ! I can copy weight and basin info of the higher neighbor to this point.
!      ! a loop that takes care of all neighbors
!!      p=xyz+(1,1,0)
!!      CALL pbc(p,chgtemp%npts)
!!      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
!!        ismax=.FALSE.
!!        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
!!        ots=6
!!        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
!!      END IF
!!      p=xyz+(1,-1,0)
!!      CALL pbc(p,chgtemp%npts)
!!      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
!!        ismax=.FALSE.
!!        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
!!        ots=6
!!        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
!!      END IF
!!      p=xyz+(-1,1,0,)
!!      p=xyz+(-1,-1,0)
!!      p=xyz+(1,0,1)
!!      p=xyz+(1,0,-1)
!!      p=xyz+(-1,0,1)
!!      p=xyz+(-1,0,-1)
!!      p=xyz+(0,1,1)
!!      p=xyz+(0,1,-1)
!!      p=xyz+(0,-1,1)
!!      p=xyz+(0,-1,-1)
!!      p=xyz+(1,1,1)
!!      p=xyz+(-1,1,1)
!!      p=xyz+(1,-1,1)
!!      p=xyz+(1,1,-1)
!!      p=xyz+(-1,-1,1)
!!      p=xyz+(1,-1,-1)
!!      p=xyz+(-1,1,-1)
!!      p=xyz+(-1,-1,-1)
!
!      ! below is the attempt of a better version, not yet successful.
!      ots=0
!!      IF (i==20) PRINT *,'point 20 with xyz',xyz 
!!      IF (sortedList(i)%pos(1)==114 .AND. &
!!          sortedList(i)%pos(2)==113 .AND. &
!!          sortedList(i)%pos(3)==38) PRINT *,'i of 114 113 38 is',i
!      isbon=.TRUE.  
!      tempvolnum=0
!      DO n1=-1,1
!        DO n2=-1,1
!          DO n3=-1,1
!            add=(/n1,n2,n3/)
!            p=xyz+add ! no problem in this
!!            IF (i==66) PRINT *,'at this point 66 to 2 has wgt', & 
!!                chgList(2581839)%volwgt(2)
!            IF (add(1)**2+add(2)**2+add(3)**2 == 0) THEN
!              ! self, pass
!              CYCLE
!            ELSE IF (add(1)**2+add(2)**2+add(3)**2 == 1) THEN
!            ! direct neighbor, calculate flux and weight
!              CALL pbc(p,chgtemp%npts)
!              ots=ots+1
!              IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
!                ismax=.FALSE.
!!                IF (i==20) PRINT *,'looking at direct p',p
!!                IF (i==20) PRINT *,'found p has higher rho'
!                cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3) 
!                CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
!              END IF
!            ELSE 
!            ! indirect neighbor, only copy info if higher in density
!              CALL pbc(p,chgtemp%npts)
!            ! also only check for indirect neighbors only when there is no known
!            ! direct neighbor with higher densities.
!              IF (rho<chgtemp%rho(p(1),p(2),p(3)) .AND. ismax) THEN
!                ismax=.FALSE.
!!                IF (i==20) PRINT *, 'looking at diagonal p',p
!!                IF (i==20) PRINT *,'found p has higher rho'
!                ! copy all infos from this basin to next
!                cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
!                sind=(sortedList(i)%pos(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+ &
!                     (sortedList(i)%pos(2)-1)*chgtemp%npts(3)+sortedList(i)%pos(3)
!                DO n4=1,bdr%nvols
!!                  sortedList(i)%volwgt(n4)=chgList(cin)%volwgt(n4)
!!                  IF (i==66 .AND. n4==2) THEN
!!                    PRINT *,'getting wgt from a neighbor p and cin',p,cin
!!                    PRINT *,'this neighbor wgt is', chgList(cin)%volwgt(n4)
!!                    PRINT *,'it substract 1 is',chgList(cin)%volwgt(n4)-1
!!                    PRINT *,'this value is 1.', chgList(cin)%volwgt(n4)==1.
!!                  END IF
!                  chgList(sind)%volwgt(n4)=chgList(cin)%volwgt(n4)
!!                  IF (i==8) THEN
!!                    PRINT *,'8 got wgt to',n4,'from cin',cin,'of &
!!                                value',chgList(cin)%volwgt(n4)
!!                  END IF
!                END DO
!              END IF
!            END IF            
!          END DO
!        END DO
!      END DO
!!      END IF
!      IF (ismax==.TRUE.) THEN
!        !print *, 'i is a maximum',i
!        cin=(xyz(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(xyz(2)-1)*chgtemp%npts(3)+xyz(3)
!        !PRINT *,'cin is',cin
!!        sortedList(i)%volwgt(nvols)=1
!        chgList(cin)%volwgt(nvols)=1.
!        nvols=nvols+1
!      END IF
!    END DO
!
  END SUBROUTINE boundries
!
!-------------------------------------------------------------------!
! flux_weight: It takes the current position and a target position as in put,
! adds the target's weight times flux to target to current position's weight.
! This subroutine needs to be called everytime a neighbor with greater density
! is called. 
!-------------------------------------------------------------------!

  SUBROUTINE flux_weight(bdr,chgtemp,chgList,walker,ions,indList)
    TYPE(weight_obj),DIMENSION(:) :: chgList
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) ::chgtemp
    TYPE(ions_obj) :: ions
    INTEGER :: n1,n2,n3,n4,walker,nein,counter
    INTEGER,DIMENSION(:,:,:) :: indList
    INTEGER,DIMENSION(3) :: p
    REAL(q2) :: denom,area,length,temp
    REAL(q2),DIMENSION(6) :: flx,nom
    denom=0
    counter=0
    temp=0
    ALLOCATE(chgList(walker)%volwgt(bdr%nvols))
    DO n1=1,bdr%nvols
      chgList(walker)%volwgt(n1)=0
    END DO
    DO n1=-1,1
      DO n2=-1,1
        DO n3=-1,1
          ! only look at direct neighbors
          ! this is first time looping over neighbors. this loop only finds
          ! denom and noms
          IF (n1**2+n2**2+n3**2 /=1) CYCLE
          counter=counter+1
          nom(counter)=0
          flx(counter)=0
          p=chgList(walker)%pos+(/n1,n2,n3/)
          CALL pbc(p,chgtemp%npts)
          IF (chgtemp%rho(p(1),p(2),p(3))>chgList(walker)%rho) THEN
            !first calculate nom and denom
            ! need to first figure out which plane this is
            IF (n3/=0) THEN
              CALL find_area(ions%lattice(1,:),ions%lattice(2,:),area)
              length=SQRT(ions%lattice(3,1)**2+ions%lattice(3,2)**2+ions%lattice(3,3)**2)
            ELSEIF (n2/=0) THEN
              CALL find_area(ions%lattice(1,:),ions%lattice(3,:),area)
              length=SQRT(ions%lattice(2,1)**2+ions%lattice(2,2)**2+ions%lattice(2,3)**2)
            ELSEIF (n1/=0) THEN
              CALL find_area(ions%lattice(2,:),ions%lattice(3,:),area)
              length=SQRT(ions%lattice(1,1)**2+ions%lattice(1,2)**2+ions%lattice(1,3)**2)
            END IF
            nom(counter)=ABS(area)*length*(chgtemp%rho(p(1),p(2),p(3))-chgList(walker)%rho)
            denom=denom+nom(counter)
          END IF 
        END DO
      END DO
    END DO
!    IF (WALKER==456497) THEN
!      PRINT *,'456497 denom',denom
!      PRINT *,'456497 nom',nom
!    END IF
    counter=0
    DO n1=-1,1
      DO n2=-1,1
        DO n3=-1,1
          IF (n1**2+n2**2+n3**2 /=1) CYCLE
          counter=counter+1
          p=chgList(walker)%pos+(/n1,n2,n3/)
          CALL pbc(p,chgtemp%npts)

          IF (chgtemp%rho(p(1),p(2),p(3))>chgList(walker)%rho) THEN
            IF (walker==456497) THEN
              PRINT *,'read in nom(conter)',nom(counter)
              PRINT *,'read in denom',denom
            END IF
            flx(counter)=nom(counter)/denom
            IF (walker==456497) THEN
              
              PRINT *,'calculated flx',flx(counter)
              PRINT *,'neighbor has volnum',bdr%volnum(p(1),p(2),p(3))
            END IF
            ! interior neighbor and boundry neighbor needs to be treated
            ! separately
            
            IF (bdr%volnum(p(1),p(2),p(3))/=0) THEN !interior neighbor
              chgList(walker)%volwgt(bdr%volnum(p(1),p(2),p(3)))= & 
                chgList(walker)%volwgt(bdr%volnum(p(1),p(2),p(3)))+ flx(counter)
            
            ELSE
              ! need to find out that point's position in chgList

              nein=indList(p(1),p(2),p(3))
              IF (WALKER==456497) THEN
                PRINT *,'neighbor nein is',nein
                PRINT *,'its wgt is',chgList(nein)%volwgt
              END IF
              DO n4=1,bdr%nvols
                chgList(walker)%volwgt(n4)=chgList(walker)%volwgt(n4)+ &
                  chgList(nein)%volwgt(n4)*flx(counter)
              END DO
            END IF
          END IF          
        END DO
      END DO
    END DO
    IF (walker==456497) THEN
      PRINT *,'456497 FLX',flx
      PRINT *,'WGT',chgList(456497)%volwgt
    END IF
    temp=0
    DO n1=1,bdr%nvols
      bdr%volchg(n1)=bdr%volchg(n1)+chgList(walker)%rho*chgList(walker)%volwgt(n1)
      temp=temp+chgList(walker)%volwgt(n1)
    END DO
!    IF ( temp /=1.) THEN
!      PRINT *,'walker',walker
!      PRINT *,'have total wgt',temp
!    END IF

            ! is that point a boundary point or interior point?
            ! an interior point has bdr%volnum
!            IF (bdr%volnum(p(1),p(2),p(3))/=0)


!    LOGICAL :: templogi
!    REAL(q2) :: temp
!      !sind should be sorted list index to CHGLIST index
!      sind=(sortedList(i)%pos(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(sortedList(i)%pos(2)-1)*chgtemp%npts(3)&
!        +sortedList(i)%pos(3)
!!      templogi=.TRUE.
!      DO n1=1,bdr%nvols
!        sortedList(i)%volwgt(n1)=sortedList(i)%volwgt(n1)+sortedList(i)%flx(ots)*chgList(cin)%volwgt(n1)
!        chgList(sind)%volwgt(n1)=sortedList(i)%volwgt(n1)   
!        x=sortedList(i)%pos(1)
!        y=sortedList(i)%pos(2)
!        z=sortedList(i)%pos(3)
!!        IF (i==8) PRINT *,'preexisting value is', chgList(sind)%volwgt(n1)
!!        IF (1.-chgList(sind)%volwgt(n1)<=10**(-10)) THEN
!!          CYCLE
!!        ELSE
!!          temp=chgList(sind)%volwgt(n1)+sortedList(i)%flx(ots)*chgList(cin)%volwgt(n1)
!!          chgList(sind)%volwgt(n1)=temp !chgList(sind)%volwgt(n1) +sortedList(i)%flx(ots)*chgList(cin)%volwgt(n1)
!!          IF (i==8)THEN
!!            PRINT *,'8 got weight to ',n1,'from cin',cin, &
!!                       'of value',chgList(sind)%volwgt(n1)
!!            PRINT *,'flx is',sortedList(i)%flx(ots)
!!            PRINT *,'the neighbor has weight',chgList(cin)%volwgt(n1)
!!          END IF
!!        END IF
!!        IF (chgList(sind)%volwgt(n1)>0 .AND. chgList(sind)%volwgt(n1)/=1.) THEN
!!        IF (i<=100 .AND. chgList(sind)%volwgt(n1)>1.) THEN
!!        IF (i==66) THEN
!!          PRINT *,'abnormal weight at point i,xyz,sind',i,x,y,z,sind
!!          PRINT *,'to basin n1',n1
!!          PRINT *,'with value',chgList(sind)%volwgt(n1)
!!          PRINT *,'this is point i',i,'. looking at neighbor cin',cin
!!          PRINT *,'neighbor has volwgt to basin',n1,' of value', &
!!                   chgList(cin)%volwgt(n1)
!!          PRINT *,''
!!        END IF
!!        IF (chgList(sind)%volwgt(n1)/=1. .AND. chgList(sind)%volwgt(n1)/=0. .AND. i<=20 ) THEN 
!!          PRINT *,'point',i,'doesn not have weight 1 or 0'
!!          PRINT *,'IT IS',chgList(sind)%volwgt(n1)
!!        END IF
!      END DO
  END SUBROUTINE flux_weight


  ! this subroutine finds a points index given its xyz position. 
  ! to make this faster, we know that the neighbor we are looking at is a
  ! boundry point. we can start from the first boundry point. 
  SUBROUTINE find_ind(p,chgList,nein,fbp)
    TYPE(weight_obj),DIMENSION(:) :: chgList
    INTEGER,DIMENSION(3) :: p
    INTEGER :: nein,walker,fbp
    DO walker=fbp,1000000
      IF (chgList(walker)%pos(1)==p(1) .AND. &
          chgList(walker)%pos(2)==p(2) .AND. &
          chgList(walker)%pos(3)==p(3) ) THEN
        nein=walker
        EXIT
      END IF
    END DO
  END SUBROUTINE find_ind
!
!
!
!  ! This subroutine gets flux as well
!
!  SUBROUTINE get_nvols (sortedList,bdr,chgtemp,ions)
!    TYPE(weight_obj),DIMENSION(:) :: sortedList
!    TYPE(bader_obj) :: bdr
!    TYPE(charge_obj) :: chgtemp
!    TYPE(ions_obj) :: ions
!    INTEGER :: i,j,n1,n2,n3
!    INTEGER,DIMENSION(3) :: xyz,p,p1,p2,p3,p4,p5,p6
!    REAL(q2),DIMENSION(3) :: v1,v2,v3
!    REAL(q2) :: area1,area2,area3,rho,denom,temp
!    REAL(q2),DIMENSION(6) :: nom
!    LOGICAL :: ismax,bein
!    DO i=1,SIZE(sortedList)
!      ismax=.TRUE.
!      xyz(1)=sortedList(i)%pos(1)
!      xyz(2)=sortedList(i)%pos(2)
!      xyz(3)=sortedList(i)%pos(3)
!      p1=xyz+(/1,0,0/)
!      p2=xyz+(/-1,0,0/)
!      p3=xyz+(/0,1,0/)
!      p4=xyz+(/0,-1,0/)
!      p5=xyz+(/0,0,1/)
!      p6=xyz+(/0,0,-1/)
!      CALL pbc(p1,chgtemp%npts)
!      CALL pbc(p2,chgtemp%npts)
!      CALL pbc(p3,chgtemp%npts)
!      CALL pbc(p4,chgtemp%npts)
!      CALL pbc(p5,chgtemp%npts)
!      CALL pbc(p6,chgtemp%npts)
!      v1=ions%lattice(1,:)
!      v2=ions%lattice(2,:)
!      v3=ions%lattice(3,:)
!      CALL find_area(v2,v3,area1)
!      CALL find_area(v1,v3,area2)
!      CALL find_area(v1,v2,area3)
!      area1=ABS(area1) ! for x directions
!      area2=ABS(area2) ! for y directions
!      area3=ABS(area3) ! for z directions
!      ! in assign weight, the access to these fluxes are in order:
!      ! -1,0,0 | 0,-1,0 | 0,0,-1
!      !  0,0,1 | 0, 1,0 | 1,0, 0
!      ! that is :
!      ! 1,2,3,3,2,1
!      ! or p2,p4,p6,p5,p3,p1
!      denom=0
!      nom=(/0,0,0,0,0,0/)
!      rho=sortedList(i)%rho
!      IF (rho<chgtemp%rho(p1(1),p1(2),p1(3))) THEN
!        ismax=.FALSE.
!        nom(1)=SQRT(v1(1)**2+v1(2)**2+v1(3)**2)* &
!         &(chgtemp%rho(p1(1),p1(2),p1(3))-rho)*area1
!        denom=denom+SQRT(v1(1)**2+v1(2)**2+v1(3)**2)* &
!         &(chgtemp%rho(p1(1),p1(2),p1(3))-rho)*area1
!      END IF
!
!      IF (rho<chgtemp%rho(p2(1),p2(2),p2(3))) THEN
!        ismax=.FALSE.
!        nom(2)=SQRT(v1(1)**2+v1(2)**2+v1(3)**2)* &
!         &(chgtemp%rho(p2(1),p2(2),p2(3))-rho)*area1
!        denom=denom+SQRT(v1(1)**2+v1(2)**2+v1(3)**2)* &
!         &(chgtemp%rho(p2(1),p2(2),p2(3))-rho)*area1
!      END IF
!
!      IF (rho<chgtemp%rho(p3(1),p3(2),p3(3))) THEN
!        ismax=.FALSE.
!        nom(3)=SQRT(v2(1)**2+v2(2)**2+v2(3)**2)* &
!         &(chgtemp%rho(p3(1),p3(2),p3(3))-rho)*area2
!        denom=denom+SQRT(v2(1)**2+v2(2)**2+v2(3)**2)* &
!         &(chgtemp%rho(p3(1),p3(2),p3(3))-rho)*area2
!      END IF
!
!      IF (rho<chgtemp%rho(p4(1),p4(2),p4(3))) THEN
!        nom(4)=SQRT(v2(1)**2+v2(2)**2+v2(3)**2)* &
!         &(chgtemp%rho(p4(1),p4(2),p4(3))-rho)*area2
!        ismax=.FALSE.
!        denom=denom+SQRT(v2(1)**2+v2(2)**2+v2(3)**2)* &
!         &(chgtemp%rho(p4(1),p4(2),p4(3))-rho)*area2
!      END IF
!
!      IF (rho<chgtemp%rho(p5(1),p5(2),p5(3))) THEN
!        ismax=.FALSE.
!        nom(5)=SQRT(v3(1)**2+v3(2)**2+v3(3)**2)* &
!         &(chgtemp%rho(p5(1),p5(2),p5(3))-rho)*area3
!        denom=denom+SQRT(v3(1)**2+v3(2)**2+v3(3)**2)* &
!         &(chgtemp%rho(p5(1),p5(2),p5(3))-rho)*area3
!      END IF
!
!      IF (rho<chgtemp%rho(p6(1),p6(2),p6(3))) THEN
!        ismax=.FALSE.
!        nom(6)=SQRT(v3(1)**2+v3(2)**2+v3(3)**2)* &
!         &(chgtemp%rho(p6(1),p6(2),p6(3))-rho)*area3
!        denom=denom+SQRT(v3(1)**2+v3(2)**2+v3(3)**2)* &
!         &(chgtemp%rho(p6(1),p6(2),p6(3))-rho)*area3
!      END IF
!
!      DO n1=-1,1
!        DO n2=-1,1
!          DO n3=-1,1
!            p=xyz+(/n1,n2,n3/)
!            CALL pbc(p,chgtemp%npts)
!            IF (n1**2+n2**2+n3**2 > 1 .AND. &
!                rho<chgtemp%rho(p(1),p(2),p(3))) THEN
!              ismax=.FALSE.
!              ! no nom and denom calc as there is no flx
!            END IF
!          END DO
!        END DO
!      END DO
!
!      ! can't know if a point is interior here because there may be a point with
!      ! flux to 2 points but to the same basin
! 
!      ! storage sequence should be:
!      ! p2,4,6,5,3,1 
!      sortedList(i)%flx(1)=nom(2)/denom
!      sortedList(i)%flx(2)=nom(4)/denom
!      sortedList(i)%flx(3)=nom(6)/denom
!      sortedList(i)%flx(4)=nom(5)/denom
!      sortedList(i)%flx(5)=nom(3)/denom
!      sortedList(i)%flx(6)=nom(1)/denom
!!      IF (i==8) THEN
!!        PRINT *,'denom is',denom
!!        temp=0
!!        DO n1=1,6
!!          PRINT *,'nom',n1,'is',nom(n1)
!!          PRINT *,'flx',n1,' is',sortedList(i)%flx(n1)
!!          temp=temp+sortedList(i)%flx(n1)
!!        END DO
!!        PRINT *,'temp is 1',temp==1.
!!      END IF
!      bein=.FALSE.
!      IF (ismax) THEN
!        ! debugging block
!        PRINT *,xyz,'is a maximum'
!        !PRINT *, 'and the rho is', rho
!        ! end block
!        bdr%nvols=bdr%nvols+1
!      END IF
!    END DO
!  END SUBROUTINE get_nvols
  !-----------------------------------------------------------------------------------!
  !Quick_sort: Modified from the original program (info given below)
  !-----------------------------------------------------------------------------------!



  RECURSIVE SUBROUTINE quick_sort(list)
  
    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified by Alan Miller to include an associated integer array which gives
    ! the positions of the elements in the original order.
    
    IMPLICIT NONE
    TYPE(weight_obj), DIMENSION (:), INTENT(IN OUT)  :: list
    ! Local variable
    INTEGER(KIND=8) :: i,length
    
    length=SIZE(list)
    CALL quick_sort_1(1, length)
    
    CONTAINS
    
      RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)
      
        INTEGER(KIND=8), INTENT(IN) :: left_end, right_end
        
        !     Local variables
        INTEGER(KIND=8)             :: i, j, itemp
        TYPE(weight_obj) :: temp
        REAL(q2)                :: reference
        INTEGER, PARAMETER  :: max_simple_sort_size = 6
        
        IF (right_end < left_end + max_simple_sort_size) THEN
          ! Use interchange sort for small lists
          CALL interchange_sort(left_end, right_end)
        
        ELSE
          ! Use partition ("quick") sort
          reference = list((left_end + right_end)/2)%rho
          i = left_end - 1; j = right_end + 1
        
          DO
            ! Scan list from left end until element >= reference is found
            DO
              i = i + 1
              IF (list(i)%rho >= reference) EXIT
            END DO
            ! Scan list from right end until element <= reference is found
            DO
              j = j - 1
              IF (list(j)%rho <= reference) EXIT
            END DO
        
        
            IF (i < j) THEN
              ! Swap two out-of-order elements
              temp = list(i); list(i) = list(j); list(j) = temp
            ELSE IF (i == j) THEN
              i = i + 1
              EXIT
            ELSE
              EXIT
            END IF
          END DO
        
          IF (left_end < j) CALL quick_sort_1(left_end, j)
          IF (i < right_end) CALL quick_sort_1(i, right_end)
        END IF
      
      END SUBROUTINE quick_sort_1
    
    
      SUBROUTINE interchange_sort(left_end, right_end)
      
        INTEGER(KIND=8), INTENT(IN) :: left_end, right_end
        
        !     Local variables
        INTEGER(KIND=8)             :: i, j, itemp
        TYPE(weight_obj)    :: temp
        
        DO i = left_end, right_end - 1
          DO j = i+1, right_end
            IF (list(i)%rho > list(j)%rho) THEN
              temp = list(i); list(i) = list(j); list(j) = temp
            END IF
          END DO
        END DO
      
      END SUBROUTINE interchange_sort
    
  END SUBROUTINE quick_sort


  
  END MODULE
