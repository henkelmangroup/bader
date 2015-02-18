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
      REAL(q2) :: rho,weight
      INTEGER,DIMENSION(3) :: pos
      REAL(q2),DIMENSION(6) :: flx
      INTEGER,DIMENSION(:),ALLOCATABLE :: volnum !  similar to bader volnum
      REAL(q2),DIMENSION(:),ALLOCATABLE :: volwgt ! weight to the bader volume
      LOGICAL :: isInterior
    END TYPE

    PUBLIC :: weight_obj
    PUBLIC :: bader_weight_calc,quick_sort

  CONTAINS 

  SUBROUTINE bader_weight_calc(bdr,ions,chgval,opts)

    TYPE(weight_obj),ALLOCATABlE,DIMENSION(:) :: chgList, sortedList,tempList
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chgval,chgtemp
    TYPE(weight_obj) :: wt
    TYPE(ions_obj) :: ions,ionstemp
    TYPE(options_obj) :: opts
    INTEGER(KIND=8) :: totalLength,temp
    INTEGER :: i, j, k, l, n1, n2, n3, walker
    REAL(q2) :: tempflx
    IF (opts%ref_flag) THEN
      CALL read_charge_ref(ionstemp,chgtemp,opts)
    ELSE
      chgtemp = chgval
    END IF
    ALLOCATE(bdr%volnum(chgtemp%npts(1),chgtemp%npts(2),chgtemp%npts(3)))
    ALLOCATE(bdr%nnion(bdr%nvols), bdr%iondist(bdr%nvols),bdr%ionchg(ions%nions))
    totalLength=chgtemp%npts(1)*chgtemp%npts(2)*chgtemp%npts(3)
    ALLOCATE (sortedList(totalLength)) ! sortedList should have the correct size
    !PRINT *,'assigned totalLength', totalLength
    temp=totalLength
    ! The length of the array to be sorted should be of size 2^n. The additional
    ! buffer slots will be empty
    totalLength=2**(CEILING(LOG(temp*1.0_q2)/LOG(2.0_q2)))
    ALLOCATE (chgList(totalLength))
    ALLOCATE (tempList(totalLength))
    walker=1  
    DO n1=1, chgtemp%npts(1)
      DO n2=1, chgtemp%npts(2)
        DO n3=1,chgtemp%npts(3)
          chgList(walker)%rho=chgtemp%rho(n1,n2,n3)
          chgList(walker)%pos(1)=n1
          chgList(walker)%pos(2)=n2
          chgList(walker)%pos(3)=n3
          walker=walker+1
          chgList(walker)%isInterior=.FALSE.
          chgList(walker)%flx=(/0,0,0,0,0,0/)
        END DO
      END DO
    END DO
    ! The empty slots are filled with junk. Lets fill it with -1
    DO i=walker,totalLength
      chgList(walker)%rho=-1
      chgList(walker)%pos(1)=-1
      chgList(walker)%pos(2)=-1
      chgList(walker)%pos(3)=-1
      walker=walker+1  
    END DO
    tempList=chgList
    PRINT *, 'sorting...'
    CALL quick_sort(tempList)
    DO i=1,SIZE(sortedList)
      sortedList(i)=tempList(totalLength-i+1)
    END DO
    CALL weight_calc(sortedList,bdr,chgtemp,chgval,chgList,ions)
    ! It appears that the weight are properly calculated. Now they need to be
    ! added together. 
!    PRINT *,'point 1 1 1 has weight', chgList(1)%volwgt

    DO n1=1,bdr%nvols
      bdr%volchg=0
    END DO
    DO n1=1,chgval%npts(1)
      DO n2=1,chgval%npts(2)
        DO n3=1,chgval%npts(3)
          i=(n1-1)*chgval%npts(2)*chgval%npts(3)+(n2-1)*chgval%npts(3)+n3
          DO walker=1,bdr%nvols
            bdr%volchg(walker)=bdr%volchg(walker)+chgList(i)%volwgt(walker)*chgval%rho(n1,n2,n3)
          END DO
        END DO
      END DO
    END DO

    bdr%volchg = bdr%volchg/REAL(chgval%nrho,q2)
    PRINT *, 'nrho is ',chgval%nrho
    PRINT *,bdr%volchg
    PRINT *, 'checking for anormalies'
    ! nothing obviously wrong with flx
    DO n1=1,SIZE(sortedList)
      tempflx=0
      DO n2=1,6
        tempflx=tempflx+sortedList(n1)%flx(n2)
      END DO
      IF (tempflx>1.0001) THEN
        PRINT *, 'PROBLEM IN ENTRY',n1
        PRINT *, 'TOTAL FLUX IS',tempflx
      END IF   
    END DO
  END SUBROUTINE bader_weight_calc

  !-----------------------------------------------------------------------------------!
  ! weight_calc : calculates the weight of a cell to a certain basin 
  !-----------------------------------------------------------------------------------!
  SUBROUTINE weight_calc(sortedList,bdr,chgtemp,chgval,chgList,ions)
    TYPE(weight_obj),DIMENSION(:) :: sortedList,chgList
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chgtemp,chgval
    TYPE(ions_obj) :: ions
    INTEGER :: length,i,n1,n2,n3,temp,tempnvol,tempvolnum,temp2,j
    REAL(q2) :: temprho,denom,rho,areasum
    REAL(q2),DIMENSION(6) :: nom,area
    INTEGER,DIMENSION(3) :: npos !neighbor position
    INTEGER,DIMENSION(3) :: xyz ! xyz indexes of current point 
    INTEGER :: nbd,nin,hdn,temppos,ots ! one to six
    LOGICAL :: fb
    bdr%nvols=0 
    CALL GET_NVOLS(sortedList,bdr,chgtemp,ions)
    ALLOCATE (bdr%volchg(bdr%nvols))
!    PRINT *, 'nvols is' , bdr%nvols
    length=SIZE(sortedList)
    nbd=0
    nin=0
    temp2=0
    temp=SIZE(sortedList)/100
    PRINT *, 'nvols is',bdr%nvols
    DO i=1, length
      ! now that all flux are calculated. It is time to do the second loop and
      ! add up all weights. 
      
      ALLOCATE(chgList(i)%volnum(bdr%nvols))
      ALLOCATE(chgList(i)%volwgt(bdr%nvols))
      ALLOCATE(sortedList(i)%volnum(bdr%nvols))
      ALLOCATE(sortedList(i)%volwgt(bdr%nvols))
      ! Initialize the arrays, fill them with zeros. 
      DO n1=1,bdr%nvols
        chgList(i)%volnum(n1)=0
        chgList(i)%volwgt(n1)=0
        sortedList(i)%volnum(n1)=0 ! is volnum array really necessary?
        sortedList(i)%volwgt(n1)=0 
      END DO
    END DO
    PRINT *, 'check bounds,size of volwgt',size(chgList(1000)%volwgt)
    PRINT *,'chgval rho has sizes',size(chgval%rho)
    PRINT *, 'chgval rho entry 9 9 9 has value', chgval%rho(8,8,8)
    CALL assign_weight(sortedList,chgList,bdr,chgtemp,chgval)
    PRINT *,'volwgt of 61,31,61', chgList(60*120*120+30*120+61)%volwgt
    PRINT *,'nin is',nin  
     ! I will not add charges up in the way below
!    ALLOCATE(bdr%volchg(bdr%nvols))
!    bdr%volchg = 0._q2
!    DO n1 = 1,chgval%npts(1)
!      DO n2 = 1,chgval%npts(2)
!        DO n3 = 1,chgval%npts(3)
!          i=(n1-1)*chgval%npts(2)*chgval%npts(3)+(n2-1)*chgval%npts(3)+n3
!          PRINT *, n1,n2,n3,'is interior',chgList(i)%isInterior
!          IF (chgList(i)%isInterior==.FALSE.) THEN
!            ! Boundary points have a 6 slot array. Some of them have information
!            ! about where the density is flowing to. Need to look at all of them
!            ! But first lets find flux
!            CYCLE
!          END IF
!
!!          PRINT *, 'bdr%volnum is',bdr%volnum(n1,n2,n3)
!          IF (bdr%volnum(n1,n2,n3) == bdr%nvols+1) CYCLE
!!          PRINT *, 'volnum of ',n1,n2,n3,'is',bdr%volnum(n1,n2,n3)
!!          PRINT *, 'calculated i is', i
!!          PRINT *, 'sortedList(i) has coordinates', chgList(i)%pos
!!          PRINT *, 'is this an interior point?', chgList(i)%isInterior
!
!          bdr%volchg(bdr%volnum(n1,n2,n3)) = &
!          &  bdr%volchg(bdr%volnum(n1,n2,n3)) + chgval%rho(n1,n2,n3)
!        END DO
!      END DO
!    END DO
!    bdr%volchg = bdr%volchg/REAL(chgval%nrho,q2)
!    PRINT *, bdr%volchg
!
    DO i=1,length
      DEALLOCATE(sortedList(i)%volnum)
      DEALLOCATE(sortedList(i)%volwgt)
    END DO
  
  END SUBROUTINE
!
  SUBROUTINE assign_weight(sortedList,chgList,bdr,chgtemp,chgval)
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj):: chgtemp,chgval
    INTEGER         :: i,ti1,ti2,nvols,cin,ots
    LOGICAL         :: ismax
    REAL(q2)        :: rho
    INTEGER,DIMENSION(3) :: xyz,p
    TYPE(weight_obj),DIMENSION(:) :: sortedList,chgList

    nvols=1
    DO i=1,SIZE(sortedList)
      xyz(1)=sortedList(i)%pos(1)
      xyz(2)=sortedList(i)%pos(2)
      xyz(3)=sortedList(i)%pos(3)
      ! First assign weights to maximas. This has to be done to both sorted list
      ! and chglist
      rho=sortedList(i)%rho
      ismax=.TRUE.
      ti1=(xyz(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(xyz(2)-1)*chgtemp%npts(3)+xyz(3)
      p=xyz+(/1,0,0/)
      CALL pbc(p,chgtemp%npts)
      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
        ismax=.FALSE.
        ots=1
        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
      END IF

      p=xyz+(/-1,0,0/)
      CALL pbc(p,chgtemp%npts)
      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
        ismax=.FALSE.
        ots=2
        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
      END IF

      p=xyz+(/0,1,0/)
      CALL pbc(p,chgtemp%npts)
      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
        ismax=.FALSE.
        ots=3
        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
      END IF

      p=xyz+(/0,-1,0/)
      CALL pbc(p,chgtemp%npts)
      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
        ismax=.FALSE.
        ots=4
        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
      END IF

      p=xyz+(/0,0,1/)
      CALL pbc(p,chgtemp%npts)
      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
        ismax=.FALSE.
        ots=5
        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
      END IF

      p=xyz+(/0,0,-1/)
      CALL pbc(p,chgtemp%npts)
      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
        ismax=.FALSE.
        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
        ots=6
        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
      END IF
      IF (ismax==.TRUE.) THEN
!        PRINT *, '---------------- A new basin found ! ------------'
        !give basin assignment and weight. The index is volnum
        !this needs to be done to the point in chgList
        cin=(xyz(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(xyz(2)-1)*chgtemp%npts(3)+xyz(3)
        sortedList(i)%volwgt(nvols)=1
        chgList(cin)%volwgt(nvols)=1
!        PRINT *,'basin 2 initialized with rho',chgtemp%rho(xyz(1),xyz(2),xyz(3))
!        PRINT *,'CHGVAL of the basin is', chgtemp%rho(xyz(1),xyz(2),xyz(3))
!        PRINT *,'index i just got', i,nvols
!        PRINT *, 'the weight is', sortedList(i)%volwgt(nvols)
        nvols=nvols+1
!        PRINT *, 'the entire weight array is',sortedList(i)%volwgt
      END IF
    END DO

  END SUBROUTINE assign_weight

!-------------------------------------------------------------------!
! flux_weight: It takes the current position and a target position as in put,
! adds the target's weight times flux to target to current position's weight.
! This subroutine needs to be called everytime a neighbor with greater density
! is called. 
!-------------------------------------------------------------------!

  SUBROUTINE flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
    TYPE(weight_obj),DIMENSION(:) :: sortedList,chgList
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) ::chgtemp,chgval
    INTEGER :: n1,ots,sind,i,cin,x,y,z
!    IF ( i <= 10 ) THEN
!      PRINT *,'----------------- New flux weight call ----------------' 
!    sind=(sortedList(i)%pos(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(sortedList(i)%pos(2)-1)*chgtemp%npts(3)&
!        +sortedList(i)%pos(3)
!    PRINT *, 'index in sortedList is',i
!    PRINT *, 'of location',sortedList(i)%pos
!    DO n1=1,bdr%nvols
!      PRINT *, '   '
!      PRINT *, 'n1 is',n1
!      PRINT *, 'pre existing weight is', sortedList(i)%volwgt(n1)
!      sortedList(i)%volwgt(n1)=sortedList(i)%volwgt(n1)+sortedList(i)%flx(ots)*chgList(cin)%volwgt(n1)
!      chgList(sind)%volwgt(n1)=sortedList(i)%volwgt(n1)  
!      PRINT *, 'given weight from index in chgList',cin
!      PRINT *, 'of location ', chgList(cin)%pos
!      PRINT *, 'it has weight to ',n1
!      PRINT *,'of value', chgList(cin)%volwgt(n1)
!
!    END DO
!
!    ELSE 
      sind=(sortedList(i)%pos(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(sortedList(i)%pos(2)-1)*chgtemp%npts(3)&
        +sortedList(i)%pos(3)
      DO n1=1,bdr%nvols
        sortedList(i)%volwgt(n1)=sortedList(i)%volwgt(n1)+sortedList(i)%flx(ots)*chgList(cin)%volwgt(n1)
        chgList(sind)%volwgt(n1)=sortedList(i)%volwgt(n1)    
!        PRINT *, 'last entry with sortedList'
        x=sortedList(i)%pos(1)
        y=sortedList(i)%pos(2)
        z=sortedList(i)%pos(3)
!        IF (n1==2 .AND. sortedList(i)%volwgt(2) /=0) THEN
!        PRINT *,'THIS IS POINT',x,y,z
!        PRINT *, 'this point has weight to 2', sortedList(i)%volwgt(2)
!        PRINT *, 'this point has rho', chgval%rho(x,y,z)
!        PRINT *, 'basin 2 volchg has changed to',bdr%volchg(n1)
!        END IF
!        IF (sortedList(i)%volwgt(n1)>1) THEN 
!          PRINT *, 'entry i in sorted List :',i
!          PRINT *, 'has acquired abnotmal weight :',sortedList(i)%volwgt(n1)
!          PRINT *, 'from n1',n1
!          PRINT *, 'the flux matrix is',sortedList(i)%flx
!        END IF
      END DO
!    END IF
  END SUBROUTINE

  SUBROUTINE  neighbors(bdr,npos,xyz,chgtemp,rho, nboundary,ninterior,i &
                ,ions,tempnvol,tempvolnum,hdn,fb,sortedList, & 
                denom,nom,ots,area,areasum )
  TYPE(weight_obj),DIMENSION(:) :: sortedList
  TYPE(bader_obj) :: bdr
  TYPE(charge_obj) :: chgtemp  
  TYPE(ions_obj) :: ions
  INTEGER,DIMENSION(3) :: npos,xyz

  INTEGER :: nboundary,ninterior,ots
  LOGICAL :: fb ! first bigger neighbor
  INTEGER :: temp,tempnvol,tempvolnum,temp2,i,hdn
  REAL(q2) :: rho,denom,areasum
  REAL(q2),DIMENSION(6) :: nom,area
    CALL pbc(npos,chgtemp%npts)
    IF (rho<chgtemp%rho(npos(1),npos(2),npos(3))) THEN
!      If (ots ==1 .OR. ots ==2) THEN
        ! a step taken in the x direction
!        CALL FIND_AREA(v2,v3,area(ots))
!        denom=denom+ SQRT(v1(1)**2+v1(2)**2+v1(3)**2)* &
!          (chgtemp%rho(npos(1),npos(2),npos(3))-rho)*area(ots)
!        nom(ots)=chgtemp%rho(npos(1),npos(2),npos(3))*area(ots)* &
!          (chgtemp%rho(npos(1),npos(2),npos(3))-rho)
!      ELSE IF ( ots == 3 .OR. ots ==4) THEN
!        CALL FIND_AREA(v1,v3,area(ots))
!        denom=denom+ SQRT(v1(1)**2+v1(2)**2+v1(3)**2)* &
!          (chgtemp%rho(npos(1),npos(2),npos(3))-rho)*area(ots)
!        nom(ots)=chgtemp%rho(npos(1),npos(2),npos(3))*area(ots)* &
!          (chgtemp%rho(npos(1),npos(2),npos(3))-rho)
!      ELSE IF ( ots == 5 .OR. ots ==6) THEN
!        CALL FIND_AREA(v1,v2,area(ots))
!        denom=denom+ SQRT(v1(1)**2+v1(2)**2+v1(3)**2)* &
!          (chgtemp%rho(npos(1),npos(2),npos(3))-rho)*area(ots)
!        nom(ots)=chgtemp%rho(npos(1),npos(2),npos(3))*area(ots)* &
 !         (chgtemp%rho(npos(1),npos(2),npos(3))-rho)
!      END IF
      areasum=areasum+area(ots)
      ! 1 neighbor with higher density
      hdn =hdn+1
      ! Get the information of that 
      ! read its volnum, if is the first, record its volnum
      IF (fb) THEN
        fb=.FALSE.
        ! The volnum of neighbors need to be kept track of
        tempnvol=1 !
        tempvolnum=bdr%volnum(npos(1),npos(2),npos(3))
      ELSE IF(bdr%volnum(npos(1),npos(2),npos(3))/=tempvolnum) THEN
        ! I'm a boundary point. 
        ! More than one neighbors with higher density has been found
        tempnvol=tempnvol+1
      END IF

!      PRINT *, 'AFTER SECOND IF'  
    END IF                    

  END SUBROUTINE


  ! This subroutine gets flux as well

  SUBROUTINE get_nvols (sortedList,bdr,chgtemp,ions)
    TYPE(weight_obj),DIMENSION(:) :: sortedList
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chgtemp
    TYPE(ions_obj) :: ions
    INTEGER :: i,j,t1,t2,t3
    INTEGER,DIMENSION(3) :: xyz,p1,p2,p3,p4,p5,p6
    REAL(q2),DIMENSION(3) :: v1,v2,v3
    REAL(q2) :: area1,area2,area3,rho,denom
    REAL(q2),DIMENSION(6) :: nom
    LOGICAL :: ismax,bein
    t2=0
    t3=0
    DO i=1,SIZE(sortedList)
      ismax=.TRUE.
      xyz(1)=sortedList(i)%pos(1)
      xyz(2)=sortedList(i)%pos(2)
      xyz(3)=sortedList(i)%pos(3)
      p1=xyz+(/1,0,0/)
      p2=xyz+(/-1,0,0/)
      p3=xyz+(/0,1,0/)
      p4=xyz+(/0,-1,0/)
      p5=xyz+(/0,0,1/)
      p6=xyz+(/0,0,-1/)
      CALL pbc(p1,chgtemp%npts)
      CALL pbc(p2,chgtemp%npts)
      CALL pbc(p3,chgtemp%npts)
      CALL pbc(p4,chgtemp%npts)
      CALL pbc(p5,chgtemp%npts)
      CALL pbc(p6,chgtemp%npts)
      v1=ions%lattice(1,:)
      v2=ions%lattice(2,:)
      v3=ions%lattice(3,:)
      CALL find_area(v2,v3,area1)
      CALL find_area(v1,v3,area2)
      CALL find_area(v1,v2,area3)
      area1=ABS(area1)
      area2=ABS(area2)
      area3=ABS(area3)
      denom=0
      nom=(/0,0,0,0,0,0/)
      rho=sortedList(i)%rho
      IF (rho<chgtemp%rho(p1(1),p1(2),p1(3))) THEN
        ismax=.FALSE.
        nom(1)=SQRT(v1(1)**2+v1(2)**2+v1(3)**2)* &
         &(chgtemp%rho(p1(1),p1(2),p1(3))-rho)*area1
        denom=denom+SQRT(v1(1)**2+v1(2)**2+v1(3)**2)* &
         &(chgtemp%rho(p1(1),p1(2),p1(3))-rho)*area1
      END IF

      IF (rho<chgtemp%rho(p2(1),p2(2),p2(3))) THEN
        ismax=.FALSE.
        nom(2)=SQRT(v1(1)**2+v1(2)**2+v1(3)**2)* &
         &(chgtemp%rho(p2(1),p2(2),p2(3))-rho)*area1
        denom=denom+SQRT(v1(1)**2+v1(2)**2+v1(3)**2)* &
         &(chgtemp%rho(p2(1),p2(2),p2(3))-rho)*area1
      END IF

      IF (rho<chgtemp%rho(p3(1),p3(2),p3(3))) THEN
        ismax=.FALSE.
        nom(3)=SQRT(v2(1)**2+v2(2)**2+v2(3)**2)* &
         &(chgtemp%rho(p3(1),p3(2),p3(3))-rho)*area2
        denom=denom+SQRT(v2(1)**2+v2(2)**2+v2(3)**2)* &
         &(chgtemp%rho(p3(1),p3(2),p3(3))-rho)*area2
      END IF

      IF (rho<chgtemp%rho(p4(1),p4(2),p4(3))) THEN
        nom(4)=SQRT(v2(1)**2+v2(2)**2+v2(3)**2)* &
         &(chgtemp%rho(p4(1),p4(2),p4(3))-rho)*area2
        ismax=.FALSE.
        denom=denom+SQRT(v2(1)**2+v2(2)**2+v2(3)**2)* &
         &(chgtemp%rho(p4(1),p4(2),p4(3))-rho)*area2
      END IF

      IF (rho<chgtemp%rho(p5(1),p5(2),p5(3))) THEN
        ismax=.FALSE.
        nom(5)=SQRT(v3(1)**2+v3(2)**2+v3(3)**2)* &
         &(chgtemp%rho(p5(1),p5(2),p5(3))-rho)*area3
        denom=denom+SQRT(v3(1)**2+v3(2)**2+v3(3)**2)* &
         &(chgtemp%rho(p5(1),p5(2),p5(3))-rho)*area3
      END IF

      IF (rho<chgtemp%rho(p6(1),p6(2),p6(3))) THEN
        ismax=.FALSE.
        nom(6)=SQRT(v3(1)**2+v3(2)**2+v3(3)**2)* &
         &(chgtemp%rho(p6(1),p6(2),p6(3))-rho)*area3
        denom=denom+SQRT(v3(1)**2+v3(2)**2+v3(3)**2)* &
         &(chgtemp%rho(p6(1),p6(2),p6(3))-rho)*area3
      END IF
      ! can't know if a point is interior here because there may be a point with
      ! flux to 2 points but to the same basin
      DO j=1,6
        sortedList(i)%flx(j)=nom(j)/denom
      END DO
      bein=.FALSE.
      t1=0
        
!      DO j=1,6
!        IF (sortedList(i)%flx(j)/=1.) t1=t1+1
!      END DO
!      IF (t1==6) THEN
!        t2=t2+1
!        PRINT *, 'i is',i
!        PRINT *, 'nom is',nom
!        PRINT *, 'denom is',denom
!        PRINT *, 'flx is',sortedList(i)%flx
!      END IF
      IF (ismax) bdr%nvols=bdr%nvols+1
    END DO
!    PRINT *, t2,t3,120*120*120-t2-t3
  END SUBROUTINE 
!
!
!  SUBROUTINE get_neighbor_weight(sortedList,chgList,npos,i,chgtemp)
!    TYPE(weight_obj),DIMENSION(:) :: sortedList,chgList
!    TYPE(charge_obj) :: chgtemp
!    INTEGER,DIMENSION(3) :: npos
!    INTEGER :: i,n1,n2,nindex
!    ! Assumes that the neighbor does not have weight to a basin that this
!    ! boundary point has no flux to
!    CALL pbc(npos,chgtemp%npts)
!    nindex=(npos(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(npos(2)-1)*chgtemp%npts(3)&
!      +npos(3)
!    DO n1=1,6
!      DO n2=1,6
!        IF (chgList(nindex)%nbrvol(n2)==sortedList(i)%nbrvol(n1)) THEN
!          sortedList(i)%volwgt(n1)=sortedList(i)%nbrflx(n1)*chgList(nindex)%volwgt(n2)
!        
!        END IF
!      END DO
!    END DO
!  END SUBROUTINE
!
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
