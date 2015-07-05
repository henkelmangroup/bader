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
    REAL(q2) :: tempflx,totalE,tempw
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
    totalE=0
    DO n1=1,bdr%nvols
      totalE=totalE+bdr%volchg(n1)
    END DO
    PRINT *, 'number of electrons:', totalE
    DO n1=1,SIZE(sortedList)
      tempflx=0
      DO n2=1,6
        tempflx=tempflx+sortedList(n1)%flx(n2)
      END DO
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
    REAL(q2) :: temprho,denom,rho,areasum,tempw
    REAL(q2),DIMENSION(6) :: nom,area
    INTEGER,DIMENSION(3) :: npos !neighbor position
    INTEGER,DIMENSION(3) :: xyz ! xyz indexes of current point 
    INTEGER :: nbd,nin,hdn,temppos,ots ! one to six
    LOGICAL :: fb
    bdr%nvols=0 
    CALL get_nvols(sortedList,bdr,chgtemp,ions)
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
    CALL assign_weight(sortedList,chgList,bdr,chgtemp,chgval)
!    DO i=100,length
!      tempw=0
!      DO n2=1,bdr%nvols
!        tempw=tempw+sortedList(i)%volwgt(n2)
!      END DO
!      IF (temp /=1) THEN 
!        PRINT *,'point',i,'in sorted list  has total weight',tempw
!        PRINT *,'they are'
!        DO n2=1,bdr%nvols
!          PRINT *,sortedList(i)%volwgt(n2)
!        END DO
!      END IF
!    END DO

    DO i=1,length
      DEALLOCATE(sortedList(i)%volnum)
      DEALLOCATE(sortedList(i)%volwgt)
    END DO
  
  END SUBROUTINE
!
  SUBROUTINE assign_weight(sortedList,chgList,bdr,chgtemp,chgval)
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj):: chgtemp,chgval
    INTEGER         :: i,nvols,cin,ots,n1,n2,n3,n4,sind
    LOGICAL         :: ismax
    REAL(q2)        :: rho,temp
    INTEGER,DIMENSION(3) :: xyz,p,add
    TYPE(weight_obj),DIMENSION(:) :: sortedList,chgList

    nvols=1
    DO i=1, SIZE(sortedList)
      xyz(1)=sortedList(i)%pos(1)
      xyz(2)=sortedList(i)%pos(2)
      xyz(3)=sortedList(i)%pos(3)
      ! First assign weights to maximas. This has to be done to both sorted list
      ! and chglist
      rho=sortedList(i)%rho
      ismax=.TRUE.
!      ! diagonal neighbors needs to be checked as well
!      p=xyz+(/1,0,0/)
!      CALL pbc(p,chgtemp%npts)
!      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
!        ismax=.FALSE.
!        ots=1
!        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
!        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
!      END IF
!
!      p=xyz+(/-1,0,0/)
!      CALL pbc(p,chgtemp%npts)
!      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
!        ismax=.FALSE.
!        ots=2
!        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
!        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
!      END IF
!
!      p=xyz+(/0,1,0/)
!      CALL pbc(p,chgtemp%npts)
!      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
!        ismax=.FALSE.
!        ots=3
!        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
!        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
!      END IF
!
!      p=xyz+(/0,-1,0/)
!      CALL pbc(p,chgtemp%npts)
!      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
!        ismax=.FALSE.
!        ots=4
!        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
!        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
!      END IF
!
!      p=xyz+(/0,0,1/)
!      CALL pbc(p,chgtemp%npts)
!      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
!        ismax=.FALSE.
!        ots=5
!        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
!        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
!      END IF
!
!      p=xyz+(/0,0,-1/)
!      CALL pbc(p,chgtemp%npts)
!      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
!        ismax=.FALSE.
!        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
!        ots=6
!        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
!      END IF
      ! diagonla neighbors needs to be treated separately. They do not get flux
      ! as there is no area between them. Dallas treat them as interior points.
      ! I can copy weight and basin info of the higher neighbor to this point.
      ! a loop that takes care of all neighbors
!      p=xyz+(1,1,0)
!      CALL pbc(p,chgtemp%npts)
!      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
!        ismax=.FALSE.
!        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
!        ots=6
!        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
!      END IF
!      p=xyz+(1,-1,0)
!      CALL pbc(p,chgtemp%npts)
!      IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
!        ismax=.FALSE.
!        cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
!        ots=6
!        CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
!      END IF
!      p=xyz+(-1,1,0,)
!      p=xyz+(-1,-1,0)
!      p=xyz+(1,0,1)
!      p=xyz+(1,0,-1)
!      p=xyz+(-1,0,1)
!      p=xyz+(-1,0,-1)
!      p=xyz+(0,1,1)
!      p=xyz+(0,1,-1)
!      p=xyz+(0,-1,1)
!      p=xyz+(0,-1,-1)
!      p=xyz+(1,1,1)
!      p=xyz+(-1,1,1)
!      p=xyz+(1,-1,1)
!      p=xyz+(1,1,-1)
!      p=xyz+(-1,-1,1)
!      p=xyz+(1,-1,-1)
!      p=xyz+(-1,1,-1)
!      p=xyz+(-1,-1,-1)

      ! below is the attempt of a better version, not yet successful.
      ots=0
!      IF (i==20) PRINT *,'point 20 with xyz',xyz 
      DO n1=-1,1
        DO n2=-1,1
          DO n3=-1,1
            add=(/n1,n2,n3/)
            p=xyz+add ! no problem in this
            IF (i==101) THEN
              PRINT *,'this is point 101 with xyz',xyz
            END IF
            IF (add(1)**2+add(2)**2+add(3)**2 == 0) THEN
              ! self, pass
              CYCLE
            ELSE IF (add(1)**2+add(2)**2+add(3)**2 == 1) THEN
            ! direct neighbor, calculate flux and weight
              CALL pbc(p,chgtemp%npts)
              ots=ots+1
              IF (rho<chgtemp%rho(p(1),p(2),p(3))) THEN
                ismax=.FALSE.
!                IF (i==20) PRINT *,'looking at direct p',p
!                IF (i==20) PRINT *,'found p has higher rho'
                cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3) 
                CALL flux_weight(sortedList,i,cin,chgList,bdr,ots,chgtemp,chgval)
                IF (i==101) THEN
                  PRINT *,'neighbor p',p,'has higher rho'
                 
                END IF
              END IF
            ELSE 
            ! indirect neighbor, only copy info if higher in density
              CALL pbc(p,chgtemp%npts)
            ! also only check for indirect neighbors only when there is no known
            ! direct neighbor with higher densities.
              IF (rho<chgtemp%rho(p(1),p(2),p(3)) .AND. ismax) THEN
                ismax=.FALSE.
!                IF (i==20) PRINT *, 'looking at diagonal p',p
!                IF (i==20) PRINT *,'found p has higher rho'
                ! copy all infos from this basin to next
                cin=(p(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(p(2)-1)*chgtemp%npts(3)+p(3)
                sind=(sortedList(i)%pos(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+ &
                     (sortedList(i)%pos(2)-1)*chgtemp%npts(3)+sortedList(i)%pos(3)
                DO n4=1,bdr%nvols
!                  sortedList(i)%volwgt(n4)=chgList(cin)%volwgt(n4)
                  chgList(sind)%volwgt(n4)=chgList(cin)%volwgt(n4)
                END DO
              END IF
            END IF            
          END DO
        END DO
      END DO
!      END IF
      IF (ismax==.TRUE.) THEN
        !print *, 'i is a maximum',i
        cin=(xyz(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(xyz(2)-1)*chgtemp%npts(3)+xyz(3)
        !PRINT *,'cin is',cin
        sortedList(i)%volwgt(nvols)=1
        chgList(cin)%volwgt(nvols)=1
        nvols=nvols+1
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
      !sind should be sorted list index to CHGLIST index
      sind=(sortedList(i)%pos(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(sortedList(i)%pos(2)-1)*chgtemp%npts(3)&
        +sortedList(i)%pos(3)
!      IF (i==20) THEN
!        PRINT *,'sind is',sind
!        PRINT *,'cin is',cin
!        PRINT *,'this cin has volwgt:'
!        DO n1=1,bdr%nvols
!          PRINT *,chgList(cin)%volwgt(n1)
!        END DO
!      END IF
      IF (i==101) THEN
        PRINT *,'looking at point cin',cin
        PRINT *,'with xyz', & 
                 chgList(cin)%pos(1),chgList(cin)%pos(2),chgList(cin)%pos(3)
        PRINT *,'this is weight array of point i'
        DO n1=1,bdr%nvols
          PRINT *, sortedList(i)%volwgt(n1)
        END DO
        PRINT *,'this is sind',sind
        PRINT *,'with xyz', chgList(sind)%pos(1), &
                chgList(sind)%pos(2),chgList(sind)%pos(3)
      END IF  
      DO n1=1,bdr%nvols
        
        sortedList(i)%volwgt(n1)=sortedList(i)%volwgt(n1)+sortedList(i)%flx(ots)*chgList(cin)%volwgt(n1)
        IF (i==101) PRINT *,chgList(cin)%volwgt(n1)
        chgList(sind)%volwgt(n1)=sortedList(i)%volwgt(n1)    
        x=sortedList(i)%pos(1)
        y=sortedList(i)%pos(2)
        z=sortedList(i)%pos(3)
        
      END DO
!      PRINT *,'i is',i,'with xyz',x,y,z
      
!      IF (i==20) THEN
!        PRINT *,'by the end this point has volwgt'
!        DO n1=1,bdr%nvols
!          PRINT *,sortedList(i)%volwgt(n1)
!        END DO
!      END IF
  END SUBROUTINE



  ! This subroutine gets flux as well

  SUBROUTINE get_nvols (sortedList,bdr,chgtemp,ions)
    TYPE(weight_obj),DIMENSION(:) :: sortedList
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chgtemp
    TYPE(ions_obj) :: ions
    INTEGER :: i,j,n1,n2,n3
    INTEGER,DIMENSION(3) :: xyz,p,p1,p2,p3,p4,p5,p6
    REAL(q2),DIMENSION(3) :: v1,v2,v3
    REAL(q2) :: area1,area2,area3,rho,denom
    REAL(q2),DIMENSION(6) :: nom
    LOGICAL :: ismax,bein
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
      area1=ABS(area1) ! for x directions
      area2=ABS(area2) ! for y directions
      area3=ABS(area3) ! for z directions
      ! in assign weight, the access to these fluxes are in order:
      ! -1,0,0 | 0,-1,0 | 0,0,-1
      !  0,0,1 | 0, 1,0 | 1,0, 0
      ! that is :
      ! 1,2,3,3,2,1
      ! or p2,p4,p6,p5,p3,p1
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

      DO n1=-1,1
        DO n2=-1,1
          DO n3=-1,1
            p=xyz+(/n1,n2,n3/)
            CALL pbc(p,chgtemp%npts)
            IF (n1**2+n2**2+n3**2 > 1 .AND. &
                rho<chgtemp%rho(p(1),p(2),p(3))) THEN
              ismax=.FALSE.
              ! no nom and denom calc as there is no flx
            END IF
          END DO
        END DO
      END DO

      ! can't know if a point is interior here because there may be a point with
      ! flux to 2 points but to the same basin
 
      ! storage sequence should be:
      ! p2,4,6,5,3,1 
!      DO j=1,6
!        sortedList(i)%flx(j)=nom(j)/denom
!      END DO
      sortedList(i)%flx(1)=nom(2)/denom
      sortedList(i)%flx(2)=nom(4)/denom
      sortedList(i)%flx(3)=nom(6)/denom
      sortedList(i)%flx(4)=nom(5)/denom
      sortedList(i)%flx(5)=nom(3)/denom
      sortedList(i)%flx(6)=nom(1)/denom
      bein=.FALSE.
      IF (ismax) THEN
        ! debugging block
        PRINT *,xyz,'is a maximum'
        !PRINT *, 'and the rho is', rho
        ! end block
        bdr%nvols=bdr%nvols+1
      END IF
    END DO
  END SUBROUTINE get_nvols
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
