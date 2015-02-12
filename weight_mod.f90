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
      ! This structure will definetly not work correctly with voronoi and
      ! irregular grid for not. For non-cubic grid, there is a chance that the
      ! area of the plane shared by two boundary cells will be not calculated
      ! correctly
      INTEGER,DIMENSION(6) :: nbrvol ! volnum of neighbors
      REAL(q2),DIMENSION(6) :: volwgt,nbrflx,nbrwgt ! weight to the basin and flux to a neighbor
      LOGICAL :: isInterior
    END TYPE

    PUBLIC :: weight_obj
    PUBLIC :: bader_weight_calc,quick_sort

  CONTAINS 

  SUBROUTINE bader_weight_calc(bdr,ions,chgval,opts)

    ! chgList has 4 slots. 1st is rho, 2 3 4 are the coordinate indices.
    TYPE(weight_obj),ALLOCATABlE,DIMENSION(:) :: chgList, sortedList,tempList
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chgval,chgtemp
    TYPE(weight_obj) :: wt
    TYPE(ions_obj) :: ions,ionstemp
    TYPE(options_obj) :: opts
    INTEGER(KIND=8) :: totalLength,temp
    INTEGER :: i, j, k, l, n1, n2, n3, walker
    
    PRINT *, 'THIS IS IN WEIGHT. OPTS%REF_FLAG IS ',opts%ref_flag   
    IF (opts%ref_flag) THEN
      CALL read_charge_ref(ionstemp,chgtemp,opts)
    ELSE
      chgtemp = chgval
    END IF
    ALLOCATE(bdr%volnum(chgtemp%npts(1),chgtemp%npts(2),chgtemp%npts(3)))
    ALLOCATE(bdr%nnion(bdr%nvols), bdr%iondist(bdr%nvols),bdr%ionchg(ions%nions))
    totalLength=chgtemp%npts(1)*chgtemp%npts(2)*chgtemp%npts(3)
    !totalLength = 8000 ! for testing purpose only
    ALLOCATE (sortedList(totalLength)) ! sortedList should have the correct size
    !PRINT *,'assigned totalLength', totalLength
    temp=totalLength
    ! The length of the array to be sorted should be of size 2^n. The additional
    ! buffer slots will be empty
    totalLength=2**(CEILING(LOG(temp*1.0_q2)/LOG(2.0_q2)))
    ! PRINT *, 'enlarged total length', totalLength
    ALLOCATE (chgList(totalLength))
    ALLOCATE (tempList(totalLength))
    ! first merge sort   
    ! the lines below are commented out so that the program does not actually
    ! run using the CHGCAR, but the mannualy put in array A
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
          chgList(walker)%nbrvol=(/0,0,0,0,0,0/)
          chgList(walker)%volwgt=(/0,0,0,0,0,0/)
          chgList(walker)%nbrflx=(/0,0,0,0,0,0/)
        END DO
      END DO
    END DO
    PRINT *,'THIS IS IN WEIGHT. CHGVAL AT 90 90 30 IS',chgval%rho(90,90,30)
    PRINT *,'THIS IS IN WEIGHT. CHGVAL AT 91 91 31 IS',chgval%rho(91,91,31)
    PRINT *,'THIS IS IN WEIGHT. CHGTEMP AT 90 90 30 IS',chgtemp%rho(90,90,30)
    PRINT *,'THIS IS IN WEIGHT. CHGTEMP AT 91 91 31 IS',chgtemp%rho(91,91,31)
    PRINT *, 'the lattice vectors are',ions%lattice(1,:)
    PRINT *, 2/2,3/2,4/2,5/2,6/2,7/2



    ! Here begins the manual array
  !  PRINT *,'PRINTING THE CONTENT OF chgList'
  !  walker=1
  !  DO n1=1, 20
  !    DO n2=1, 20
  !      DO n3=1, 20
  !      chgList(walker)%rho = n1*n2*n3
  !      chgList(walker)%x = n1
  !      chgList(walker)%y = n2
  !      chgList(walker)%z = n3
  !      walker=walker+1
  !      END DO
  !    END DO
  !  END DO
  !  PRINT *, 'walker is ', walker
    ! after the last loop is finished, 
    ! The empty slots are filled with junk. Lets fill it with -1
    DO i=walker,totalLength
      chgList(walker)%rho=-1
      chgList(walker)%pos(1)=-1
      chgList(walker)%pos(2)=-1
      chgList(walker)%pos(3)=-1
      walker=walker+1  
    END DO
    tempList=chgList
    CALL quick_sort(tempList)
    DO i=1,SIZE(sortedList)
      sortedList(i)=tempList(totalLength-i+1)
    END DO
 
    PRINT *, 'going to call weight_calc'
    CALL weight_calc(sortedList,bdr,chgtemp,chgval,chgList,ions)
 
  END SUBROUTINE bader_weight_calc

  ! lets say the sortedList is correct. Now let's start to walk the list from top
  ! to bottom. Relavent variablees: bdr%volnum. bdr%volpos
  
  
  !-----------------------------------------------------------------------------------!
  ! weight calc 
  !-----------------------------------------------------------------------------------!
  SUBROUTINE weight_calc(sortedList,bdr,chgtemp,chgval,chgList,ions)
      ! apply periodic boudary condition first. ?????
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
    PRINT *, 'CALLED'
    bdr%nvols=0 
    length=SIZE(sortedList)
    nbd=0
    nin=0
    temp2=0
    temp=SIZE(sortedList)/100
    DO i=1, length
      !PRINT *, 'ENTERED LOOP'
      ! count number of neighbors with higher density
!      i=9
      hdn=0
      xyz(1)=sortedList(i)%pos(1)
      xyz(2)=sortedList(i)%pos(2)
      xyz(3)=sortedList(i)%pos(3)

      npos=xyz+(/1,0,0/)
      denom=0
      areasum=0
      DO j=1,6
        nom(j) = 0
        sortedList(i)%nbrvol(j)=0
      END DO
      rho=sortedList(i)%rho
      tempnvol=0
      fb=.TRUE.
      ! %nbr(1),%vonwgt(1)
      ots=1  
      CALL neighbors(bdr,npos,xyz,chgtemp,rho,nbd,nin,i,ions, & 
tempnvol,tempvolnum,hdn,fb,sortedList,denom,nom,ots,area,areasum)
      npos=xyz+(/-1,0,0/)
      ! %nbr(2),%vonwgt(2)
      ots=2  
      CALL neighbors(bdr,npos,xyz,chgtemp,rho,nbd,nin,i,ions, & 
tempnvol,tempvolnum,hdn,fb,sortedList,denom,nom,ots,area,areasum)
      npos=xyz+(/0,1,0/)
      ! %nbr(3),%vonwgt(3)
      ots=3  
      CALL neighbors(bdr,npos,xyz,chgtemp,rho,nbd,nin,i,ions, & 
tempnvol,tempvolnum,hdn,fb,sortedList,denom,nom,ots,area,areasum)
      npos=xyz+(/0,-1,0/)
      ! %nbr(4),%vonwgt(4)
      ots=4  
      CALL neighbors(bdr,npos,xyz,chgtemp,rho,nbd,nin,i,ions, &
tempnvol,tempvolnum,hdn,fb,sortedList,denom,nom,ots,area,areasum)
      npos=xyz+(/0,0,1/)
      ! %nbr(5),%vonwgt(5)
      ots=5  
      CALL neighbors(bdr,npos,xyz,chgtemp,rho,nbd,nin,i,ions, &
tempnvol,tempvolnum,hdn,fb,sortedList,denom,nom,ots,area,areasum)
      npos=xyz+(/0,0,-1/)
      ! %nbr(6),%vonwgt(6)
      ots=6  
      CALL neighbors(bdr,npos,xyz,chgtemp,rho,nbd,nin,i,ions, & 
tempnvol,tempvolnum,hdn,fb,sortedList,denom,nom,ots,area,areasum)
      ! This finds the coresponding entry in chgList before sorting
!      IF (xyz(1)==1 .AND. xyz(2)==1 .AND. xyz(3)==1) THEN
!        PRINT *,'------------------- NEW ENTRY ------------------'
!        PRINT *, 'i is,',i
!        PRINT *,'XYZ are',xyz(1),xyz(2),xyz(3)
        temppos=(xyz(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+ & 
          (xyz(2)-1)*chgtemp%npts(3)+ xyz(3)
!        PRINT *, 'Position',xyz(1),xyz(2),xyz(3),'translates to'
!        PRINT *, 'temp index', temppos
!        PRINT *, 'which corresponds to position',chgList(temppos)%pos
!      END IF
!      IF (xyz(1)==1 .AND. xyz(2)==2 .AND. xyz(3)==3) THEN
!        PRINT *,'tempnvol is', tempnvol
!        PRINT *, 'the calculated index in chglist is', temppos
!        PRINT *, 'correspond to this point in chglist', chgList(temppos)%pos
!        PRINT *, 'tempvolnum is', tempvolnum
!      END IF
      IF (tempnvol >1 .OR. tempvolnum==0 ) THEN
        ! is a boundary point
        nbd=nbd+1
        sortedList(i)%isInterior = .FALSE.
        chgList(temppos)%isInterior=.FALSE.
        ! Need to calculate the weight for each of those neighbors. denom should
        ! be working. Need to find area, cross ptoducts of vectors.
        
      END IF
      IF (tempnvol==1 .AND. tempvolnum/=0 ) THEN
        ! is an interior point. 
!        PRINT *, 'tempnvol is 1 for',xyz(1),xyz(2),xyz(3)
        bdr%volnum(xyz(1),xyz(2),xyz(3))=tempvolnum
!        nin=nin+1
        sortedList(i)%isInterior= .TRUE.
        chgList(temppos)%isInterior= .TRUE.
        sortedList(i)%nbrvol(1)=tempvolnum
        sortedList(i)%volwgt(1)=1
      END IF  

      IF (sortedList(i)%isInterior) THEN
        nin=nin+1
      END IF 
      ! I need to take care of boundary points right here. After ots has went
      ! from 1 to 6, all flux can be calculated.



!      IF (xyz(1)==1 .AND. xyz(2)==1 .AND. xyz(3)==1) THEN
!        PRINT *, 'tempnvol is ',tempnvol
!        PRINT *, 'tempvolnum is ',tempvolnum
!        PRINT *, 'volnum is', bdr%volnum(1,1,1)
!        PRINT *, 'xyz are',xyz(1),xyz(2),xyz(3)
!        PRINT *, chgList(temppos)%pos,'is interior',chgList(temppos)%isInterior
!        PRINT *, 'this point has volnum', &
!          bdr%volnum(chgList(temppos)%pos(1),chgList(temppos)%pos(2),chgList(temppos)%pos(3))
!        PRINT *, 'this point index in sortedList is', i
!      END IF

      DO j=1,6
        sortedList(i)%nbrflx(j)=nom(j)/denom
      END DO
      ! after all flux is calculated, look at each neighbor's basin assignment
      ! and weight. Need to look into either chgList or sortedList
      ! USE some subroutine to do this.
!      CALL get_neighbor_weight(sortedList,chgList,npos,i,chgtemp)

      IF (i>temp) THEN
        temp2=temp2+1
        PRINT *, temp2*20 ,'percent done'
        temp=temp+SIZE(sortedList)/5
      END IF
       

      IF (hdn==0) THEN
!          PRINT *,'FOUND A LOCAL MAXIMUM'
        sortedList(i)%isInterior=.TRUE.
        bdr%nvols=bdr%nvols+1
!        PRINT *, 'nvols is now', bdr%nvols
        bdr%volnum(xyz(1),xyz(2),xyz(3))=bdr%nvols
!          PRINT *,'TO POINT',xyz(1),xyz(2),xyz(3)
!          PRINT *, 'it now has volnum',bdr%volnum(xyz(1),xyz(2),xyz(3))
!          PRINT *, 'GIVEN volnum'
!          PRINT *, '91,91,31 has volnum',  bdr%volnum(91,91,31)
      END IF 
!      PRINT *, sortedList(i)%isInterior
!      IF (.NOT. sortedList(i)%isInterior .AND. tempnvol<=1) THEN
!        PRINT *, sortedList(i)%pos
!      END IF
!      IF (bdr%volnum(xyz(1),xyz(2),xyz(3))==0) THEN
!        PRINT *, xyz
!      END IF

    END DO
    ALLOCATE(bdr%volchg(bdr%nvols))
    bdr%volchg = 0._q2
!    PRINT *, 'Adding up charges'
!    PRINT *, 'bdr%nvols is',bdr%nvols
    ! how to go from x y z to sorted i?

    ! This is flux to each neighbor


    DO n1 = 1,chgval%npts(1)
      DO n2 = 1,chgval%npts(2)
        DO n3 = 1,chgval%npts(3)
          i=(n1-1)*chgval%npts(2)*chgval%npts(3)+(n2-1)*chgval%npts(3)+n3
!          PRINT *, n1,n2,n3,'is interior',chgList(i)%isInterior
          IF (chgList(i)%isInterior==.FALSE.) THEN
            ! Boundary points have a 6 slot array. Some of them have information
            ! about where the density is flowing to. Need to look at all of them
            ! But first lets find flux
            CYCLE
          END IF

!          PRINT *, 'bdr%volnum is',bdr%volnum(n1,n2,n3)
          IF (bdr%volnum(n1,n2,n3) == bdr%nvols+1) CYCLE
!          PRINT *, 'volnum of ',n1,n2,n3,'is',bdr%volnum(n1,n2,n3)
!          PRINT *, 'calculated i is', i
!          PRINT *, 'sortedList(i) has coordinates', chgList(i)%pos
!          PRINT *, 'is this an interior point?', chgList(i)%isInterior

          bdr%volchg(bdr%volnum(n1,n2,n3)) = &
          &  bdr%volchg(bdr%volnum(n1,n2,n3)) + chgval%rho(n1,n2,n3)
        END DO
      END DO
    END DO
    bdr%volchg = bdr%volchg/REAL(chgval%nrho,q2)
    PRINT *, bdr%volchg

    PRINT *, 'chgval%nrho is', chgval%nrho
    PRINT *, 'the number of boundary points are ', nbd 
    PRINT *, 'the number of interior points are', nin
    PRINT *, 'the total number of point is', 120*120*120
    PRINT *, 'the number of points not identified is', &
      120*120*120-nbd-nin
!    PRINT *, bdr%volnum(10,10,10)
    !PRINT *, sortedList(990)%pos +(/3,3,3/)
  END SUBROUTINE

  SUBROUTINE  neighbors(bdr,npos,xyz,chgtemp,rho, nboundary,ninterior,i &
                ,ions,tempnvol,tempvolnum,hdn,fb,sortedList, & 
                denom,nom,ots,area,areasum )
  TYPE(weight_obj),DIMENSION(:) :: sortedList
  TYPE(bader_obj) :: bdr
  TYPE(charge_obj) :: chgtemp  
  TYPE(ions_obj) :: ions
  INTEGER,DIMENSION(3) :: npos,xyz
  REAL(q2),DIMENSION(3) ::  v1,v2,v3
  INTEGER :: nboundary,ninterior,ots
  LOGICAL :: fb ! first bigger neighbor
  INTEGER :: temp,tempnvol,tempvolnum,temp2,i,hdn
  REAL(q2) :: rho,denom,areasum
  REAL(q2),DIMENSION(6) :: nom,area
    CALL pbc(npos,chgtemp%npts)
    IF (rho<chgtemp%rho(npos(1),npos(2),npos(3))) THEN
      v1=ions%lattice(1,:)
      v2=ions%lattice(2,:)
      v3=ions%lattice(3,:)
      If (ots ==1 .OR. ots ==2) THEN
        ! a step taken in the x direction
        CALL FIND_AREA(v2,v3,area(ots))
        denom=denom+ SQRT(v1(1)**2+v1(2)**2+v1(3)**2)* &
          (chgtemp%rho(npos(1),npos(2),npos(3))-rho)*area(ots)
        nom(ots)=chgtemp%rho(npos(1),npos(2),npos(3))*area(ots)* &
          (chgtemp%rho(npos(1),npos(2),npos(3))-rho)
      ELSE IF ( ots == 3 .OR. ots ==4) THEN
        CALL FIND_AREA(v1,v3,area(ots))
        denom=denom+ SQRT(v1(1)**2+v1(2)**2+v1(3)**2)* &
          (chgtemp%rho(npos(1),npos(2),npos(3))-rho)*area(ots)
        nom(ots)=chgtemp%rho(npos(1),npos(2),npos(3))*area(ots)* &
          (chgtemp%rho(npos(1),npos(2),npos(3))-rho)
      ELSE IF ( ots == 5 .OR. ots ==6) THEN
        CALL FIND_AREA(v1,v2,area(ots))
        denom=denom+ SQRT(v1(1)**2+v1(2)**2+v1(3)**2)* &
          (chgtemp%rho(npos(1),npos(2),npos(3))-rho)*area(ots)
        nom(ots)=chgtemp%rho(npos(1),npos(2),npos(3))*area(ots)* &
          (chgtemp%rho(npos(1),npos(2),npos(3))-rho)
      END IF
      areasum=areasum+area(ots)
      ! 1 neighbor with higher density
      hdn =hdn+1
      sortedList(i)%nbrvol(ots)=bdr%volnum(npos(1),npos(2),npos(3))
      ! Thus far I have the information of each cell's flux towards neighbor
      ! cells. I still need to look up neighbor's weight. It is easy when
      ! neighbors are interior points. But if neighbor is a boundary point it
      ! self, then this flux goes to different basis. Need to look up neighbors'
      ! weight towards each basin. 


      ! read its volnum, if is the first, record its volnum
      IF (fb) THEN
        fb=.FALSE.
        ! The volnum of neighbors need to be kept track of
        tempnvol=1 !
        tempvolnum=bdr%volnum(npos(1),npos(2),npos(3))
!        IF (xyz(1)== 1 .AND. xyz(2)== 1 .AND. xyz(3)== 1 ) THEN
!          PRINT *,'1 1 1 in basin', bdr%volnum(1,1,1)
!          PRINT *,'looking into basin ',bdr%volnum(npos(1),npos(2),npos(3))
!          PRINT *, 'temp volnum is', tempvolnum
!        END IF
      ELSE IF(bdr%volnum(npos(1),npos(2),npos(3))/=tempvolnum) THEN
        ! I'm a boundary point. 
        ! More than one neighbors with higher density has been found
        tempnvol=tempnvol+1
        sortedList(i)%nbrvol(ots)=bdr%volnum(npos(1),npos(2),npos(3))
      END IF

!      PRINT *, 'AFTER SECOND IF'  
    END IF                    
    IF (xyz(1)==31 .AND. xyz(2)==31 .AND. xyz(3)==66) THEN
      PRINT *, 'SITTING ON 1 1 1, SHOULD BE A BOUNDARY POINT'
      PRINT *, 'LOOKING AT',npos(1),npos(2),npos(3)
      PRINT *, 'IT HAS bdr volnum',bdr%volnum(npos(1),npos(2),npos(3)) 
      PRINT *, 'TEMPNVOL IS', tempnvol       
      PRINT *, 'tempvolnum is', tempvolnum
      PRINT *, 'nbrvol of this point is', sortedList(i)%nbrvol
    END IF
      

  END SUBROUTINE


  SUBROUTINE get_neighbor_weight(sortedList,chgList,npos,i,chgtemp)
    TYPE(weight_obj),DIMENSION(:) :: sortedList,chgList
    TYPE(charge_obj) :: chgtemp
    INTEGER,DIMENSION(3) :: npos
    INTEGER :: i,n1,n2,nindex
    ! Assumes that the neighbor does not have weight to a basin that this
    ! boundary point has no flux to
    CALL pbc(npos,chgtemp%npts)
    nindex=(npos(1)-1)*chgtemp%npts(2)*chgtemp%npts(3)+(npos(2)-1)*chgtemp%npts(3)&
      +npos(3)
    DO n1=1,6
      DO n2=1,6
        IF (chgList(nindex)%nbrvol(n2)==sortedList(i)%nbrvol(n1)) THEN
          sortedList(i)%volwgt(n1)=sortedList(i)%nbrflx(n1)*chgList(nindex)%volwgt(n2)
        
        END IF
      END DO
    END DO
  END SUBROUTINE

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
