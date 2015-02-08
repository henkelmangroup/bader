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
      INTEGER :: hdn ! number of higher density neighbor cells
      INTEGER,DIMENSION(3) :: pos
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
        END DO
      END DO
    END DO
    PRINT *,'THIS IS IN WEIGHT. CHGVAL AT 90 90 30 IS',chgval%rho(90,90,30)
    PRINT *,'THIS IS IN WEIGHT. CHGVAL AT 91 91 31 IS',chgval%rho(91,91,31)
    PRINT *,'THIS IS IN WEIGHT. CHGTEMP AT 90 90 30 IS',chgtemp%rho(90,90,30)
    PRINT *,'THIS IS IN WEIGHT. CHGTEMP AT 91 91 31 IS',chgtemp%rho(91,91,31)



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
    CALL weight_calc(sortedList,bdr,chgtemp)
 
  END SUBROUTINE bader_weight_calc

  ! lets say the sortedList is correct. Now let's start to walk the list from top
  ! to bottom. Relavent variablees: bdr%volnum. bdr%volpos
  
  
  !-----------------------------------------------------------------------------------!
  ! weight calc 
  !-----------------------------------------------------------------------------------!
  SUBROUTINE weight_calc(sortedList,bdr,chgtemp)
      ! apply periodic boudary condition first. 
    TYPE(weight_obj),DIMENSION(:) :: sortedList
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chgtemp
    INTEGER :: length,i,n1,n2,n3,temp,tempnvol,tempvolnum,temp2
    REAL(q2) :: temprho
    INTEGER,DIMENSION(3) :: npos !neighbor position
    INTEGER,DIMENSION(3) :: xyz ! xyz indexes of current point 
    INTEGER :: nboundary,ninterior
    PRINT *, 'CALLED'
    bdr%nvols=0 
    length=SIZE(sortedList)
    nboundary=0
    ninterior=0
    temp2=0
    temp=SIZE(sortedList)/100
    DO i=1, length
      !PRINT *, 'ENTERED LOOP'
      ! count number of neighbors with higher density
!      i=9
      sortedList(i)%hdn=0
      xyz(1)=sortedList(i)%pos(1)
      xyz(2)=sortedList(i)%pos(2)
      xyz(3)=sortedList(i)%pos(3)
      npos=xyz+(/1,0,0/)
      CALL  &
neighbors(bdr,npos,xyz,chgtemp,sortedList,nboundary,ninterior,i,tempnvol)
      npos=xyz+(/-1,0,0/)
      CALL  &
neighbors(bdr,npos,xyz,chgtemp,sortedList,nboundary,ninterior,i,tempnvol)
      npos=xyz+(/0,1,0/)
      CALL  &
neighbors(bdr,npos,xyz,chgtemp,sortedList,nboundary,ninterior,i,tempnvol)
      npos=xyz+(/0,-1,0/)
      CALL  &
neighbors(bdr,npos,xyz,chgtemp,sortedList,nboundary,ninterior,i,tempnvol)
      npos=xyz+(/0,0,1/)
      CALL  &
neighbors(bdr,npos,xyz,chgtemp,sortedList,nboundary,ninterior,i,tempnvol)
      npos=xyz+(/0,0,-1/)
      CALL  &
neighbors(bdr,npos,xyz,chgtemp,sortedList,nboundary,ninterior,i,tempnvol)

    IF (tempnvol >1 ) THEN
        ! is a boundary point
      nboundary=nboundary+1
    END IF
    IF (tempnvol==1) THEN
      ! is an interior point. 
      bdr%volnum(xyz(1),xyz(2),xyz(3))=tempvolnum
      ninterior=ninterior+1
    END IF  
      IF (i>temp) THEN
        temp2=temp2+1
        PRINT *, temp2 ,'percent done'
        temp=temp+SIZE(sortedList)/100
      END IF


    END DO
    PRINT *, 'the number of boundary points are ', nboundary 
    PRINT *, 'the number of interior points are', ninterior
    PRINT *, 'the total number of point is', 120*120*120
    PRINT *, 'the number of points not identified is', &
      120*120*120-nboundary-ninterior
!    PRINT *, bdr%volnum(10,10,10)
    !PRINT *, sortedList(990)%pos +(/3,3,3/)
  END SUBROUTINE

  SUBROUTINE  neighbors(bdr,npos,xyz,chgtemp,sortedList,&
                nboundary,ninterior,i,tempnvol)
  TYPE(weight_obj),DIMENSION(:) :: sortedList
  TYPE(bader_obj) :: bdr
  TYPE(charge_obj) :: chgtemp  
  INTEGER,DIMENSION(3) :: npos,xyz
  INTEGER :: nboundary,ninterior
  LOGICAL :: fb ! first bigger neighbor
  INTEGER :: temp,tempnvol,tempvolnum,temp2,i
    fb=.TRUE.
    CALL pbc(npos,chgtemp%npts)
!    PRINT *,'BEFORE FIRST IF'
!    PRINT *, 'looking at point',npos(1),npos(2),npos(3)
!    PRINT *, 'value is',chgtemp%rho(npos(1),npos(2),npos(3))
!    PRINT *, 'We are sitting on',sortedList(i)%pos(1),&
!      sortedList(i)%pos(2),sortedList(i)%pos(3)
    IF (sortedList(i)%rho<chgtemp%rho(npos(1),npos(2),npos(3))) THEN
      ! 1 neighbor with higher density
      sortedList(i)%hdn = sortedList(i)%hdn+1
      ! read its volnum, if is the first, record its volnum
!      tempvolnum=bdr%volnum(npos(1),npos(2),npos(3))
!      PRINT *, 'BEFORE SECOND IF'
      IF (fb) THEN
        fb=.FALSE.
        ! The volnum of neighbors need to be kept track of
        tempnvol=1 !
        tempvolnum=bdr%volnum(npos(1),npos(2),npos(3))
        IF (xyz(1)==61 .AND. xyz(2)==61 .AND. xyz(3)==91) THEN
          PRINT *,'IF'
        END IF
!        PRINT *, 'WE ARE AT POINT',xyz(1),xyz(2),xyz(3)
!        PRINT *,'WE ARE LOOKING AT POINT',npos(1),npos(2),npos(3)
!        PRINT *,'bdr%volnum(npos(1),npos(2),npos(3)) IS'
!        PRINT *,bdr%volnum(npos(1),npos(2),npos(3))
!        PRINT *, 'and tempvolnum is'
!        PRINT *, tempvolnum
!        PRINT *,'91 91 31 has volnum', bdr%volnum(91,91,31)
!        PRINT *, 'tempvolnum is', tempvolnum
        ELSE IF(bdr%volnum(npos(1),npos(2),npos(3))/=tempvolnum) THEN
          ! More than one neighbors with higher density has been found
          IF (xyz(1)==61 .AND. xyz(2)==61 .AND. xyz(3)==91) THEN
            PRINT *, 'ELSE IF'
          END IF
          tempnvol=tempnvol+1
!        ELSE 
          ! All bigger neighbors has the same volnum
!          bdr%volnum(npos(1),npos(2),npos(3))=tempvolnum
!          PRINT *, 'ELSE'
      END IF
      IF (xyz(1)==61 .AND. xyz(2)==61 .AND. xyz(3)==91) THEN
        PRINT *, 'SITTING ON 61 61 91, SHOULD BE A BOUNDARY POINT'
        PRINT *, 'LOOKING AT',npos(1),npos(2),npos(3)
        PRINT *, 'IT HAS bdr volnum',bdr%volnum(npos(1),npos(2),npos(3)) 
        PRINT *, 'TEMPNVOL IS', tempnvol       
        PRINT *, 'tempvolnum is', tempvolnum
      END IF
!      PRINT *, 'AFTER SECOND IF'  
    END IF                    
    !PRINT *, 'AFTER FIRST IF'
    !PRINT *,'N1 N2 N3 ARE '
    !PRINT *, n1,n2,n3            
!      PRINT *,'OUT OF DO'
!      PRINT *, 'tempnvol is ', tempnvol
    IF (sortedList(i)%hdn==0) THEN
!        PRINT *,'FOUND A LOCAL MAXIMUM'
      bdr%nvols=bdr%nvols+1
!        PRINT *, 'nvols is now', bdr%nvols
      bdr%volnum(xyz(1),xyz(2),xyz(3))=bdr%nvols
!        PRINT *,'TO POINT',xyz(1),xyz(2),xyz(3)
!        PRINT *, 'it now has volnum',bdr%volnum(xyz(1),xyz(2),xyz(3))
!        PRINT *, 'GIVEN volnum'
!        PRINT *, '91,91,31 has volnum',  bdr%volnum(91,91,31)
    END IF 
!      PRINT *,'tempnvol is ' ,tempnvol 



!    PRINT *, 'i is ',i
!    IF (I>=10) EXIT
      

  END SUBROUTINE



  ! This subroutine is not needed. The needed information is in chgtemp
  SUBROUTINE get_unordered_rho(chgList,p1,p2,p3,temprho,chgtemp)
    INTEGER :: p1,p2,p3,origin_index
    REAL(q2) :: temprho
    TYPE(weight_obj),DIMENSION(:) :: chgList
    TYPE(charge_obj) :: chgtemp
    ! find the entry in chgList using coordinates. 
    origin_index=(p1-1)*chgtemp%npts(2)*chgtemp%npts(3)+chgtemp%npts(3)*(p2-1)+p3
    temprho=chgList(origin_index)%rho
    RETURN
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
