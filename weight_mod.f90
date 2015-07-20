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
    INTEGER :: t1,t2,cr,cm
    LOGICAL :: nbb ! never boundry before
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: indList
    nbb=.TRUE.
    tempw=0
    totalE=0
    bdr%nvols=0
    CALL SYSTEM_CLOCK(t1,cr,cm)
    IF (opts%ref_flag) THEN
      CALL read_charge_ref(ionstemp,chgtemp,opts)
    ELSE
      chgtemp = chgval
    END IF
    totalLength=chgtemp%npts(1)*chgtemp%npts(2)*chgtemp%npts(3)
    ALLOCATE (sortedList(totalLength)) ! sortedList should have the correct size
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
    i=0
    DO walker=1,totalLength
      indList(chgList(walker)%pos(1), &
              chgList(walker)%pos(2), &
              chgList(walker)%pos(3) &
             )=walker
      CALL interiors(bdr,chgtemp,chgList(walker),walker,indList)
      
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
    CALL SYSTEM_CLOCK(t2,cr,cm)
    PRINT *,'Time elapsed:',(t2-t1)/REAL(cr,q2),'seconds'
  END SUBROUTINE bader_weight_calc

  !-----------------------------------------------------------------------------------!
  ! interiors: this subroutine loops through all grid points, finding out which
  ! are interiors and which are boundries, and number of basins.
  !-----------------------------------------------------------------------------------!
  SUBROUTINE interiors(bdr,chgtemp,cLW,walker,indList) ! chgListWalker
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: indList
    TYPE(weight_obj) :: cLW
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chgtemp,chgval
    INTEGER :: length,i,n1,n2,n3,temp,tempnvol,tempvolnum,temp2,j,walker
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
      cLW%volnum=0
      bdr%volnum(cLW%pos(1),cLW%pos(2),cLW%pos(3))=0
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
    TYPE(weight_obj),DIMENSION(:) :: chgList
    IF (chgList(walker)%isInterior) THEN
      bdr%volchg(chgList(walker)%volnum)=bdr%volchg(chgList(walker)%volnum)+chgList(walker)%rho
    ELSE
      CALL flux_weight(bdr,chgtemp,chgList,walker,ions,indList)
    END IF
  END SUBROUTINE boundries
!
!-------------------------------------------------------------------!
! flux_weight: It takes the current position which is a boundry pooint, looks at
! its neighbors to find flx and weight to all basins.
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
    counter=0
    DO n1=-1,1
      DO n2=-1,1
        DO n3=-1,1
          IF (n1**2+n2**2+n3**2 /=1) CYCLE
          counter=counter+1
          p=chgList(walker)%pos+(/n1,n2,n3/)
          CALL pbc(p,chgtemp%npts)

          IF (chgtemp%rho(p(1),p(2),p(3))>chgList(walker)%rho) THEN
            flx(counter)=nom(counter)/denom
            ! interior neighbor and boundry neighbor needs to be treated
            ! separately
            IF (bdr%volnum(p(1),p(2),p(3))/=0) THEN !interior neighbor
              chgList(walker)%volwgt(bdr%volnum(p(1),p(2),p(3)))= & 
                chgList(walker)%volwgt(bdr%volnum(p(1),p(2),p(3)))+ flx(counter)
            ELSE
              ! need to find out that point's position in chgList
              nein=indList(p(1),p(2),p(3))
              DO n4=1,bdr%nvols
                chgList(walker)%volwgt(n4)=chgList(walker)%volwgt(n4)+ &
                  chgList(nein)%volwgt(n4)*flx(counter)
              END DO
            END IF
          END IF          
        END DO
      END DO
    END DO
    DO n1=1,bdr%nvols
      bdr%volchg(n1)=bdr%volchg(n1)+chgList(walker)%rho*chgList(walker)%volwgt(n1)
      temp=temp+chgList(walker)%volwgt(n1)
    END DO
  END SUBROUTINE flux_weight


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
