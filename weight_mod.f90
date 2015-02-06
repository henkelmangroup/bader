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
    IMPLICIT NONE
    PRIVATE

    TYPE weight_obj
      REAL(q2) :: rho
      INTEGER(KIND=8) :: x, y, z
    END TYPE

    PUBLIC :: weight_obj
    PUBLIC :: bader_weight_calc

  CONTAINS 

  SUBROUTINE bader_weight_calc(bdr,ions,chg,opts)

    ! chgList has 4 slots. 1st is rho, 2 3 4 are the coordinate indices.
    TYPE(weight_obj),ALLOCATABlE,DIMENSION(:) :: chgList, sortedList
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg,chgtemp,chgval
    TYPE(weight_obj) :: wt
    TYPE(ions_obj) :: ions,ionstemp
    TYPE(options_obj) :: opts
    INTEGER(KIND=8) :: totalLength,temp
    INTEGER :: i, j, k, l, n1, n2, n3, walker
    
    
    IF (opts%ref_flag) THEN
      CALL read_charge_ref(ionstemp,chgtemp,opts)
    ELSE
      chgtemp = chgval
    END IF
     
    totalLength=chg%npts(1)*chg%npts(2)*chg%npts(3)
!    totalLength = 8000 ! for testing purpose only
    ALLOCATE (sortedList(totalLength)) ! sortedList should have the correct size
!    PRINT *,'assigned totalLength', totalLength
    temp=totalLength
    ! The length of the array to be sorted should be of size 2^n. The additional
    ! buffer slots will be empty
    totalLength=2**(CEILING(LOG(temp*1.0_q2)/LOG(2.0_q2)))
!    PRINT *, 'enlarged total length', totalLength
    ALLOCATE (chgList(totalLength))
    ! first merge sort   
    ! the lines below are commented out so that the program does not actually
    ! run using the CHGCAR, but the mannualy put in array A
    walker=1  
    DO n1=1, chgtemp%npts(1)
      DO n2=1, chgtemp%npts(2)
        DO n3=1,chgtemp%npts(3)
          chgList(walker)%rho=rho_val(chgtemp,n1,n2,n2)
          chgList(walker)%x=n1
          chgList(walker)%y=n2
          chgList(walker)%z=n3
          walker=walker+1
        END DO
      END DO
    END DO
    ! Here begins the manual array
!    PRINT *,'PRINTING THE CONTENT OF chgList'
!    walker=1
!    DO n1=1, 20
!      DO n2=1, 20
!        DO n3=1, 20
!        chgList(walker)%rho = n1*n2*n3
!        chgList(walker)%x = n1
!        chgList(walker)%y = n2
!        chgList(walker)%z = n3
!        walker=walker+1
!        END DO
!      END DO
!    END DO
!    PRINT *, 'walker is ', walker
    ! after the last loop is finished, 
    ! The empty slots are filled with junk. Lets fill it with -1
    DO i=walker,totalLength
      chgList(walker)%rho=-1
      chgList(walker)%x = -1
      chgList(walker)%y = -1
      chgList(walker)%z = -1
      walker=walker+1  
    END DO

!    PRINT *, chgList
!    PRINT *, 'flag a'
!    PRINT *,'THIS IS CHGLIST'
!    DO i=1, SIZE(chgList)
!      PRINT *,'entry',i,'value is ', chgList(i)%rho
!    END DO
    CALL quick_sort(chgList)
!    PRINT *, 'THIS IS SORTED LIST'
!    DO i=1,SIZE(chgList)
!      PRINT *,'entry ',i,'value is', chgList(i)%rho
!    END DO
    DO i=1,SIZE(sortedList)
      sortedList(i)=chgList(totalLength-i+1)
      PRINT *, sortedList(i)%rho
    END DO
!    CALL bubble(chgList,totalLength) ! This is even worse
!    DO n1=1,totalLength
!      PRINT *, 'n1', n1
!      PRINT *, sortedList(n1)%rho
!      PRINT *, sortedList(n1)%x,sortedList(n1)%y,sortedList(n1)%z
!      PRINT *, ' '
!    END DO
  END SUBROUTINE bader_weight_calc


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
