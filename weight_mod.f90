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
    TYPE(charge_obj) :: chg
    TYPE(weight_obj) :: wt
    TYPE(ions_obj) :: ions
    TYPE(options_obj) :: opts
    INTEGER(KIND=8) :: totalLength
    INTEGER :: i, j, k, l, n1, n2, n3, walker
        
!    totalLength=chg%npts(1)*chg%npts(2)*chg%npts(3)
    totalLength = 64 ! for testing purpose only
    PRINT *, totalLength
    ALLOCATE (chgList(totalLength))
    ALLOCATE (sortedList(totalLength))
    ! first merge sort   
    ! the lines below are commented out so that the program does not actually
    ! run using the CHGCAR, but the mannualy put in array A
!    walker=1  
!    DO n1=1, chg%npts(1)
!      DO n2=1, chg%npts(2)
!        DO n3=1,chg%npts(3)
!          chgList(walker)%rho=rho_val(chg,n1,n2,n2)
!          chgList(walker)%x=n1
!          chgList(walker)%y=n2
!          chgList(walker)%z=n3
!          walker=walker+1
!        END DO
!      END DO
!    END DO
    ! Here begins the manual array
    PRINT *,'PRINTING THE CONTENT OF chgList'
    walker=0
    DO n1=1, 4
      DO n2=1,4
        DO n3=1,4
        walker=walker+1
        chgList(walker)%rho = n1*n2*n3
        chgList(walker)%x = n1
        chgList(walker)%y = n2
        chgList(walker)%z = n3
        PRINT *, 'chgList(',walker,') is, ',chgList(walker)%rho
        PRINT *, chgList(walker)%x, chgList(walker)%y, chgList(walker)%z
        END DO
      END DO
    END DO
    walker = 0
    CALL merge_sort(chgList,sortedList)
    DO n1=1,totalLength
      PRINT *, 'n1', n1
      PRINT *, sortedList(n1)%rho
      PRINT *, sortedList(n1)%x,sortedList(n1)%y,sortedList(n1)%z
      PRINT *, ' '
    END DO
    END SUBROUTINE bader_weight_calc


!------------------------------------------------------------------------------------!
! merge_sort: sort the array using the 1st element of it. In desending order
!------------------------------------------------------------------------------------!

    SUBROUTINE merge_sort(A,B)
      TYPE(weight_obj),DIMENSION(:) :: A
      TYPE(weight_obj),INTENT(OUT),DIMENSION(:) :: B
      TYPE(weight_obj) :: tempw
      REAL(q2) :: totalStep, tempq2
      INTEGER :: i, j, k, l, walker, n1, n2, n3, length, tempint
      INTEGER :: loopStep ! the steps it takes to go through A once
      totalStep = CEILING(LOG(SIZE(A)*1.0_q2)/LOG(2.0_q2))
      walker = 0
      DO i=1, totalStep
        PRINT *, 'i is ', i
        length=2**(i-1)
        loopStep=CEILING(SIZE(A)/(2.*length))
        DO j=1,loopStep
          PRINT *, 'j is ',j
          PRINT *, 'ACCESSING',(j-1)* length*2 ,' to ',j*length*2
          CALL sort(A,B,length,i,j)
        END DO
      END DO
      B=A
      PRINT *, 'ALL OVER'

      ! Each time only 2 arrays will be operated on. The length of the arrays
      ! are 2**n, where n is the number of steps. The first loop over the entire
      ! array to sort, the length is 2. The second loop, length is 4. The last
      ! loop, each working aray is of half the length of the entire array to be
      ! sorted. There are 2 senarios, the array may be of odd length or even
      ! length.   

    END SUBROUTINE
 
    SUBROUTINE sort(A,B,length,i,j)
    ! A is the original chgList, B is partially sorted chgList, it will be
    ! returned to merge_sort. After merge_sort loop finishes, both B will be
    ! sortedList. length is the length of array that partially sorts A. The
    ! length varies every time i changes in merge_sort
      TYPE(weight_obj),DIMENSION(:) :: A
      TYPE(weight_obj),INTENT(OUT),DIMENSION(:) :: B
      INTEGER,INTENT(IN) :: length,i,j
      TYPE(weight_obj),DIMENSION(length) :: C,D
      TYPE(weight_obj) :: tempw
      INTEGER(KIND=8) :: walker,n1,n2,n3
      walker=0
      DO n1=1,length
        PRINT *, 'getting from A:', (j-1)*length*2+n1,'and',(2*j-1)*length+n1
        C(n1)=A((j-1)*length*2+n1)
        D(n1)=A((2*j-1)*length+n1)
!        PRINT *, 'C(',n1,') is ', C(n1)%rho
!        PRINT *, C(n1)%x,C(n1)%y,C(n1)%z
!        PRINT *, 'D(',n1,') is ', D(n1)%rho
!        PRINT *, D(n1)%x,D(n1)%y,D(n1)%z
      END DO
      PRINT *, 'FINISHED READING'
      n1=1
      n2=1
      walker=1
      DO WHILE(n1<=length .AND. n2<=length)
        IF (C(n1)%rho>=D(n2)%rho) THEN
          PRINT *, 'C',n1,'is greater than D',n2
          A((j-1)*length*2+walker)=C(n1)
          walker=walker+1
          ! Existing order in A is correct. Check the next element in C agains
          ! the same element in D
          IF (n1==length) THEN ! The last element in C is bigger than D(n2). 
          ! The rest of spaces in A should be written with the rest of the terms
          ! in D 
            DO n3=walker,2*length
              A((j-1)*length*2+walker)=D(n2)
              n2=n2+1
              walker=walker+1
            END DO
            STOP
          END IF
          n1=n1+1
        ELSE
          PRINT *, 'D',n2,'is greater than C',n2
          ! D(n2) is written at the place where C(n1) sits 
          A((j-1)*length*2+walker)=D(n2)
          walker=walker+1
          IF (n2==length) THEN ! The last element in D is bigger than C(n1)
          ! The rest of spaces in A should be weitten with the rest of the terms
          ! in C
            DO n3=n1+n2,length*2
              PRINT *, 'ENTERED INNER LOOP'
              A((j-1)*length*2+walker)=C(n1)
              n1=n1+1
              walker=walker+1
            END DO
          END IF
          n2=n2+1
        END IF
        PRINT *, 'LOOPED ONCE'
      END DO
      B = A

    END SUBROUTINE

!-----------------------------------------------------------------------------------!
! Weight method
!-----------------------------------------------------------------------------------!

!  SUBROUTINE refine_weights(chgval, bdr, p)
!
!    TYPE(bader_obj) :: bdr
!    TYPE(charge_obj) :: chgval
!    INTEGER :: num_edge, n1, n2, n3, d1, d2, d3
!    INTEGER :: i, iter, num_change, mycount
!    INTEGER,DIMENSION(3) :: p, pn
!    REAL(q2) :: sum_top, sum_bottom, length, facet_a, R_
!    REAL(q2) :: new_weight, current_weight, wn
!
!!    write(*,*), ' bnum',bdr%bnum
!
!    DO i = 1,bdr%bnum
!    ! i is the current basin
!
!      num_edge = 0
!
!      ! loop through grid points and assign initial weights
!      DO n1 = 1,chgval%npts(1)
!        DO n2 = 1,chgval%npts(2)
!          DO n3 = 1,chgval%npts(3)
!
!            p = (/n1,n2,n3/)
!            chgval%weight(p(1),p(2),p(3))%w(i) = 0
!
!            ! if p is a vacuum point, skip it
!            !IF (bdr%volnum(n1,n2,n3) == bdr%bnum+1) CYCLE
!
!            !IF ((.NOT. is_vol_edge(bdr,chgval,p)) .AND. &
!            !(bdr%volnum(p(1),p(2),p(3))==i)) THEN
!
!            IF (bdr%volnum(p(1),p(2),p(3)) == i) THEN
!              chgval%weight(p(1),p(2),p(3))%w(i) = 1
!            END IF
!
!            ! count the number of edge points
!            IF (is_vol_edge(bdr, chgval, p) .AND. &
!              &  ((bdr%volnum(p(1),p(2),p(3)) == i) .OR. is_vol_neighbor(bdr, chgval, p, i))) THEN
!            !IF (is_vol_edge(bdr, chgval, p) .AND. &
!            !  &  (bdr%volnum(p(1),p(2),p(3)) == i)) THEN
!
!            !IF (is_vol_edge(bdr, chgval, p)) THEN
!              chgval%weight(p(1),p(2),p(3))%w(i)= 0
!              num_edge = num_edge+1
!            END IF


  END MODULE
