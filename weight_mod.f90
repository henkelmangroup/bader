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
      INTEGER(KIND=8) :: x,y,z
    END TYPE

    PUBLIC :: weight_obj
    PUBLIC :: bader_weight_calc

  CONTAINS 

  SUBROUTINE bader_weight_calc(bdr,chg)

    ! chgList has 4 slots. 1st is rho, 2 3 4 are the coordinate indices.
    TYPE(weight_obj),ALLOCATABlE,DIMENSION(:) :: chgList, sortedList
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    TYPE(weight_obj) :: wt
    INTEGER(KIND=8) :: totalLength, i, j, k, l, n1, n2, n3

    totalLength = chg%npts(1)*chg%npts(2)*chg%npts(3)
    PRINT *, totalLength
    ALLOCATE (chgList(totalLength))
    ALLOCATE (sortedList(totalLength))
    ! first merge sort
    DO n1=1,chg%npts(1)
      DO n2=1,chg%npts(2)
        DO n3=1,chg%npts(3)
          chgList(n1+n2+n3-2)%rho = chg%
        END DO
      END DO
    END DO
    CALL merge_sort(chgList,sortedList)

    END SUBROUTINE bader_bader_calc

!------------------------------------------------------------------------------------!
! merge_sort: sort the array using the 1st element of it. In desending order
!------------------------------------------------------------------------------------!

    SUBROUTINE merge_sort(A,B)

      TYPE(weight_obj),DIMENSION(:) :: A,B
      TYPE(weight_obj),ALLOCATABLE,DIMENSION(:) :: C,D,E,F
      REAL(q2) :: totalStep
      INTEGER :: i,j,k,l,walker

      totalStep = ceiling(Log(Size(A)*1.0_q2)/Log(2.0_q2))
!      PRINT *, totalStep
!      PRINT *, ceiling(totalStep)
      walker = 0
      DO i=1, totalStep
        ALLOCATE (C(2**i))
        ALLOCATE (D(2**i))
        DO j=1,i*2
          C(j) = A(walker+j)
        END DO
        DEALLOCATE (C,D)
      END DO

      ! Each time only 2 arrays will be operated on. The length of the arrays
      ! are 2**n, where n is the number of steps. The first loop over the entire
      ! array to sort, the length is 2. The second loop, length is 4. The last
      ! loop, each working aray is of half the length of the entire array to be
      ! sorted. There are 2 senarios, the array may be of odd length or even
      ! length.   

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
