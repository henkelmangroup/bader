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
      INTEGER, DIMENSION(3) :: pos
      INTEGER ::  volnum !  similar to bader volnum
    END TYPE

    TYPE rvert_obj
      REAL(q2), DIMENSION(3) :: ot3
      REAL(q2) :: angle
    END TYPE

    PUBLIC :: weight_obj
    PUBLIC :: bader_weight_calc, quick_sort

  CONTAINS

  SUBROUTINE bader_weight_calc(bdr, ions, chgval, opts)

    TYPE(weight_obj), ALLOCATABlE, DIMENSION(:) :: chgList
    TYPE(weight_obj) :: tempwobj
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chgval, chgtemp
    TYPE(ions_obj) :: ions, ionstemp
    TYPE(options_obj) :: opts
    INTEGER :: totalLength
    INTEGER :: i, n1, n2, n3, walker, nnvect
    REAL(q2) :: vol, tsum, tw, temp
    INTEGER :: t1, t2, cr, cm, nabove, tbasin, m
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: indList
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: vect, neigh
    INTEGER, ALLOCATABLE, DIMENSION(:) :: numbelow, basin, above
    INTEGER, DIMENSION(3) :: p
    REAL(q2), ALLOCATABLE, DIMENSION(:) :: alpha
    REAL(q2), ALLOCATABLE, DIMENSION(:) :: t, w
    REAL(q2), ALLOCATABLE, DIMENSION(:,:) :: prob
    LOGICAL :: boundary

    ! above : an array of idexes of cells with higher rho
    ! tbasin : a temperary basin entry
    ! basin : an array used to store basin info of all neighbors with higher
    ! rho
    ! numbelow : 

    bdr%nvols = 0

    CALL ws_voronoi(ions, nnvect, vect, alpha)
    CALL SYSTEM_CLOCK(t1, cr, cm)

    IF (opts%ref_flag) THEN
      CALL read_charge_ref(ionstemp, chgtemp, opts)
      ! Assert that chgval and chgtemp are of the same size
      IF ((chgval%npts(1) /= chgtemp%npts(1)) .OR. &
          (chgval%npts(2) /= chgtemp%npts(2)) .OR. &
          (chgval%npts(3) /= chgtemp%npts(3))) THEN
         WRITE(*,'(/,2x,A,/)') 'The dimensions of the primary and reference charge densities must be the same, stopping.'
         STOP
      END IF
    ELSE
      chgtemp = chgval
    END IF


    totalLength = chgtemp%npts(1)*chgtemp%npts(2)*chgtemp%npts(3)
    bdr%tol = opts%badertol
    ALLOCATE (numbelow(totalLength))
    ALLOCATE (w(totalLength))
    ALLOCATE (neigh(totalLength,nnvect))
    ALLOCATE (prob(totalLength,nnvect))
    ALLOCATE (basin(totalLength))
    ALLOCATE (chgList(totalLength))
    ALLOCATE (bdr%volnum(chgtemp%npts(1), chgtemp%npts(2), chgtemp%npts(3)))
    ALLOCATE (indList(chgtemp%npts(1), chgtemp%npts(2), chgtemp%npts(3)))
    bdr%bdim = 64 ! will expand as needed
    ALLOCATE (bdr%volpos_lat(bdr%bdim, 3)) ! will be expanded as needed
    ! Find vacuum points
    IF (opts%vac_flag) THEN
      DO n1=1,chgval%npts(1)
        DO n2=1,chgval%npts(2)
          DO n3=1,chgval%npts(3)
            IF (ABS(rho_val(chgval,n1,n2,n3)/vol) <= opts%vacval) THEN
               bdr%volnum(n1,n2,n3) = -1
               bdr%vacchg = bdr%vacchg + chgval%rho(n1,n2,n3)
               bdr%vacvol = bdr%vacvol + 1
            END IF
          END DO
        END DO
      END DO
    END IF
    bdr%vacchg = bdr%vacchg/REAL(chgval%nrho,q2)
    bdr%vacvol = bdr%vacvol*vol/chgval%nrho
    walker = 1
    DO n1=1, chgtemp%npts(1)
      DO n2=1, chgtemp%npts(2)
        DO n3=1, chgtemp%npts(3)
          chgList(walker)%rho = chgtemp%rho(n1,n2,n3)
          chgList(walker)%pos = (/n1,n2,n3/)
          chgList(walker)%volnum = 0
          walker = walker + 1
        END DO
      END DO
    END DO
    WRITE(*,'(/,2x,A,$)') 'SORTING CHARGE VALUES ... '
    CALL quick_sort(chgList)
    ! this is an ascending loop. so future loop starts from the end
    ! or flipping it should be faster
    DO walker = 1, totalLength/2
      tempwobj = chgList(totalLength + 1 - walker)
      chgList(totalLength + 1 - walker) = chgList(walker)
      chgList(walker) = tempwobj
    END DO
    DO walker = 1, totalLength
      indList(chgList(walker)%pos(1), &
              chgList(walker)%pos(2), &
              chgList(walker)%pos(3)) = walker
    END DO
    WRITE(*,'(A)'), 'DONE'

    ! first loop, deal with all interior points
    WRITE(*,'(2x,A,$)') 'CALCULATING FLUX ... '
    DO n1 = 1, totalLength
      basin(n1) = 0
      numbelow(n1) = 0
      nabove = 0
      tsum = 0
      ALLOCATE (t(nnvect))
      ALLOCATE (above(nnvect))
      DO n2 = 1, nnvect
        p = chgList(n1)%pos + vect(n2,:)
        CALL pbc(p, chgtemp%npts)
        m = indList(p(1),p(2),p(3))
       IF (m < n1 ) THEN ! point p has higher rho
         nabove = nabove + 1
         above(nabove) = m 
         t(nabove) = alpha(n2)*(chgList(m)%rho-chgList(n1)%rho)
         tsum = tsum + t(nabove)
       END IF
      END DO
      IF (nabove == 0) THEN ! maxima
        bdr%bnum = bdr%bnum + 1
        bdr%nvols = bdr%nvols + 1
        basin(n1) = bdr%nvols
        bdr%volnum(chgList(n1)%pos(1),chgList(n1)%pos(2),chgList(n1)%pos(3)) = bdr%nvols 
        IF (bdr%bnum >= bdr%bdim) THEN
          CALL reallocate_volpos(bdr,bdr%bdim*2)
        END  IF
        DEALLOCATE(t)
        DEALLOCATE(above)
        bdr%volpos_lat(bdr%bnum,:) = REAL(p,q2)
        CYCLE
      END IF
      tbasin = basin(above(1))
      boundary = .FALSE.
      DO n2=1,nabove
        IF (basin(above(n2))/=tbasin .OR. tbasin==0) THEN 
          boundary = .TRUE.
        END IF
      END DO 
      IF (boundary) THEN ! boundary
        basin(n1) = 0
        temp = 0
        DO n2 = 1, nabove
          m = above(n2)
          numbelow(m) = numbelow(m)+1
          neigh(m,numbelow(m)) = n1
          prob(m,numbelow(m)) = t(n2) / tsum
          IF (prob(m,numbelow(m)) > temp) THEN
            temp = prob(m,numbelow(m))
            bdr%volnum( &
              chgList(n1)%pos(1),chgList(n1)%pos(2),chgList(n1)%pos(3)) = &
              bdr%volnum( chgList(m)%pos(1),chgList(m)%pos(2),chgList(m)%pos(3) )
          END IF
        END DO
      ELSE ! interior
        basin(n1) = tbasin
        bdr%volnum(chgList(n1)%pos(1),chgList(n1)%pos(2),chgList(n1)%pos(3)) = tbasin
      END IF
      DEALLOCATE(t)
      DEALLOCATE(above)
    END DO
    ! restore chglist rho to values from chgval
    DO walker = 1, totalLength
      n1 = chgList(walker)%pos(1)
      n2 = chgList(walker)%pos(2)
      n3 = chgList(walker)%pos(3)
      chgList(walker)%rho=chgval%rho(n1,n2,n3)
    END DO
    WRITE(*,'(A)'), 'DONE'

    WRITE(*,'(2x,A,$)') 'INTEGRATING CHARGES ... '
    ALLOCATE (bdr%volchg(bdr%nvols))
    ALLOCATE (bdr%ionvol(bdr%nvols))
    DO n1 = 1,bdr%nvols
      bdr%volchg(n1) = 0
      bdr%ionvol(n1) = 0
    END DO
    ! bdr%volnum is written here during integration. so that each cell is
    ! assigned to the basin where it has most of the weight to. This should not
    ! affect the result of the integration.
    temp = 0
    DO n1 = 1, bdr%nvols
      DO n2 = 1, totalLength
        IF (basin(n2) == n1) THEN
          w(n2) = 1
        ELSE
          w(n2) = 0
        END IF
      END DO
      DO n2 = 1,totalLength
        tw = w(n2)
        IF (tw /= 0) THEN
          DO walker = 1, numbelow(n2)
            w(neigh(n2, walker)) = w(neigh(n2, walker)) + prob(n2, walker) * tw
          END DO
          bdr%volchg(n1) = bdr%volchg(n1) + tw * chgList(n2)%rho
          bdr%ionvol(n1) = bdr%ionvol(n1) + tw
        END IF
     
      END DO
    END DO
    bdr%volchg = bdr%volchg / REAL(chgval%nrho,q2)
    WRITE(*,'(A)'), 'DONE'

    vol = matrix_volume(ions%lattice)
    vol = vol / chgtemp%nrho
    bdr%ionvol = bdr%ionvol * vol

    CALL SYSTEM_CLOCK(t2, cr, cm)
    WRITE(*,'(/,1A12,1F10.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

    DO walker = 1, totalLength
      IF (bdr%volnum(chgList(walker)%pos(1), chgList(walker)%pos(2), &
          chgList(walker)%pos(3)) == 0) THEN
        PRINT *,'still zero'
      END IF
    END DO
    ALLOCATE (bdr%nnion(bdr%nvols), bdr%iondist(bdr%nvols), bdr%ionchg(ions%nions))
    ALLOCATE (bdr%volpos_dir(bdr%nvols, 3))
    ALLOCATE (bdr%volpos_car(bdr%nvols, 3))
    DO i = 1, bdr%nvols
      bdr%volpos_dir(i,:) = lat2dir(chgtemp, bdr%volpos_lat(i,:))
      bdr%volpos_car(i,:) = lat2car(chgtemp, bdr%volpos_lat(i,:))
    END DO
    CALL assign_chg2atom(bdr, ions, chgval)

    DEALLOCATE (numbelow)
    DEALLOCATE (w)
    DEALLOCATE (neigh)
    DEALLOCATE (prob)
    DEALLOCATE (chgList)
    DEALLOCATE (indList)
    DEALLOCATE (basin)

  END SUBROUTINE bader_weight_calc


  ! modified version to sort rvert
  RECURSIVE SUBROUTINE rvert_sort(list,length)

    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified by Alan Miller to include an associated integer array which gives
    ! the positions of the elements in the original order.

    IMPLICIT NONE
    TYPE(rvert_obj), DIMENSION (:), INTENT(IN OUT) :: list
    ! Local variable
    INTEGER :: i, length

    CALL rvert_sort_1(1, length)

    CONTAINS

      RECURSIVE SUBROUTINE rvert_sort_1(left_end, right_end)

        INTEGER, INTENT(IN) :: left_end, right_end
        !  Local variables
        INTEGER :: i, j, itemp
        TYPE(rvert_obj) :: temp
        REAL(q2) :: reference
        INTEGER, PARAMETER :: max_simple_sort_size = 6

        IF (right_end < left_end + max_simple_sort_size) THEN
          ! Use interchange sort for small lists
          CALL rinterchange_sort(left_end, right_end)
        ELSE
          ! Use partition ("quick") sort
          reference = list((left_end + right_end)/2)%angle
          i = left_end - 1
          j = right_end + 1

          DO
            ! Scan list from left end until element >= reference is found
            DO
              i = i + 1
              IF (list(i)%angle >= reference) EXIT
            END DO
            ! Scan list from right end until element <= reference is found
            DO
              j = j - 1
              IF (list(j)%angle <= reference) EXIT
            END DO

            IF (i < j) THEN
              ! Swap two out-of-order elements
              temp = list(i)
              list(i) = list(j)
              list(j) = temp
            ELSE IF (i == j) THEN
              i = i + 1
              EXIT
            ELSE
              EXIT
            END IF
          END DO

          IF (left_end < j) CALL rvert_sort_1(left_end, j)
          IF (i < right_end) CALL rvert_sort_1(i, right_end)
        END IF
      
      END SUBROUTINE rvert_sort_1

      SUBROUTINE rinterchange_sort(left_end, right_end)

        INTEGER, INTENT(IN) :: left_end, right_end
        !     Local variables
        INTEGER :: i, j, itemp
        TYPE(rvert_obj) :: temp

        DO i = left_end, right_end - 1
          DO j = i + 1, right_end
            IF (list(i)%angle > list(j)%angle) THEN
              temp = list(i)
              list(i) = list(j)
              list(j) = temp
            END IF
          END DO
        END DO
      
      END SUBROUTINE rinterchange_sort
    
  END SUBROUTINE rvert_sort

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
    ! Local variables
    INTEGER :: i, length

    length = SIZE(list)
    CALL quick_sort_1(1, length)

    CONTAINS

      RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

        INTEGER, INTENT(IN) :: left_end, right_end

        !     Local variables
        INTEGER :: i, j, itemp
        TYPE(weight_obj) :: temp
        REAL(q2) :: reference
        INTEGER, PARAMETER :: max_simple_sort_size = 6

        IF (right_end < left_end + max_simple_sort_size) THEN
          ! Use interchange sort for small lists
          CALL interchange_sort(left_end, right_end)

        ELSE
          ! Use partition ("quick") sort
          reference = list((left_end + right_end)/2)%rho
          i = left_end - 1
          j = right_end + 1

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
              temp = list(i)
              list(i) = list(j)
              list(j) = temp
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

        INTEGER, INTENT(IN) :: left_end, right_end
        ! Local variables
        INTEGER :: i, j, itemp
        TYPE(weight_obj) :: temp

        DO i = left_end, right_end - 1
          DO j = i + 1, right_end
            IF (list(i)%rho > list(j)%rho) THEN
              temp = list(i)
              list(i) = list(j)
              list(j) = temp
            END IF
          END DO
        END DO

      END SUBROUTINE interchange_sort

  END SUBROUTINE quick_sort

!---------------------------
! code from Dallas Trinkle
!---------------------------
  SUBROUTINE ws_voronoi(ions, nnvect, vect, alpha)
    TYPE(ions_obj) :: ions
    TYPE(rvert_obj),ALLOCATABLE,DIMENSION(:) :: rvert
    REAL(q2), ALLOCATABLE,DIMENSION(:,:) :: R
    REAL(q2), ALLOCATABLE,DIMENSION(:) :: alph, alpha
    REAL(q2), DIMENSION(3,3) :: Rdot,Rinv
    REAL(q2), DIMENSION(3) :: nv, tempR, R2, rx, ry, tempR2, tempR3
    REAL(q2) :: detR, tol, temp, rdRn
    INTEGER, ALLOCATABLE,DIMENSION(:,:) :: nvect, vect
    INTEGER :: nv1, nv2, nv3, i, n
    INTEGER :: Nrange, Nneigh, maxvert, nnvect, nvert
    INTEGER :: prunN, zeroarea, nvi
    LOGICAL :: ic
    tol = 1E-8
    Nrange = 3
    Nneigh = (2*Nrange + 1)**3 - 1
    ALLOCATE (nvect(Nneigh,3))
    ALLOCATE (R(Nneigh,5))
    Nneigh = 0
    DO nv1 = -Nrange, Nrange
      DO nv2 = -Nrange,Nrange
        DO nv3 = -Nrange,Nrange
          IF (nv1==0 .AND. nv2==0 .AND. nv3==0) THEN
            CYCLE
          ELSE
            Nneigh = Nneigh + 1
            nvect(Nneigh,1) = nv1
            nvect(Nneigh,2) = nv2
            nvect(Nneigh,3) = nv3
            nv = (/nv1,nv2,nv3/)
            !CALL matrix_vector(ions%lattice,nv,tempR)
            CALL mult_vect(ions%lattice, nv, tempR)
            DO i = 1, 3
              R(Nneigh,i) = tempR(i)
            END DO
!            R(Nneigh,1) = temp
!            R(Nneigh,1) = R(Nneigh,2)
!            R(Nneigh,2) = temp
            R(Nneigh,4) = 0.5*(tempR(1)**2 + tempR(2)**2 + tempR(3)**2)
            R(Nneigh,5) = 0
          END IF
        END DO
      END DO
    END DO
    DO nv1 = 1, Nneigh
      DO nv2 = 1, 3
        tempR(nv2) = R(Nneigh - nv1 + 1, nv2)*0.5
      END DO
      CALL INCELL(tempR, Nneigh, R, ic)
      IF (.NOT. ic) THEN
        R(nv1,5) = 1
      END IF
    END DO
    CALL TRIMARRAY(R, nvect, Nneigh)
    ! next step is to find vortex points
    maxvert = (Nneigh-2)*(Nneigh-4)
    ALLOCATE (rvert(maxvert))
    ALLOCATE (alph(Nneigh))
    nnvect = Nneigh
    DO n = 1, Nneigh
      nvert = 1
      DO nv1 = 1, 3
        Rdot(1,nv1) = R(n,nv1)
      END DO
      R2(1) = R(N,4)
      DO nv1 = 1, Nneigh
        DO i = 1, 3
          Rdot(2,i) = R(nv1,i)
        END DO
        R2(2) = R(nv1,4)
        DO nv2 = nv1 + 1, Nneigh
          DO i = 1, 3
            Rdot(3,i) = R(nv2,i)
          END DO
          R2(3) = R(nv2,4)
          CALL det(Rdot,detR)
          IF (ABS(detR) >= tol) THEN
            CALL matrix_3x3_inverse(Rdot, Rinv)
            Rinv = Rinv*detR
            tempR = (/0,0,0/)
            CALL matrix_vector(Rinv, R2, tempR)
            DO i = 1, 3
              rvert(nvert)%ot3(i) = tempR(i)/detR
            END DO
            CALL INCELL(tempR/detR, Nneigh, R, ic)
            IF (ic) nvert = nvert + 1
          END IF
        END DO
      END DO
      zeroarea = 0
      DO nvi = 1, nvert
        IF (ABS(rvert(nvi)%ot3(1)**2 + rvert(nvi)%ot3(2)**2 + rvert(nvi)%ot3(3)**2 - 0.5*R(n,4)) < tol) & 
          zeroarea = 1
      END DO
      IF (zeroarea == 1 .OR. nvert == 0) THEN
        alph(n) = 0
        nnvect = nnvect - 1
        CYCLE
      END IF
      DO i = 1, 3
        rx(i) = rvert(1)%ot3(i)
      END DO
      rdRn= (rx(1)*R(n,1) + rx(2)*R(n,2) + rx(3)*R(n,3) ) /  & 
            (R(n,1)**2 + R(n,2)**2 + R(n,3)**2)
      DO i = 1, 3
        rx(i) = rx(i) - rdRn*R(n,i)
      END DO
      rdRn = SQRT(rx(1)**2 + rx(2)**2 + rx(3)**2)
      DO i = 1, 3
        rx(i) = rx(i)*1/rdRn
        tempR(i) = R(n,i)
      END DO
      CALL cross_product(tempR, rx, ry)
      rdRn = SQRT(ry(1)**2 + ry(2)**2 + ry(3)**2)
      DO i = 1, 3
        ry(i) = ry(i)*1/rdRn
      END DO
      DO nvi = 1, nvert - 1
        rvert(nvi)%angle = ATAN2( (rvert(nvi)%ot3(1)*ry(1) + rvert(nvi)%ot3(2)*ry(2) + rvert(nvi)%ot3(3)*ry(3)), &
                       (rvert(nvi)%ot3(1)*rx(1) + rvert(nvi)%ot3(2)*rx(2) + rvert(nvi)%ot3(3)*rx(3)) )
      END DO
      CALL rvert_sort(rvert, nvert - 1)
      alph(n) = 0
      DO nvi = 1, nvert - 1
        DO i = 1, 3
          tempR(i) = rvert(nvi)%ot3(i)
          tempR2(i) = rvert(MOD(nvi, (nvert - 1)) + 1)%ot3(i)
          tempR3(i) = R(n, i)
        END DO
        CALL triple_product(tempR, tempR2, tempR3, temp)
        alph(n) = alph(n) + temp
      END DO
      alph(n) = alph(n)*0.25/R(n,4)
      IF (ABS(alph(n)) < tol) THEN
        alph(n) = 0
        nnvect = nnvect - 1
      END IF
    END DO 
    DO n = 1, Nneigh
    END DO
    ALLOCATE (vect(nnvect,3))
    ALLOCATE (alpha(nnvect))
    nvi = 1
    DO n = 1, Nneigh
      IF (alph(n) /=0 ) THEN
        DO i = 1, 3
          vect(nvi,i) = nvect(n,i)
        END DO
        alpha(nvi) = alph(n)
        nvi = nvi + 1
      END IF
    END DO 

  END SUBROUTINE ws_voronoi

  ! this subroutine takes a given array and trims its size, deleting unwanted entries
  SUBROUTINE trimarray(R, nvect, Nneigh)
    REAL(q2), ALLOCATABLE, DIMENSION(:,:) :: R, tR
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nvect, tnvect
    INTEGER :: i, Nneigh, j, nNneigh
    INTEGER :: n1
    nNneigh = 0
    DO i = 1, Nneigh
      IF ( R(i,5) /= 1 ) THEN
        nNneigh = nNneigh + 1
      END IF
    END DO
    ALLOCATE (tR(nNneigh,5))
    ALLOCATE (tnvect(nNneigh,3))
    j = 0
    DO i = 1, Nneigh
      IF (R(i,5) == 1) THEN
        CYCLE
      ELSE
        j = j + 1
        DO n1 = 1, 3
          tR(j,n1) = R(i,n1)
          tnvect(j,n1) = nvect(i,n1)
        END DO 
        tR(j,4) = R(i,4)
      END IF
    END DO
    DEALLOCATE (R)
    DEALLOCATE (nvect)
    ALLOCATE (R(nNneigh,5))
    ALLOCATE (nvect(nNneigh,3))
    DO i = 1, nNneigh
      DO j = 1, 3
        R(i,j) = tR(i,j)
        nvect(i,j) = tnvect(i,j)
      END DO
      R(i,4) = tR(i,4)
      R(i,5) = 0

    END DO
    DEALLOCATE (tR)
    DEALLOCATE (tnvect)
    Nneigh = nNneigh
  END SUBROUTINE TRIMARRAY

  SUBROUTINE  INCELL(r,Nneigh,Ra,ic)
    REAL(q2), DIMENSION(:,:) :: Ra
    REAL(q2), DIMENSION(3) :: r
    INTEGER :: Nneigh, i
    LOGICAL :: ic
    REAL(q2) :: tol, dot
    tol = 1E-8
    ic = .TRUE.
    DO i = 1, Nneigh 
      dot = Ra(i,1)*r(1) + Ra(i,2)*r(2) + Ra(i,3)*r(3)
      IF (dot > Ra(i,4) + tol) THEN 
        ic = .FALSE. 
        RETURN
      END IF
    END DO
  END SUBROUTINE INCELL 

  END MODULE
