  !-----------------------------------------------------------------------------------!
  ! Bader charge density analysis program
  !  Module implementing the weight method by Yu and Trinkle [JCP 134, 064111 (2011)]
  !-----------------------------------------------------------------------------------!

!GH todo
! - check bader tol
! - delete incell
! - replace quicksort
! - fix ws routine

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
      INTEGER ::  volnum  ! similar to bader volnum
    END TYPE

    TYPE rvert_obj
      REAL(q2), DIMENSION(3) :: ot3
      REAL(q2) :: angle
    END TYPE

    PUBLIC :: weight_obj
    PUBLIC :: bader_weight_calc

  CONTAINS

  SUBROUTINE bader_weight_calc(bdr, ions, chgval, opts)

    TYPE(weight_obj), ALLOCATABlE, DIMENSION(:) :: chgList
    TYPE(weight_obj) :: tempwobj
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chgval, chgtemp
    TYPE(ions_obj) :: ions, ionstemp
    TYPE(options_obj) :: opts
    INTEGER :: nPts, i, n1, n2, n3, walker, numVect
    INTEGER :: t1, t2, cr, cm, nabove, tbasin, m
    REAL(q2) :: vol, tsum, tw, temp
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: indList
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: vect, neigh
    INTEGER, ALLOCATABLE, DIMENSION(:) :: numbelow, basin, above
    INTEGER, DIMENSION(3) :: p
    REAL(q2), ALLOCATABLE, DIMENSION(:) :: alpha
    REAL(q2), ALLOCATABLE, DIMENSION(:) :: t, w
    REAL(q2), ALLOCATABLE, DIMENSION(:,:) :: prob
    REAL(q2), DIMENSION(3,3) :: cell
    LOGICAL :: boundary

    ! above : an array of idecies of cells with higher rho
    ! tbasin : a temporary basin entry
    ! basin : an array used to store basin info of all neighbors with higher rho
    ! numbelow : 

    bdr%nvols = 0

    DO i=1,3
      cell(i,:) = ions%lattice(i,:)/chgval%npts(i)
    END DO

    CALL ws_voronoi(cell, numVect, vect, alpha)
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

    WRITE(*,*) "HERE1"
    WRITE(*,*) "numVect: ", numVect

    nPts = chgtemp%npts(1)*chgtemp%npts(2)*chgtemp%npts(3)
!    bdr%tol = opts%badertol  ! delete this?
    ALLOCATE (numbelow(nPts))
    ALLOCATE (w(nPts))
    ALLOCATE (neigh(nPts, numVect))
    ALLOCATE (prob(nPts, numVect))
    ALLOCATE (basin(nPts))
    ALLOCATE (chgList(nPts))
    ALLOCATE (bdr%volnum(chgtemp%npts(1), chgtemp%npts(2), chgtemp%npts(3)))
    ALLOCATE (indList(chgtemp%npts(1), chgtemp%npts(2), chgtemp%npts(3)))
    bdr%bdim = 64 ! will expand as needed
    ALLOCATE (bdr%volpos_lat(bdr%bdim, 3))

    WRITE(*,*) "HERE1.5"

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

    WRITE(*,*) "HERE2"

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
!    CALL quick_sort(chgList)
    CALL sort_weight(nPts, chgList)

    ! GH: fix this nonsense

    ! this is an ascending loop. so future loop starts from the end
    ! or flipping it should be faster
    DO walker = 1, nPts/2
      tempwobj = chgList(nPts + 1 - walker)
      chgList(nPts + 1 - walker) = chgList(walker)
      chgList(walker) = tempwobj
    END DO
    DO walker = 1, nPts
      indList(chgList(walker)%pos(1), &
              chgList(walker)%pos(2), &
              chgList(walker)%pos(3)) = walker
    END DO
    WRITE(*,'(A)'), 'DONE'

    ! first loop, deal with all interior points
    WRITE(*,'(2x,A,$)') 'CALCULATING FLUX ... '
    DO n1 = 1, nPts
      basin(n1) = 0
      numbelow(n1) = 0
      nabove = 0
      tsum = 0
      ALLOCATE (t(numVect))
      ALLOCATE (above(numVect))
      DO n2 = 1, numVect
        p = chgList(n1)%pos + vect(n2,:)
        CALL pbc(p, chgtemp%npts)
        m = indList(p(1),p(2),p(3))
       IF (m < n1 ) THEN ! point p has higher rho
         nabove = nabove + 1
         above(nabove) = m 
         t(nabove) = alpha(n2)*(chgList(m)%rho - chgList(n1)%rho)
         tsum = tsum + t(nabove)
       END IF
      END DO
      IF (nabove == 0) THEN ! maxima
        bdr%bnum = bdr%bnum + 1
        bdr%nvols = bdr%nvols + 1
        basin(n1) = bdr%nvols
        bdr%volnum(chgList(n1)%pos(1), chgList(n1)%pos(2), chgList(n1)%pos(3)) = bdr%nvols 
        IF (bdr%bnum >= bdr%bdim) THEN
          CALL reallocate_volpos(bdr, bdr%bdim*2)
        END IF
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
    DO walker = 1, nPts
      n1 = chgList(walker)%pos(1)
      n2 = chgList(walker)%pos(2)
      n3 = chgList(walker)%pos(3)
      chgList(walker)%rho = chgval%rho(n1,n2,n3)
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
      DO n2 = 1, nPts
        IF (basin(n2) == n1) THEN
          w(n2) = 1
        ELSE
          w(n2) = 0
        END IF
      END DO
      DO n2 = 1, nPts
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
    WRITE(*,'(/,1A12,1F10.2,1A8)') 'RUN TIME: ', (t2-t1)/REAL(cr,q2), ' SECONDS'

    DO walker = 1, nPts
      IF (bdr%volnum(chgList(walker)%pos(1), chgList(walker)%pos(2), &
          chgList(walker)%pos(3)) == 0) THEN
        PRINT *,'still zero'
      END IF
    END DO

    ALLOCATE (bdr%nnion(bdr%nvols))
    ALLOCATE (bdr%iondist(bdr%nvols))
    ALLOCATE (bdr%ionchg(ions%nions))
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


  !-----------------------------------------------------------------------------------!
  !  Source:  adapted from ws_voronoi.H
  !  Author:  D. Trinkle
  !  Date:    2010 December 27
  !  Purpose: Determines the prefactors for computation of flux in Wigner-Seitz
  !           grid cells, based on the Voronoi decomposition.
  !-----------------------------------------------------------------------------------!

  SUBROUTINE ws_voronoi(cell, numVect, vect, alpha)

! Construct a list of the neighboring vectors that define the Wigner-Seitz cell
! and compute the "alpha" factors needed for flux; you multiply the difference
! in densities by alpha, and use this to compute the transition probabilities.

    REAL(q2), DIMENSION(3,3) :: cell
    INTEGER, INTENT(INOUT) :: numVect
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: vect
    REAL(q2), ALLOCATABLE, DIMENSION(:) :: alpha

    TYPE(rvert_obj), ALLOCATABLE, DIMENSION(:) :: rVert
    REAL(q2), ALLOCATABLE, DIMENSION(:,:) :: R, Rtmp
    REAL(q2), ALLOCATABLE, DIMENSION(:) :: alphtmp
    REAL(q2), DIMENSION(3,3) :: Rdot, Rinv
    REAL(q2), DIMENSION(3) :: neighR, R2, Rx, Ry, nv
    REAL(q2) :: detR, tol, temp, rdRn
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nVect
    INTEGER :: nv1, nv2, nv3, nA, nB, i, n, nvi, zeroarea
    INTEGER :: nRange, numVert, maxVert, neigh, numNeigh, maxNeigh
    LOGICAL :: incell

! Generate a list of neighboring vectors that bound the Wigner-Seitz cell.
! Note for future: should precompute this by making a sphere of radius
! with the largest length vector multiplied by, say, 2.

    tol = 1E-8
    nRange = 3
    maxNeigh = (2*nRange + 1)**3 - 1

    ALLOCATE (nVect(maxNeigh,3))
    ALLOCATE (R(maxNeigh,3))
    ALLOCATE (Rtmp(maxNeigh,3))

    neigh = 0
    DO nv1 = -nRange, nRange
      DO nv2 = -nRange,nRange
        DO nv3 = -nRange,nRange
          nv = (/nv1,nv2,nv3/)
          IF (ALL(nv == 0)) CYCLE
          neigh = neigh + 1
          nVect(neigh,:) = nv
!          CALL mult_vect(cell, nv, R(neigh,:))
          R(neigh,:) = MATMUL(cell, nv)
        END DO
      END DO
    END DO

    write(*,*) "maxNeigh: ", maxNeigh

    ! find the number of neighbor vector in the WS cell, numNeigh
    numNeigh = 0
    DO neigh = 1, maxNeigh

      incell = .TRUE.
      ! loop over all neighbors to see if this is in the WS cell
      DO n = 1, maxNeigh
        IF(SUM(R(n,:)*R(neigh,:)) > (SUM(R(n,:)**2) + tol)) THEN
          ! this vector is outside the WS cell
          incell = .FALSE.
          EXIT
        END IF
      END DO
      IF(incell) THEN
        ! this vector is in the WS cell
        numNeigh = numNeigh + 1
        Rtmp(numNeigh,:) = R(neigh,:)
!        WRITE(*,*) "HIT, numNeigh: ",numNeigh
      END IF
    END DO

    write(*,*) "numNeigh: ", numNeigh

    DEALLOCATE(R)
    ALLOCATE(R(numNeigh,3))
    R(1:numNeigh,:) = Rtmp(1:numNeigh,:)
    DEALLOCATE(Rtmp)
           
    ! next step is to find all of the vortex points
    maxVert = (numNeigh-2)*(numNeigh-4)
    ALLOCATE (rVert(maxVert))
    ALLOCATE (alphtmp(numNeigh))
    numVect = numNeigh
    DO neigh = 1, numNeigh
      numVert = 1
      Rdot(1,:) = R(neigh,:)
      R2(1) = SUM(R(neigh,:)**2)
      DO nA = 1, numNeigh
        Rdot(2,:) = R(nA,:)
        R2(2) = SUM(R(nA,:)**2)
        DO nB = nA + 1, numNeigh
          Rdot(3,:) = R(nB,:)
          R2(3) = SUM(R(nB,:)**2)
!          CALL det(Rdot, detR)
          detR = determinant(Rdot)
          IF (ABS(detR) >= tol) THEN
            Rinv = matrix_3x3_inverse(Rdot)*detR
            rVert(numVert)%ot3 = MATMUL(Rinv, R2)/detR

            ! check if this vertex is in the WS cell
            incell = .TRUE.
            DO n = 1, numNeigh
              IF(SUM(rVert(numVert)%ot3(:)*R(n,:)) > SUM(R(n,:)**2) + tol) THEN
                ! outside the cell
                incell = .FALSE.
                EXIT
              END IF
            END DO
            IF (incell) THEN
                ! inside the cell
              numVert = numVert + 1
            END IF

          END IF
        END DO
      END DO

      zeroarea = 0
      DO n = 1, numVert
        IF (ABS(SUM(rVert(n)%ot3(:)**2)) < 0.5*SUM(R(neigh,:)**2) + tol) zeroarea = 1
      END DO
      IF (zeroarea == 1 .OR. numVert == 0) THEN
        alphtmp(neigh) = 0
        numVect = numVect - 1
        CYCLE
      END IF
      Rx = rvert(1)%ot3
      rdRn = SUM(rx(:)*R(neigh,:)) / SUM(R(neigh,:)**2)
      rx(:) = rx(:) - rdRn*R(neigh,:)
      rdRn = SQRT(SUM(rx(:)**2))
      Rx = Rx/rdRn
      Ry = cross_product(R(neigh,:), Rx)
      rdRn = SQRT(SUM(Ry(:)**2))
      Ry = Ry/rdRn
      DO nvi = 1, numVert - 1
        rvert(nvi)%angle = ATAN2(SUM(rvert(nvi)%ot3(:)*ry(:)), SUM(rvert(nvi)%ot3(:)*rx(:)))
      END DO
!      CALL rvert_sort(rvert, numVert - 1)
      CALL sort_vert(numVert-1, rvert)
      alphtmp(neigh) = 0
      DO nvi = 1, numVert - 1
!        CALL triple_product(rvert(nvi)%ot3, rvert(MOD(nvi, (numVert-1))+1)%ot3, R(neigh,:), temp)
!        alphtmp(neigh) = alphtmp(neigh) + temp
        alphtmp(neigh) = alphtmp(neigh) + &
                         triple_product(rvert(nvi)%ot3, rvert(MOD(nvi, (numVert-1))+1)%ot3, R(neigh,:))
      END DO
      alphtmp(neigh) = alphtmp(neigh)*0.25/SUM(R(neigh,:)**2)
      IF (ABS(alphtmp(neigh)) < tol) THEN
        alphtmp(neigh) = 0
        numVect = numVect - 1
      END IF
    END DO 

    ! assign the vertex array with the known number of verticies
    ALLOCATE (vect(numVect,3))
    ALLOCATE (alpha(numVect))
    nvi = 1
    DO n = 1, numNeigh
      IF (alphtmp(n) /= 0 ) THEN
        vect(nvi,:) = nvect(n,:)
        alpha(nvi) = alphtmp(n)
        nvi = nvi + 1
      END IF
    END DO 

    write(*,*) "numVect: ",numVect
    DO nvi = 1, numVect
      write(*,*) vect(nvi,1), " ", vect(nvi,2), " ", vect(nvi,3), " ", alpha(nvi)
    END DO
  END SUBROUTINE ws_voronoi


  !---------------------------------------------------------------
  ! Sort a charge list
  !---------------------------------------------------------------

  SUBROUTINE sort_weight(array_size, weightList)

    INTEGER, INTENT(IN) :: array_size
    TYPE(weight_obj), INTENT(INOUT), DIMENSION(array_size) :: weightList
    INCLUDE "qsort_inline.inc"
  CONTAINS
    SUBROUTINE init()
    END SUBROUTINE init
    LOGICAL &

    FUNCTION less_than(a,b)
      INTEGER, INTENT(IN) :: a,b
      IF ( weightList(a)%rho == weightList(b)%rho ) then
        less_than = a < b
      ELSE
        less_than = weightList(a)%rho == weightList(b)%rho
      END IF
    END FUNCTION less_than

    SUBROUTINE swap(a,b)
      INTEGER, INTENT(IN) :: a,b
      TYPE(weight_obj) :: hold
      hold = weightList(a)
      weightList(a) = weightList(b)
      weightList(b) = hold
    END SUBROUTINE swap

  ! circular shift-right by one:
    SUBROUTINE rshift(left,right)
      INTEGER, INTENT(in) :: left, right
      TYPE(weight_obj) :: hold
      hold = weightList(right)
      weightList(left+1:right) = weightList(left:right-1)
      weightList(left) = hold
    END SUBROUTINE rshift
  END SUBROUTINE sort_weight

  !---------------------------------------------------------------
  ! Sort a vertex list
  !---------------------------------------------------------------

  SUBROUTINE sort_vert(array_size, vertList)

    INTEGER, INTENT(IN) :: array_size
    TYPE(rvert_obj), INTENT(INOUT), DIMENSION(array_size) :: vertList
    INCLUDE "qsort_inline.inc"
  CONTAINS
    SUBROUTINE init()
    END SUBROUTINE init
    LOGICAL &

    FUNCTION less_than(a,b)
      INTEGER, INTENT(IN) :: a,b
      IF ( vertList(a)%angle == vertList(b)%angle ) then
        less_than = a < b
      ELSE
        less_than = vertList(a)%angle == vertList(b)%angle
      END IF
    END FUNCTION less_than

    SUBROUTINE swap(a,b)
      INTEGER, INTENT(IN) :: a,b
      TYPE(rvert_obj) :: hold
      hold = vertList(a)
      vertList(a) = vertList(b)
      vertList(b) = hold
    END SUBROUTINE swap

  ! circular shift-right by one:
    SUBROUTINE rshift(left,right)
      INTEGER, INTENT(IN) :: left, right
      TYPE(rvert_obj) :: hold
      hold = vertList(right)
      vertList(left+1:right) = vertList(left:right-1)
      vertList(left) = hold
    END SUBROUTINE rshift
  END SUBROUTINE sort_vert

  !---------------------------------------------------------------

  END MODULE
