  MODULE critpoints_mod
    USE kind_mod
    USE matrix_mod
    USE bader_mod
    USE charge_mod 
    USE options_mod
    USE ions_mod
    USE io_mod
    USE ions_mod
    USE dsyevj3_mod
    IMPLICIT NONE

    PRIVATE 
    PUBLIC :: critpoint_find

    TYPE hessian

      ! du dv dw are derivatives of the three original lattice vectors read from CHGCAR
      REAL(q2),DIMENSION(3) ::  du, dv, dw
      REAL(q2) :: dudu, dvdv, dwdw, dudv, dudw, dvdw
      ! eigval and eigvec are eigenvalues and eigenvectors of hessian matrix
    END TYPE
    CONTAINS

!-----------------------------------------------------------------------------------!
!critpoint_find: find critical points on the boundary of the Bader volumes
!NOTE: this subroutine should be called after refine_edge
!      in order to restrict the calculation to edge points
!-----------------------------------------------------------------------------------!
  SUBROUTINE critpoint_find(bdr,chg,opts,ions)

! Candidate CP which may be excluded due to numerical error.
    TYPE cpc ! stands for critical point candidate
      INTEGER,DIMENSION(3) :: ind ! indices of the cp
      REAL(q2),DIMENSION(3) :: force
      REAL(q2),DIMENSION(3) :: tempcart, tempind
      REAL(q2),DIMENSION(3,3,3) :: dx, dy, dz ! first derivatives of neighbors
!      REAL(q2),DIMENSION(3,3,3) :: du, dv, dw
      REAL(q2),DIMENSION(3) :: du, dv, dw ! 1 2 3 are backward, current, forward
      REAL(q2),DIMENSION(3) :: eigvals, r
      REAL(q2),DIMENSION(8,3) :: gnn ! gradients of neighbors
      REAL(q2),DIMENSION(3,3) :: eigvecs
      LOGICAL :: proxy, keep
    END TYPE

    TYPE(hessian) :: hes
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    TYPE(options_obj) :: opts
    TYPE(ions_obj) :: ions
    TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: cpl, cplt ! critical point list and a temporary copy
    INTEGER,DIMENSION(3) :: p, pt, ptt, ptx1, ptx2, pty1, pty2, ptz1, ptz2 ! points 1 and 2 are +1, -1

    INTEGER :: n1, n2, n3, d1, d2, d3, cptnum, ucptnum, i, j, debugnum
    REAL(q2),DIMENSION(3) :: eigvec1, eigvec2, eigvec3, tempVec
    REAL(q2),DIMENSION(3,3) ::  hessianMatrix, bkhessianMatrix
    ! vectors orthogonal to eigenvectors
    REAL(q2),DIMENSION(3) :: tem, tem2a, tem2b, tem2c, force, eigvals, carts
    REAL(q2) :: umag, vmag, wmag, threshhold, minmag
    REAL(q2),DIMENSION(3,3) :: eigvecs, A, inverseHessian
    ! linearized approximated derivatives for critical point screening
    REAL(q2),DIMENSION(3,3) :: transformationmatrix ! normalized lattice vectors
    LOGICAL :: trilinear ! determines if hessian calculation uses interpolation
    REAL(q2) :: dx0,dx1,dy0,dy1,dz0,dz1 ! outputs from interpolated gradients
    REAL(q2),DIMENSION(6,3) :: intcarts ! positions in cartesian coordinates of 6 interpolated points
    ! row 1 2 are + and 1 x, then + and - y, then + and - z
    REAL(q2),DIMENSION(6,3) :: intgrads ! gradients of interpolated points
    REAL(q2),DIMENSION(6) :: intrhos ! rhos of interpolated points 
    REAL(q2),DIMENSION(6,3) :: intinds ! fraction indicies for interpolated
    REAL(q2) :: rhocur ! rho of current point
    REAL(q2),DIMENSION(3) :: preal
    INTEGER,DIMENSION(8,3) :: nn ! alternative trilinear approx points
    LOGICAL,DIMENSION(3) :: cartcoor ! check if axis are cartesian
    WRITE(*,'(A)')  'FINDING CRITICAL POINTS'

    trilinear = .FALSE. ! do not interpolate, as it can get messy
    ucptnum = 0
    cptnum = 0
    PRINT * , "These code requires -vac auto or -vac #"
    PRINT *, '-----------                  WARNING                -----------'
    PRINT *, ' An analysis of valence charges will not yield sensible results'
    PRINT *, ' Use the total charge density for finding CPs                  '
    PRINT *, '_______________________________________________________________'
    OPEN(97,FILE='CPF.dat',STATUS='REPLACE',ACTION='WRITE')
    OPEN(98,FILE='CPFU.dat',STATUS='REPLACE',ACTION='WRITE')
    debugnum = 0
    umag = SQRT(ions%lattice(1,1)**2 + ions%lattice(1,2)**2 + &
      ions%lattice(1,3)**2) / REAL(chg%npts(1),q2)
    vmag = SQRT(ions%lattice(2,1)**2 + ions%lattice(2,2)**2 + &
      ions%lattice(2,3)**2) / REAL(chg%npts(2),q2)
    wmag = SQRT(ions%lattice(3,1)**2 + ions%lattice(3,2)**2 + &
      ions%lattice(3,3)**2) / REAL(chg%npts(3),q2)
    minmag = MIN(umag,vmag,wmag)
    ! check if axis are cartesian
!    cartcoor = coorcheck(ions%lattice)
    IF ( ALL(cartcoor) ) trilinear = .FALSE.
    DO n1 = 1,chg%npts(1)
      DO n2 = 1,chg%npts(2)
        DO n3 = 1,chg%npts(3)
            ! check to see if this point is in the vacuum
            IF (bdr%volnum(n1,n2,n3) == bdr%bnum + 1) THEN
              debugnum = debugnum + 1
              CYCLE
            END IF
            p = (/n1,n2,n3/)
            preal = p
            trilinear = .FALSE.
            IF (trilinear) THEN
              ! get cartesian coordinates of current point
              carts = getcart(p,chg%lat2car)
              ! get positions of where interpolated points should be
              intcarts = getintcarts(carts,umag,vmag,wmag)
              ! get indices of the interpolated points
              intinds = getinds(chg%car2lat,intcarts)
              ! get gradients of interpolated points
              DO i = 1, 6
                ! indices 1 to 6 are + - x, + - y, + - z
                intgrads(i,:) = rho_grad(chg,intinds(i,:),intrhos(i))
!                PRINT *, 'old grads', intgrads(1,:)
                ! alternative : find closest neighbors to do trilinear
                ! approximation
!                nn = findnn(p,intcarts(i,:),chg)
!                intgrads(i,:) = nn_grad(chg,intinds(i,:),intrhos(i),nn)
!                PRINT *, 'new grads', intgrads(1,:)
              END DO
              ! get second derivatives
              ! dxdx
              hessianMatrix(1,1) = (intgrads(1,1) - intgrads(2,1)) / (2 * minmag)
              ! dydy
              hessianMatrix(2,2) = (intgrads(3,2) - intgrads(4,2)) / (2 * minmag)
              ! dzdz
              hessianMatrix(3,3) = (intgrads(5,3) - intgrads(6,3)) / (2 * minmag)
              ! dxdy
              hessianMatrix(1,2) = (intgrads(1,2) - intgrads(2,2)) / (2 * minmag)
              hessianMatrix(2,1) = hessianMatrix(1,2)
              ! dxdz
              hessianMatrix(1,3) = (intgrads(1,3) - intgrads(2,3)) / (2 * minmag)
              hessianMatrix(3,1) = hessianMatrix(1,3)
              ! dydz
              hessianMatrix(2,3) = (intgrads(3,3) - intgrads(4,3)) / (2 * minmag)
              hessianMatrix(3,2) = hessianMatrix(3,2)
              ! force
              force = rho_grad(chg,preal,rhocur)
            ELSE


  !-----------------------------------------------------------------------------------!
  ! after finding candidate critical edge points, now  find the hessian
  !-----------------------------------------------------------------------------------!
               CALL getgradhes(p,chg,hes,force)
               hessianMatrix(1,1) = hes%dudu
               hessianMatrix(1,2) = hes%dudv
               hessianMatrix(1,3) = hes%dudw
               hessianMatrix(2,1) = hes%dudv
               hessianMatrix(2,2) = hes%dvdv
               hessianMatrix(2,3) = hes%dvdw
               hessianMatrix(3,1) = hes%dudw
               hessianMatrix(3,2) = hes%dvdw
               hessianMatrix(3,3) = hes%dwdw
               hessianMatrix = MATMUL(chg%car2lat,hessianMatrix)
               hessianMatrix = MATMUL(hessianMatrix,TRANSPOSE(chg%car2lat))

               force = MATMUL(chg%car2lat,force)
               ! now everything is cartesian
             END IF
!             tem2 = tem2 * umag 
!             tem2a = (/ ions%lattice(1,1),ions%lattice(2,1), & 
!               ions%lattice(3,1) /) / chg%npts(1)
             inverseHessian = inverse(hessianMatrix)
             tem = MATMUL(inverseHessian,force)
             ! convert from cartesian to lattice
             tem = MATMUL(chg%car2lat,tem)
             IF (ABS(tem(1)) <= 0.5) THEN
               IF (ABS(tem(2)) <= 0.5) THEN
                 IF (ABS(tem(3)) <= 0.5) THEN
                   cptnum = cptnum + 1
                   WRITE(97,*) '*********** A NEW ENTRY *************'
                   bkhessianMatrix = hessianMatrix
                   CALL DSYEVJ3(hessianMatrix,eigvecs,eigvals)
                   WRITE(97,*) 'Critical point number: ', cptnum
                   WRITE(97,*) "Indices are"
                   WRITE(97,*) p(1),p(2),p(3)
                   WRITE(97,*) "Density at this point is" 
                   WRITE(97,*) rho_val(chg,p(1),p(2),p(3))
                   WRITE(97,*) 'Eigenvalues: '
                   WRITE(97,*) eigvals
                   WRITE(97,'(3(1X,E18.11))') 
                   WRITE(97,*) 'Eigenvectors:'
                   WRITE(97,*) eigvecs(1,:)
                   WRITE(97,*) eigvecs(2,:)
                   WRITE(97,*) eigvecs(3,:)
                   WRITE(97,*) 'tem', tem
                   !WRITE(97,*) 'force is'
                   !WRITE(97,*)  force
                   !WRITE(97,*) 'Hessian is'
                   !WRITE(97,*)  bkhessianMatrix
                   IF (cptnum == 1)  THEN
                     ALLOCATE(cpl(1))
                     cpl(1)%du = hes%du
                     cpl(1)%dv = hes%dv
                     cpl(1)%dw = hes%dw
                     cpl(1)%ind(1) = n1
                     cpl(1)%ind(2) = n2
                     cpl(1)%ind(3) = n3
                     cpl(1)%force = force
                     cpl(1)%proxy = .FALSE.
                     cpl(1)%keep = .TRUE.
                     cpl(1)%eigvals = eigvals
                     cpl(1)%eigvecs = eigvecs
                     cpl(1)%r = tem
                     cpl(1)%tempcart = MATMUL(tem + p,chg%car2lat)
                   ELSE
                     ALLOCATE(cplt(cptnum))
                     DO i = 1, cptnum -1
                       cplt(i) = cpl(i)
                     END DO
                     DEALLOCATE(cpl)
                     ALLOCATE(cpl(cptnum))
                     DO i = 1, cptnum - 1
                       cpl(i) = cplt(i)
                     END DO
                     DEALLOCATE(cplt)
                     cpl(cptnum)%du = hes%du
                     cpl(cptnum)%dv = hes%dv
                     cpl(cptnum)%dw = hes%dw
                     cpl(cptnum)%ind(1) = n1
                     cpl(cptnum)%ind(2) = n2
                     cpl(cptnum)%ind(3) = n3
                     cpl(cptnum)%force = force
                     cpl(cptnum)%proxy = .FALSE.
                     cpl(cptnum)%keep = .TRUE.
                     cpl(cptnum)%eigvals = eigvals
                     cpl(cptnum)%eigvecs = eigvecs
                     cpl(cptnum)%r = tem
                     cpl(cptnum)%tempcart = MATMUL(tem + p, chg%car2lat)
                   END IF
                 END IF
               END IF
!               ELSE 
             END IF
        END DO
      END DO
    END DO
    PRINT *, "CRITICAL POINTS FOUND: ", cptnum 
    PRINT *,'CRITICAL POINT INFO WRITEN TO CPF.dat'
!    PRINT *, 'FINDING UNIQUE CRITICAL POINTS...'
!!------------------ Attempt A -------------------------------!
!!    DO i = 1, cptnum
!!      ! first filter out critical points that are less than 0.1
!!      IF (rho_val(chg,cpl(i)%ind(1),cpl(i)%ind(2),cpl(i)%ind(3)) <= 10 ) THEN
!!        CYCLE
!!      END IF
!!      IF (cpl(i)%keep == .FALSE.) THEN
!!        CYCLE
!!      END IF
!!      DO j = 1, cptnum
!!        IF (cpl(j)%keep == .FALSE. .OR. i == j) THEN
!!          CYCLE
!!        END IF
!!        IF (ABS(cpl(j)%ind(1) - cpl(i)%ind(1)) <= 1) THEN
!!          IF (ABS(cpl(j)%ind(2) - cpl(i)%ind(2)) <= 1) THEN
!!            IF (ABS(cpl(j)%ind(3) - cpl(i)%ind(3)) <= 1) THEN
!!              IF ((cpl(j)%tem(1)**2 + cpl(j)%tem(2)**2 + & 
!!                cpl(j)%tem(3)**2) <= (cpl(i)%tem(1)**2 + & 
!!                cpl(i)%tem(2)**2 + cpl(i)%tem(3)**2)) THEN
!!                cpl(i)%keep = .FALSE.
!!              ELSE 
!!                cpl(j)%keep = .FALSE.
!!              END IF
!!            END IF
!!          END IF
!!        END IF
!!      END DO
!!      IF (cpl(i)%keep == .TRUE.) THEN
!!        WRITE(98,*) '***** A NEW UNIQUE ENTRY *****'
!!        WRITE(98,*) 'Critical point number'
!!        WRITE(98,*) i
!!        ucptnum = ucptnum + 1
!!        WRITE(98,*) 'Unique Critical point number'
!!        WRITE(98,*) ucptnum
!!        WRITE(98,*) 'Density at this point is'
!!        WRITE(98,*) rho_val(chg,cpl(i)%ind(1),cpl(i)%ind(2),cpl(i)%ind(3))
!!        WRITE(98,*) 'Indicies are '
!!        WRITE(98,*) cpl(i)%ind
!!        WRITE(98,*) 'Eigenvalues are'
!!        WRITE(98,*) cpl(i)%eigvals
!!        WRITE(98,*) 'Eigenvectors are'
!!        WRITE(98,*) cpl(i)%eigvecs(1,:)
!!        WRITE(98,*) cpl(i)%eigvecs(2,:)
!!        WRITE(98,*) cpl(i)%eigvecs(2,:)
!!      END IF
!!    END DO
!!!----------------- END ATTEMPT A ------------------------------
!
!!    PRINT *, 'Found unique critical points: ', ucptnum
!!    DO n1=1,chg%npts(1)
!!      DO n2=1,chg%npts(2)
!!        DO n3=1,chg%npts(3)
!!          IF (n1 == n2 .AND. n2==n3) THEN
!!            CYCLE
!!          END IF
!!        END DO
!!      END DO
!!    END DO
!!    PRINT *, cpl
!!****************************************************************
!! ************* DONT DELETE THESE CODE **************************
!!    ! now go through candidates, for each candidate, see if there is another
!!    ! within 1 distance on u v w direction.
!!    DO i = 1 , cptnum - 1
!!      DO j = i + 1 , cptnum
!!        IF (i < j) THEN
!!          IF (ABS(cpl(i)%ind(1) - cpl(j)%ind(1)) <= 1 .AND. &
!!             ABS(cpl(i)%ind(2) - cpl(j)%ind(2)) <= 1 .AND. & 
!!             ABS(cpl(i)%ind(3) - cpl(j)%ind(3)) <= 1) THEN
!!             cpl(i)%proxy = .TRUE.
!!             cpl(j)%proxy = .TRUE.
!!             ! we know the cp is being shared. follow r to see
!!             ! if we go out of cell. 
!!             ! first step is to linearize dx
!!             ldu = cpl(i)%du(2) - cpl(j)%du(2)
!!             ldv = cpl(i)%dv(2) - cpl(j)%dv(2)
!!             ldw = cpl(i)%dw(2) - cpl(j)%dw(2) 
!!!             END IF
!!             WRITE(97,*)  i, j   
!!          END IF
!!        END IF
!!      END DO
!!    END DO
!!    PRINT *, 'number of unique cp :',debugnum
!!*******************************************************************
!!*******************************************************************
!    ! To find critical points (unique), start with a cell that contains a
!    ! critical point and its hessian and force. Use Newton's method to make a
!    ! move. Once moved, get the new force through trilinear interpolation, and
!    ! get the new hessian which will be a matrix of constants, make moves until
!    ! r is zero. get the coordinates of the new true critical point. If this
!    ! point is within half lattice to another, do not record this new point.
!    DO i = 1, cptnum
!      ! move to r     
!      PRINT *, 'CP is right now at'
!      PRINT *, cpl(i)%tempcart
!      ! find neighbors to interpolate points
!    END DO    
    DEALLOCATE (cpl)
    CLOSE(97)
    CLOSE(98)
    END SUBROUTINE critpoint_find


    REAL(q2) FUNCTION projection(r,rp)
      ! r is the vector pointing towards a CP
      ! rp is a cell vector to be projected onto
      REAL(q2),DIMENSION(3) :: r, rp
      projection = DOT_PRODUCT(r,rp)/sqrt((rp(1)**2 + &
      rp(2)**2 + rp(3)**2))
      RETURN
    END FUNCTION

    ! get cartesian coordinates of a point
    FUNCTION getcart(ind,lat2car)
      !ind is indicies of the current point
      !cart is the cartisian coordinates of the current point
      INTEGER,DIMENSION(3),INTENT(IN) :: ind
      REAL(q2),DIMENSION(3) :: getcart
      REAL(q2),DIMENSION(3,3) :: lat2car
      getcart(1) = ind(1) * lat2car(1,1) + & 
        ind(2) * lat2car(1,2) + & 
        ind(3) *  lat2car(1,3)
      getcart(2) = ind(1) * lat2car(2,1) + &
        ind(2) * lat2car(2,2) + &
        ind(3) *  lat2car(2,3)
      getcart(3) = ind(1) * lat2car(3,1) + &
        ind(2) * lat2car(3,2) + &
        ind(3) *  lat2car(3,3)

      RETURN
    END FUNCTION

    ! get cartesian coordinates of interpolated points
    FUNCTION  getintcarts(carts,umag,vmag,wmag)
      REAL(q2) :: minmag,umag,vmag,wmag
      REAL(q2),DIMENSION(3) :: carts
      REAL(q2),DIMENSION(6,3) :: getintcarts
      INTEGER :: i,j,k
      ! first find which one is the smallest mag
      minmag = MIN(umag,vmag,wmag)
      DO i = 1, 6
        getintcarts(i,:) = carts
      END DO
      getintcarts(1,1) = getintcarts(1,1) + minmag
      getintcarts(2,1) = getintcarts(2,1) - minmag
      getintcarts(3,2) = getintcarts(3,2) + minmag
      getintcarts(4,2) = getintcarts(4,2) - minmag
      getintcarts(5,3) = getintcarts(5,3) + minmag
      getintcarts(6,3) = getintcarts(6,3) - minmag
      RETURN
    END FUNCTION
    
    ! get indices from cartesian coordinates
    FUNCTION  getinds(car2lat,intcarts) 
      REAL(q2),DIMENSION(6,3) :: getinds, intcarts
      REAL(q2),DIMENSION(3,3) :: car2lat,W
      REAL(q2),DIMENSION(3) :: Z
      INTEGER :: i,j,k
      DO i = 1, 6 
        getinds(i,:) = MATMUL(car2lat,intcarts(i,:))
      END DO
      RETURN
    END FUNCTION

    ! finds nearest grid point to a interpolated point
    ! used for doing trilinear interpolation
    FUNCTION findnn(p,intcart,chg)
      TYPE(charge_obj) :: chg
      INTEGER,DIMENSION(8,3) :: findnn,tfindnn
      REAL(q2),DIMENSION(3) :: intcart,p2cart
      INTEGER,DIMENSION(3) :: p,maxs
      REAL(q2),DIMENSION(27) :: dist
      INTEGER,DIMENSION(27) :: scores
      INTEGER,DIMENSION(27,3) :: p2
      INTEGER :: i,j,k,counter
      ! since only closest neigherbos of current grid point
      ! will be used for interpolation
      counter = 1
      DO i = -1,1
        DO j = -1,1
          DO k = -1,1
            p2(counter,:) = (/i,j,k/)
            CALL pbc(p2,chg%npts)
            p2cart = getcart(p2,chg%lat2car)
            dist(counter) = SQRT( &
              (intcart(1) - p2cart(1))**2 + & 
              (intcart(2) - p2cart(2))**2 + &
              (intcart(3) - p2cart(3))**2 &
              )
            counter = counter + 1
          END DO
        END DO
      END DO
      ! now find the top 4 smallest array
      ! compare each element to all
      ! if it is smaller, it gets score
      ! keep the ones with highest scores
      ! there should not be points with equal scores
      ! Each point should have a score ranging from 0 to 26
      DO i = 1, 27
        scores(i) = 0
        DO j = 1, 27
          IF (i == j) THEN
            CYCLE
          END IF
          IF (dist(i) <= dist(j)) THEN
            scores(i) = scores(i) + 1
          END IF
          IF (scores(i) == 26) THEN
            tfindnn(1,:) = p2(i,:)
            CYCLE
          END IF
          IF (scores(i) == 25) THEN
            tfindnn(2,:) = p2(i,:)
            CYCLE
          END IF
          IF (scores(i) == 24) THEN
            tfindnn(3,:) = p2(i,:)
            CYCLE
          END IF
          IF (scores(i) == 23) THEN
            tfindnn(4,:) = p2(i,:)
            CYCLE
          END IF
          IF (scores(i) == 22) THEN
            tfindnn(5,:) = p2(i,:)
            CYCLE
          END IF
          IF (scores(i) == 21) THEN
            tfindnn(6,:) = p2(i,:)
            CYCLE
          END IF
          IF (scores(i) == 20) THEN
            tfindnn(7,:) = p2(i,:)
            CYCLE
          END IF
          IF (scores(i) == 19) THEN
            tfindnn(8,:) = p2(i,:)
            CYCLE
          END IF
        END DO
      END DO     
      ! then these needs to be rearranged 
      ! so that is follows this format:
      !  1   2   3   4   5   6   7   8
      ! 000 001 010 100 011 101 110 111
      ! indices should have only 2 values
      ! on each axis. min gives 0 on that 
      ! axis and max gives 1. 
      maxs = (/0,0,0/)
      DO i = 1, 8
        DO j = 1,3
          IF (tfindnn(i,j)>= maxs(j) ) THEN
            maxs(j) = tfindnn(i,j)
          END IF
        END DO
      END DO
      DO i = 1, 8
        IF (ALL(tfindnn(i,:) == maxs)) THEN
          findnn(8,:) = tfindnn(i,:)
          CYCLE
        END IF
        IF (ALL(tfindnn(i,:) == maxs - (/1,1,1/))) THEN
          findnn(1,:) = tfindnn(i,:)
          CYCLE
        END IF
        IF (ALL(tfindnn(i,:) == maxs - (/1,1,0/))) THEN
          findnn(2,:) = tfindnn(i,:)
          CYCLE
        END IF
        IF (ALL(tfindnn(i,:) == maxs - (/1,0,1/))) THEN
          findnn(3,:) = tfindnn(i,:)
          CYCLE
        END IF
        IF (ALL(tfindnn(i,:) == maxs - (/0,1,1/))) THEN
          findnn(4,:) = tfindnn(i,:)
          CYCLE
        END IF
        IF (ALL(tfindnn(i,:) == maxs - (/1,0,0/))) THEN
          findnn(5,:) = tfindnn(i,:)
          CYCLE
        END IF
        IF (ALL(tfindnn(i,:) == maxs - (/0,1,0/))) THEN
          findnn(6,:) = tfindnn(i,:)
          CYCLE
        END IF
        IF (ALL(tfindnn(i,:) == maxs - (/0,0,1/))) THEN
          findnn(7,:) = tfindnn(i,:)
          CYCLE
        END IF
      END DO 
      RETURN
    END FUNCTION

    ! This function takes in a list of nearest neighbors
    FUNCTION nn_grad(chg,r,rho,nn)
      TYPE(charge_obj) :: chg
      REAL(q2),DIMENSION(3),INTENT(IN) :: r
      REAL(q2),INTENT(OUT) :: rho
      REAL(q2),DIMENSION(3) :: nn_grad
      INTEGER :: p1, p2, p3
      REAL(q2),DIMENSION(3) :: rho_grad_lat
      REAL(q2) :: f1, f2, f3, g1, g2, g3
      REAL(q2) :: rho000, rho001, rho010, rho100, rho011, rho101, rho110, rho111
      REAL(q2) :: rho00_, rho01_, rho10_, rho11_
      REAL(q2) :: rho0__, rho1__, rho_0_, rho_1_, rho__0, rho__1
      REAL(q2) :: rho_00, rho_01, rho_10, rho_11
      INTEGER,DIMENSION(8,3) :: nn
      p1 = FLOOR(r(1))
      p2 = FLOOR(r(2))
      p3 = FLOOR(r(3))
      f1 = r(1) - REAL(p1,q2)
      f2 = r(2) - REAL(p2,q2)
      f3 = r(3) - REAL(p3,q2)
      g1 = 1._q2-f1
      g2 = 1._q2-f2
      g3 = 1._q2-f3
      rho000 = rho_val(chg,nn(1,1),nn(1,2),nn(1,3))
      rho001 = rho_val(chg,nn(2,1),nn(2,2),nn(2,3))
      rho010 = rho_val(chg,nn(3,1),nn(3,2),nn(3,3))
      rho100 = rho_val(chg,nn(4,1),nn(4,2),nn(4,3))
      rho011 = rho_val(chg,nn(5,1),nn(5,2),nn(5,3))
      rho101 = rho_val(chg,nn(6,1),nn(6,2),nn(6,3))
      rho110 = rho_val(chg,nn(7,1),nn(7,2),nn(7,3))
      rho111 = rho_val(chg,nn(8,1),nn(8,2),nn(8,3))
      rho00_ = rho000*g3 + rho001*f3
      rho01_ = rho010*g3 + rho011*f3
      rho10_ = rho100*g3 + rho101*f3
      rho11_ = rho110*g3 + rho111*f3
      rho0__ = rho00_*g2 + rho01_*f2
      rho1__ = rho10_*g2 + rho11_*f2
      rho = rho0__*g1 + rho1__*f1
  ! More work for gradients
      rho_0_ = rho00_*g1 + rho10_*f1
      rho_1_ = rho01_*g1 + rho11_*f1
      rho_00 = rho000*g1 + rho100*f1
      rho_01 = rho001*g1 + rho101*f1
      rho_10 = rho010*g1 + rho110*f1
      rho_11 = rho011*g1 + rho111*f1
      rho__0 = rho_00*g2 + rho_10*f2
      rho__1 = rho_01*g2 + rho_11*f2
      rho_grad_lat(1) = rho1__ - rho0__
      rho_grad_lat(2) = rho_1_ - rho_0_
      rho_grad_lat(3) = rho__1 - rho__0
  !    CALL vector_matrix(rho_grad_lat, chg%car2lat, rho_grad)
      nn_grad = MATMUL(rho_grad_lat, chg%car2lat)
    RETURN
    END FUNCTION nn_grad

    FUNCTION coorcheck(lattice)
      LOGICAL,DIMENSION(3) :: coorcheck
      REAL(q2) :: prod1,prod2,prod3
      REAL(q2),DIMENSION(3,3) :: lattice
      INTEGER :: i
      ! elements 1, 2, 3 are for whether coordinates align with x, y, z axis.
      coorcheck = (/.TRUE.,.TRUE.,.TRUE./)
      IF (lattice(1,2) .NE. 0 .OR. lattice(1,3) .NE. 0) THEN
        coorcheck(1)=.FALSE.
      END IF
      IF (lattice(2,1) .NE. 0 .OR. lattice(2,3) .NE. 0) THEN
        coorcheck(2)=.FALSE.
      END IF
      IF (lattice(3,2) .NE. 0 .OR. lattice(3,1) .NE. 0) THEN
        coorcheck(3)=.FALSE.
      END IF
    END FUNCTION
    
    SUBROUTINE getgradhes(p,chg,hes,force)
    TYPE(hessian) :: hes
    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3) :: ptxy1, ptxy2, ptxy3, ptxy4, ptxz1, ptxz2, ptxz3 
    INTEGER,DIMENSION(3) :: ptxz4, ptyz1, ptyz2, ptyz3, ptyz4
    REAL(q2),DIMENSION(3) :: force
    INTEGER,DIMENSION(3) :: p, pt, ptt, ptx1, ptx2, pty1, pty2, ptz1, ptz2
    ! calculate hessian matrix, the second derivative at a point
    ! is calculated in the following way:
    ! first order derivatives are calculated between the point
    ! and its first neighbor point to get a central difference derivative
    ! the forward and backward central derivatives are used to calculate
    ! the second order derivatives
    ptx1 = p + (/1,0,0/)
    ptx2 = p + (/-1,0,0/)
    pty1 = p + (/0,1,0/)
    pty2 = p + (/0,-1,0/)
    ptz1 = p + (/0,0,1/)
    ptz2 = p + (/0,0,-1/)
    CALL pbc(ptx1,chg%npts)
    CALL pbc(ptx2,chg%npts)
    CALL pbc(pty1,chg%npts)
    CALL pbc(pty2,chg%npts)
    CALL pbc(ptz1,chg%npts)
    CALL pbc(ptz2,chg%npts)
    ! these are the first order derivatives
!    hes%du(1) = 0.5_q2 / (umag / 2) * & 
    hes%du(1) = &
      (rho_val(chg,p(1),p(2),p(3)) - &
       rho_val(chg,ptx2(1),ptx2(2),ptx2(3)))
    hes%du(2) = 0.5 * & 
      (rho_val(chg,ptx1(1),ptx1(2),ptx1(3)) - &
       rho_val(chg,ptx2(1),ptx2(2),ptx2(3)))
    hes%du(3) = & 
      (rho_val(chg,ptx1(1),ptx1(2),ptx1(3)) - &
       rho_val(chg,p(1),p(2),p(3)))
    hes%dv(1) = & 
      (rho_val(chg,p(1),p(2),p(3)) - &
       rho_val(chg,pty2(1),pty2(2),pty2(3)))
    hes%dv(2) = 0.5 * & 
      (rho_val(chg,pty1(1),pty1(2),pty1(3)) - &
       rho_val(chg,pty2(1),pty2(2),pty2(3)))
    hes%dv(3) = & 
      (rho_val(chg,pty1(1),pty1(2),pty1(3)) - &
       rho_val(chg,p(1),p(2),p(3)))
    hes%dw(1) = & 
      (rho_val(chg,p(1),p(2),p(3)) - &
       rho_val(chg,ptz2(1),ptz2(2),ptz2(3)))
    hes%dw(2) = 0.5 * & 
      (rho_val(chg,ptz1(1),ptz1(2),ptz1(3)) - &
       rho_val(chg,ptz2(1),ptz2(2),ptz2(3)))
    hes%dw(3) = & 
      (rho_val(chg,ptz1(1),ptz1(2),ptz1(3)) - &
       rho_val(chg,p(1),p(2),p(3)))
    ! these are the second order derivatives
!    hes%dudu = 1 / umag * (hes%du(3) - hes%du(1))
    hes%dudu = (hes%du(3) - hes%du(1)) 
    hes%dvdv = (hes%dv(3) - hes%dv(1)) 
    hes%dwdw = (hes%dw(3) - hes%dw(1)) 
    ! for dudv, calculate dv of half points on u direction, using
    ! extrapolated points, such as - 0.5 u, +- 1 v
    ptxy1 = p + (/-1,-1,0/)
    ptxy2 = p + (/-1,+1,0/)
    ptxy3 = p + (/+1,+1,0/)
    ptxy4 = p + (/+1,-1,0/)
    ptxz1 = p + (/-1,0,-1/)
    ptxz2 = p + (/-1,0,+1/)
    ptxz3 = p + (/+1,0,+1/)
    ptxz4 = p + (/+1,0,-1/)
    ptyz1 = p + (/0,-1,-1/)
    ptyz2 = p + (/0,-1,+1/)
    ptyz3 = p + (/0,+1,+1/)
    ptyz4 = p + (/0,+1,-1/)
    CALL pbc(ptxy1,chg%npts)
    CALL pbc(ptxy2,chg%npts)
    CALL pbc(ptxy3,chg%npts)
    CALL pbc(ptxy4,chg%npts)
    CALL pbc(ptxz1,chg%npts)
    CALL pbc(ptxz2,chg%npts)
    CALL pbc(ptxz3,chg%npts)
    CALL pbc(ptxz4,chg%npts)
    CALL pbc(ptyz1,chg%npts)
    CALL pbc(ptyz2,chg%npts)
    CALL pbc(ptyz3,chg%npts)
    CALL pbc(ptyz4,chg%npts)
!    hes%dudv = 1 / umag * & 
    hes%dudv = &
      ( &
      ! this is the backward dv
!      - 0.5_q2 / vmag * &
      - 0.5_q2 * & 
      ((rho_val(chg,ptxy2(1),ptxy2(2),ptxy2(3)) + &
          rho_val(chg,pty1(1),pty1(2),pty1(3))) / 2 - &
        (rho_val(chg,ptxy1(1),ptxy1(2),ptxy1(3)) + &
          rho_val(chg,pty2(1),pty2(2),pty2(3))) / 2 ) & 
      ! this is the forward dv
      + 0.5_q2 * & 
      ((rho_val(chg,ptxy3(1),ptxy3(2),ptxy3(3)) + &
          rho_val(chg,pty1(1),pty1(2),pty1(3))) / 2 - &
        (rho_val(chg,ptxy4(1),ptxy4(2),ptxy4(3)) + &
          rho_val(chg,pty2(1),pty2(2),pty2(3))) / 2 ) &
      )
    hes%dudw = &
      ( &
      ! this is the bacward dw
      - 0.5_q2 * &
      ((rho_val(chg,ptxz2(1),ptxz2(2),ptxz2(3)) + &
          rho_val(chg,ptz1(1),ptz1(2),ptz1(3))) / 2 - &
        (rho_val(chg,ptxz1(1),ptxz1(2),ptxz1(3)) + &
          rho_val(chg,ptz2(1),ptz2(2),ptz2(3))) / 2 )  &
      ! this is the forward dw
      + 0.5_q2 * & 
      ((rho_val(chg,ptxz3(1),ptxz3(2),ptxz3(3)) + &
          rho_val(chg,ptz1(1),ptz1(2),ptz1(3))) / 2 - &
        (rho_val(chg,ptxz4(1),ptxz4(2),ptxz4(3)) + &
          rho_val(chg,ptz2(1),ptz2(2),ptz2(3))) / 2 ) &
      ) 
  
    hes%dvdw = &
      ( &
      ! this is the bacward dw
      - 0.5_q2 * &
      ((rho_val(chg,ptyz2(1),ptyz2(2),ptyz2(3)) + &
          rho_val(chg,ptz1(1),ptz1(2),ptz1(3))) / 2 - &
        (rho_val(chg,ptyz1(1),ptyz1(2),ptyz1(3)) + &
          rho_val(chg,ptz2(1),ptz2(2),ptz2(3))) / 2 ) & 
      ! this is the forward dw
      + 0.5_q2 * & 
      ((rho_val(chg,ptyz3(1),ptyz3(2),ptyz3(3)) + &
          rho_val(chg,ptz1(1),ptz1(2),ptz1(3))) / 2 - &
        (rho_val(chg,ptyz4(1),ptyz4(2),ptyz4(3)) + &
          rho_val(chg,ptz2(1),ptz2(2),ptz2(3))) / 2 ) &
      )

      force(1) = hes%du(2)
      force(2) = hes%dv(2)
      force(3) = hes%dw(2)
    END SUBROUTINE

  END MODULE

