  MODULE critpoints_mod
    USE kind_mod
    USE matrix_mod
    USE bader_mod
!    USE linsolve_mod
    USE charge_mod 
    USE options_mod
    USE ions_mod
    USE io_mod
    USE ions_mod
    USE weight_mod
    USE dsyevj3_mod
    IMPLICIT NONE

    PRIVATE 
    PUBLIC :: critpoint_find

    TYPE hessian

      ! du dv dw are derivatives of the three original lattice vectors read from
      ! CHGCAR
      REAL(q2),DIMENSION(3) ::  du, dv, dw
      REAL(q2) :: dudu, dvdv, dwdw, dudv, dudw, dvdw
      ! eigval and eigvec are eigenvalues and eigvectors of hessian matrix
    END TYPE
    CONTAINS

!-----------------------------------------------------------------------------------!
!critpoint_find: find critical points on the edge of the Bader volumes
!NOTE: this subroutine should be called after refine_edge
!      in order to restrict the calculation to edge points
!-----------------------------------------------------------------------------------!
  SUBROUTINE critpoint_find(bdr,chg,opts,ions,stat)
! These are for screening CP due to numerical error. 
    TYPE cpc ! stands for critical point candidate
      INTEGER,DIMENSION(3) :: ind  ! these are the indices of the cp
      REAL(q2),DIMENSION(3) :: trueind, truer
      REAL(q2),DIMENSION(3) :: force 
      REAL(q2),DIMENSION(3) :: tempcart, tempind
      REAL(q2),DIMENSION(3,3,3) :: dx, dy, dz ! first derivatives of neighbors
!      REAL(q2),DIMENSION(3,3,3) :: du, dv, dw
      REAL(q2),DIMENSION(3) :: du, dv, dw ! 1 2 3 are backward, present, forward
      REAL(q2),DIMENSION(3) :: eigvals, r, cocart, colat, tempr
      INTEGER,DIMENSION(8,3) :: nnind ! indices of neighbors 
      ! indices from 0 to 8 are zyx 000 001 010 011 100 101 110 111
!      REAL(q2),DIMENSION(8,3) :: nngrad ! gradients of nn mentioned above.
!      REAL(q2),DIMENSION(6,3) :: intnngrad ! gradients of interpolated neighbors
      ! used to find interpolated hessians.
      REAL(q2),DIMENSION(3,3) :: eigvecs
      INTEGER :: negcount
      LOGICAL :: proxy, isunique
    END TYPE

    TYPE(hessian) :: hes
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    TYPE(options_obj) :: opts
    TYPE(ions_obj) :: ions
    TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: cpl, cplt ! critical point list and a temporary
! copy
! for points, 1 and 2 are +1, -1
    INTEGER :: stat
    INTEGER,DIMENSION(3) :: p, pt, ptt, ptx1, ptx2, pty1, pty2, ptz1, ptz2
    INTEGER,DIMENSION(3) :: tempind
    INTEGER :: n1, n2, n3, d1, d2, d3, cptnum, ucptnum, i, j, k, debugnum
    INTEGER :: negcount, bondcount, ringcount, maxcount, cagecount
    INTEGER :: ubondcount, uringcount, umaxcount, ucagecount
    REAL(q2),DIMENSION(3) :: eigvec1, eigvec2, eigvec3, tempVec
    REAL(q2),DIMENSION(3) :: tempforce, truer
    INTEGER, DIMENSION(3) :: tempr
    REAL(q2), DIMENSION(3) :: temprealr, tempreal3d
    ! to be used in newton method in finding unique critical points.
    REAL(q2),DIMENSION(3,3) ::  hessianMatrix, bkhessianMatrix
    ! these are vectors orthogonal to eigenvectors
    REAL(q2),DIMENSION(3) :: tem, tem2a,tem2b,tem2c, force,eigvals,carts
    REAL(q2) :: umag, vmag, wmag, threshhold, minmag, tempreal
    REAL(q2),DIMENSION(3,3) :: eigvecs, inverseHessian
    ! linearized approximated derivatives for proxy critical screening
    REAL(q2) :: dx0,dx1,dy0,dy1,dz0,dz1 ! outputs from interpolated gradients
    REAL(q2),DIMENSION(6,3) :: intcarts ! positions in cart of 6 interpolated points
    ! row 1 2 are + and 1 x, then + and - y, then + and - z
    REAL(q2),DIMENSION(6,3) :: intgrads ! gradients of interpolated points
    REAL(q2),DIMENSION(6) :: intrhos ! rhos of interpolated points 
    REAL(q2),DIMENSION(6,3) :: intinds ! fraction indicies for interpolated
    REAL(q2) :: rhocur ! rho of current point
    REAL(q2) :: stepsize
    REAL(q2),DIMENSION(3) :: distance ! vector to 000 in trilinear
    REAL(q2),DIMENSION(3) :: preal
    INTEGER,DIMENSION(8,3) :: nn ! alternative trilinear approx.
    REAL(q2),DIMENSION(8) :: vals
    ! points
    LOGICAL,DIMENSION(3) :: cartcoor ! check if axis are alone cartesian.
    ! The followings are for finding unique critical points
    REAL(q2),DIMENSION(8,3) :: nngrad
    REAL(q2),DIMENSION(8,3,3) :: nnhes !hessian of 8 nn
    INTEGER, DIMENSION(8,3) :: nnind
    REAL(q2),DIMENSION(6,3) :: intnngrad
    REAL(q2),DIMENSION(3,3) :: temphessian
    REAL(q2),DIMENSION(3) :: nexttem, previoustem, averager
    INTEGER :: averagecount
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: nucleiInd
    INTEGER :: stepcount
    ! The following are for least square calculations
    INTEGER, DIMENSION(3,26) :: vi, matw
    INTEGER, DIMENSION(26,3) :: vit
    REAL(q2), DIMENSION(3,3) :: ggrid
    REAL(q2), DIMENSION(26) :: wi
    REAL(q2), DIMENSION(3,13) :: matwprime
    REAL(q2), DIMENSION(3,3) :: matm, outerproduct
    ! below are variables for least sqaures gradient
    stat = 0 ! 0 means nothing
    PRINT *, ''//achar(27)//'[31m Finding Critical points'//achar(27)//'[0m'
    PRINT *, ''//achar(27)//'[91m Interrogation of the soul:'//achar(27)//'[0m'
    PRINT *, ''//achar(27)//'[33m Did I turn on vacuum ?'//achar(27)//'[0m'
    PRINT *, ''//achar(27)//'[92m Did I tell if this is a crystall or molecule?'//achar(27)//'[0m'
    PRINT *, ''//achar(27)//'[36m Did I use the CHGCAR_sum ?'//achar(27)//'[0m'
    PRINT *, ''//achar(27)//'[34m Did I use the least square flag -ls ? '//achar(27)//'[0m'
    PRINT *, ''//achar(27)//'[95m Critical point is like a box of chocolates. &
              You never know what you are gonna get.'//achar(27)//'[0m'
    !WRITE(*,'(A)')  'FINDING CRITICAL POINTS'
    IF (opts%leastsquare_flag .eqv. .TRUE. ) THEN
      PRINT *, 'Using least square gradient'
    END IF
    ! get expected nucleus indices
!    WRITE (97,*) , 'expecting nuclei at:'
    ALLOCATE(nucleiInd(ions%nions,3))
    DO d1 = 1, ions%nions
      nucleiInd(d1,1) = NINT(ions%r_lat(d1,1))
      nucleiInd(d1,2) = NINT(ions%r_lat(d1,2))
      nucleiInd(d1,3) = NINT(ions%r_lat(d1,3))
!      WRITE (97,*), nucleiInd(d1,:)
    END DO
    bondcount = 0
    ubondcount = 0
    ringcount = 0
    uringcount = 0
    maxcount = ions%nions
    ucagecount = 0
    cagecount = 0
    ucptnum = ions%nions
    cptnum = 0

!    PRINT * , "These code requires -vac auto or -vac #"
!    PRINT *, '-----------                  WARNING             -----------'
!    PRINT *, 'Using valence charge may yield useless and confusing results'
!    PRINT *, '    It is recommended to use total charge for finding CPs   '
!    PRINT *, '____________________________________________________________'
!    OPEN(97,FILE='CPF.dat',STATUS='REPLACE',ACTION='WRITE')
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
    IF ( opts%leastsquare_flag .eqv. .TRUE. )THEN
      vi = makevi()
      vit = TRANSPOSE(vi)
      ggrid = makeggrid(chg)
      DO i = 1,26
        wi(i) = 1/DOT_PRODUCT(MATMUL(vit(i,:),ggrid),vi(:,i))
      END DO
      matm = 0.0_q2
      DO i = 1,13
        matwprime(:,i) = vi(:,i) * wi(i)
      END DO
      DO i = 1, 26
        outerproduct(1,1)= vi(1,i) * vit(i,1)
        outerproduct(1,2)= vi(1,i) * vit(i,2)
        outerproduct(1,3)= vi(1,i) * vit(i,3)
        outerproduct(2,1)= vi(2,i) * vit(i,1)
        outerproduct(2,2)= vi(2,i) * vit(i,2)
        outerproduct(2,3)= vi(2,i) * vit(i,3)
        outerproduct(3,1)= vi(3,i) * vit(i,1)
        outerproduct(3,2)= vi(3,i) * vit(i,2)
        outerproduct(3,3)= vi(3,i) * vit(i,3)
        matm = matm + wi(i) * outerproduct
      END DO
    END IF



!    truer = (/10.192,3.778,7.001/)
!    DO n1 = 1,50
!        nnind(1,:) = (/floor(truer(1)),floor(truer(2)),    &
!                floor(truer(3))/)
!        nnind(2,:) = (/ceiling(truer(1)),floor(truer(2)),  &
!                floor(truer(3))/)
!        nnind(3,:) = (/floor(truer(1)),ceiling(truer(2)),  &
!                floor(truer(3))/)
!        nnind(4,:) = (/ceiling(truer(1)),ceiling(truer(2)),&
!                floor(truer(3))/)
!        nnind(5,:) = (/floor(truer(1)),floor(truer(2)),    &
!                ceiling(truer(3))/)
!        nnind(6,:) = (/ceiling(truer(1)),floor(truer(2)),  &
!                ceiling(truer(3))/)
!        nnind(7,:) = (/floor(truer(1)),ceiling(truer(2)),  &
!                ceiling(truer(3))/)
!        nnind(8,:) = (/ceiling(truer(1)),ceiling(truer(2)),&
!                ceiling(truer(3))/)
!        distance = truer - nnind(1,:)
!        do j = 1,8
!          if (opts%leastsquare_flag == .true.) then
!            nngrad(j,:) = lsg(nnind(j,:),chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
!            hessianmatrix = lsh(nnind(j,:),chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
!            nnhes(j,1,:) = hessianmatrix(1,:)
!            nnhes(j,2,:) = hessianmatrix(2,:)
!            nnhes(j,3,:) = hessianmatrix(3,:)
!          else 
!            call getgradhes(nnind(j,:),chg,hes,nngrad(j,:))
!            nnhes(j,1,:) = hessianmatrix(1,:)
!            nnhes(j,2,:) = hessianmatrix(2,:)
!            nnhes(j,3,:) = hessianmatrix(3,:)
!          end if
!        end do
!        ! row find the nearest neighbors at this new locaiton
!        ! first update critical point location
!        ! the next big step is to interpolate the force at predicted critical
!        ! point.
!        tempforce = trilinear_interpol_grad(nngrad,distance) ! val r interpol
!        temphessian = trilinear_interpol_hes(nnhes,distance)
!        nexttem = - matmul(inverse(temphessian),tempforce)
!        nexttem = matmul(nexttem,chg%car2lat)
!      !tempforce = rho_grad(chg,truer,tempreal)
!      !PRINT *, tempforce(1) + tempforce(2) + tempforce(3)
!      tempreal3d = (/0.1,0.1,0.1/)
!      truer = truer + tempreal3d
!    END DO
!    STOP



    ! First start with nuclear critical points
    ALLOCATE (cpl(10000)) ! start with 10000 capacity
    DO n1 = 1, ions%nions
      tempr(1) = NINT(ions%r_dir(n1,1) * chg%npts(1)) 
      tempr(2) = NINT(ions%r_dir(n1,2) * chg%npts(2)) 
      tempr(3) = NINT(ions%r_dir(n1,3) * chg%npts(3)) 
      temprealr = lsgascension(tempr,chg,matm,matwprime, &
                 wi,vi,vit,ggrid,outerproduct,opts)
      cptnum = cptnum + 1
      cpl(cptnum)%trueind = temprealr
      WRITE(98,*) '_________________________________________'
      WRITE(98,*) 'Nucleus critical point found at'
      WRITE(98,*) temprealr
      WRITE(98,*) 'Coordinates in cartesian are'
      WRITE(98,*) MATMUL(temprealr,chg%lat2car)
      WRITE(98,*) 'Direct coordinates are'
      WRITE(98,*)  temprealr(1)/chg%npts(1),temprealr(2)/chg%npts(2), &
                   temprealr(3)/chg%npts(3)
      WRITE(98,*) '_________________________________________'
      WRITE(98,*) ' ' 
    END DO
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
  !-----------------------------------------------------------------------------------!
  ! now that this subroutine can find the correct amount of edge points, lets have
  ! it find the hessian
  !-----------------------------------------------------------------------------------!
            IF (opts%leastsquare_flag .eqv. .TRUE.) THEN
              force = lsg(p,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
              hessianMatrix = &
                lsh(p,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
              tem = - MATMUL(INVERSE(hessianMatrix),force)
              ! tem is now in cartesian. convert it back to lattice
              tem = MATMUL(tem, chg%car2lat)
            ELSE 
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
              inverseHessian = inverse(hessianMatrix)
              tem = - MATMUL(inverseHessian,force)
              hessianMatrix = MATMUL(chg%car2lat,hessianMatrix)
              hessianMatrix = MATMUL(hessianMatrix,TRANSPOSE(chg%car2lat))
              force = MATMUL(chg%car2lat,force)
              ! Now everything is cartesian.
            END IF
            ! convert from cartesian to lattice
            IF ( ABS(tem(1)) <= 1 + opts%knob_tem) THEN
              IF (ABS(tem(2)) <= 1 + opts%knob_tem) THEN
                IF (ABS(tem(3)) <= 1 + opts%knob_tem) THEN              
                  cptnum = cptnum + 1
!                  WRITE(97,*) '*********** A NEW ENTRY *************'
                  bkhessianMatrix = hessianMatrix
!                  CALL DSYEVJ3(hessianMatrix,eigvecs,eigvals)
!                  WRITE(97,*) 'Critical candidate point number: ', cptnum
!                  WRITE(97,*) "Indices are"
!                  WRITE(97,*) p(1),p(2),p(3)
!                  WRITE(97,*) "Density at this point is" 
!                  WRITE(97,*) rho_val(chg,p(1),p(2),p(3))
!                  WRITE(97,'(3(1X,E18.11))') 
!                  WRITE(97,*) 'tem', tem
!                  WRITE(97,*) 'in cartesian units, force is'
!                  WRITE(97,*)  MATMUL(force,chg%lat2car)
                  IF (cptnum < SIZE(cpl) ) THEN
                    cpl(cptnum)%du = hes%du
                    cpl(cptnum)%dv = hes%dv  
                    cpl(cptnum)%dw = hes%dw
                    cpl(cptnum)%ind(1) = n1
                    cpl(cptnum)%ind(2) = n2
                    cpl(cptnum)%ind(3) = n3
                    cpl(cptnum)%force = force
                    cpl(cptnum)%proxy = .FALSE.
                    cpl(cptnum)%r = tem
                    cpl(cptnum)%tempcart = MATMUL(tem + p,chg%car2lat)
                  ELSE 
                    ALLOCATE(cplt(cptnum))
                    DO i = 1, cptnum -1
                      cplt(i) = cpl(i)
                    END DO
                    DEALLOCATE(cpl)
                    ALLOCATE(cpl(cptnum*10))
                    DO i = 1, cptnum - 1
                      cpl(i)=cplt(i)
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
                    cpl(cptnum)%r = tem
                    cpl(cptnum)%tempcart = MATMUL(tem + p, chg%car2lat)
                END IF
              END IF
            END IF
            ELSE 
              DO i = 1, ions%nions
                IF (n1 == nucleiInd(i,1) .AND. n2 == nucleiInd(i,2) &
                    .AND. n3 == nucleiInd(i,3)) THEN
!                  WRITE (97,*), '****************** WARNING ******************'
!                  WRITE (97,*), 'Expected Nucleus Critical Point Nout Found at:'
!                  WRITE (97,*), '           ', p
!                  WRITE (97,*), '****************** WARNING ******************'
                END IF
              END DO
            END IF
        END DO
      END DO
    END DO
    PRINT *,'CRITICAL POINTS INFO WRITEN TO CPF.dat'
    PRINT *, "CRITICAL POINTS FOUND: ", cptnum 
    PRINT *, 'FINDING UNIQUE CRITICAL POINTS...'
!!*******************************************************************
    ! To find critical points (unique), start with a cell that contains a
    ! critical point and its hessian and force. Use Newton's method to make a
    ! move. Interpolate the force inside the voxel. 
    ! Once moved, get the new force through trilinear interpolation, and
    ! get the new hessian which will be a matrix of constants, make moves until
    ! r is zero. get the coordinates of the new true critical point. If this
    ! point is within half lattice to another, do not record this new point.
    DO i = ions%nions + 1, cptnum
      cpl(i)%isunique = .FALSE.
      stepcount = 0
      averagecount = 0
      ! move to r in lattice units
      ! find the nearest neighbors to this point.
      ! first find a nearby grid point.
      ! get cartesians of the current location
      ! this is not needed if we don't want to guarentee that that NN are 
      ! actually the closest. 
      ! find nearest neighbors
      ! the following code uses only one side of the criticalpoints
      ! x-1, y-1, z-1 000
      cpl(i)%nnind(1,:) = cpl(i)%ind + & 
        (/FLOOR(cpl(i)%r(1)),FLOOR(cpl(i)%r(2)),FLOOR(cpl(i)%r(3))/)
      ! x+1, y-1, z-1 100
      cpl(i)%nnind(2,:) = cpl(i)%ind + &
        (/CEILING(cpl(i)%r(1)),FLOOR(cpl(i)%r(2)),FLOOR(cpl(i)%r(3))/)
      ! x-1, y+1, z-1 010
      cpl(i)%nnind(3,:) = cpl(i)%ind + &
        (/FLOOR(cpl(i)%r(1)),CEILING(cpl(i)%r(2)),FLOOR(cpl(i)%r(3))/)
      ! x+1, y+1, z-1 110
      cpl(i)%nnind(4,:) = cpl(i)%ind + &
        (/CEILING(cpl(i)%r(1)),CEILING(cpl(i)%r(2)),FLOOR(cpl(i)%r(3))/)
      ! x-1, y-1, z+1 001
      cpl(i)%nnind(5,:) = cpl(i)%ind + &
        (/FLOOR(cpl(i)%r(1)),FLOOR(cpl(i)%r(2)),CEILING(cpl(i)%r(3))/)
      ! x+1, y-1, z+1 101
      cpl(i)%nnind(6,:) = cpl(i)%ind + &
        (/CEILING(cpl(i)%r(1)),FLOOR(cpl(i)%r(2)),CEILING(cpl(i)%r(3))/)
      ! x-1, y+1, z+1 011
      cpl(i)%nnind(7,:) = cpl(i)%ind + &
        (/FLOOR(cpl(i)%r(1)),CEILING(cpl(i)%r(2)),CEILING(cpl(i)%r(3))/)
      ! x            111
      cpl(i)%nnind(8,:) = cpl(i)%ind + &
        (/CEILING(cpl(i)%r(1)),CEILING(cpl(i)%r(2)),CEILING(cpl(i)%r(3))/)
      ! get gradients "force" of all nearest neighbors
      ! These codes uses a larger box.
      DO j = 1, 8
        ! getting nn indices are checked. pass. 
        IF (opts%leastsquare_flag .eqv. .TRUE.) THEN
          nngrad(j,:) = lsg(cpl(i)%nnind(j,:),chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
        ELSE 
          CALL getgradhes(tempind,chg,hes,nngrad(j,:))
        END IF
        ! this force is not adjusted with lattice
        !cpl(i)%nngrad(j,:) = MATMUL(chg%car2lat,cpl(i)%nngrad(j,:))
        !PRINT *, cpl(i)%nngrad(j,:)
        ! now things are cartesian. forces are checked to be correct.
        ! correct: if I make things cartesian, forces are checked out to be
        ! correct. However I still need to find hessian in lattice first. 
        ! So I do not convert force to cartesian yet.
        ! force is checked to be correct.
      END DO
      ! Now start newton method iterations
      truer =  cpl(i)%ind + cpl(i)%r
      CALL pbc_r_lat(truer,chg%npts)
      previoustem = cpl(i)%r
      averager = (/0.,0.,0./)
      DO stepcount = 1,400
        
        CALL pbc_r_lat(truer,chg%npts)
        nnind(1,:) = (/floor(truer(1)),floor(truer(2)),    &
                floor(truer(3))/)
        nnind(2,:) = (/ceiling(truer(1)),floor(truer(2)),  &
                floor(truer(3))/)
        nnind(3,:) = (/floor(truer(1)),ceiling(truer(2)),  &
                floor(truer(3))/)
        nnind(4,:) = (/ceiling(truer(1)),ceiling(truer(2)),&
                floor(truer(3))/)
        nnind(5,:) = (/floor(truer(1)),floor(truer(2)),    &
                ceiling(truer(3))/)
        nnind(6,:) = (/ceiling(truer(1)),floor(truer(2)),  &
                ceiling(truer(3))/)
        nnind(7,:) = (/floor(truer(1)),ceiling(truer(2)),  &
                ceiling(truer(3))/)
        nnind(8,:) = (/ceiling(truer(1)),ceiling(truer(2)),&
                ceiling(truer(3))/)
        distance = truer - nnind(1,:)
        do j = 1,8
          if (opts%leastsquare_flag .eqv. .true.) then
            nngrad(j,:) = lsg(nnind(j,:),chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
            hessianmatrix = lsh(nnind(j,:),chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
            nnhes(j,1,:) = hessianmatrix(1,:)
            nnhes(j,2,:) = hessianmatrix(2,:)
            nnhes(j,3,:) = hessianmatrix(3,:)
          else 
            call getgradhes(nnind(j,:),chg,hes,nngrad(j,:))
            nnhes(j,1,:) = hessianmatrix(1,:)
            nnhes(j,2,:) = hessianmatrix(2,:)
            nnhes(j,3,:) = hessianmatrix(3,:)
          end if
        end do
        ! row find the nearest neighbors at this new locaiton
        ! first update critical point location
        ! the next big step is to interpolate the force at predicted critical
        ! point.
        tempforce = trilinear_interpol_grad(nngrad,distance) ! val r interpol
        temphessian = trilinear_interpol_hes(nnhes,distance)
        nexttem = - matmul(inverse(temphessian),tempforce)
        nexttem = matmul(nexttem,chg%car2lat)
!        nexttem = 0.8*nexttem
        ! nexttem should not be too large
        nexttem(1) = MIN(nexttem(1),1.)
        nexttem(2) = MIN(nexttem(2),1.)
        nexttem(3) = MIN(nexttem(3),1.)
        
        tempr(1) = NINT(truer(1))
        tempr(2) = NINT(truer(2))
        tempr(3) = NINT(truer(3))
        CALL pbc(tempr,chg%npts)
        IF (stepcount >= 1900) THEN
          ! take the average over 100 steps as the critical point position.
          averager = averager + truer
          averagecount = averagecount + 1
        END IF
        IF (stepcount == 2000) THEN
          WRITE(98,*) 'Special treatment done at'
          WRITE(98,*) cpl(i)%ind
          averager = averager / averagecount
          truer = averager
          nexttem = 0.
        END IF
        IF (bdr%volnum(tempr(1), &
            tempr(2),tempr(3)) == bdr%bnum + 1) THEN
          ! We are heading into the vacuum space, cosmonaughts! 
          cpl(i)%isunique = .FALSE.
          EXIT
        END IF
!        IF ( truer(1) == truer(1) + nexttem(1) .AND. &
!             truer(2) == truer(2) + nexttem(2) .AND. &
!             truer(3) == truer(3) + nexttem(3) )  THEN
        IF (cpl(i)%ind(1) == 1 .AND. cpl(i)%ind(2) == 6 &
            .AND. cpl(i)%ind(3) == 18 ) THEN
          WRITE(*,*) 'truer is', truer
        END IF
!        IF (cpl(i)%ind(1) == 6 .AND. & cpl(i)%ind(2) == 10 &
!            .AND. cpl(i)%ind(3) == 19 ) THEN
!          PRINT *, 'truer is'
!          PRINT *, truer
!        END IF
        IF ( ABS(nexttem(1)) <= 0.00001 + opts%knob_newtonr .AND. &
             ABS(nexttem(2)) <= 0.00001 + opts%knob_newtonr .AND. &
             ABS(nexttem(3)) <= 0.00001 + opts%knob_newtonr ) THEN
          cpl(i)%trueind = truer 
          cpl(i)%isunique = .TRUE.
          DO j = 1, i - 1 
            IF (ABS(truer(1) - cpl(j)%trueind(1) )<= 0.5 + &
                opts%knob_distance .AND. &
                ABS(truer(2) - cpl(j)%trueind(2) )<= 0.5 + &
                opts%knob_distance .AND. &
                ABS(truer(3) - cpl(j)%trueind(3) )<= 0.5 + &
                opts%knob_distance) THEN
              cpl(i)%isunique = .FALSE.
              EXIT
            END IF
          END DO
          IF (cpl(i)%isunique .eqv. .TRUE.) THEN
!            WRITE (98,*) 'Inspecting critical point number: ', i 
!            WRITE (98,*) 'Indicies of this point is'
!            WRITE (98,*) cpl(i)%ind
            WRITE (98,*) '_______________________________________'
            WRITE (98,*) 'Critical point is found at indices'
            WRITE (98,*) truer
            WRITE (98,*) 'Coordinates in cartesian are'
            WRITE (98,*) MATMUL(truer,chg%lat2car)
            ucptnum = ucptnum + 1
            WRITE (98,*) 'Direct coordinates are'
            WRITE (98,*) truer(1)/chg%npts(1),truer(2)/chg%npts(2), &
                          truer(3)/chg%npts(3)
            WRITE (98,*) 'Gradiant is'
            WRITE (98,*) tempforce
            !PRINT *, 'original critical point gradient is'
            !PRINT *, cpl(i)%force
            !PRINT *, 'original critical point gradient magnitude is'
            !PRINT *, SQRT(cpl(i)%force(1)**2 + cpl(i)%force(2)**2 + &
            !         cpl(i)%force(3)**2)
            WRITE (98,*) 'Hessian is'
            WRITE (98,*) temphessian(1,:)
            WRITE (98,*) temphessian(2,:)
            WRITE (98,*) temphessian(3,:)
            CALL DSYEVJ3(temphessian,eigvecs,eigvals)

            WRITE (98,*) 'Eigenvalues are'
            WRITE (98,*) eigvals
            WRITE (98,*) 'Eigenvectors are'
            WRITE (98,*) eigvecs
            negcount = 0
            DO d1 = 1,3
              IF (eigvals(d1) <0) THEN
                negcount = negcount + 1
              END IF
            END DO
            IF (negcount == 0) THEN
              cagecount = cagecount + 1
              WRITE(*,*) 'Found a unique cage critical point'
              WRITE(98,*) 'This is a cage critical point'
              WRITE(98,*) ' '
            END IF
            IF (negcount == 2) THEN
              bondcount = bondcount + 1
              WRITE(*,*) 'Found a unique bond critical point'
              WRITE(98,*) 'This is a bond critical point'
              WRITE(98,*) ' '
            ELSEIF(negcount == 1) THEN
              ringcount = ringcount + 1
              WRITE(*,*) 'Found a unique ring critical point'
              WRITE(98,*) 'This is a ring critical point'
              WRITE(98,*) ' '
            ELSEIF(negcount == 3) THEN
              maxcount = maxcount + 1
              WRITE(*,*) 'Found a unique nuclear critical point'
              WRITE(98,*) 'This is a nuclear critical point'
              WRITE(98,*) ' ' 
              IF (maxcount > ions%nions) THEN
                WRITE(*,*) 'WARNING 1: Finding more nucleus points than  &
                          number of atoms!'
              END IF
            END IF
            WRITE(98,*) '________________________________________'
          END IF
          EXIT
        ELSE 
          previoustem = nexttem
        END IF
        IF (stepcount == 1000 ) THEN
          WRITE(98,*) 'Inspecting critical point number: ', i
          WRITE(98,*) ' ******* 1000 steps not enough **********'
          WRITE(*,*) 'WARNING 2: Fail to exit Newton wrapping'
        END IF
          ! for the line above:
          ! note here that cpl(i)%r is the vector tem starting from cpl(i)%ind
          ! tem is in lattice unit.
        ! should expect gradient of neighbors with different diretions
        ! nngrad forces are in lattice. so these forces are also in lattice. 
        ! assuming it is correct.
        truer = truer + nexttem ! this keeps track the total movement
        CALL pbc_r_lat(truer,chg%npts)
      END DO
    END DO    
    WRITE(*,*) 'Unique critical point count: ', ucptnum
    WRITE(*,*) 'Unique bond critical point count: ', bondcount
    WRITE(*,*) 'Unique ring critical point count: ', ringcount
    WRITE(*,*) 'Unique cage critical point count: ', cagecount
    WRITE(*,*) 'Unique nucl critical point count: ', maxcount
    WRITE(98,*) 'Unique critical point count: ', ucptnum
    WRITE(98,*) 'Unique bond critical point count: ', bondcount
    WRITE(98,*) 'Unique ring critical point count: ', ringcount
    WRITE(98,*) 'Unique cage critical point count: ', cagecount
    WRITE(98,*) 'Unique nucl critical point count: ', maxcount
    IF (opts%ismolecule) THEN
      IF (maxcount - bondcount + ringcount - cagecount == 1) THEN
        WRITE(*,*) 'Satisfies The Poincare-Hopf rule for a molecule'
        WRITE(98,*) 'Satisfies The Poincare-Hopf rule for a molecule'
      ELSE
        WRITE(*,*) 'WARNING 3: Fails The poincare-Hopf rule for a molecule'
        WRITE(98,*)  'WARNING 3: Fails The poincare-Hopf rule for a molecule'
      END IF
    ELSE IF (opts%iscrystal) THEN
      IF (maxcount - bondcount + ringcount - cagecount == 0) THEN
        WRITE(*,*) 'Satisfies The Poincare-Hopf rule for a crystal'
        WRITE(98,*) 'Satisfies The Poincare-Hopf rule for a crystal'
      ELSE
        WRITE(*,*) 'WARNING 4: Fails The poincare-Hopf rule for a crystal'
        WRITE(98,*)  'WARNING 4: Fails The poincare-Hopf rule for a crystal'

      END IF
    END IF
    DEALLOCATE (cpl)
!    CLOSE(97)
    CLOSE(98)
    END SUBROUTINE critpoint_find

    ! get cartesian coordinates of a point
    FUNCTION getcart(ind,lat2car)
      !ind is indicies of the current point
      !cart is the cartisian coordinates of the current point
      INTEGER,DIMENSION(3),INTENT(IN) :: ind
      REAL(q2),DIMENSION(3) :: getcart
      REAL(q2),DIMENSION(3,3) :: lat2car
      getcart(1) = ind(1) * lat2car(1,1) + & 
        ind(2) * lat2car(1,2) + & 
        ind(3) * lat2car(1,3)
      getcart(2) = ind(1) * lat2car(2,1) + &
        ind(2) * lat2car(2,2) + &
        ind(3) * lat2car(2,3)
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
    
    ! get indices from cartesian coordinates.
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

    ! find neares on grid points for a interpolated point
    ! to do trilinear interpolation
    ! p is the centered on grid point. intcart is a nearby point to be
    ! interpolated in cartesian coordinates.
    FUNCTION findnn(p,intcart,chg) ! THIS FUNCTION IS NOT STABLE
      TYPE(weight_obj),DIMENSION(27) :: nndist ! record distances of a point to all nearest on
      !grid points
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
      !CALL pbc(p,chg%npts)
      !to calculate distance it is not necessary to run pbc
      !infact pbc should be avoided at this stage
      DO i = -1,1
        DO j = -1,1
          DO k = -1,1
            p2(counter,:) = (/i,j,k/) + p
            !CALL pbc(p2(counter,:),chg%npts)
            PRINT *,'p2 is', p2(counter,:)
            !p2cart = getcart(p2,chg%lat2car)
            p2cart = MATMUL(p2(counter,:),chg%lat2car)
!            PRINT *, 'p2 cart is', p2cart
!            PRINT *,'counter is', counter
            nndist(counter)%rho = SQRT( &
              (intcart(1) - p2cart(1))**2 + & 
              (intcart(2) - p2cart(2))**2 + &
              (intcart(3) - p2cart(3))**2 &
              )
            WRITE(*,*) 'dist is',dist(counter)
            counter = counter + 1
          END DO
        END DO
      END DO
      ! now find the top 4 smallest array
      ! compare each element to all
      ! if it is smaller, it gets score
      ! keep the ones with highest scores
      ! there should not be points with equal 
      ! scores. Each point should have a score
      ! ranging from 0 to 26.
      DO i = 1, 27
        scores(i) = 0
        DO j = 1, 27
          IF (i == j) THEN
            CYCLE
          END IF
          IF (dist(i) <= dist(j)) THEN
            scores(i) = scores(i) + 1
          END IF
        END DO
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
      WRITE(*,*) 'scores', scores
      WRITE(*,*) 'tfindnn'
      DO i = 1, 8
        WRITE(*,*) tfindnn(i,:)
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
      PRINT *, 'maxs is', maxs
      DO i = 1, 8
        PRINT *, 'tfindnn here is'
        PRINT * ,tfindnn(i,:)
        IF (ALL(tfindnn(i,:) == maxs)) THEN
          PRINT *, '8'
          findnn(8,:) = tfindnn(i,:)
          CYCLE
        END IF
        IF (ALL(tfindnn(i,:) == maxs - (/1,1,1/))) THEN
          findnn(1,:) = tfindnn(i,:)
          PRINT *, '1'
          CYCLE
        END IF
        IF (ALL(tfindnn(i,:) == maxs - (/1,1,0/))) THEN
          findnn(2,:) = tfindnn(i,:)
          PRINT *, '2'
          CYCLE
        END IF
        IF (ALL(tfindnn(i,:) == maxs - (/1,0,1/))) THEN
          findnn(3,:) = tfindnn(i,:)
          PRINT *, '3'
          CYCLE
        END IF
        IF (ALL(tfindnn(i,:) == maxs - (/0,1,1/))) THEN
          findnn(4,:) = tfindnn(i,:)
          PRINT *, '4'
          CYCLE
        END IF
        IF (ALL(tfindnn(i,:) == maxs - (/1,0,0/))) THEN
          findnn(5,:) = tfindnn(i,:)
          PRINT *, '5'
          CYCLE
        END IF
        IF (ALL(tfindnn(i,:) == maxs - (/0,1,0/))) THEN
          findnn(6,:) = tfindnn(i,:)
          PRINT *, '6'
          CYCLE
        END IF
        IF (ALL(tfindnn(i,:) == maxs - (/0,0,1/))) THEN
          findnn(7,:) = tfindnn(i,:)
          PRINT *, '7'
          CYCLE
        END IF
      END DO 
      RETURN
    END FUNCTION

    ! This function should take in a list of nearest neighbors predetermined.
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
  !   CALL vector_matrix(rho_grad_lat, chg%car2lat, rho_grad)
      nn_grad = MATMUL(rho_grad_lat, chg%car2lat)
    RETURN
    END FUNCTION nn_grad

    ! this funciton takes in 8 values, return a
    ! trilinear interpolated gradient of the values.
    ! the 8 value list value order is 
    ! 000 001 010 100 011 101 110 111
    ! Note 02042019: the above order is what I wrote previously
    ! Note 02042019: I believe the actuall order is the following
    ! 000 100 010 110 001 101 011 111
    !  1   2   3   4   5   6   7   8
    ! r is the indice of the predicted critical point
    ! The interpolation result is checked to be OK by mathematica
    FUNCTION trilinear_interpol_grad(vals,r)
      ! varls come nngrad
      ! could r be the problem?
      ! r needs to be a vector where each component 
      ! number is between 0 and 1
      REAL(q2),DIMENSION(3),INTENT(IN) :: r
      REAL(q2),DIMENSION(3) :: trilinear_interpol_grad
      INTEGER :: p1, p2, p3
      REAL(q2) :: f1, f2, f3, g1, g2, g3
      REAL(q2),DIMENSION(8,3) :: vals
      !REAL(q2) :: rho000, rho001, rho010, rho100, rho011, rho101, rho110, rho111
      REAL(q2),DIMENSION(3) :: val00_, val01_, val10_, val11_
      REAL(q2),DIMENSION(3) :: val0__, val1__, val_0_, val_1_, val__0, val__1
      REAL(q2),DIMENSION(3) :: val_00, val_01, val_10, val_11
 
      p1 = FLOOR(r(1))
      p2 = FLOOR(r(2))
      p3 = FLOOR(r(3))
      f1 = r(1) - REAL(p1,q2)
      f2 = r(2) - REAL(p2,q2)
      f3 = r(3) - REAL(p3,q2)
      ! f1 f2 f3 are checked to be correct. 
      ! they should equal to tem for the first step
      g1 = 1._q2-f1
      g2 = 1._q2-f2
      g3 = 1._q2-f3
      val00_ = vals(1,:)*g3 + vals(5,:)*f3
      val01_ = vals(3,:)*g3 + vals(7,:)*f3
      val10_ = vals(2,:)*g3 + vals(6,:)*f3
      val11_ = vals(4,:)*g3 + vals(8,:)*f3
      val0__ = val00_*g2 + val01_*f2
      val1__ = val10_*g2 + val11_*f2
      trilinear_interpol_grad = val0__*g1 + val1__*f1
    RETURN
    END FUNCTION trilinear_interpol_grad
 
    FUNCTION trilinear_interpol_hes(vals,r)
      ! varls come nnhes
      ! could r be the problem?
      ! r needs to be a vector where each component
      ! number is between 0 and 1
      REAL(q2),DIMENSION(3),INTENT(IN) :: r
      REAL(q2),DIMENSION(3,3) :: trilinear_interpol_hes
      INTEGER :: p1, p2, p3
      REAL(q2) :: f1, f2, f3, g1, g2, g3
      REAL(q2),DIMENSION(8,3,3) :: vals
      !REAL(q2) :: rho000, rho001, rho010, rho100, rho011, rho101, rho110,
      !rho111
      REAL(q2),DIMENSION(3,3) :: val00_, val01_, val10_, val11_
      REAL(q2),DIMENSION(3,3) :: val0__, val1__, val_0_, val_1_, val__0, val__1
      REAL(q2),DIMENSION(3,3) :: val_00, val_01, val_10, val_11

      p1 = FLOOR(r(1))
      p2 = FLOOR(r(2))
      p3 = FLOOR(r(3))
      f1 = r(1) - REAL(p1,q2)
      f2 = r(2) - REAL(p2,q2)
      f3 = r(3) - REAL(p3,q2)
      ! f1 f2 f3 are checked to be correct.
      ! they should equal to tem for the first step
      g1 = 1._q2-f1
      g2 = 1._q2-f2
      g3 = 1._q2-f3
      val00_(1,:) = vals(1,1,:)*g3 + vals(5,1,:)*f3
      val00_(2,:) = vals(1,2,:)*g3 + vals(5,2,:)*f3
      val00_(3,:) = vals(1,3,:)*g3 + vals(5,3,:)*f3
      val01_(1,:) = vals(3,1,:)*g3 + vals(7,1,:)*f3
      val01_(2,:) = vals(3,2,:)*g3 + vals(7,2,:)*f3
      val01_(3,:) = vals(3,3,:)*g3 + vals(7,3,:)*f3
      val10_(1,:) = vals(2,1,:)*g3 + vals(6,1,:)*f3
      val10_(2,:) = vals(2,2,:)*g3 + vals(6,2,:)*f3
      val10_(3,:) = vals(2,3,:)*g3 + vals(6,3,:)*f3
      val11_(1,:) = vals(4,1,:)*g3 + vals(8,1,:)*f3
      val11_(2,:) = vals(4,2,:)*g3 + vals(8,2,:)*f3
      val11_(3,:) = vals(4,3,:)*g3 + vals(8,3,:)*f3
      val0__ = val00_*g2 + val01_*f2
      val1__ = val10_*g2 + val11_*f2
      trilinear_interpol_hes = val0__*g1 + val1__*f1
      
    END FUNCTION

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
    
    ! The following subroutine gets gradient using central difference
    FUNCTION cdgrad(p,chg)
      TYPE(charge_obj) :: chg
      INTEGER, DIMENSION(3) :: p
      INTEGER, DIMENSION(3) :: pzm,pzp,pxm,pxp,pym,pyp
      REAL(q2), DIMENSION(3) :: cdgrad
      pzm = p + (/0,0,-1/)
      pzp = p + (/0,0,1/)
      pxm = p + (/-1,0,0/)
      pxp = p + (/1,0,0/)
      pym = p + (/0,-1,0/)
      pyp = p + (/0,1,0/)
      CALL pbc(pxm,chg%npts)
      CALL pbc(pym,chg%npts)
      CALL pbc(pzm,chg%npts)
      CALL pbc(pxp,chg%npts)
      CALL pbc(pyp,chg%npts)
      CALL pbc(pzp,chg%npts)
      cdgrad(3) = rho_val(chg,pzp(1),pzp(2),pzp(3)) - &
                  rho_val(chg,pzm(1),pzm(2),pzm(3))
      cdgrad(2) = rho_val(chg,pyp(1),pyp(2),pyp(3)) - &
                  rho_val(chg,pym(1),pym(2),pym(3))
      cdgrad(1) = rho_val(chg,pxp(1),pxp(2),pxp(3)) - &
                  rho_val(chg,pxm(1),pxm(2),pxm(3))
      cdgrad = MATMUL(cdgrad,chg%car2lat)
      RETURN
      ! now the gradient should be in cartesian
    END FUNCTION
    
    ! the following subroutine gets hes and force in lattice units
    SUBROUTINE getgradhes(p,chg,hes,force)
    TYPE(hessian) :: hes
    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3) :: ptxy1, ptxy2, ptxy3, ptxy4, ptxz1, ptxz2, ptxz3 
    INTEGER,DIMENSION(3) :: ptxz4, ptyz1, ptyz2, ptyz3, ptyz4
    REAL(q2),DIMENSION(3) :: force
    INTEGER,DIMENSION(3) :: p, pt, ptt, ptx1, ptx2, pty1, pty2, ptz1, ptz2
    ! to calculate hessian matrix, second derivatives at a point

    ! is calculated in the following way:
    ! first order derivatives are calculated between the point
    ! and its first neighbor point to get the first derivative at
    ! half point.
    ! the forward and backward half point first derivatives are
    ! used to calculate second order derivatives at the point.
    CALL pbc(p,chg%npts)
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
      - 0.25_q2 * & 
      ((rho_val(chg,ptxy2(1),ptxy2(2),ptxy2(3)) + &
          rho_val(chg,pty1(1),pty1(2),pty1(3)))  - &
        (rho_val(chg,ptxy1(1),ptxy1(2),ptxy1(3)) + &
          rho_val(chg,pty2(1),pty2(2),pty2(3)))  ) & 
      ! this is the forward dv
      + 0.25_q2 * & 
      ((rho_val(chg,ptxy3(1),ptxy3(2),ptxy3(3)) + &
          rho_val(chg,pty1(1),pty1(2),pty1(3)))  - &
        (rho_val(chg,ptxy4(1),ptxy4(2),ptxy4(3)) + &
          rho_val(chg,pty2(1),pty2(2),pty2(3)))  ) &
      )
    hes%dudw = &
      ( &
      ! this is the bacward dw
      - 0.25_q2 * &
      ((rho_val(chg,ptxz2(1),ptxz2(2),ptxz2(3)) + &
          rho_val(chg,ptz1(1),ptz1(2),ptz1(3)))   - &
        (rho_val(chg,ptxz1(1),ptxz1(2),ptxz1(3)) + &
          rho_val(chg,ptz2(1),ptz2(2),ptz2(3)))   )  &
      ! this is the forward dw
      + 0.25_q2 * & 
      ((rho_val(chg,ptxz3(1),ptxz3(2),ptxz3(3)) + &
          rho_val(chg,ptz1(1),ptz1(2),ptz1(3)))   - &
        (rho_val(chg,ptxz4(1),ptxz4(2),ptxz4(3)) + &
          rho_val(chg,ptz2(1),ptz2(2),ptz2(3)))   ) &
      ) 
  
    hes%dvdw = &
      ( &
      ! this is the bacward dw
      - 0.25_q2 * &
      ((rho_val(chg,ptyz2(1),ptyz2(2),ptyz2(3)) + &
          rho_val(chg,ptz1(1),ptz1(2),ptz1(3)))  - &
        (rho_val(chg,ptyz1(1),ptyz1(2),ptyz1(3)) + &
          rho_val(chg,ptz2(1),ptz2(2),ptz2(3)))  ) & 
      ! this is the forward dw
      + 0.25_q2 * & 
      ((rho_val(chg,ptyz3(1),ptyz3(2),ptyz3(3)) + &
          rho_val(chg,ptz1(1),ptz1(2),ptz1(3)))  - &
        (rho_val(chg,ptyz4(1),ptyz4(2),ptyz4(3)) + &
          rho_val(chg,ptz2(1),ptz2(2),ptz2(3)))  ) &
      )
    
      force(1) = hes%du(2)
      force(2) = hes%dv(2)
      force(3) = hes%dw(2)
    END SUBROUTINE
    

    ! This function takes in current position and grid point position in
    ! lattice, finds the distance this point is to the nearest cell, take half
    ! the distance as step size. 
    FUNCTION findstepsize(r)
      REAL(q2) :: findstepsize    
      REAL(q2),DIMENSION(3) :: r
      REAL(q2) :: f1,f2,f3,c1,c2,c3
      
      f1 = MIN(ABS(r(1)),ABS(0.5-ABS(r(1)))) 
      f2 = MIN(ABS(r(2)),ABS(0.5-ABS(r(2))))
      f3 = MIN(ABS(r(3)),ABS(0.5-ABS(r(3))))
      findstepsize = MIN(MIN(f1,f2),MIN(f1,f3))/2
    RETURN
    END FUNCTION findstepsize

    ! this funciton finds hessian of a interpolated point using interpolated
    ! nearest neighbor gradiants. Also gradiants taken in here should be in
    ! lattice.
    FUNCTION inthessian(grad,stepsize)
      REAL(q2),DIMENSION(6,3) :: grad
      REAL(q2),DIMENSION(3,3) :: inthessian
      REAL(q2) :: stepsize
      ! again, intnngrad is following this order:
      ! +x -x +y -y +z -z
      inthessian(1,1) = (grad(1,1)-grad(2,1))*0.5_q2/stepsize
      inthessian(2,2) = (grad(3,2)-grad(4,2))*0.5_q2/stepsize
      inthessian(3,3) = (grad(5,3)-grad(6,3))*0.5_q2/stepsize
      inthessian(1,2) = (grad(3,1)-grad(4,1))*0.5_q2/stepsize
      inthessian(2,1) = inthessian(1,2)
      inthessian(1,3) = (grad(5,1)-grad(6,1))*0.5_q2/stepsize
      inthessian(3,1) = inthessian(1,3)
      inthessian(2,3) = (grad(5,2)-grad(6,2))*0.5_q2/stepsize
      inthessian(3,2) = inthessian(2,3)
    ! assuming that this function is fine
    RETURN
    END FUNCTION inthessian

    FUNCTION newtonstep(car2lat,lat2car,hessian,force)
      REAL(q2),DIMENSION(3,3) :: car2lat,lat2car,hessian,inversehessian
      REAL(q2),DIMENSION(3) :: force,newtonstep
!      PRINT *, 'in newtonstep, force is'
!      PRINT *, force
!      PRINT *, 'hessian is'
!      PRINT *, hessian
      PRINT *, 'in newtonstep'
      inversehessian = INVERSE(hessian)
      PRINT *, 'still in newtonstep'
      newtonstep = MATMUL( inversehessian,force)
      newtonstep = MATMUL(car2lat,newtonstep)
    RETURN
    END FUNCTION newtonstep
   
    ! below is the function for least sqare gradient
    FUNCTION lsg(r0,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct) 
      TYPE(charge_obj) :: chg
      REAL(q2), DIMENSION(3) :: lsg
      !INTEGER, DIMENSION(3) :: p
      ! input may either be integers or not. May need to change later on.
      INTEGER, DIMENSION(3) :: r0
      ! r0 is the position of the current grid point
      !INTEGER, DIMENSION(26,3) :: vi , nbp
      INTEGER,DIMENSION(3,26) :: vi, nbp
      ! v is a column of coordinates. 
      ! nbp is indecies of neighbors, to be used for pbc
      INTEGER, DIMENSION(26,3) :: vit
      ! vit is the transpose of vi
      ! vi are the vectors from point r0 to ri, where ri are all neighbors of r0
      REAL(q2), DIMENSION(3) :: frakturavi
      REAL(q2), DIMENSION(26) :: w, di 
      REAL(q2), DIMENSION(3,26) :: matw
      !REAL(q2), DIMENSION(26,3) :: matw
      REAL(q2), DIMENSION(3,13) :: matwprime
      !REAL(q2), DIMENSION(13,3) :: matwprime
      
      REAL(q2), DIMENSION(3,3) ::  matm
      REAL(q2), DIMENSION(26) :: wi
      ! wi are the weights of each neighbor
      REAL(q2), DIMENSION(3,26) :: bi
      ! bi is weight times direction times difference towards each neighbor
      REAL(q2), DIMENSION(3) :: vecb
      REAL(q2), DIMENSION(13) :: deltarho
      ! deltarho is the difference between +ri and -ri
      REAL(q2), DIMENSION(3,3) :: ggrid
      ! ggrid is the metric tensor, from grid basis to dual basis
      INTEGER :: i, j, k

      REAL(q2), DIMENSION(3,3) :: testermatrix

      REAL(q2), DIMENSION(3,3) :: outerproduct
      ! ggrid is 
!      ggrid = makeggrid(chg)
      ! first get all vi.
!      vi = makevi()
!      vit = TRANSPOSE(vi)
!      DO i = 1,26
!        wi(i) = 1/DOT_PRODUCT(MATMUL(vit(i,:),ggrid),vi(:,i))
!      END DO
        ! If locations are written in columns, do lat2car needs transpose too?
        ! Or simply the output in cartesian is also in column order?
      ! get the value differences to all neighbors
      CALL pbc(r0,chg%npts)
      DO i = 1, 26
        nbp(:,i) = vi(:,i) + r0
        CALL pbc(nbp(:,i),chg%npts)
      END DO
!      DO i = 1, 13
!        matwprime(:,i) = vi(:,i) * wi(i)
!      END DO
!      DO i = 1, 26
!        matw(:,i) = wi(i) * vi(:,i)
!      END DO
!      vecb = MATMUL(matw,di)
      DO i = 1 , 13
        deltarho(i) = rho_val(chg,nbp(1,i),nbp(2,i),nbp(3,i)) - &
          rho_val(chg,nbp(1,i+13),nbp(2,i+13),nbp(3,i+13))
      END DO
!      matm = 0.0_q2
!      DO i = 1, 26
!        matm = matm + wi(i)*MATMUL(vi,vit)
!      END DO
!      DO i = 1, 26
!        outerproduct(1,1)= vi(1,i) * vit(i,1)
!        outerproduct(1,2)= vi(1,i) * vit(i,2)
!        outerproduct(1,3)= vi(1,i) * vit(i,3)
!        outerproduct(2,1)= vi(2,i) * vit(i,1)
!        outerproduct(2,2)= vi(2,i) * vit(i,2)
!        outerproduct(2,3)= vi(2,i) * vit(i,3)
!        outerproduct(3,1)= vi(3,i) * vit(i,1)
!        outerproduct(3,2)= vi(3,i) * vit(i,2)
!        outerproduct(3,3)= vi(3,i) * vit(i,3)
!        matm = matm + wi(i) * outerproduct
!      END DO
      lsg = 0._q2
      lsg = &
        MATMUL( &
          MATMUL( & 
            MATMUL(INVERSE(ggrid),INVERSE(matm)), &
            matwprime) & 
        , deltarho)
      ! lsg up till now seems to be larger than force in cartesian units by car2lat
      lsg = MATMUL(lsg,chg%lat2car) ! now force should be in cartesian units

      RETURN
    END FUNCTION lsg

    ! This function finds hessian by 
    ! finding lsg of gradients found with lsg
    FUNCTION lsh(r0,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
      TYPE(charge_obj) :: chg
      TYPE(hessian) :: hes
      REAL(q2), DIMENSION(3,3) :: lsh
      INTEGER, DIMENSION(3) :: r0
      ! to store neighbor lsg's
      REAL(q2), DIMENSION(3,26) :: lsg_val
      ! vi, as usual
      INTEGER, DIMENSION(3,26) :: vi
      INTEGER, DIMENSION(26,3) :: vit
      ! ggrid, as usual
      REAL(q2), DIMENSION(3,3) :: ggrid
      ! one time use for neighbor positions
      INTEGER, DIMENSION(3) :: nbp
      REAL(q2), DIMENSION(26) :: wi
      REAL(q2), DIMENSION(3,13) :: matwprime
 
      ! differences in gradient components
      REAL(q2), DIMENSION(13) :: deltadx, deltady, deltadz
!      REAL(q2), DIMENSION(13) :: deltadxdx, deltadxdy, deltadxdz
!      REAL(q2), DIMENSION(13) :: deltadydx, deltadydy, deltadydz
!      REAL(q2), DIMENSION(13) :: deltadzdx, deltadzdy, deltadzdz
      REAL(q2), DIMENSION(3,3) :: matm
      REAL(q2), DIMENSION(3,3) :: outerproduct
      REAL(q2) :: average
      INTEGER :: i, j, k
!      vi = makevi()
!      vit = TRANSPOSE(vi)
!      ggrid = makeggrid(chg)
      ! find weight of all neighbors
!      DO i = 1, 26
!        wi(i) = 1/DOT_PRODUCT(MATMUL(vit(i,:),ggrid),vi(:,i))
        ! wi matches with lsg
!      END DO
      ! find lsg of all neighbors
!      matm = 0.0_q2
      DO i = 1, 26
        nbp = vi(:,i) + r0
        CALL pbc(nbp,chg%npts) ! nbps are fine
!        CALL getgradhes(nbp,chg,hes,lsg_val(:,i))
!        PRINT *, 'i is', i
!        PRINT *, 'nbp is ', nbp
!        PRINT *, 'default force is ', MATMUL(lsg_val(:,i),chg%car2lat)
        lsg_val(:,i) = lsg(nbp,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
!        PRINT *, 'its lsg is', lsg_val(:,i)
!        matm = matm + wi(i) * MATMUL(vi,vit)

        ! matm matches with lsg
      END DO
!      matm = 0.0_q2

!      DO i = 1, 26
!        outerproduct(1,1)= vi(1,i) * vit(i,1)
!        outerproduct(1,2)= vi(1,i) * vit(i,2)
!        outerproduct(1,3)= vi(1,i) * vit(i,3)
!        outerproduct(2,1)= vi(2,i) * vit(i,1)
!        outerproduct(2,2)= vi(2,i) * vit(i,2)
!        outerproduct(2,3)= vi(2,i) * vit(i,3)
!        outerproduct(3,1)= vi(3,i) * vit(i,1)
!        outerproduct(3,2)= vi(3,i) * vit(i,2)
!        outerproduct(3,3)= vi(3,i) * vit(i,3)
!        matm = matm + wi(i) * outerproduct
!      END DO
      DO i = 1, 13
!        matwprime(:,i) = vi(:,i) * wi(i)
        ! matwprime matches with lsg
        deltadx(i) = lsg_val(1,i) - lsg_val(1,i+13)
        deltady(i) = lsg_val(2,i) - lsg_val(2,i+13)
        deltadz(i) = lsg_val(3,i) - lsg_val(3,i+13)
      END DO
      lsh(1,:) = &
        MATMUL( &
          MATMUL( &
            MATMUL(INVERSE(ggrid),INVERSE(matm)), &
            matwprime &  
          ), &        
        deltadx &
        )
      lsh(2,:) = &
        MATMUL( &
          MATMUL( &
            MATMUL(INVERSE(ggrid),INVERSE(matm)), &
            matwprime &  
          ), &        
        deltady &
        )
      lsh(3,:) = &
        MATMUL( &
          MATMUL( &
            MATMUL(INVERSE(ggrid),INVERSE(matm)), &
            matwprime &  
          ), &        
        deltadz &
        )

      ! the hessian up till now seems larger than cartesian by car2lat
      lsh(1,:) = MATMUL(lsh(1,:),chg%lat2car)
      lsh(2,:) = MATMUL(lsh(2,:),chg%lat2car)
      lsh(3,:) = MATMUL(lsh(3,:),chg%lat2car)

      average = (lsh(1,2) + lsh(2,1))/2
      lsh(1,2) = average
      lsh(2,1) = average
      average = (lsh(1,3) + lsh(3,1))/2
      lsh(1,3) = average
      lsh(3,1) = average
      average = (lsh(2,3) + lsh(3,2)) /2
      lsh(2,3) = average
      lsh(3,2) = average
      RETURN
    END FUNCTION lsh

    FUNCTION makevi()
      INTEGER, DIMENSION(3,26) :: makevi
      makevi(:,1)=(/-1,-1,-1/)
      makevi(:,2)=(/-1,-1,0/)
      makevi(:,3)=(/-1,-1,1/)
      makevi(:,4)=(/-1,0,-1/)
      makevi(:,5)=(/-1,0,0/)
      makevi(:,6)=(/-1,0,1/)
      makevi(:,7)=(/-1,1,-1/)
      makevi(:,8)=(/-1,1,0/)
      makevi(:,9)=(/-1,1,1/)
      makevi(:,10)=(/0,-1,-1/)
      makevi(:,11)=(/0,-1,0/)
      makevi(:,12)=(/0,-1,1/)
      makevi(:,13)=(/0,0,-1/)
      makevi(:,14)=(/1,1,1/)
      makevi(:,15)=(/1,1,0/)
      makevi(:,16)=(/1,1,-1/)
      makevi(:,17)=(/1,0,1/)
      makevi(:,18)=(/1,0,0/)
      makevi(:,19)=(/1,0,-1/)
      makevi(:,20)=(/1,-1,1/)
      makevi(:,21)=(/1,-1,0/)
      makevi(:,22)=(/1,-1,-1/)
      makevi(:,23)=(/0,1,1/)
      makevi(:,24)=(/0,1,0/)
      makevi(:,25)=(/0,1,-1/)
      makevi(:,26)=(/0,0,1/)
      RETURN
    END FUNCTION

    FUNCTION makeggrid(chg)
      TYPE(charge_obj) :: chg
      REAL(q2),DIMENSION(3,3) :: makeggrid
      makeggrid(1,1) = DOT_PRODUCT(chg%lat2car(1,:),chg%lat2car(1,:))
      makeggrid(1,2) = DOT_PRODUCT(chg%lat2car(1,:),chg%lat2car(2,:))
      makeggrid(1,3) = DOT_PRODUCT(chg%lat2car(1,:),chg%lat2car(3,:))
      makeggrid(2,1) = DOT_PRODUCT(chg%lat2car(2,:),chg%lat2car(1,:))
      makeggrid(2,2) = DOT_PRODUCT(chg%lat2car(2,:),chg%lat2car(2,:))
      makeggrid(2,3) = DOT_PRODUCT(chg%lat2car(2,:),chg%lat2car(3,:))
      makeggrid(3,1) = DOT_PRODUCT(chg%lat2car(3,:),chg%lat2car(1,:))
      makeggrid(3,2) = DOT_PRODUCT(chg%lat2car(3,:),chg%lat2car(2,:))
      makeggrid(3,3) = DOT_PRODUCT(chg%lat2car(3,:),chg%lat2car(3,:))
      RETURN
    END FUNCTION

    FUNCTION lsgascension(ind,chg,matm,matwprime,wi,vi,vit, & 
                       ggrid,outerproduct,opts)
      ! this function finds nucleus critical points. 
      INTEGER, DIMENSION(3) :: ind
      TYPE(charge_obj) :: chg
      TYPE(options_obj) :: opts
      REAL(q2),DIMENSION(3) :: lsgascension, distance
      INTEGER, DIMENSION(8,3) :: nnind
      REAL(q2), DIMENSION(8,3) :: nngrad
      INTEGER :: j, stepcount
      REAL(q2), DIMENSION(3) :: tempr, rn, rnm1 ! rn minus 1
      REAL(q2), DIMENSION(3) :: grad, stepsize, gradnm1
      INTEGER, DIMENSION(3,26) :: vi, matw
      INTEGER, DIMENSION(26,3) :: vit
      REAL(q2), DIMENSION(3,3) :: ggrid
      REAL(q2), DIMENSION(26) :: wi
      REAL(q2), DIMENSION(3,13) :: matwprime
      REAL(q2), DIMENSION(3,3) :: matm, outerproduct
      stepcount = 0
      ! initialize the process
      CALL pbc(ind,chg%npts)
      rn(1) = REAL(ind(1),q2)
      rn(2) = REAL(ind(2),q2)
      rn(3) = REAL(ind(3),q2)
      stepsize(1) = 0.5
      stepsize(2) = 0.5
      stepsize(3) = 0.5
      !grad = lsg(ind,chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
      grad = cdgrad(ind,chg)
      ! this gradient is in cartesian. convert it to lattice
      grad = MATMUL(grad,chg%lat2car)
      DO WHILE (stepsize(1) >= 0.1 .AND. &
                stepsize(2) >= 0.1 .AND. &
                stepsize(3) >= 0.1 )
        ! the gradient is in cartesian. 
        gradnm1 = grad
        rnm1 = rn
        ! determine where to go, a unit vector
        tempr = grad / SQRT(SUM(grad*grad))

        tempr(1) = stepsize(1) * tempr(1)
        tempr(2) = stepsize(2) * tempr(2)
        tempr(3) = stepsize(3) * tempr(3)
        rn = rn + tempr
        CALL pbc_r_lat(rn,chg%npts)
        nnind(1,:) = (/FLOOR(rn(1)),FLOOR(rn(2)),    &
                FLOOR(rn(3))/)
        nnind(2,:) = (/CEILING(rn(1)),FLOOR(rn(2)),  &
                FLOOR(rn(3))/)
        nnind(3,:) = (/FLOOR(rn(1)),CEILING(rn(2)),  &
                FLOOR(rn(3))/)
        nnind(4,:) = (/CEILING(rn(1)),CEILING(rn(2)),&
                FLOOR(rn(3))/)
        nnind(5,:) = (/FLOOR(rn(1)),FLOOR(rn(2)),    &
                CEILING(rn(3))/)
        nnind(6,:) = (/CEILING(rn(1)),FLOOR(rn(2)),  &
                CEILING(rn(3))/)
        nnind(7,:) = (/FLOOR(rn(1)),CEILING(rn(2)),  &
                CEILING(rn(3))/)
        nnind(8,:) = (/CEILING(rn(1)),CEILING(rn(2)),&
                CEILING(rn(3))/)
        DO j = 1,8
          CALL pbc(nnind(j,:),chg%npts)
          !IF (opts%leastsquare_flag == .TRUE.) THEN
          IF ( .TRUE. ) THEN
            !nngrad(j,:) = lsg(nnind(j,:),chg,matm,matwprime,wi,vi,vit,ggrid,outerproduct)
            nngrad(j,:) = cdgrad(nnind(j,:),chg)
          ELSE 
          END IF
        END DO
        distance = rn - nnind(1,:)
        ! Row find the nearest neighbors at this new locaiton
        ! First update critical point location
        ! The next big step is to interpolate the force at predicted critical
        ! point.
        grad = trilinear_interpol_grad(nngrad,distance) ! val r interpol
        ! this grad is in cartesian. convert it to lattice
        grad = MATMUL(grad,chg%lat2car)
        ! if grad points backwards, reduce stepsize
        IF (grad(1) * gradnm1(1) < 0) THEN
          stepsize(1) = stepsize(1)/2
        END IF
        IF (grad(2) * gradnm1(2) < 0) THEN
          stepsize(2) = stepsize(2)/2
        END IF
        IF (grad(3) * gradnm1(3) < 0) THEN
          stepsize(3) = stepsize(3)/2
        END IF
        stepcount = stepcount + 1
      END DO 
      lsgascension = rn
      RETURN 


      ! first calculate the 
    END FUNCTION
    

    
  END MODULE

