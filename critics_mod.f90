  MODULE critics_mod
    USE kind_mod
    USE matrix_mod
    USE bader_mod
    USE charge_mod 
    USE options_mod
    USE ions_mod
    USE io_mod
    IMPLICIT NONE

    PRIVATE 
    PUBLIC :: critical_find
    CONTAINS

!    !----------------------------------------------------------------------
!    ! find_vector : find the eigenvector
!    !----------------------------------------------------------------------
!    SUBROUTINE find_vector(yita2,iDM,dM,s1,s2,v1,v2,v3)
!      REAL(q2),INTENT(IN):: yita2 
!      REAL(q2),INTENT(IN),DIMENSION(3)::s1,s2,v1
!      REAL(q2),INTENT(IN),DIMENSION(3,3)::iDM,dM
!      REAL(q2),INTENT(OUT),DIMENSION(3)::v2,v3
!      REAL(q2):: temp,norm
!      REAL(q2),DIMENSION(3)::tempVec,u1,u2,w1
!      REAL(q2),DIMENSION(3,3)::tempMat
!      CALL scalar_matrix(yita2,iDM,tempMat)
!      CALL matrix_substraction(dM,tempMat,tempMat)
!      CALL matrix_vector(tempMat,s1,u1)
!      CALL matrix_vector(tempMat,s2,u2)
!      norm=1.0_q2/SQRT(u1(1)**2+u1(2)**2+u1(3)**2)
!      CALL scalar_vector(norm,u1,w1)
!      CALL cross_product(w1,v1,v2)
!      CALL cross_product(v1,v2,v3)
!      RETURN
!      
!
!    END SUBROUTINE 

!-----------------------------------------------------------------------------------!
!critical_find: find critical points on the edge
!NOTE: this subroutine should be called after refine_edge
!      inorder to restrict the calculation on edge points only
!      this subroutine needs to use the edge points found after refine edge
!-----------------------------------------------------------------------------------!
  SUBROUTINE critical_find(bdr,chg,opts)


    TYPE hessian
      REAL(q2),ALLOCATABLE,DIMENSION(:,:,:) :: rho,dx,dy,dz
      REAL(q2) :: dxdx,dydy,dzdz,dxdy,dxdz,dydz,r1,r2,r3,eigval1,eigval2,eigval3
      ! eigval and eigvec are eigenvalues and eigvectors of hessian matrix.
    END TYPE


    TYPE(hessian) :: hes
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    TYPE(options_obj) :: opts
    TYPE(charge_obj) :: chgtemp
    TYPE(ions_obj) :: ionstemp
    TYPE(charge_obj) :: chgval
! ptoo is another point, for the sake of running pbc to make sure that we do not
! select a point outside of the box. usage is same as pt.


! for points, 1 and 2 are +1, -1
    INTEGER,DIMENSION(3) :: p,pt,ptt,ptx1,ptx2,pty1,pty2,ptz1,ptz2
    REAL(q2),DIMENSION(3,3) :: dM ! the deviatoric matrix
    REAL(q2),DIMENSION(3,3) :: dMSQ,dMVSS ! dM squared and dM in vss basis
    REAL(q2):: j2,j3 ! second and third invariant of dM
    INTEGER :: n1,n2,n3,path_volnum,bvolnum,i,i2
    INTEGER :: num_edge,num_reassign,num_check
    INTEGER :: d1,d2,d3,dummy,switch
    INTEGER :: lc1,lc2,lc3
    REAL(q2)::x,y,z,xx,yy,zz,xy,xz,yz,denom,nomx,& 
nomy,nomz,ODS,trace,num1,num2,num3,phi,PI
    REAL(q2):: traceOver3,temp,alpha
    ! variables for degenerate eigenvalues.
    REAL(q2):: yita1,yita2,yita3
    REAL(q2):: MDE !Most Distinct Eigenvalue
    ! ODS and trace are the off diagonal sum and trace of the hessian matrix.
    REAL(q2),DIMENSION(3) ::eigvec1,eigvec2,eigvec3,cartX,cartY,cartZ,tempVec
    REAL(q2),DIMENSION(3,3) :: nIdentity,identityM,devSubNIden
    ! these are vectors orthogonal to eigenvectors
    REAL(q2),DIMENSION(3) :: orthoR1,orthoR2,orthoR3,S1,S2,orT2,orT3
    REAL(q2):: norm,s1r2,s1r3


    PRINT *, ' '//achar(27)//'[35;40;1m Doing Honest &
Critical Finding Buisiness'//achar(27)//'[0m.'
    ALLOCATE (hes%rho(chg%npts(1),chg%npts(2),chg%npts(3)))
    ALLOCATE (hes%dx(chg%npts(1),chg%npts(2),chg%npts(3)))
    ALLOCATE (hes%dy(chg%npts(1),chg%npts(2),chg%npts(3)))
    ALLOCATE (hes%dz(chg%npts(1),chg%npts(2),chg%npts(3)))
    IF (opts%ref_flag) THEN
      CALL read_charge_ref(ionstemp,chgtemp,opts)
    ELSE
      chgtemp = chgval
    END IF
    ! constants section
    PI=4.0*ATAN(1.0)
    switch=0
    dummy=0
    DO n1=1,3
      DO n2=1,3
        identityM(n1,n2)=0.0_q2
      END DO
    END DO
    identityM(1,1)=1.0_q2
    identityM(2,2)=1.0_q2
    identityM(3,3)=1.0_q2
    DO n1=1,3
      cartX(n1)=0.0_q2
      cartY(n1)=0.0_q2
      cartZ(n1)=0.0_q2
    END DO
    cartX(1)=1.0_q2
    cartY(2)=1.0_q2
    cartZ(3)=1.0_q2

!****************************************

    DO n1 = 1,chg%npts(1)
      DO n2 = 1,chg%npts(2)
        DO n3 = 1,chg%npts(3)
! the problem why this subroutine doesnt work previously is that the object p
! was not properly declared
! holy cow now, segmentation fault starts showing up. could it be the fact that
! this p object is global?

            p=(/n1,n2,n3/)
            If (rho_val(chg,p(1),p(2),p(3))<=opts%vacval) THEN
              CYCLE
            END IF
!            IF (is_vol_edge(bdr,chg,p) .AND. (.NOT.is_max(chg,p))) THEN

!-----------------------------------------------------------------------------------!
! now that this subroutine can find the correct amount of edge points, lets have
! it find the hessian
!-----------------------------------------------------------------------------------!

! this loop finds all the neighboring points that will be needed and stores its
! rho
             DO d1=-2,2
               DO d2=-2,2
                 DO d3=-2,2
                   pt=p+(/d1,d2,d3/)
                   CALL pbc(pt,chg%npts)
                   ptt=p+(/d1,d2,d3/)
                   CALL pbc(ptt,chg%npts)
                   hes%rho(pt(1),pt(2),pt(3))= rho_val(chg,pt(1),pt(2),pt(3))
                 END DO
               END DO
             END DO


! now there is a problem.  ptt(1)+1 may get out of boundary. How to do that
! boundary check? --- by creating points for every point that could possibily be
! needed and run pbc on it. Freezing is because of points out of boundary.
! However, the compiler flags seems to not catch this problem.



             DO d1=-1,1
               DO d2=-1,1
                 DO d3=-1,1
                   ptt=p+(/d1,d2,d3/)
                   CALL pbc(ptt,chg%npts)
                   ptx1=ptt+(/1,0,0/)
                   ptx2=ptt+(/-1,0,0/)
                   pty1=ptt+(/0,1,0/)
                   pty2=ptt+(/0,-1,0/)
                   ptz1=ptt+(/0,0,1/)
                   ptz2=ptt+(/0,0,-1/)
                   CALL pbc(ptx1,chg%npts)
                   CALL pbc(ptx2,chg%npts)
                   CALL pbc(pty1,chg%npts)
                   CALL pbc(pty2,chg%npts)
                   CALL pbc(ptz1,chg%npts)
                   CALL pbc(ptz2,chg%npts)
                   hes%dx(ptt(1),ptt(2),ptt(3))= 0.5 *&
 (hes%rho(ptx1(1),ptt(2),ptt(3)) - &
hes%rho(ptx2(1),ptt(2),ptt(3)) )
                   hes%dy(ptt(1),ptt(2),ptt(3))= 0.5 *&
 (hes%rho(ptt(1),pty1(2),ptt(3)) - &
hes%rho(ptt(1),pty2(2),ptt(3)) )
                   hes%dz(ptt(1),ptt(2),ptt(3))= 0.5 *&
(hes%rho(ptx1(1),ptt(2),ptz1(3)) - &
hes%rho(ptt(1),ptt(2),ptz2(3)) )
                 END DO
               END DO
             END DO
             ptx1=p+(/1,0,0/)
             ptx2=p+(/-1,0,0/)
             pty1=p+(/0,1,0/)
             pty2=p+(/0,-1,0/)
             ptz1=p+(/0,0,1/)
             ptz2=p+(/0,0,-1/)
             CALL pbc(ptx1,chg%npts)
             CALL pbc(ptx2,chg%npts)
             CALL pbc(pty1,chg%npts)
             CALL pbc(pty2,chg%npts)
             CALL pbc(ptz1,chg%npts)
             CALL pbc(ptz2,chg%npts)

             hes%dx(p(1),p(2),p(3))= 0.5 * ( hes%rho(ptx1(1),p(2),p(3)) - &
 hes%rho(ptx2(1),p(2),p(3)) )
             hes%dy(p(1),p(2),p(3))= 0.5 * ( hes%rho(p(1),pty1(2),p(3)) - &
 hes%rho(p(1),pty2(2),p(3)) )
             hes%dz(p(1),p(2),p(3))= 0.5 * ( hes%rho(p(1),p(2),ptz1(3)) - &
 hes%rho(p(1),p(2),ptz2(3)) )
! now calculate the acceleration on the point p
             hes%dxdx=0.5_q2 *&
(hes%dx(ptx1(1),p(2),p(3))-hes%dx(ptx2(1),p(2),p(3)))
             hes%dydy=0.5_q2 *&
(hes%dy(p(1),pty1(2),p(3))-hes%dy(p(1),pty2(2),p(3)))
             hes%dzdz=0.5_q2 *&
(hes%dz(p(1),p(2),ptz1(3))-hes%dz(p(1),p(2),ptz2(3)))
             hes%dxdy=0.25_q2 *  (hes%rho(ptx1(1),pty1(2),p(3))-&
hes%rho(ptx2(1),pty1(2),p(3))-hes%rho(ptx1(1),pty2(2),p(3))+hes%rho(ptx2(1),pty2(2),p(3)))
             hes%dxdz=0.25_q2 * (hes%rho(ptx1(1),p(2),ptz1(3))-&
hes%rho(ptx2(1),p(2),ptz1(3))-hes%rho(ptx1(1),p(2),ptz2(3))+hes%rho(ptx2(1),p(2),ptz2(3)))
             hes%dydz=0.25_q2 * (hes%rho(p(1),pty1(2),ptz1(3))-&
hes%rho(p(1),pty2(2),ptz1(3))-hes%rho(p(1),pty1(2),ptz2(3))+hes%rho(p(1),pty2(2),ptz2(3)))
!solutions for the vector
!x
!(-dxdz dy dydz + dx dydz^2 + dxdz dydy dz -  dxdy dydz dz + dxdy dy dzdz - dx
!dydy dzdz)/(dxdz^2 dydy - 2 dxdy dxdz dydz + dxdx dydz^2 + dxdy^2 dzdz -  dxdx
!dydy dzdz)
!y
!(dxdz^2 dy - dx dxdz dydz - dxdy dxdz dz + dxdx dydz dz +  dx dxdy dzdz - dxdx
!dy dzdz)/(dxdz^2 dydy - 2 dxdy dxdz dydz + dxdx dydz^2 + dxdy^2 dzdz -  dxdx
!dydy dzdz)
!z
!(-dxdy dxdz dy + dx dxdz dydy - dx dxdy dydz + dxdx dy dydz +  dxdy^2 dz - dxdx
!dydy dz)/(dxdz^2 dydy - 2 dxdy dxdz dydz + dxdx dydz^2 + dxdy^2 dzdz -  dxdx
!dydy dzdz)
               x=hes%dx(p(1),p(2),p(3))
               y=hes%dy(p(1),p(2),p(3))
               z=hes%dz(p(1),p(2),p(3))
               xx=hes%dxdx
               yy=hes%dydy
               zz=hes%dzdz
               xy=hes%dxdy
               xz=hes%dxdz
               yz=hes%dydz
               denom= xz*xz*yy - 2.0_q2*xy*xz*yz + xx*yz*yz + xy*xy*zz -&
xx*yy*zz
               nomx=-xz*y*yz + x*yz*yz + xz*yy*z - xy*yz*z + xy*y*zz - x*yy*zz
               nomy=xz*xz*y - x*xz*yz - xy*xz*z + xx*yz*z + x*xy*zz - xx*y*zz
               nomz=-xy*xz*y + x*xz*yy - x*xy*yz + xx*y*yz + xy*xy*z - xx*yy*z
               hes%r1=nomx/denom
               hes%r2=nomy/denom
               hes%r3=nomz/denom

               IF (ABS(hes%r1)<=0.5_q2*bdr%stepsize) THEN
                 IF (ABS(hes%r2)<=0.5_q2*bdr%stepsize) THEN
                   IF (ABS(hes%r3)<=0.5_q2*bdr%stepsize) THEN
                     dummy=dummy+1
                     ! below are eigenstuff
                     trace=hes%dxdx+hes%dydy+hes%dzdz
                       traceOver3=(hes%dxdx+hes%dydy+hes%dzdz)/3.0_q2
                       dM(1,1)=hes%dxdx-traceOver3
                       dM(2,2)=hes%dydy-traceOver3
                       dM(3,3)=hes%dzdz-traceOver3
                       dM(1,2)=hes%dxdy
                       dM(1,3)=hes%dxdz
                       dM(2,1)=hes%dxdy
                       dM(2,3)=hes%dydz
                       dM(3,1)=hes%dxdz
                       dM(3,2)=hes%dydz
                       ! Now a loop to calculate the dMSQ
                       temp=0.
                       ! this loop is functional
                       CALL matrix_mult(dM,dM,dMSQ)
                       ! j2 is 1/2 tr(dM dot dM)
                       j2=0.5_q2*(dMSQ(1,1)+dMSQ(2,2)+dMSQ(3,3))
                       ! j3 is det(dM)
                       CALL det(dM,j3)
                       alpha=&
ACOS(j3/2.0_q2*(3.0_q2/j2)**(3.0_q2/2.0_q2))/3.0_q2
                       yita1=2.0_q2*SQRT(j2/3.0_q2)*COS(alpha)
                       ! Find the most distinct eigenvalue first
                       ! every eigenvalue are equaly distinct when alpha = PI/6
                       ! can use the code when alpha < PI/6
                       IF (alpha>PI/6.0_q2) THEN
                           ! Start with eigval 3 first
                           temp=hes%eigval1
                           hes%eigval1=hes%eigval3
                           hes%eigval3=temp
                           temp=0
                       END IF
                       ! this is the default when yita1 is the most distinct
                       CALL scalar_matrix(yita1,identityM,nIdentity)
                       CALL matrix_substraction(dM,nIdentity,devSubNIden)
                       CALL matrix_vector(devSubNIden,cartX,orthoR1)
                       CALL matrix_vector(devSubNIden,cartY,orthoR2)
                       CALL matrix_vector(devSubNIden,cartZ,orthoR3)
                       ! assume r1 is the largest, normalize all orthoR
                       ! vectors
                       norm=1.0_q2/SQRT(orthoR1(1)**2+orthoR1(2)**2+orthoR1(3)**2)
                       CAll scalar_vector(norm,orthoR1,S1)
                       CALL scalar_vector(DOT_PRODUCT(S1,orthoR2),S1,tempVec)
                       CALL vector_substraction(orthoR2,tempVec,orT2)
                       CALL scalar_vector(DOT_PRODUCT(S1,orthoR3),S1,tempVec)
                       CALL vector_substraction(orthoR3,tempVec,orT3)
                       ! assume t2 is the larger one, normalize it
                       norm=1.0_q2/SQRT(orT2(1)**2+orT2(2)**2+orT2(3)**2)
                       CALL scalar_vector(norm,orT2,S2)
                       ! eigenvector 1 is the cross product of s1 s2
                       CALL cross_product(S1,S2,eigvec1)
                       ! now that I have a eigenvalue and a eigenvector,
                       ! write the deviatoric matrix using v1,s1,s2 basis
                       ! hes%eigval1 |         0        |      0
                       ! 0           | s1 dot dM dot s1 | s1 dot dM dot s2
                       ! 0           | s2 dot dM dot s1 | s2 dot dM dot s2
                       DO lc1=1,3
                         DO lc2=1,3
                           dMVSS(lc1,lc2)=0
                         END DO
                       END DO
                       dMVSS(1,1)=yita1
                       CALL vector_matrix(S1,dM,tempVec)
                       dMVSS(2,2)=DOT_PRODUCT(tempVec,S1)
                       dMVSS(2,3)=DOT_PRODUCT(tempVec,S2)
                       CALL vector_matrix(S2,dM,tempVec)
                       dMVSS(3,2)=DOT_PRODUCT(tempVec,S1)
                       dMVSS(3,3)=DOT_PRODUCT(tempVec,S2)
                       ! This is a sign function right below
                       IF (dMVSS(2,2)-dMVSS(3,3)<0) THEN
                         temp=-1
                         ELSE IF (dMVSS(2,2)-dMVSS(3,3)>0) THEN
                         temp=1
                         ELSE IF (dMVSS(2,2)-dMVSS(3,3)==0) THEN
                         temp=0
                       END IF
                         yita2=(dMVSS(2,2)+dMVSS(3,3))/2.0_q2-1.0_q2/2.0_q2*temp*&
SQRT((dMVSS(2,2)-dMVSS(3,3))**2+4*dMVSS(2,3)*dMVSS(3,2))
                       yita3=dMVSS(2,2)+dMVSS(3,3)-yita2
                       ! these eigenvalues are shifted by traceOver3
                       hes%eigval1=yita1+traceOver3
                       hes%eigval2=yita2+traceOver3
                       hes%eigval3=yita3+traceOver3

                       CALL &
find_vector(yita2,identityM,dM,s1,s2,eigvec1,eigvec2,eigvec3)

                     OPEN(97,FILE="critics")
                     WRITE(97,*)'*********** A NEW ENTRY *************'
                     WRITE(97,*),dummy
                     WRITE(97,*), p(1),p(2),p(3)
                     WRITE(97,*), 'The treshhold is ', 0.5_q2*bdr%stepsize
                     WRITE(97,*), 'r1' , hes%r1
                     WRITE(97,*), 'r2' , hes%r2
                     WRITE(97,*), 'r3' , hes%r3
                     WRITE(97,'(3(1X,E18.11))'),0,&
hes%rho(ptz1(1),ptz1(2),ptz1(3)),hes%rho(ptx2(1),ptx2(2),ptx2(3))
                     WRITE(97,'(3(1X,E18.11))'),hes%rho(pty2(1),pty2(2),pty2(3)),hes%rho(p(1),p(2),p(3)),hes%rho(pty1(1),pty1(2),pty1(3))
                     WRITE(97,'(3(1X,E18.11))'),hes%rho(ptx1(1),ptx1(2),ptx1(3)),hes%rho(ptz2(1),ptz2(2),ptz2(3)),0
                     WRITE(97,*),'---------------------------'
                     WRITE(97,*),'The hessian matrix is'
                     WRITE(97,'(3(1X,E18.11))'),hes%dxdx,hes%dxdy,hes%dxdz
                     WRITE(97,'(3(1X,E18.11))'),hes%dxdy,hes%dydy,hes%dydz
                     WRITE(97,'(3(1X,E18.11))'),hes%dxdz,hes%dydz,hes%dzdz
                     WRITE(97,*),'the eigenvalues are '
                     WRITE(97,'(3(1X,E18.11))'),hes%eigval1,hes%eigval2,hes%eigval3
                     WRITE(97,*),'eigenvector 1'
                     WRITE(97,*),eigvec1
                     WRITE(97,*),'eigenvector 2'
                     WRITE(97,*),eigvec2
                     WRITE(97,*),'eigvenctor 3'
                     WRITE(97,*),eigvec3
                     CLOSE(97)
                   END IF
                 END IF
               END IF

        END DO
      END DO
    END DO

    PRINT *, "FOUND STATIONARY POINT NUM: ", dummy

    DEALLOCATE (hes%rho)
    DEALLOCATE (hes%dx)
    DEALLOCATE (hes%dy)
    DEALLOCATE (hes%dz)
    END SUBROUTINE critical_find





  END MODULE

    
