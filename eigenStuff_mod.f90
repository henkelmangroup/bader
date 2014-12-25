  MODULE eigenStuff_mod
    USE kind_mod
    USE matrix_mod
    
    IMPLICIT NONE
    PRIVATE 
    PUBLIC :: find_vector

    CONTAINS

    !----------------------------------------------------------------------
    ! find_vector : find the eigenvector
    !----------------------------------------------------------------------
    SUBROUTINE find_vector(yita2,iDM,dM,s1,s2,v1,v2,v3)
      REAL(q2),INTENT(IN):: yita2 
      REAL(q2),INTENT(IN),DIMENSION(3)::s1,s2,v1
      REAL(q2),INTENT(IN),DIMENSION(3,3)::iDM,dM
      REAL(q2),INTENT(OUT),DIMENSION(3)::v2,v3
      REAL(q2):: temp,norm
      REAL(q2),DIMENSION(3)::tempVec,u1,u2,w1
      REAL(q2),DIMENSION(3,3)::tempMat
      CALL scalar_matrix(yita2,iDM,tempMat)
      CALL matrix_substraction(dM,tempMat,tempMat)
      CALL matrix_vector(tempMat,s1,u1)
      CALL matrix_vector(tempMat,s2,u2)
      norm=1.0_q2/SQRT(u1(1)**2+u1(2)**2+u1(3)**2)
      CALL scalar_vector(norm,u1,w1)
      CALL cross_product(w1,v1,v2)
      CALL cross_product(v1,v2,v3)
      RETURN
      

    END SUBROUTINE 


  END MODULE

    
