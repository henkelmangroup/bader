!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module containing matrix functions
!
! By Andri Arnaldsson and Graeme Henkelman
! Last modified by GH on Feb 14 2003
!-----------------------------------------------------------------------------------!
MODULE matrix_mod
  USE vars_mod , ONLY : q2
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: matrix_multip,matrix_vector,transpose_matrix,volume
!  PUBLIC :: q1,q2,pi,matrix_multip,matrix_vector,transpose_matrix
!  INTEGER,PARAMETER :: q1=SELECTED_REAL_KIND(6,30)         ! single precision
!  INTEGER,PARAMETER :: q2=SELECTED_REAL_KIND(15,305)       ! double precision
!  REAL(q2),PARAMETER :: pi=3.141592653589793238462643_q2

  CONTAINS
!-----------------------------------------------------------------------------------!
! matrix_multip:  Multiply the matricies A and B and return the product C
!-----------------------------------------------------------------------------------!
  SUBROUTINE matrix_multip(A,C,B,l,m,n)
    INTEGER,INTENT(IN) :: l,m,n
    REAL(q2),INTENT(IN),DIMENSION(l,m) :: A
    REAL(q2),INTENT(IN),DIMENSION(m,n) :: C
    REAL(q2),INTENT(OUT),DIMENSION(l,n) :: B

    INTEGER :: j,k

!    WRITE(*,*) A
!    WRITE(*,*) C
!    WRITE(*,*) B
!    pause

    DO j=1,n
      B(:,j)=0.0_q2
      DO k=1,m
        B(:,j)=B(:,j)+c(k,j)*A(:,k)
      END DO
    END DO

  RETURN
  END SUBROUTINE matrix_multip

!-----------------------------------------------------------------------------------!
! matrix_vector:  Multiply the matrix A with the vector c and return the product b
!-----------------------------------------------------------------------------------!
  SUBROUTINE matrix_vector(A,v,b,n,m)
    INTEGER,INTENT(IN) :: n,m
    REAL(q2),INTENT(IN),DIMENSION(n,m) :: A
    REAL(q2),INTENT(IN),DIMENSION(m) :: v
    REAL(q2),INTENT(OUT),DIMENSION(n) :: b

    INTEGER :: i

    b=0.0_q2
    DO i=1,m
      b=b+v(i)*A(:,i)
    END DO

  RETURN
  END SUBROUTINE matrix_vector

!-----------------------------------------------------------------------------------!
! transpose_matrix:  Set matrix B to be the transpose of A
!-----------------------------------------------------------------------------------!
  SUBROUTINE transpose_matrix(A,B,n,m)
    INTEGER,INTENT(IN) :: n,m
    REAL(q2),INTENT(IN),DIMENSION(n,m) :: A
    REAL(q2),INTENT(INOUT),DIMENSION(m,n) :: B

    INTEGER :: i,j

    DO i=1,n
      DO j=1,m
        B(j,i)=A(i,j)
      END DO
    END DO

  RETURN
  END SUBROUTINE transpose_matrix

!------------------------------------------------------------------------------------!
! volume: Function returning the triple product of the lattice vectors.
!------------------------------------------------------------------------------------!
    
  FUNCTION volume(h)
    REAL(q2),INTENT(IN),DIMENSION(3,3) :: h
    REAL(q2) :: volume
    
    volume = h(1,1)*(h(2,2)*h(3,3)-h(2,3)*h(3,2))                                    &
  &       -h(1,2)*(h(2,1)*h(3,3)-h(3,1)*h(2,3))                                      &
  &       +h(1,3)*(h(2,1)*h(3,2)-h(3,1)*h(2,2))

  RETURN
  END FUNCTION volume

!-----------------------------------------------------------------------------------!

END MODULE matrix_mod

