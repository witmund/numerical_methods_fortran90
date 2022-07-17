! System of linear equations solver powered by näive Gauss elimination
! Double precision is being used

! There's no pivoting in this näive method so mind your matrix
! This code won't work if during the elimination and
! backward substitution phases division by zero occurs.

PROGRAM naivegaussel

IMPLICIT NONE
INTEGER, PARAMETER :: dp=KIND(1.0D0)
INTEGER :: i, j, k, n, m
REAL(dp), ALLOCATABLE :: aug(:,:), a(:,:), b(:), x(:), s(:)

! Number of unknowns
WRITE(*,*) 'Enter the number of unknowns'
READ(*,*) n

! Number of column of the augmented matrix
m = n + 1

ALLOCATE(aug(1:n, 1:m), a(1:n, 1:n), b(1:n), x(1:n), s(1:n))

! Initialization
DO i = 1, n
  b(i) = 0.0_dp
  x(i) = 0.0_dp
  s(i) = 0.0_dp
END DO

DO i = 1, n
  DO j = 1, n
    a(i,j) = 0.0_dp
  END DO
END DO

DO i = 1, n
  DO j = 1, m
    aug(i,j) = 0.0_dp
  END DO
END DO

! Read the (external input) matrix
OPEN(UNIT=10,FILE="matrix.txt",STATUS='OLD',ACTION='READ')
! Write your augmented matrix in a file named "matrix.txt" (without the quotation mark)
! in the same folder (directory) as this .f90 file. You can also place it or name it whatever
! you like if you modify the line above.

! Read the augmented matrix
DO i = 1, n
    READ(10,*,END=100) (aug(i,j), j = 1, m)
END DO

100 CLOSE(10)

! Read the matrix of coefficients
DO i = 1, n
  DO j = 1, n
    a(i,j) = aug(i,j)
  END DO
END DO

! Read the vector of constant terms
DO i = 1, n
    b(i) = aug(i,m)
END DO

DEALLOCATE(aug)

! Print the matrices
WRITE(*,*) "Matrix of coefficients"
DO i = 1, n
  WRITE(*,*) (a(i,j),j=1,n)
END DO

WRITE(*,*) "Vector of constant terms"
DO i = 1, n
  WRITE(*,*) b(i)
END DO

! Elimination
DO k = 1, (n-1)
  DO i = (k+1), n
    b(i) = b(i) - (b(k)*(a(i,k)/a(k,k)))
    DO j = k+1,n
      a(i,j) = a(i,j) - (a(k,j)*(a(i,k)/a(k,k)))
    END DO
  END DO
END DO

! Backward substitution
DO i = n, 1, -1
  DO k = 1, n
    s(i) = (a(i,k)*x(k)) + s(i)
  END DO
  x(i)=(b(i)-s(i))/a(i,i)
END DO

! Solution of the SLE
WRITE(*,*) 'Solutions'
DO i = 1, n
  WRITE(*,*) 'x',i,'=', x(i)
END DO

! Exporting solutions
!OPEN(UNIT=20,FILE="solutions.txt",STATUS='NEW',ACTION='WRITE')
!WRITE(20,*) 'Solutions'
!DO i = 1, n
!  WRITE(20,*) 'x',i,'=', x(i)
!END DO
!CLOSE(20)

DEALLOCATE(a, b, x, s)

STOP

END PROGRAM naivegaussel
