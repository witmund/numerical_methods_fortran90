! System of linear equations solver powered by Gauss elimination
! Double precision is being used

! There's pivoting in this code so you don't have to worry about your matrix

PROGRAM gaussel

IMPLICIT NONE
INTEGER, PARAMETER :: dp=KIND(1.0D0)
INTEGER :: i, j, k, n, m, p, q
REAL(dp), ALLOCATABLE :: a(:,:), b(:,:), x(:,:), a0(:,:), b0(:,:), aug(:,:)
REAL(dp) :: sum

! Number of unknowns
WRITE(*,*) 'Enter the number of unknowns'
READ(*,*) n

! This program supports multiple systems of linear equations
! i.e. systems with same coefficient matrix of one another,
! but with different vector of constant terms
! Change the line below to adjust it (default = 1)

! Number of systems of linear equations
WRITE(*,*) 'Enter the number of systems of linear equations'
!READ(*,*) m
m = 1   ! default
WRITE(*,*) m

! Total column of augmented matrix
p = n + m

ALLOCATE(a(1:n, 1:n), b(1:n,1:m), a0(1:n,1:n), b0(1:n,1:m), x(1:n,1:m), aug(1:n,1:p))

! Initialization
DO i = 1, n
  DO j = 1, n
    a(i,j) = 0.0_dp
    a0(i,j) = 0.0_dp
  END DO
END DO

DO i = 1, n
  DO j = 1, m
    b(i,j) = 0.0_dp
    b0(i,j) = 0.0_dp
    x(i,j) = 0.0_dp
  END DO
END DO

DO i = 1, n
  DO j = 1, p
    aug(i,j) = 0.0_dp
  END DO
END DO

! Read the (external input) matrix
OPEN(UNIT=10,FILE="matrix.txt",STATUS='OLD',ACTION='READ')
! Write your augmented matrix in a file named "matrix.txt" (without the quotation mark)
! in the same folder (directory) as this .f90 file. You can also place it or name it whatever
! you like if you modify the line above.
! Structure: <A|b>
! where A is your coefficients matrix and b is your vector(s) of constant terms

! Read the augmented matrix
DO i = 1, n
    READ(10,*,END=100) (aug(i,j), j = 1, p)
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
  DO j = 1, m
    DO k = (n+1), p
      b(i,j) = aug(i,k)
    END DO
  END DO
END DO

DEALLOCATE(aug)

! Print the matrices
WRITE(*,*) "Matrix of coefficients"
DO i = 1, n
  WRITE(*,*) (a(i,j),j=1,n)
END DO

WRITE(*,*) "Vector of constant terms"
DO i = 1, n
  DO j = 1, m
    WRITE(*,*) b(i,j)
  END DO
END DO

! Elimination
DO k = 1, (n-1)
  CALL swap(a,b,n,m)
  a0 = a
  b0 = b
  DO i = k+1, n
    DO j = k, n
    a(i,j) = a0(i,j) - a0(k,j)*(a0(i,k)/a0(k,k))
    END DO
    DO q = 1,m
      b(i,q) = b0(i,q) - (b0(k,q)*(a0(i,k)/a0(k,k)))
    END DO
  END DO
END DO

DEALLOCATE(a0,b0)

! Backward substitution
DO i = 1, n
  DO j = 1, m
    x(i,j) = b(i,j)/a(i,i)
  END DO
END DO

DO j = 1,n-1
  DO q = 1, m
    sum = 0.0_dp
    DO k = (n-j+1),n
      sum = (a(n-j,k)*x(k,q)) + sum
    END DO
  x(n-j,q) = (b(n-j,q) - sum)/a(n-j,n-j)
  END DO
END DO

DEALLOCATE(a,b)

! Solutions
WRITE(*,*) 'Solutions'
DO i = 1, n
  DO j = 1,m
  WRITE(*,*) 'x',i,j,'=', x(i,j)
  END DO
END DO

! Exporting solutions
!OPEN(UNIT=20,FILE="solutions.txt",STATUS='NEW',ACTION='WRITE')
!WRITE(20,*) 'Solutions'
!DO i = 1, n
!  DO j = 1,m
!    WRITE(20,*) 'x',i,j,'=', x(i,j)
!  END DO
!END DO
!CLOSE(20)

DEALLOCATE(x)

STOP

CONTAINS

! Subroutine for row swap
SUBROUTINE swap(a, b, n, m)
IMPLICIT NONE
INTEGER, INTENT(IN) :: n, m
REAL(dp), DIMENSION(n,n) :: a, a0
REAL(dp), DIMENSION (n,m) :: b, b0
INTEGER :: k

! Row swap?
88 DO k = 1, n
      a0 = a
      b0 = b
      IF (a(k,k) == 0) THEN
        IF (a(k+1, k) .NE. 0.0_dp .AND. k < n) THEN
          a(k,:) = a0(k+1,:)
          b(k,:) = b0(k+1,:)
          a(k+1,:) = a0(k,:)
          b(k+1,:) = b0(k,:)
        ELSE IF (a(k,k+1) .NE. 0.0_dp .AND. k < n) THEN
          a(:,k) = a0(:,k+1)
          b(k,:) = b0(k+1,:)
          a(k+1,:) = a0(k,:)
          b(k+1,:) = b0(k,:)
        END IF
        GOTO 88
      END IF
    END DO

RETURN
END SUBROUTINE swap

END PROGRAM gaussel
