! System of linear equations solver powered by LU decomposition
! Double precision is being used

! There's pivoting in this code so you don't have to worry about your matrix

PROGRAM ludecompose

IMPLICIT NONE
INTEGER, PARAMETER :: dp=KIND(1.0D0)
INTEGER :: i, j, k, m, n, p, q
REAL(dp), ALLOCATABLE :: a(:,:), b(:,:), L(:,:), U(:,:), x(:,:), y(:,:)
REAL(dp), ALLOCATABLE :: aug(:,:), L0(:,:), U0(:,:), a0(:,:), b0(:,:)
REAL(dp) :: sum

! Number of unknowns
WRITE(*,*) 'Enter the number of unknowns'
READ(*,*) n

! This program supports multiple systems of linear equations
! i.e. systems with same coefficient matrix of one another,
! but with different vector of constant terms
! Change the line below to adjust it (default = 1)

! Number of systems of linear equations
!WRITE(*,*) 'Enter the number of systems of linear equations'
!READ(*,*) m
m = 1   ! default

! Total column of augmented matrix
p = n + m

ALLOCATE(a(1:n,1:n), b(1:n,1:m), L(1:n, 1:n), U(1:n, 1:n), x(1:n,1:m), y(1:n,1:m),&
& L0(1:n, 1:n), U0(1:n, 1:n), a0(1:n,1:n), b0(1:n,1:m), aug(1:n,1:p))

! Initialization
DO i = 1, n
  DO j = 1, n
    a(i,j) = 0.0_dp
    L(i,j) = 0.0_dp
    U(i,j) = 0.0_dp
    L0(i,j) = 0.0_dp
    U0(i,j) = 0.0_dp
    a0(i,j) = 0.0_dp
  END DO
END DO

DO i = 1, n
  DO j = 1, m
    b(i,j) = 0.0_dp
    x(i,j) = 0.0_dp
    y(i,j) = 0.0_dp
    b0(i,j) = 0.0_dp
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

! Initializing L matrix
DO i = 1, n
  L(i,1) = a(i,1)
END DO

! Initializing U matrix
DO j = 1, n
  U(1,j) = a(1,j)/L(1,1)
END DO

! Completing the L matrix
DO j = 2, n             
77  DO i = j, n        
      sum = 0.0_dp
      DO k = 1, j-1
        sum = L(i,k)*U(k,j) + sum
      END DO
      L(i,j) = a(i,j) - sum

! Row swap?
      IF (L(j,j)== 0.0_dp .AND. j < n) THEN
        L0 = L
        a0 = a
        b0 = b
        a(j,:) = a0(j+1,:)
        b(j,:) = b0(j+1,:)
        L(j,:) = L0(j+1,:)
        a(j+1,:) = a0(j,:)
        b(j+1,:) = b0(j,:)
        L(j+1,:) = L0(j,:)
        GOTO 77       ! Periksa lagi bila terjadi pertukaran
      END IF

! Completing the U matrix
      sum = 0.0_dp
        DO k = 1, i-1
          sum = L(j,k)*U(k,i) + sum
        END DO
      U(j,i) = (a(j,i) - sum)/L(j,j)
    END DO
END DO

! Print L & U matrices
! L matrix
WRITE(*,*) 'L ='
DO i = 1, n
  WRITE(*,*) (L(i,j), j = 1,n)
END DO

! U matrix
WRITE(*,*) 'U ='
DO i = 1, n
  WRITE(*,*) (U(i,j), j = 1,n)
END DO

! Forward substitution
y(1,:) = b(1,:)/L(1,1)
DO q = 1, m
  DO i = 2, n
    sum = 0.0_dp
    DO j = 1, (i-1)
      sum = L(i,j)*y(j,q) + sum
    END DO
    y(i,q) = (b(i,q) - sum)/L(i,i)
  END DO
END DO

! Backward substitution
x(n,:) = y(n,:)
DO q = 1, m
  DO i = 1, n-1
    sum = 0.0_dp
    DO j = (n-i+1), n
      sum = U(n-i,j)*x(j,q) + sum
    END DO
    x(n-i,q) = y(n-i,q) - sum
  END DO
END DO

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

DEALLOCATE(a,b,L,U,x,y,L0,U0,a0,b0)

STOP

END PROGRAM ludecompose
