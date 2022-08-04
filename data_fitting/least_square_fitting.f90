! Least square fitting program powered by Gauss elimination
! Double precision is being used

PROGRAM lesquare

IMPLICIT NONE
INTEGER, PARAMETER :: dp=KIND(1.0D0)  ! double precision parameter
INTEGER(dp) :: i, j, n, m      
REAL(dp), ALLOCATABLE :: x(:), f(:), d(:,:), sol(:)

! Number of data
WRITE(*,*) 'Enter the number of data that you want to fit'
READ(*,*) n

! Order of best fit polynomial
WRITE(*,*) 'Enter the order of best fit polynomial'
READ(*,*) m

ALLOCATE(d(1:n, 1:2), x(1:n), f(1:n))

! Read the data
OPEN(UNIT=10,FILE="data.txt",STATUS='OLD',ACTION='READ')
! Write your data on a file named "data.txt"
! The first column contains the x components
! The second column contains the y components

DO i = 1, n
    READ(10,*,END=100) (d(i,j), j = 1, 2)
END DO

100 CLOSE(10)

! Read the x
DO i = 1, n
  x(i) = d(i,1)
END DO

! Read the y
DO i = 1, n
  f(i) = d(i,2)
END DO

DEALLOCATE(d)

ALLOCATE(sol(0:m))

sol = fit(x,f,n,m)

DEALLOCATE(x,f,sol)

STOP

CONTAINS

! Subprogram for least square fitting
FUNCTION fit(x,f,n,m) RESULT(a)
IMPLICIT NONE
REAL(dp), DIMENSION(1:n), INTENT(IN) :: x, f
REAL(dp), DIMENSION(0:m) :: a, b, b0
REAL(dp), DIMENSION(0:m,0:m) :: c, c0
INTEGER(dp), INTENT(IN) :: n, m
INTEGER(dp) :: k, p      
REAL(dp) :: sum

! Matrix of coefficients for Gauss elimination
DO k = 0, m
  DO j = 0, m
    c(k,j) = 0.0_dp
    DO i = 1, n
      c(k,j) = x(i)**(j+k) + c(k,j)
    END DO
  END DO 
END DO

! Vector of constant terms for Gauss elimination
DO k = 0, m
  b(k) = 0.0_dp
  DO i = 1, n
    b(k) = f(i)*(x(i)**k) + b(k)
  END DO
END DO

! Elimination
DO k = 0, m
  CALL swap(c,b,m)
  c0 = c
  b0 = b
  DO i = k+1, m
    DO j = k, m
      c(i,j) = c0(i,j) - c0(k,j)*(c0(i,k)/c0(k,k))
    END DO
    b(i) = b0(i) - (b0(k))*(c0(i,k)/c0(k,k))
  END DO
END DO

! Backward substitution
DO i = 0, m
  a(i) = b(i)/c(i,i)
END DO

DO j = 0,m
  sum = 0.0_dp
    DO k = (m-j+1),m
      sum = (c(m-j,k)*a(k)) + sum
    END DO
  a(m-j) = (b(m-j) - sum)/c(m-j,m-j)
END DO

! Polynomial coefficients
! Form: y = a(0) + a(1)*x + a(2)*x^2 + ... + a(n)*x^n
WRITE(*,*) 'Polynomial coefficients:'
DO i = 0, m
  WRITE(*,*) 'a',i,'=', a(i)
END DO

! Exporting result
!OPEN(UNIT=30,FILE="fitting.txt",STATUS='NEW',ACTION='WRITE')
!WRITE(30,*) 'Polynomial coefficients:'
!DO i = 0, m
!  WRITE(30,*) 'a',i,'=', a(i)
!END DO
!CLOSE(30)

RETURN
END FUNCTION fit

! Subroutine for row swap
SUBROUTINE swap(a, b, n)
IMPLICIT NONE
INTEGER(dp), INTENT(IN) :: n
REAL(dp), DIMENSION(n,n) :: a, a0
REAL(dp), DIMENSION (n) :: b, b0
INTEGER :: k

! Row swap?
20 DO k = 1, n
      a0 = a
      b0 = b
      IF (a(k,k) == 0) THEN
        IF (a(k+1, k) .NE. 0.0_dp .AND. k < n) THEN
          a(k,:) = a0(k+1,:)
          b(k) = b0(k+1)
          a(k+1,:) = a0(k,:)
          b(k+1) = b0(k)
        ELSE IF (a(k,k+1) .NE. 0.0_dp .AND. k < n) THEN
          a(:,k) = a0(:,k+1)
          b(k) = b0(k+1)
          a(k+1,:) = a0(k,:)
          b(k+1) = b0(k)
        END IF
        GOTO 20
      END IF
    END DO

RETURN
END SUBROUTINE swap

END PROGRAM lesquare
