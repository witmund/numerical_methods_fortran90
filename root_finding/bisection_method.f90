! Root finding program using bisection method (Bolzano's method)
! Double precession is used 

PROGRAM bisection

IMPLICIT NONE
INTEGER, PARAMETER :: dp = KIND(1.0D0)      ! double precession parameter
REAL(dp), PARAMETER :: es = 1.0D-5          ! stopping criterion
REAL(dp) :: xl, xu, xr, rootold, ea, root
INTEGER :: iter = 0

5   WRITE(*,*) 'Enter two initial guesses in which the root is expected to be found'
    READ(*,*) xl, xu

IF (func(xl)*func(xu) == 0.0_dp) THEN
    IF ((func(xl) == 0.0_dp) .AND. (func(xu) == 0.0_dp)) THEN
        WRITE(*,*) 'root =', xl
        WRITE(*,*) 'root =', xu
    ELSE IF (func(xl) == 0.0_dp) THEN
        WRITE(*,*) 'root =', xl
    ELSE IF (func(xu) == 0.0_dp) THEN
        WRITE(*,*) 'root =', xu
    END IF
ELSE IF (func(xl)*func(xu) > 0.0_dp) THEN
    WRITE(*,*) "Error!"
    WRITE(*,*) "The root does not lie between", xl, " and ", xu
    WRITE(*,*) "Reinitializing..."
    GOTO 5
ELSE IF (func(xl)*func(xu) < 0.0_dp) THEN
    DO
        xr = (xl + xu)/2.0_dp
        rootold = xr
        IF (func(xl)*func(xr) < 0.0_dp) THEN
            xu = xr
        ELSE
            xl = xr
        END IF

        xr = (xl + xu)/2.0_dp
        root = xr
        ea = DABS((rootold - root)/root)

        IF(ea < es) THEN
            EXIT
        END IF
        rootold = root
        iter = iter + 1     ! number of iteration
    END DO
END IF

WRITE(*,*) 'root =', root
WRITE(*,*) 'total iterations =', iter
WRITE(*,*) 'approximate relative error  =', ea

STOP

CONTAINS

! Subprogram for the function you wish to find its root
FUNCTION func(x) RESULT(y)
IMPLICIT NONE
REAL(dp), INTENT(IN) :: x
REAL(dp) :: y

y = (x**3.0_dp) + (5.0_dp*x) + 7.0_dp   ! the function you wish to find its root

RETURN
END FUNCTION func

END PROGRAM bisection
