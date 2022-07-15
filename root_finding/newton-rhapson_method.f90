! Root finding program using Newton-Rhapson method
! Double precession is being used

PROGRAM newrhap

IMPLICIT NONE
INTEGER, PARAMETER :: dp = KIND(1.0D0)      ! double precession parameter
REAL(dp), PARAMETER :: es = 1.0D-5          ! stopping criterion
REAL(dp) :: xi, rootold, ea, root
INTEGER :: iter = 0

WRITE(*,*) 'Enter the initial guess'
READ(*,*) xi

DO
	rootold = xi
	root = rootold - (func(rootold)/funcder(rootold))

	ea = DABS((rootold - root)/root)

	IF (ea < es) THEN
		EXIT
	END IF
	xi = root
	iter = iter + 1
END DO

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

! Subprogram for the function first derivative
FUNCTION funcder(x) RESULT(y)
IMPLICIT NONE
REAL(dp), INTENT(IN) :: x
REAL(dp) :: y

y = (3.0_dp*x**2.0_dp) + 5.0_dp       ! first derivative

RETURN
END FUNCTION funcder

END PROGRAM newrhap
