! Root finding program using secant method
! Double precision is being used

PROGRAM secant

IMPLICIT NONE
INTEGER, PARAMETER :: dp = KIND(1.0D0)      ! double precession parameter
REAL(dp), PARAMETER :: es = 1.0D-5          ! stopping criterion
REAL(dp) :: xi1, xi2, rootold1, rootold2, rootold, ea, root
INTEGER :: iter = 0

WRITE(*,*) 'Enter two initial estimates'
WRITE(*,*) "Reminder: This is not a bracketing method!"
READ(*,*) xi1, xi2

DO
	rootold1 = xi2
	rootold2 = xi1
	rootold = ((rootold2*func(rootold1))-(rootold1*func(rootold2)))/(func(rootold1)-func(rootold2))

	ea = DABS((rootold - rootold1)/rootold1)
	root = rootold

	IF (ea < es) THEN
		EXIT
	END IF
	xi2 = rootold
	xi1 = rootold1
	iter = iter + 1		! number of iteration
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

END PROGRAM secant
