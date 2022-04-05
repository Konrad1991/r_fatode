SUBROUTINE solver()

	use ISO_C_BINDING, only : C_DOUBLE, C_PTR
	IMPLICIT NONE
	integer var

	INTERFACE
		SUBROUTINE Rcppfct(inp) BIND(C, name ='Rcppfct_')
			integer inp
		END SUBROUTINE Rcppfct
	END INTERFACE

	var = 23

	call Rcppfct(var)

END SUBROUTINE
