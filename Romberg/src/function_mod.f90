! in order to caculate different funciont f(x),
! this module is used to define f(x) 
! you can also add other f(x) 
MODULE function_mod
    IMPLICIT NONE
    PRIVATE
    PUBLIC::fun
CONTAINS 
    SUBROUTINE fun(k,x,y)
        REAL(8),INTENT(IN):: x
        REAL(8),INTENT(INOUT) ::y
        INTEGER(4)::k ! k is used to select the function
!-------------------------------------------------------
    SELECT CASE (k)
        CASE(1)
            y=1/(1+x)
        CASE(2)
            y=log(1+x)/(1+x*x)
        CASE(3)
            y=log(1+x)/x
        CASE(4)
            y=sin(x)/x
    END SELECT
!-------------------------------------------------------
    END SUBROUTINE
END MODULE