!in order to try to use type,and complete the work of jisuanfangfa
!so ,I write this code ,though it's complex and unnecessary
MODULE tscr_mod
    USE function_mod
    IMPLICIT NONE
    PRIVATE 
    PUBLIC :: jifen_type
    PUBLIC :: tixing,simpson,cotes,romberg
    TYPE :: jifen_type
    REAL(8) ::T
    REAL(8) ::S
    REAL(8) ::C
    REAL(8) ::R
    ENDTYPE
CONTAINS
    SUBROUTINE tixing(kk,a,b,jifen1,k,jifen2)
        CLASS(jifen_type),INTENT(IN)   ::jifen1
        CLASS(jifen_type),INTENT(INOUT)::jifen2

        REAL(8),INTENT(IN) ::   a,b
        INTEGER(4),INTENT(IN) ::k,kk
        REAL(8)               ::sum,sum1
        INTEGER(4)            ::i
        jifen2%t=0
        sum=0
        sum1=0
        DO i=1,2**k
            CALL fun(kk,(a+(2*i-1)*(b-a)/2**(k+1)),sum1)
            sum=sum+sum1
        END DO
        jifen2%t=0.5*jifen1%t+(b-a)*sum/2**(k+1)
    END SUBROUTINE tixing
    !-------------------------------------
    SUBROUTINE simpson(jifen1,jifen2)
        CLASS(jifen_type),INTENT(IN)::jifen1
        CLASS(jifen_type),INTENT(INOUT)::jifen2
        jifen2%s=jifen2%t+(jifen2%t-jifen1%t)/3.0
    END SUBROUTINE simpson
    !-------------------------------------
    SUBROUTINE cotes(jifen1,jifen2)
        CLASS(jifen_type),INTENT(IN)::jifen1
        CLASS(jifen_type),INTENT(INOUT)::jifen2
        jifen2%c=jifen2%s+(jifen2%s-jifen1%s)/15.0
    END SUBROUTINE cotes    
    !-------------------------------------
    SUBROUTINE romberg(jifen1,jifen2)
        CLASS(jifen_type),INTENT(IN)::jifen1
        CLASS(jifen_type),INTENT(INOUT)::jifen2
        jifen2%r=jifen2%c+(jifen2%c-jifen1%c)/63.0
    END SUBROUTINE romberg      
END MODULE tscr_mod