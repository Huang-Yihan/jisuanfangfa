!this module is use to solve the three diag matrix Ax=b
MODULE threediag_mod
    IMPLICIT NONE
    PUBLIC threediag
CONTAINS
    SUBROUTINE threediag(n,a,b,c,d,x)
        INTEGER(4),INTENT(IN)::n
        REAL(8),DIMENSION(:),INTENT(IN)::a(n),b(n),c(n),d(n)
        REAL(8),DIMENSION(:),INTENT(INOUT)::x(n)
        REAL(8),ALLOCATABLE::y(:),u(:),l(:)
        INTEGER(4)::i
        ALLOCATE(y(n))
        ALLOCATE(u(n))
        ALLOCATE(l(n))
        x=0
        y=0
        u=0
        l=0
        u(1)=b(1)
        y(1)=d(1)
        DO i=2,n
            l(i)=a(i)/u(i-1)
            u(i)=b(i)-l(i)*c(i-1)
            y(i)=d(i)-l(i)*y(i-1)
        END DO
        x(n)=y(n)/u(n)
        DO i=n-1,1,-1
            x(i)=(y(i)-c(i)*x(i+1))/u(i)
        END DO
        DEALLOCATE(y)
        DEALLOCATE(u)
        DEALLOCATE(l)
    END SUBROUTINE
END MODULE