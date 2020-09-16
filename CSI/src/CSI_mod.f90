!this module is used to calculate the M and h in the CSI 
!this module have a function to calculate f[x1,x2,x3]
MODULE csi_mod
    USE threediag_mod
    IMPLICIT NONE
    PRIVATE
    PUBLIC::CSI
CONTAINS
    SUBROUTINE CSI(n_in,x,y,condition,a_in,b_in,M,h)
        INTEGER(4),INTENT(IN):: n_in
        REAL(8),DIMENSION(:),INTENT(IN)::x(n_in),y(n_in)
        CHARACTER(10),INTENT(IN)::condition
        REAL(8),INTENT(IN)::a_in,b_in
        REAL(8),ALLOCATABLE::miu(:),lamda(:),d(:)
        REAL(8),DIMENSION(:),INTENT(OUT)::M(n_in),h(n_in-1)
        REAL(8),ALLOCATABLE::a_three(:),b_three(:),c_three(:),d_three(:),M_three(:)
        INTEGER(4)::i
        M=0
        h=0
        ALLOCATE(miu(n_in-2))
        ALLOCATE(lamda(n_in-2))
        ALLOCATE(d(n_in))
        miu=0
        lamda=0
        d=0
        DO i=2,n_in
            h(i-1)=x(i)-x(i-1) !calculate h
        END DO
        DO i=2,n_in-1
            d(i)=6*diff3(x(i-1),x(i),x(i+1),y(i-1),y(i),y(i+1))
        END DO
        DO i=1,n_in-2
            miu(i)=h(i)/(h(i)+h(i+1))
            lamda(i)=1-miu(i)
        END DO
        SELECT CASE(condition)
            CASE('first') !the first condition
                M(1)=a_in
                M(n_in)=b_in
                d(2)=d(2)-miu(1)*M(1)
                d(n_in-1)=d(n_in-1)-lamda(n_in-1)*M(n_in)
                ALLOCATE(a_three(n_in-2))
                ALLOCATE(b_three(n_in-2))
                ALLOCATE(c_three(n_in-2))
                ALLOCATE(d_three(n_in-2))
                ALLOCATE(M_three(n_in-2)) !define a,b,c,d,m to use the threediag Ax=b
                a_three=0
                b_three=2     !b=2
                c_three=0
                d_three=0
                M_three=0
                c_three(1)=lamda(1)
                d_three(1)=d(2)
                a_three(n_in-2)=miu(n_in-2)
                d_three(n_in-2)=d(n_in-1)
                DO i=2,n_in-3
                    a_three(i)=miu(i)
                    c_three(i)=lamda(i)
                    d_three(i)=d(i+1)
                END DO
                CALL threediag(n_in-2,a_three,b_three,c_three,d_three,&
                                M_three)
                DO i=1,n_in-2
                    M(i+1)=M_three(i) !get the M
                END DO
                DEALLOCATE(a_three)
                DEALLOCATE(b_three)
                DEALLOCATE(c_three)
                DEALLOCATE(d_three)
                DEALLOCATE(M_three)
            CASE('second') !the second condition
                d(1)=(6/h(1))*((y(2)-y(1))/h(1)-a_in)
                d(n_in)=(6/h(n_in-1))*(b_in-(y(n_in)-y(n_in-1))/h(n_in-1))
                ALLOCATE(a_three(n_in))
                ALLOCATE(b_three(n_in))
                ALLOCATE(c_three(n_in))!define a,b,c to use the threediag Ax=b
                a_three=0
                b_three=2     !b=2
                c_three=0
                c_three(1)=1
                a_three(n_in)=1
                DO i=2,n_in-1
                    a_three(i)=miu(i-1)
                    c_three(i)=lamda(i-1)
                END DO
                CALL threediag(n_in,a_three,b_three,c_three,d,&
                                M)                
                DEALLOCATE(a_three)
                DEALLOCATE(b_three)
                DEALLOCATE(c_three)
            CASE('third')
                ALLOCATE(a_three(n_in-1))
                ALLOCATE(b_three(n_in-1))
                ALLOCATE(c_three(n_in-1))!define a,b,c,d,m to use the threediag Ax=b
                b_three=2
                a_three(1:n_in-2)=miu(1:n_in-2)
                c_three(1:n_in-2)=lamda(1:n_in-2)
                a_three(n_in-1)=h(n_in-1)/(h(n_in-1)+h(1))
                c_three(n_in-1)=1-a_three(n_in-1)
                d(n_in)=(6/(h(n_in-1)+h(1)))*((y(2)-y(n_in))/h(1)-(y(n_in)-y(n_in-1))/h(n_in-1))
                ALLOCATE(d_three(n_in-1))
                ALLOCATE(M_three(n_in-1))
                d_three(1:n_in-1)=d(2:n_in)
                CALL boundary3(n_in-1,a_three,b_three,c_three,&
                                d_three,M_three)
                M(2:n_in)=M_three(1:n_in-1)
                M(1)=M_three(n_in-1)
                DEALLOCATE(a_three)
                DEALLOCATE(b_three)
                DEALLOCATE(c_three)
                DEALLOCATE(d_three)
                DEALLOCATE(M_three)
            END SELECT
            DEALLOCATE(miu)
            DEALLOCATE(lamda)
            DEALLOCATE(d)
    END SUBROUTINE
    SUBROUTINE boundary3(n,a,b,c,d,x)
        INTEGER(4),INTENT(IN)::n
        REAL(8),DIMENSION(:),INTENT(IN)::a(n),b(n),c(n),d(n)
        REAL(8),DIMENSION(:),INTENT(INOUT)::x(n)
        REAL(8),ALLOCATABLE::u(:),l(:),y(:)
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
        l(1)=a(1)
        u(n)=c(n)/u(1)
        DO i=2,n
            l(i)=a(i)/u(i-1)
            u(i)=b(i)-l(i)*c(i-1)
            y(i)=d(i)-l(i)*y(i-1)
        END DO
        y(n)=y(n)-u(n)*l(1)
        x(n)=y(n)/u(n)
        DO i=n-1,2,-1
            x(i)=(y(i)-c(i)*x(i+1))/u(i)
        END DO
        x(1)=(y(1)-c(2)*x(2)-l(1)*x(n))/u(1)
        DEALLOCATE(y)
        DEALLOCATE(u)
        DEALLOCATE(l)
    END SUBROUTINE



    FUNCTION diff3(x1,x2,x3,y1,y2,y3)
        REAL(8),INTENT(IN)::x1,x2,x3
        REAL(8),INTENT(IN)::y1,y2,y3
        REAL(8)::f1,f2,diff3
        f1=(y2-y1)/(x2-x1)
        f2=(y3-y2)/(x3-x2)
        diff3=(f2-f1)/(x3-x1)
    END FUNCTION
END MODULE