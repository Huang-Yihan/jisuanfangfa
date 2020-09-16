PROGRAM MAIN
    USE csi_mod
    INTEGER(4)::i
    CHARACTER(10)::condition
    REAL(8)::a_condition,b_condition !a and b is the condition of boundary
    REAL(8),ALLOCATABLE::x(:),y(:)
    REAL(8),ALLOCATABLE::M(:),h(:)
    REAL(8),ALLOCATABLE::a(:),b(:),c(:),d(:)
    OPEN(10,file='input')
    READ(10,*) number        !read the total number
    READ(10,'(A10)') condition
    SELECT CASE (condition)
        CASE('first')
            READ(10,*) a_condition,b_condition
        CASE('second')
            READ(10,*) a_condition,b_condition
        CASE('third')
            READ(10,*) ! the third boundary inputfile should have a blank row
            a_condition=0
            b_condition=0
        CASE DEFAULT
            WRITE(*,*) 'ERROR:This problem have no boundary,please choose one boundary'
            STOP
    END SELECT
    ALLOCATE(x(number))
    ALLOCATE(y(number))
    DO i=1,number
        READ(10,*) x(i),y(i) !read the point (x(i),y(i))
    END DO
    CLOSE(10)
    IF (condition=='third'.and.y(1)/=y(number)) THEN
        WRITE(*,*) 'ERROR: the third boundary requires y(0)==y(n),please check out!'
        STOP
    END IF
    ALLOCATE(M(number))
    ALLOCATE(h(number-1))
    CALL CSI(number,x,y,condition,a_condition,b_condition,M,h)
    ALLOCATE(a(number-1))
    ALLOCATE(b(number-1))
    ALLOCATE(c(number-1))
    ALLOCATE(d(number-1))
    OPEN(81,FILE='output')
    DO i=1,number-1
        CALL SIMPLI(x(i),y(i),x(i+1),y(i+1),h(i),&
                    M(i),M(i+1),a(i),b(i),c(i),d(i))
        WRITE(*, "('S(x)=',F9.4'x^3',F9.4,'x^2',F9.4,'x',F9.4,' ,',F5.2'<=x<=',F5.2)") &
                a(i),b(i),c(i),d(i),x(i),x(i+1)
        WRITE(81,"('S(x)=',F9.4'x^3',F9.4,'x^2',F9.4,'x',F9.4,' ,',F5.2'<=x<=',F5.2)") &
                a(i),b(i),c(i),d(i),x(i),x(i+1)
    END DO
    CLOSE(81)
    DEALLOCATE(a)
    DEALLOCATE(b)
    DEALLOCATE(c)
    DEALLOCATE(d)
    DEALLOCATE(h)
    DEALLOCATE(M)
    DEALLOCATE(x)
    DEALLOCATE(y)
CONTAINS
    SUBROUTINE SIMPLI(x1,y1,x2,y2,h,M1,M2,a,b,c,d)
        REAL(8),INTENT(IN)::x1,y1,x2,y2,h,M1,M2
        REAL(8),INTENT(INOUT)::a,b,c,d
        a=0
        b=0
        c=0
        d=0
        a=(M2-M1)/(6*h)
        b=(0-x1*M2+x2*M1)/(2*h)
        c=(x1*x1*M2-x2*x2*M1)/(2*h)-(y1-h*h*M1/6)/h+(y2-h*h*M2/6)/h
        d=(0-x1*x1*x1*M2+x2*x2*x2*M1)/(6*h) & 
           +(y1-h*h*M1/6)*x2/h-(y2-h*h*M2/6)*x1/h
    END SUBROUTINE
END PROGRAM