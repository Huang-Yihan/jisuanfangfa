PROGRAM MAIN
    USE tscr_mod
    USE function_mod
    IMPLICIT NONE
    REAL(8) :: a,b
    INTEGER::k
    CLASS(jifen_type),ALLOCATABLE:: jifen(:)
    INTEGER(4),PARAMETER :: number=50
    REAL(8),PARAMETER :: EPS = 1.0E-8
    INTEGER(4)           :: i,number_fun
    ALLOCATE(jifen(number))
    DO i=1,number
        jifen(i)%t=0
        jifen(i)%s=0
        jifen(i)%c=0
        jifen(i)%r=0
    END DO
    WRITE(*,*)'choose a and b ([a,b]):'
    READ(*,*)a,b
    WRITE(*,*)'choose the function number of the problem'
    WRITE(*,*)'if you first caculate this function,please modify the file:function_mod.f90'
    WRITE(*,*)'then you can define a function number'
    READ(*,*) number_fun
    CALL fun(number_fun,(b-a)/2*(a+b),jifen(1)%t)
    OPEN(81,FILE='output')
    WRITE(*,*)&
    ' k        T             S             C             R'
    WRITE(81,*)&
    ' k        T             S             C             R'
    DO i=1,number
        IF (i>=2) CALL tixing(number_fun,a,b,jifen(i-1),i,jifen(i))
        IF (i>=2) CALL simpson(jifen(i-1),jifen(i))
        IF (i>=3) CALL cotes(jifen(i-1),jifen(i))
        IF (i<=3) WRITE(81,"(I3,2X,F12.9,2x,F12.9,2X,F12.9,2X,F12.9)") i-1,jifen(i)%t,jifen(i)%S,jifen(i)%C,jifen(i)%R
        IF (i<=3) WRITE(*, "(I3,2X,F12.9,2x,F12.9,2X,F12.9,2X,F12.9)") i-1,jifen(i)%t,jifen(i)%S,jifen(i)%C,jifen(i)%R
        IF (i>=4) THEN
            CALL romberg(jifen(i-1),jifen(i))
            WRITE(*, "(I3,2X,F12.9,2x,F12.9,2X,F12.9,2X,F12.9)") i-1,jifen(i)%t,jifen(i)%S,jifen(i)%C,jifen(i)%R
            WRITE(81,"(I3,2X,F12.9,2x,F12.9,2X,F12.9,2X,F12.9)") i-1,jifen(i)%t,jifen(i)%S,jifen(i)%C,jifen(i)%R
            IF (abs(jifen(i)%r-jifen(i-1)%r)<=EPS) EXIT
        END IF
    END DO
    CLOSE(81)
    DEALLOCATE(jifen)
    WRITE(*,*)'This code was written by Huang Yihan'
END PROGRAM
