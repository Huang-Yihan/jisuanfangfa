PROGRAM MAIN
    USE wrtmat_mod
    USE cgm_mod
    IMPLICIT NONE
    INTEGER(8)::m,n,i
    CHARACTER(80):: words
    REAL(8),ALLOCATABLE::A(:,:)
    REAL(8),ALLOCATABLE::b(:,:),x(:,:)
    !----------------------------------------
    ! read the input file and make the matrix
    OPEN(99,FILE='input')
    m=0
    n=0
    DO 
        READ(99,'(A80)',END=100) words
        IF(words=='m =') READ(99,*) m
        IF(words=='n =') READ(99,*) n
        IF(m<20) THEN
            IF(words=='A =') then
                ALLOCATE(A(m,n))
                DO i=1,m
                    READ(99,*) A(i,1:n)
                END DO
            END IF
            IF(words=='bt=') then
                ALLOCATE(b(m,1))
                READ(99,*) b(1:m,1)
            END IF
        ELSE
    !------------------------------------------
    !define A and b according to the  P113-3.2        
            ALLOCATE(A(m,m))
            A=0
            ALLOCATE(b(m,1))
            b=0

            A(1,1)  =-2
            A(1,2)  = 1
            A(m,m)  =-2
            A(m,m-1)= 1
            b(1,1)  =-1
            b(m,1)  =-1

            DO i=2,m-1
                A(i,i)=-2
                A(i,i+1)=1
                A(i,i-1)=1
            END DO
            EXIT
        END IF       
    END DO
100 CLOSE(99)
    !----------------------------------------
    !----------------------------------------
    !do CGM,and write A,b,x in the output file
    OPEN(100,file='output')
    WRITE(*,*)'A='
    WRITE(100,*)'A='
    CALL wrtmat(A)
    WRITE(*,*)'b='
    WRITE(100,*)'b='
    CALL wrtmat(b)
    ALLOCATE(x(m,1))
    x=0
    CALL CGM(A,b,x)

    WRITE(*,*)'x='
    WRITE(100,*)'x='
    CALL wrtmat(x)
    CLOSE(100)
    !--------------------------------------- 
    DEALLOCATE(b)
    DEALLOCATE(A)
    DEALLOCATE(x)
END PROGRAM