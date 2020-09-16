PROGRAM MAIN
    USE QR_mod
    USE wrtmat_mod
    USE upmatrix_mod
    IMPLICIT NONE
    INTEGER(4)::m,n,i
    CHARACTER(20):: words
    REAL(4),ALLOCATABLE::A(:,:),Q(:,:),R(:,:)
    REAL(4),ALLOCATABLE::b(:,:)
    !----------------------------------------
    ! read the input file and make the matrix
    OPEN(99,FILE='input')
    DO 
        READ(99,'(A80)',END=100) words
        IF(words=='m =') READ(99,*) m
        IF(words=='n =') READ(99,*) n
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
    END DO
100 CLOSE(99)
    !----------------------------------------
    !----------------------------------------
    !do A=QR,then b=Qt*b,next x=b*R-1, print A,b,Q,R,x in the output file
    OPEN(100,file='output')
    WRITE(*,*)'A='
    WRITE(100,*)'A='
    CALL wrtmat(A)
    WRITE(*,*)'b='
    WRITE(100,*)'b='
    CALL wrtmat(b)
    ALLOCATE(Q(m,m))
    ALLOCATE(R(m,n))
    CALL makeQR(m,n,A,Q,R) !make Q,R,st A=QR
    WRITE(*,*)'R='
    WRITE(100,*)'R='
    CALL wrtmat(R)
    WRITE(*,*)'Q='
    WRITE(100,*)'Q='
    CALL wrtmat(Q)
    b=MATMUL(transpose(Q),b) !caculate b=Qt*b
    CALL upmatrix(R,b)       !CALCULATE THE Rx=b ,then b=x
    WRITE(*,*)'x='
    WRITE(100,*)'x='
    CALL wrtmat(b)
    CLOSE(100)
    !--------------------------------------- 
    DEALLOCATE(b)
    DEALLOCATE(A)
    DEALLOCATE(Q)
    DEALLOCATE(R)
END PROGRAM