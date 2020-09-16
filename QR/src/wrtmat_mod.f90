MODULE wrtmat_mod
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: wrtmat
CONTAINS
    SUBROUTINE wrtmat(A)
        REAL(4),DIMENSION(:,:),INTENT(IN):: A
        REAL(4),ALLOCATABLE::B(:,:)
        INTEGER(4)::i,j,irow,icol
        irow = size(a,1)
        icol = size(a,2)
        ALLOCATE(B(irow,icol))
        B=A

        DO i=1,irow
            DO j=1,icol
                IF(abs(A(i,j))<1.0E-4) B(i,j)=0
                WRITE(*,'(2X,F12.4,$)') B(i,j)
                WRITE(100,'(2X,F12.4,$)') B(i,j)
            END DO
            WRITE(*,*)' '
            WRITE(100,*) ' '
        END DO
        DEALLOCATE(B)
    END SUBROUTINE wrtmat
END MODULE wrtmat_mod