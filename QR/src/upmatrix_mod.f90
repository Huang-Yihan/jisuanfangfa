!this mod is use to find  x such as Ux=b(LU) or Rx=b(QR)
!you should make sure that R b x have the same row and R is a upmatrix(n)
MODULE upmatrix_mod
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: upmatrix
CONTAINS 
    SUBROUTINE upmatrix(R,b)
        INTEGER(4)::m,n,i,ii
        REAL(4)::sum
        REAL(4),DIMENSION(:,:),INTENT(IN)::R
        REAL(4),DIMENSION(:,:),INTENT(INOUT)::b
        m=size(R,1)
        n=size(R,2)
        DO i=n,1,-1
            IF (i==n) b(i,1)=b(i,1)/R(i,i)
            IF (i/=n) THEN
                sum=0
                DO ii=i+1,n
                    sum=sum+R(i,ii)*b(ii,1)
                END DO
                b(i,1)=(b(i,1)-sum)/R(i,i)
            END IF
        END DO
     END SUBROUTINE upmatrix
END MODULE
