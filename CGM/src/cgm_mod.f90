! this mod is used to find the x in Ax=b
! you should make sure that A is a stmmetrical matrix
MODULE cgm_mod
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: CGM
CONTAINS
    SUBROUTINE CGM(A,b,x)
        INTEGER(8)::n,i,ii
        REAL(8),PARAMETER ::delta = 1.0E-5
        REAL(8)::alpha,oldnorm,beta
        REAL(8),DIMENSION(:,:),INTENT(IN):: A
        REAL(8),DIMENSION(:,:),INTENT(IN):: b
        REAL(8),DIMENSION(:,:),INTENT(INOUT)::x
        REAL(8),ALLOCATABLE::r(:,:)
        REAL(8),ALLOCATABLE::d(:,:)
        REAL(8),ALLOCATABLE::r1(:)
        REAL(8)::rtr(1,1),dtad(1,1)
        

        n=size(A,1)
        ALLOCATE(r(n,1))
        ALLOCATE(d(n,1))
        ALLOCATE(r1(n))

        !step 0:
        r=b-MATMUL(A,x)
        d=r
        r1(1:n)=r(1:n,1)
        DO i=0,n-1
            rtr=MATMUL(transpose(r),r)
            dtad=MATMUL(transpose(d),MATMUL(A,d))
            alpha=rtr(1,1)/dtad(1,1)
            x=x+alpha*d
            oldnorm=NORM2(r1)
            r=b-MATMUL(A,x)
            r1(1:n)=r(1:n,1)
            IF(NORM2(r1)<=delta .OR. i+1==n) EXIT
            beta=NORM2(r1)*NORM2(r1)/(oldnorm*oldnorm)
            d=r+beta*d
        END DO
        DEALLOCATE(r)
        DEALLOCATE(d)
        DEALLOCATE(r1)
    END SUBROUTINE CGM
END MODULE