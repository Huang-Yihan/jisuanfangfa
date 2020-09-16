MODULE QR_mod
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: makeQR
CONTAINS
    SUBROUTINE makeQR(m,n,A,Q,R)
        INTEGER(4)::m,n
        REAL(4),DIMENSION(:,:),INTENT(IN):: A
        REAL(4),DIMENSION(:,:),INTENT(OUT)::Q(m,m),R(m,n)
        REAL(4),ALLOCATABLE::H(:,:)
        REAL(4),ALLOCATABLE::Hp(:,:)
        INTEGER(4),ALLOCATABLE::diag(:,:)
        REAL(4),ALLOCATABLE::e(:)
        REAL(4),ALLOCATABLE::x(:)
        REAL(4),ALLOCATABLE::sigma(:)
        INTEGER(4):: i,ii,ij
        REAL(4),ALLOCATABLE::w(:,:)

        ALLOCATE(sigma(min(m-1,n)))
        ALLOCATE(H(m,m))
        R=A
        Q=0
        DO i=1,m
            Q(i,i)=1
        END DO
        DO i=1,min(m-1,n)
            IF (i==1) THEN      !the first step
                ALLOCATE(e(m))
                ALLOCATE(x(m))
                ALLOCATE(w(m,1))
                ALLOCATE(diag(m,m))
                diag=0
                DO ii=1,m
                    diag(ii,ii)=1
                END DO
                e=0
                e(1)=1
                x=A(1:m,1)
                sigma(1)=-sgn(R(1,1))*norm2(x)       
                if(R(1,1)>=0) sigma(1)=-sigma(1) 
                w(:,1)=x-sigma(1)*e
                H=diag-MATMUL(w,transpose(w))/(sigma(1)*(sigma(1)-R(1,1)))
                R=MATMUL(H,R)
            
                Q=MATMUL(Q,transpose(H))
                DEALLOCATE(e)
                DEALLOCATE(x)
                DEALLOCATE(w)
                DEALLOCATE(diag)
            END IF

            IF(i/=1) THEN      ! the Ith step
                ALLOCATE(e(m-i+1))
                ALLOCATE(x(m-i+1))
                ALLOCATE(w(m-i+1,1))
                ALLOCATE(diag(m-i+1,m-i+1))
                ALLOCATE(Hp(m-i+1,m-i+1))
                e=0
                e(1)=1
                x=R(i:m,i)
                sigma(i)=-sgn(R(i,i))*norm2(x)
                if(R(i,i)>=0) sigma(i)=-sigma(i) 
                w(:,1)=x-sigma(i)*e
                diag=0
                DO ii=1,m-i+1
                    diag(ii,ii)=1
                END DO
                Hp=diag-MATMUL(w,transpose(w))/(sigma(i)*(sigma(i)-R(i,i)))
                H=0
                DO ii=1,i-1
                    H(ii,ii)=1
                END DO
                DO ii=i,m
                    DO ij=i,m
                        H(ii,ij)=Hp(ii-i+1,ij-i+1)
                    END DO
                END DO
                R=MATMUL(H,R)
                Q=MATMUL(Q,transpose(H))
                DEALLOCATE(e)
                DEALLOCATE(x)
                DEALLOCATE(w)
                DEALLOCATE(diag)
                DEALLOCATE(Hp)
            END IF
        END DO
    END SUBROUTINE makeQR
    FUNCTION sgn(x)
        REAL(4)::sgn
        REAL(4)::x
        if (x==0) sgn = 0
        if (x>0)  sgn = 1
        if (x<0)  sgn = -1
    END FUNCTION
END MODULE

