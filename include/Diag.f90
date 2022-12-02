module Diag
	implicit none
	contains

	subroutine Diagonalization(Matrix,Eigenvalues,Eigenvectors,N)
       implicit none
        real*8,dimension(N,N),intent(in) :: Matrix 
        real*8,dimension (N,N),intent(inout) :: Eigenvectors
        real*8,dimension(:,:),allocatable :: Matrix_tmp
        real*8,dimension(N),intent(inout) :: Eigenvalues
        real*8,dimension(:),allocatable :: WORK
        integer,intent(in) :: N
        integer :: INFO, LWORK
        allocate (Matrix_tmp(N,N))
        allocate (work(1))
        Matrix_tmp (:,:)=Matrix (:,:)
        call DSYEV('V','U',N,Matrix_tmp,N,Eigenvalues,WORK,-1,INFO)
        if (INFO .NE. 0) stop
        LWORK=WORK(1)
        deallocate (WORK)
        allocate (WORK(LWORK))
        Matrix_tmp (:,:)=Matrix (:,:) 
        call DSYEV('V','U',N,Matrix_tmp,N,Eigenvalues,WORK,LWORK,INFO)
        If (INFO .LT. 0) print*, "diagonalization failure: wrong argument"
        If (INFO .GT. 0) print*, "diagonalization failure and  convergence not reached"
        if (INFO .NE. 0) stop
        EigenVectors(:,:)=Matrix_tmp(:,:)
        deallocate (work,matrix_tmp)
	endsubroutine
	
end module