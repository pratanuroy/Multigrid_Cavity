! This subroutine exchange boundary data of each block with the neighbouring block
subroutine exchange_lsdata(dummy_var, chunk_size, start_index, end_index)
use parallel_variables
use kind_parameters
implicit none
include "mpif.h"
integer:: status_array(MPI_STATUS_SIZE,4)
integer:: chunk_size, start_index, end_index, chunk_size1
double precision:: dummy_var(-1:(chunk_size-2), start_index-3:end_index+3)
!write(*,*) 'chunk size', chunk_size, -1+chunk_size
!write(*,*) 'dummy var size = ', size(dummy_var)
!write(*,*) 'start_x =', start_x
chunk_size1 = 3*chunk_size
call MPI_IRECV(dummy_var(-1,start_index-3), chunk_size1, MPI_DOUBLE_PRECISION, iprevid, 1, MPI_COMM_WORLD, request(1), ierr)
call MPI_IRECV(dummy_var(-1,end_index+1), chunk_size1, MPI_DOUBLE_PRECISION, inextid, 0, MPI_COMM_WORLD, request(2), ierr)
call MPI_ISEND(dummy_var(-1,start_index), chunk_size1, MPI_DOUBLE_PRECISION, iprevid, 0, MPI_COMM_WORLD, request(3), ierr)
call MPI_ISEND(dummy_var(-1,end_index-2), chunk_size1, MPI_DOUBLE_PRECISION, inextid, 1, MPI_COMM_WORLD, request(4), ierr)

call MPI_Waitall(4, request, status_array, ierr)

end subroutine
