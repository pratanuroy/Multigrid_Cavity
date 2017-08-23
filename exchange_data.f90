! This subroutine exchange boundary data of each block with the neighbouring block
subroutine exchange_data(dummy_var, chunk_size, start_index, end_index)
use parallel_variables
use kind_parameters
implicit none
include "mpif.h"
integer:: status_array(MPI_STATUS_SIZE,4)
integer:: chunk_size, start_index, end_index
double precision:: dummy_var(chunk_size, start_index-1:end_index+1)

call MPI_IRECV(dummy_var(1,start_index-1), chunk_size, MPI_DOUBLE_PRECISION, iprevid, 1, MPI_COMM_WORLD, request(1), ierr)
call MPI_IRECV(dummy_var(1,end_index+1), chunk_size, MPI_DOUBLE_PRECISION, inextid, 0, MPI_COMM_WORLD, request(2), ierr)
call MPI_ISEND(dummy_var(1,start_index), chunk_size, MPI_DOUBLE_PRECISION, iprevid, 0, MPI_COMM_WORLD, request(3), ierr)
call MPI_ISEND(dummy_var(1,end_index), chunk_size, MPI_DOUBLE_PRECISION, inextid, 1, MPI_COMM_WORLD, request(4), ierr)

call MPI_Waitall(4, request, status_array, ierr)

end subroutine
