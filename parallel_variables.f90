module parallel_variables
use kind_parameters
       integer:: ierr, taskid, numtasks, iproc, root, block_size, dest, source, status, tag, inextid, iprevid, request(4), start_x, end_x
       integer, allocatable:: displs_x(:), rcounts_x(:), displs_ls(:), rcounts_ls(:)
!       real(kind=dp):: local_sum, global_sum
        real(kind=dp):: communication_time, computation_time
end module
