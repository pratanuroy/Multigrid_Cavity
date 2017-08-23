subroutine map_task(ilevel)
USE Multigrid_Levels
USE Parallel_Variables
use kind_parameters
implicit none
include "mpif.h"

! Local variables
integer :: imax, jmax
integer, intent(in) :: ilevel
real(kind=dp):: res

!write(*,*) 'num_levels =', num_level

 imax = mg_level(ilevel)%nx
 jmax = mg_level(ilevel)%ny

block_size = (jmax-2)/numtasks
res = mod((jmax-2), numtasks)

if (taskid < res) then

 block_size = block_size + 1
 start_x = 2 + taskid*block_size

else

 start_x = 2 + res + taskid*block_size

end if

end_x = min(start_x + block_size - 1, jmax-1)

mg_level(ilevel)%xi = start_x
mg_level(ilevel)%xf = end_x

end subroutine
