subroutine postprocess
    use global_variables
    use parallel_variables
use kind_parameters
    implicit none
include "mpif.h"

    !real(kind=dp):: interface_location_x, interface_location_y
    !real(kind=dp):: minloc_from_centroid_x, minloc_from_centroid_y
    integer:: fileid, i, j

    fileid = 10

    fileid = fileid+1

    if (taskid .eq. root) then

            filename = 'uvelocity_centerline.curve'

            open (unit = fileid, file = filename)

            write(fileid,*) 'yposition', '   ', 'uvelocity'

           do i = 1,imax
            write(fileid,*) yp(i), u(i)
           end do
    end if

end subroutine
