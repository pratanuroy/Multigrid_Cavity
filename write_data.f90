subroutine write_data(ilevel)
    use global_variables
    use parallel_variables
    use parameters
    use multigrid_levels
    use VTR
use kind_parameters
    implicit none
include "mpif.h"

! Note: Currently the VTK module (binary format) for 2D case is not working.
! Need to fix the bug

    type(VTR_file_handle):: fdin
    integer*4, intent(in) :: ilevel
    integer*4:: restart, i, imax, jmax, xi, xf
    logical:: file_exists
    character*20:: snapshot_tcounter, snapshot_proc, snapshot_numtasks, snapshot_time

    tcounter = 0

    imax = mg_level(ilevel)%nx
    jmax = mg_level(ilevel)%ny

    xi = mg_level(ilevel)%xi
    xf = mg_level(ilevel)%xf
    
    u  => mg_level(ilevel)%u
    v  => mg_level(ilevel)%v

    xp  => mg_level(ilevel)%xp
    yp  => mg_level(ilevel)%yp
   
    if (parallel_debug .eq. 2) then
    write(*,*) imax, jmax
    write(*,*) xi, xf 
    
    if (taskid .eq. root) write(*,*) xp
    if (taskid .eq. root) write(*,*) yp

    if (taskid .eq. root) write(*,*) u(imax/2,:)
    if (taskid .eq. root) write(*,*) v(imax/2,:)

    write(*,*) shape(u), shape(v)

    end if

    restart = 0

    if (write_VTK) then

        if (parallel_debug)  write(*,*) 'Started writing VTR file from processor = ', taskid
        
        call VTR_open_file(PREFIX=prog_name, proc_rank=taskid, num_procs=numtasks, restart=0, FD = fdin)
        !call VTK_write_mesh(FD =fdin, X = xp(1:imax), Y = yp(xi-1:xf+1), Z = xp(1:imax))
        call VTR_write_mesh(FD =fdin, X = xp(1:imax), Y = yp(xi-1:xf+1))
        !call VTK_write_var(FD = fdin, NAME="Velocity", VX = u(1:imax,xi-1:xf+1), VY = v(1:imax,xi-1:xf+1), VZ = v(1:imax,1))
        call VTR_write_var(FD = fdin, NAME="Velocity", VX = u(1:imax,xi-1:xf+1), VY = v(1:imax,xi-1:xf+1))

        call VTR_close_file(FD = fdin)

        if (parallel_debug) write(*,*) 'Finished writing VTR file from processor = ', taskid

        if (taskid .eq. root) then
            ! Write .visit file
            write(snapshot_tcounter,'(i4)') tcounter/output_step+1
            write(snapshot_numtasks,'(i8)') numtasks
            write(snapshot_time,'(e13.5)') time
            filename = trim(prog_name)//'.visit'
            inquire(file=filename, exist=file_exists)
            if(file_exists) then

                if (tcounter .eq. 0) then

                    open(unit=9, file=filename, status='replace', action='write')
                    write(9,'(a)')adjustl('!NBLOCKS ')//trim(adjustl(snapshot_numtasks))
                    !write(9,'(a)')adjustl('!TIME ')//trim(adjustl(snapshot_time))
                    do i = 0, numtasks-1
                        write(snapshot_proc,'(i8)') i
                        write(9,'(a,a)')trim(adjustl(prog_name))//'_'//trim(adjustl(snapshot_proc))//'_',trim(adjustl(snapshot_tcounter))//'.vtr'
                    end do
                    close(9)

                else

                    open(unit=9, file=filename, status='old',position='append', action='write')
                    !write(9,'(a)')adjustl('!TIME ')//trim(adjustl(snapshot_time))
                    do i = 0, numtasks-1
                        write(snapshot_proc,'(i8)') i
                        write(9,'(a,a)')trim(adjustl(prog_name))//'_'//trim(adjustl(snapshot_proc))//'_',trim(adjustl(snapshot_tcounter))//'.vtr'
                    end do
                    close(9)

                end if

            else

                open(unit=9, file=filename, status='replace', action='write')
                write(9,'(a)')adjustl('!NBLOCKS ')//trim(adjustl(snapshot_numtasks))
                !write(9,'(a)')adjustl('!TIME ')//trim(adjustl(snapshot_time))
                do i = 0, numtasks-1
                    write(snapshot_proc,'(i8)') i
                    write(9,'(a,a)')trim(adjustl(prog_name))//'_'//trim(adjustl(snapshot_proc))//'_',trim(adjustl(snapshot_tcounter))//'.vtr'
                end do
                close(9)

            end if

        end if ! taskid .eq. root or not

    end if

end subroutine
