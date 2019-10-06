program main
use global_variables
use input_output
use initialize_update
use compute_energy
implicit none

!##########Data Dictionary############!
  integer :: i
  real*8  :: EE        ! Total Energy before move
  real*8  :: EE1       ! Total Energy after move
  real*8  :: DeltaE    ! Energy difference
  !   real*8  :: wn        ! Rosenbluth factor of new configuration
  !   real*8  :: wo        ! Rosenbluth factor of old configuration
!#####################################!

!#############Initialize##############!
  call cpu_time(started)
  call random_seed()
  !
  !input and initialize system, timing and histogram parameters.
  call initialize_parameters
  !
  !
  if (restart_or_continue /= 1 ) then
    !
    !initialize position
    call Initialize_position
    call write_pos
    call write_pos1(1)
    !
    !initialize energy and parameters of potential
    call initialize_energy_parameters
    !
    !Compute total energy
    call total_energy(EE)
    i=1
  else
    !
    !read position and histogram data
    call continue_read_data(i)
    !
    !initialize energy and parameters of potential
    call initialize_energy_parameters
    !
    !Compute total energy
    call total_energy(EE)
  end if
!#####################################!


!##############Preheation#############!
  if ( i <= StepNum0 ) then
    do step = i, StepNum0
      call Monte_Carlo_Move(EE, DeltaE)
!       write(*,*) EE,DeltaE
      if ( mod(step,DeltaStep1) == 0 ) then
        call compute_physical_quantities
        call total_energy(EE1)
        call write_physical_quantities( step, EE, EE1, DeltaE )
        call adjust_move_distance
      end if
      call update_verlet_list
      if ( mod(step,DeltaStep2) == 0 ) then
        call write_pos1(step)
      end if
    end do
    i = step
  end if
!#####################################!
  
  call total_energy(EE)
!###############Running###############!
  do step=i, StepNum+StepNum0
    call Monte_Carlo_Move(EE, DeltaE)
    if ( mod(step,DeltaStep1) == 0 ) then 
      call compute_physical_quantities
      call compute_radial_distribution_function
      call total_energy(EE1)
      call write_physical_quantities( step, EE, EE1, DeltaE )
    end if
    call update_verlet_list
    if ( mod(step, DeltaStep2) == 0 ) then
      call write_pos1(step)
    end if
  end do
!#####################################!

!###############Finished##############!
  call cpu_time(finished)
  total_time=finished-started+total_time
  call write_time(total_time)
  write(*,*) 'finished!'
!#####################################!

end program main








