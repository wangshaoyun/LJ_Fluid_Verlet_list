module initialize_update
  !------------------------------------!
  !Use system parameters to initialize
  !the positions and update the positions.
  !------------------------------------!
  implicit none

  contains


subroutine Initialize_position
  !--------------------------------------!
  !Initialize linear chains.
  !
  !Input
  !  pos
  !Output
  !  pos
  !External Variables
  !  
  !Routine Referenced:
  !   rij_and_rr, period_condition_rij
  !Reference:
  !   Spherical coordinate on Wiki to generate uniform distribution
  !   on the sphere surface.
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, m
  real*8 :: rnd1, rnd2, rnd3, rsqr
  real*8, dimension(3) :: rij

  do i=1, NN
    m=1
    do while (m==1)
      m=0
      call random_number(rnd1)
      call random_number(rnd2)
      call random_number(rnd3)
      pos(i,1)=rnd1 * Lx - Lx/2
      pos(i,2)=rnd2 * Ly - Ly/2
      pos(i,3)=rnd3 * Lz - Lz/2
      !
      !periodic condition
      call periodic_condition(pos(i,1:3))
      !
      !Judge whether the particle is close the former paritcle
      !too much.
      do j=1,i-1
        call rij_and_rr(rij,rsqr,i,j)
        if (rsqr<0.81) then
          m=1
          cycle
        end if
      end do
    end do
  end do

end subroutine Initialize_position


subroutine Monte_Carlo_Move( EE, DeltaE )
  !------------------------------------!
  !
  !   
  !Input
  !   
  !Output
  !   
  !External Variables
  !   
  !Routine Referenced:
  !1.
  !------------------------------------!
  use global_variables
  use compute_energy
  implicit none
  real*8, intent(inout) :: EE
  real*8, intent(out)   :: DeltaE
  integer :: j
  real*8 :: EE1, EE2

  do j = 1, NN
!       call total_energy(EE1)

    call Choose_Particle
    call New_Position
    call Delta_Energy(DeltaE)
    call Move_or_not(EE, DeltaE)
    !
    !test EE2-EE1 = DeltaE
!       call total_energy(EE2)
!       write(*,*) EE2 - EE1, DeltaE, EE2, EE1  
  end do

end subroutine Monte_Carlo_Move


subroutine choose_particle
  !------------------------------------!
  !This subroutine is used to choose a particle ip to move.
  !   
  !Input
  !   
  !Output
  !   ip
  !External Variables
  !   NN, Nm, Npe, ip
  !Routine Referenced:
  !1.
  !------------------------------------!
  use global_variables
  implicit none
  real*8 :: rnd                 

  call random_number(rnd)
  ip = int(rnd*NN) + 1

end subroutine choose_particle


subroutine New_Position
  !--------------------------------------!
  !This program is used to generate new position.
  !   
  !Input
  !   ip
  !Output
  !   pos1
  !External Variables
  !   pos, pos1, dr
  !Routine Referenced:
  !1. Periodic_condition( rr(2) )
  !--------------------------------------!
  use global_variables
  implicit none
  real*8 :: rnd(3)

  pos_ip0 = pos(ip,:)
  call random_number(rnd)
  pos_ip1(1:3) = pos_ip0(1:3) + (2*rnd - 1) * dr
  call periodic_condition( pos_ip1 )

end subroutine New_Position


subroutine Move_or_not(EE, DeltaE)
  !--------------------------------------!
  !
  !   
  !Input
  !   
  !Output
  !   
  !External Variables
  !   pos, pos_ip0, pos_ip1, ip, Beta
  !Routine Referenced:
  !1.
  !Reference:
  !
  !--------------------------------------!
  use global_variables
  use compute_energy
  implicit none
  real*8,  intent(in)   :: DeltaE
  real*8,  intent(inout) :: EE
  real*8  :: rnd
  !
  !Judge whether move or not
  if ( DeltaE < 0 ) then
    pos(ip,1:3) = pos_ip1(1:3)
    EE = EE + DeltaE
    accpt_num = accpt_num + 1
  else 
    call random_number(rnd)
    if ( rnd < Exp(-Beta*DeltaE) ) then
      pos(ip,1:3) = pos_ip1(1:3)
      EE = EE + DeltaE
      accpt_num = accpt_num + 1
    endif
  endif
  total_num = total_num + 1
end subroutine Move_or_not


subroutine adjust_move_distance
  !--------------------------------------!
  !
  !   
  !Input
  !   
  !Output
  !   
  !External Variables
  !   pos, pos_ip0, pos_ip1, ip, Beta
  !Routine Referenced:
  !1.
  !Reference:
  !
  !--------------------------------------!
  use global_variables
  implicit none
  real*8 :: delta_accpt_ratio

  accpt_ratio = accpt_num / total_num
  delta_accpt_ratio = accpt_ratio - best_accpt_ratio 
  dr = dr + delta_dr * delta_accpt_ratio / abs(delta_accpt_ratio)

  accpt_num = 0
  total_num = 0

end subroutine adjust_move_distance


end module initialize_update



