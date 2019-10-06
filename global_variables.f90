module global_variables

  implicit none
  save

!########################constants#########################!
  real*8, parameter :: pi=3.141592653589793D0     !Circumference ratio pi
  real*8, parameter :: gamma=.5772156649015329D0  !Euler constants gamma
!########################constants#########################!

!####################systems coefficient###################!
  integer :: NN       !Total particles in the system
  real*8  :: rho      !Polymer monomer number density 
  real*8  :: Lx       !Length of cell in x direction
  real*8  :: Ly       !Length of cell in y direction
  real*8  :: Lz       !Length of cell in z direction
  real*8  :: Beta     !Beta=1/(kB*T), T is temperature, 
                      !kB is Boltzmann constant
!##################end systems coefficient#################!


!##################running and Histogram###################!
  integer :: restart_or_continue  !restart or continue after breaking off 
  integer :: StepNum0             !steps of preheating
  integer :: StepNum              !steps of running
  integer :: DeltaStep1           !step inteval, physical quantities
  integer :: DeltaStep2           !step inteval, write data
  integer :: step                 !steps of calculate the physical quantities
  real*8  :: dr                   !length of each moving
  real*8  :: total_num = 0        !Total choose number
  real*8  :: accpt_num = 0        !accepted number
  real*8  :: accpt_ratio          !accepted ratio
  real*8  :: best_accpt_ratio     !best accepted ratio
  real*8  :: delta_dr             !adjust move distance
  !
  !timing
  real*8  :: started    = 0       !time at starting
  real*8  :: finished   = 0       !time at finishing
  real*8  :: total_time = 0       !total time of the simulation
  !
  !histogram
  integer :: SizeHist             !number of histogram which is equally divided
!################end running and Histogram#################!

!##########################arrays##########################!
 real*8, allocatable, dimension(:,:) :: pos     !old position array
 real*8,  dimension(4)   :: pos_ip0
 real*8,  dimension(4)   :: pos_ip1
!real*8, allocatable, dimension(:,:) :: pos_ip0 !old position of part of chains
!real*8, allocatable, dimension(:,:) :: pos_ip1 !new position of part of chains
 integer :: ip                                  !The chain that is choosed
!  integer :: im                                  !The monomer that is choosed
!########################end arrays########################!

  contains 

subroutine periodic_condition(rr)
  !--------------------------------------!
  !3D Peridodic condition of position vector
  !   
  !Input
  !   rr
  !Output
  !   rr
  !External Variables
  !   Lx, Ly, Lz
  !Routine Referenced:
  !1.
  !--------------------------------------!
  implicit none
  real*8, intent(inout) :: rr(3)

  if ( rr(1) > Lx/2 ) then
    rr(1) = rr(1) - Lx
  elseif( rr(1) <= -Lx/2 ) then
    rr(1) = rr(1) + Lx
  end if
  if ( rr(2) > Ly/2 ) then
    rr(2) = rr(2) - Ly
  elseif( rr(2) <= -Ly/2 ) then
    rr(2) = rr(2) + Ly
  end if
  if ( rr(3) > Lz/2 ) then
    rr(3) = rr(3) - Lz
  elseif( rr(3) <= -Lz/2 ) then
    rr(3) = rr(3) + Lz
  end if

end subroutine periodic_condition


subroutine rij_and_rr(rij, rsqr, i, j)
  !-----------------------------------------!
  !compute displacement vector and displacement of two particles
  !input:
  !  i, j(particle number) 
  !output:
  !  rij(displacement vecter), rr(square of displacement)
  !External Variant:
  !  Lx,Ly,Lz(used in period condition)
  !  pos
  !note:
  !  including period condition
  !-----------------------------------------!
  implicit none
  real*8, dimension(3), intent(out) :: rij
  real*8, intent(out) :: rsqr
  integer, intent(in) :: i
  integer, intent(in) :: j

  rij = pos(i,1:3) - pos(j,1:3)

  ! Periodic Condition
  if ( rij(1) > Lx/2 ) then
    rij(1) = rij(1) - Lx
  elseif( rij(1) <= -Lx/2 ) then
    rij(1) = rij(1) + Lx
  end if
  if ( rij(2) > Ly/2 ) then
    rij(2) = rij(2) - Ly
  elseif( rij(2) <= -Ly/2 ) then
    rij(2) = rij(2) + Ly
  end if
  if ( rij(3) > Lz/2 ) then
    rij(3) = rij(3) - Lz
  elseif( rij(3) <= -Lz/2 ) then
    rij(3) = rij(3) + Lz
  end if

  rsqr = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)

end subroutine rij_and_rr

end module global_variables

