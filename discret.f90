
!       xxxxxx               xxxxxx                       
!       x    x               x    x
!       x    x     xxxxx     x    x
!    iy x O- x  jy x R x     x O+ x                     
!       x    x     xxxxx     x    x                  
!       x    x      jx       x    x         
!       xxxxxx  d            xxxxxx          
!         ix

!------------------------------------------------------------------------------
subroutine make_borders(ix, iy, jx, jy, d, h, nodes, edges) 
!------------------------------------------------------------------------------

  implicit none

  double precision, intent(in) :: ix, iy, jx, jy, d
  double precision, intent(in) :: h !discretization step
  double precision, allocatable, intent(out) :: nodes(:,:), edges(:,:)

  integer :: node_count
  
  node_count = 4 * ceiling(ix / h) + 4 * ceiling(iy / h) &
             + 2 * ceiling(jx / h) + 2 * ceiling(jy / h)
 
  write (*,*) 'Discretization size is ', node_count 
end subroutine make_borders
