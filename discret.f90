!******************************************************************************
          module discretization
!******************************************************************************

contains
!******************************************************************************

!------------------------------------------------------------------------------
integer function nodeCount(ix, iy, jx, jy, h) 
!------------------------------------------------------------------------------

  implicit none

  real, intent(in) :: ix, iy, jx, jy, h
  integer :: nodeCount

 nodeCount = 4 * ceiling(ix / h) + 4 * ceiling(iy / h) &
             + 2 * ceiling(jx / h) + 2 * ceiling(jy / h)

end function nodeCount

!------------------------------------------------------------------------------
!       +----+               +----+                       
!       |    |               |    |
!       |    |     +---+     |    |
!    iy | O- |  jy | R |     | O+ |                     
!       |    |     +---+     |    |                  
!       |    |<---> jx       |    |         
!       +----+  d            +----+          
!         ix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine make_borders(ix, iy, jx, jy, d, h, nodes, edges, area_index_bounds) 
!------------------------------------------------------------------------------

  implicit none

  real, intent(in) :: ix, iy, jx, jy, d
  real, intent(in) :: h !discretization step
  real, intent(out) :: nodes(:,:)
  integer, intent(out) :: edges(:,:)
  integer, intent(out) :: area_index_bounds(3,2)

  integer :: counter, i

  counter = 1

  call discrete_rect( (/ -ix - d - jx / 2.0, -iy / 2.0 /), (/ - d - jx / 2.0, iy / 2.0 /), &
                      h, nodes, edges, counter)
  area_index_bounds(1,:) = (/ 1, counter /)

  counter = counter + 1
  area_index_bounds(2,1) = counter
  
  call discrete_rect( (/ -jx / 2.0, -jy / 2.0 /), (/ jx / 2.0, jy / 2.0 /), &
                      h, nodes, edges, counter)

  area_index_bounds(2,2) = counter
  counter = counter + 1
  area_index_bounds(3,1) = counter
  call discrete_rect( (/ + d + jx / 2.0, -iy / 2.0 /), (/ jx / 2.0 + d + ix, iy / 2.0 /),  &
                      h, nodes, edges, counter)
  area_index_bounds(3,2) = counter

end subroutine make_borders

!------------------------------------------------------------------------------
subroutine discrete_rect(bl, tr, h, nodes, edges, counter)
!------------------------------------------------------------------------------
  implicit none

  real, intent(in), dimension(2) :: bl, tr
  real, intent(out) :: nodes(:,:)
  real, intent(in) :: h
  integer, intent(out) :: edges(:,:)
  integer, intent(inout) :: counter
  integer :: start_index, i
  integer :: no_parts_x, no_parts_y

  start_index = counter 
  
  nodes(counter, :) = bl(:) 

  no_parts_x = ceiling((tr(1) - bl(1)) / h)
  no_parts_y = ceiling((tr(2) - bl(2)) / h)
  
  i = 0
  do
    nodes(counter, :) = (/ bl(1) + float(i) / float(no_parts_x) * (tr(1) - bl(1)), bl(2)/)
    edges(counter, :) = (/ counter, counter + 1/)
    if (i == no_parts_x) then
      exit
    end if
    counter = counter + 1
    i = i + 1
  end do

  i = 0
  do
    nodes(counter, :) = (/ tr(1), bl(2) + float(i) / float(no_parts_y) * (tr(2) - bl(2)) /)
    edges(counter, :) = (/ counter, counter + 1/)
    if (i == no_parts_y) then
      exit
    end if
    counter = counter + 1
    i = i + 1
  end do

  i = 0
  do
    nodes(counter, :) = (/ tr(1) - float(i) / float(no_parts_x) * (tr(1) - bl(1)), tr(2)/)
    edges(counter, :) = (/ counter, counter + 1/)
    if (i == no_parts_x) then
      exit
    end if
    counter = counter + 1
    i = i + 1
  end do

  i = 0
  do
    nodes(counter, :) = (/ bl(1), tr(2) - float(i) / float(no_parts_y) * (tr(2) - bl(2)) /)
    edges(counter, :) = (/ counter, counter + 1/)
    if (i == no_parts_y - 1) then
      exit
    end if
    counter = counter + 1
    i = i + 1
  end do

  edges(counter, :) =  (/ counter, start_index /)

end subroutine discrete_rect

end module discretization
!******************************************************************************
