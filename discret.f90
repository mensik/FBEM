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
subroutine make_borders(ix, iy, jx, jy, d, h, nodes, edges) 
!------------------------------------------------------------------------------

  implicit none

  real, intent(in) :: ix, iy, jx, jy, d
  real, intent(in) :: h !discretization step
  real, intent(out) :: nodes(:,:)
  integer, intent(out) :: edges(:,:)

  integer :: counter, i

  counter = 1

  call discrete_rect( (/ -ix - d - jx / 2.0, -iy / 2.0 /), (/ - d - jx / 2.0, iy / 2.0 /), &
                      h, nodes, edges, counter)

  counter = counter + 1

  call discrete_rect( (/ -jx / 2.0, -jy / 2.0 /), (/ jx / 2.0, jy / 2.0 /), &
                      h, nodes, edges, counter)

  counter = counter + 1
 
  call discrete_rect( (/ + d + jx / 2.0, -iy / 2.0 /), (/ jx / 2.0 + d + ix, iy / 2.0 /),  &
                      h, nodes, edges, counter)

  do i = 1, counter
    write (*,*) i, nodes(i,:)
  end do

  do i = 1, counter
    write (*,*) i, edges(i,:)
  end do

end subroutine make_borders

!------------------------------------------------------------------------------
subroutine discrete_rect(bl, tr, h, nodes, edges, counter)
!------------------------------------------------------------------------------
  implicit none

  real, intent(in), dimension(1,2) :: bl, tr
  real, intent(out) :: nodes(:,:)
  real, intent(in) :: h
  integer, intent(out) :: edges(:,:)
  integer, intent(inout) :: counter
  integer :: start_index

  start_index = counter 
  
  nodes(counter, :) = bl(1,:) 

  do
    counter = counter + 1
   
    nodes(counter, :) = (/ nodes(counter - 1, 1) + h, nodes(counter - 1, 2) /)
    edges(counter, :) = (/ counter - 1, counter /)

    if (nodes(counter, 1) >= tr(1,1)) then
      nodes(counter, 1) = tr(1,1)
      exit
    end if
  end do

  do
    counter = counter + 1
   
    nodes(counter, :) = (/ nodes(counter - 1, 1), nodes(counter - 1, 2) + h /)
    edges(counter, :) = (/ counter - 1, counter /)
    
    if (nodes(counter, 2) >= tr(1,2)) then
      nodes(counter, 2) = tr(1,2)
      exit
    end if
  end do

  do
    counter = counter + 1
   
    nodes(counter, :) = (/ nodes(counter - 1, 1) - h, nodes(counter - 1, 2) /)
    edges(counter, :) = (/ counter - 1, counter /)
    
    if (nodes(counter, 1) <= bl(1,1)) then
      nodes(counter, 1) = bl(1,1)
      exit
    end if
  end do

  do
    if (nodes(counter, 2) - h <= bl(1,2)) then
      exit
    end if
    
    counter = counter + 1
   
    nodes(counter, :) = (/ nodes(counter - 1, 1), nodes(counter - 1, 2) - h/)
    edges(counter, :) = (/ counter - 1, counter /)
    
    if (nodes(counter, 2) <= bl(1,2)) then
      nodes(counter, 2) = tr(1,2)
      exit
    end if
  end do

  edges(counter, :) = (/ counter, start_index /)

end subroutine discrete_rect

end module discretization
!******************************************************************************
