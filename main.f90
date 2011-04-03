program fbem

  use discretization
  use bem

  real , parameter :: h = 3
  real, allocatable :: nodes(:,:)
  real, allocatable :: V(:,:)
  integer, allocatable :: edges(:,:)
  integer :: bounds(3,2)
  integer :: node_count
  integer :: v_size, i, j
  real :: a(2), b(2), x(2)

  node_count = nodeCount(2.0,6.0,2.0,2.0,h)

  write (*,*) node_count

  allocate(nodes(node_count, 2), edges(node_count,2))

  call make_borders(2.0,6.0,2.0,2.0,2.0,h, nodes, edges, bounds)

  v_size = bounds(1,2) - bounds(1,1) + 1
  allocate(V(v_size, v_size))
  
  call assembleV(bounds(1,1), bounds(1,2), bounds(1,1), bounds(1,2), nodes, edges, V)

  do i = bounds(1,1), bounds(1,2)
    write (*,*) V(i, :)
  end do
  
  deallocate(nodes, edges, V)
end program fbem
