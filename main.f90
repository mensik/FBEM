program fbem

  use discretization

  real :: ix = 2.0
  real, allocatable :: nodes(:,:)
  integer, allocatable :: edges(:,:)
  integer :: node_count

  node_count = nodeCount(2.0,6.0,2.0,2.0,0.25)

  allocate(nodes(node_count, 2), edges(node_count,2))

  call make_borders(2.0,6.0,2.0,2.0,2.0,0.25, nodes, edges) 

  deallocate(nodes, edges)

end program fbem
