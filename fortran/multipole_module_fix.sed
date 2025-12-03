# Fix pointer issues
s/sys%particles => temp_particles/call move_alloc(temp_particles, sys%particles)/g
s/cell%particle_ids => temp_ids/call move_alloc(temp_ids, cell%particle_ids)/g
s/type(cell_type), pointer :: cell$/type(cell_type) :: cell/g
s/cell => sys%cells/cell = sys%cells/g
s/type(cell_type), pointer :: cell1, cell2$/type(cell_type) :: cell1, cell2/g
s/cell1 => sys%cells/cell1 = sys%cells/g
s/cell2 => sys%cells/cell2 = sys%cells/g
