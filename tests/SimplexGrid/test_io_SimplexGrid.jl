function import_dummy_simplex_grid_file()
  try
    dummy_sg = Grid.SimplexGrid("./data/unavailable_simplexgrid.sg")
  catch import_dummy_sg_error
    println(import_dummy_sg_error)
    return true
  end
  return false
end

function import_1d_simplex_grid_file()
  try
    pn_grid = Grid.SimplexGrid("./SimplexGrid/data/pn_1d.sg")
    print(pn_grid; full = true)
    print(pn_grid)
  catch import_1d_sg_error
    @error "Importing the testing simplex grid file is failed!"
    throw(import_1d_sg_error)
  end
  return true
end

@test import_dummy_simplex_grid_file() == true
@test import_1d_simplex_grid_file() == true
