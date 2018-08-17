function import_dummy_simplex_grid_file()
  try
    dummy_sg = Grid.SimplexGrid("./data/unavailable_simplexgrid.sg")
  catch import_dummy_sg_error
    println(import_dummy_sg_error)
    throw(import_dummy_sg_error)
  end
  return false
end

function import_1d_simplex_grid_file()
  #= try =#
  #=   pn_grid = Grid.SimplexGrid("./SimplexGrid/data/pn_1d.sg") =#
  #= catch import_1d_sg_error =#
  #=   println(import_1d_sg_error) =#
  #=   throw(import_1d_sg_error) =#
  #= end =#
  pn_grid = Grid.SimplexGrid("./SimplexGrid/data/pn_1d.sg")
  return true
end

#= @test_throws ErrorException import_dummy_simplex_grid_file() =#
import_1d_simplex_grid_file()
