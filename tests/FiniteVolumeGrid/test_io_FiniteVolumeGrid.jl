function convert_simplexgrid_to_finitevolumegrid_1d()
  sg = Grid.SimplexGrid("./SimplexGrid/data/pn_1d.sg")
  try
    fvg = Grid.FiniteVolumeGrid(sg)
    #= println(fvg.faceMarkers) =#
    println(fvg; full = true)
  catch convert_error
    println(convert_error)
    return false
  end
  #= print(fvg) =#
  return true
end

@test convert_simplexgrid_to_finitevolumegrid_1d() == true
