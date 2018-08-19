""" 
    FiniteVolumeGrid(sg :: SimplexGrid)

A composed type containing information on a finite volume grid and extra information.

# Data fields of `FiniteVolumeGrid`

## Dimensional invariants
* `dimentions :: Int64` -  the dimension ``d`` of grid.
* `numberOfVerticesPerCell :: Int64` -  the number of vertices of a ``d``-simplex
    which is ``d + 1``.
* `numberOfEdgesPerCell :: Int64` -  the number of edges (1-simplex) of a ``d``-simplex
    which is ``\\frac{(d + 1)d}{2}``.
* `numberOfFacetsPerCell :: Int64` -  the number of ``(d - 1)``-simplices of a ``d``-simplex
    which is ``d + 1``.
* `numberOfVerticesPerFacet :: Int64` -  the number of vertices of a ``(d - 1)``-simplex
    which is ``d``.

## Fundamental data fields
* `numberOfNodes :: Int64` -  the number of nodes ``N_{\\text{nodes}}`` of grid.
* `nodes :: Array{Float64, 2}` -  the coordinates of nodes given as a 2D array size of
    ``N_{\\text{nodes}} \\times d`` in which `nodes[i, :]` is ``d``-coordinate of the `i`-th
    node.
* `numberOfSimplices :: Int64` -  the number of ``d``-simplices ``N_{\\text{simplices}}`` of grid.
* `simplices :: Array{Int64, 2}` -  the nodal indices of simplices' vertices given as a 2D array size of
    ``N_{\\text{simplices}} \\times (d + 1)`` in which `simplices[i, :]` is nodal indices of the 
    `i`-th simplex's vertices.
    Simplices are sorted increasingly with respect to their markers (see `simplexMarkers` and 
    `simplexMarkerPointers`).
* `numberOfFaces :: Int64` -  the number of boundary facets ``N_{\\text{faces}}`` of grid.
* `faces :: Array{Int64, 2}` -  the nodal indices of boundary facets' vertices given as a 2D array
    size of ``N_{\\text{faces}} \\times d`` in which `faces[i, :]` is nodal indices of
    the `i`-th face's vertices.
    Faces are sorted increasingly with respect to their markers (see `faceMarkers` and
    `faceMarkerPointers`)
* `numberOfCellMarkers :: Int64` -  the number of markers ``N_{\\text{simplex-markers}}``
    assigned to all simplices. They are going to be refered to materials assigned to simplices.
* `simplexMarkers :: Array{Int64, 1}` - the simplex markers given as a 1D array size of
    ``N_{\\text{simplices}}`` in which `simplexMarkers[i]` is marker assigned to `i`-th simplex.
    These markers must be strictly positive.
* `numberOfFaceMarkers :: Int64` -  the number of markers ``N_{\\text{face-markers}}``
    assigned to all boundary facet. They are going to be refered to boundary conditions.
* `faceMakers :: Array{Int64, 1}` -  the face markers given as a 1D array size of
    ``N_{\\text{faces}}`` in which `faceMarkers[i]` is marker assigned to `i`-th boundary facet.
    These markers must be strictly positive.
* `simplexMarkerPointers :: Array{Int64, 1}` -  the 1D array size of ``N_{\\text{simplex-markers}} + 1`` 
    in which indices of simplices having marker `i` are numbered from `simplexMarkerPointers[i]` to
    `simplexMarkerOffet[i + 1] - 1`.
    All simplices must be ordered such that simplices having the same marker are numbered consecutively. 
* `faceMarkerPointers :: Array{Int64, 1}` -  the 1D array size of ``N_{\\text{face-markers}} + 1``
    in which indices of boundary facets having marker `i` are numbered from `faceMarkerPointers[i]`
    to `faceMarkerPointers[i + 1] - 1`.

## Extra data fields
* `localEdge2Vertices :: Array{Int64, 2}` -  a customized definition of local edges of a simplex.
    Each edge is defined by local indices of two vertices of simplex:
    
    * For 1D simplex, there is only one edge and hence `localEdge2Vertices = [1 2]`.

    * For 2D simplex, there are three edges and hence `localEdge2Vertices = [2 3; 3 1; 1 2]`.

* `nodeFactors :: Array{Float64, 2}`

* `edgeFactors :: Array{Float64, 2}`

* `faceFactors :: Array{Float64, 2}`

"""
mutable struct FiniteVolumeGrid
  dimension                  :: Int64
  numberOfVerticesPerSimplex :: Int64
  numberOfEdgesPerSimplex    :: Int64
  numberOfFacetsPerSimplex   :: Int64
  numberOfVerticesPerFacet   :: Int64

  numberOfNodes     :: Int64
  nodes             :: Array{Float64, 2}
  numberOfSimplices :: Int64
  simplices         :: Array{Int64, 2}
  numberOfFaces     :: Int64
  faces             :: Array{Int64, 2}

  numberOfSimplexMarkers :: Int64
  simplexMarkerPointers  :: Array{Int64, 1}
  numberOfFaceMarkers    :: Int64
  faceMarkerPointers     :: Array{Int64, 1}
  simplexMarkers         :: Array{Int64, 1}
  faceMarkers            :: Array{Int64, 1}


  localEdge2Vertices :: Array{Int64, 2}
  nodeFactors        :: Array{Float64, 2}
  edgeFactors        :: Array{Float64, 2}
  faceFactors        :: Array{Float64, 2}

  function FiniteVolumeGrid()
    empty_fvg = new()
    empty_fvg.dimension         = 0
    empty_fvg.numberOfNodes     = 0
    empty_fvg.numberOfSimplices = 0
    empty_fvg.numberOfFaces     = 0
    return empty_fvg
  end

  function FiniteVolumeGrid(sg :: SimplexGrid)
    fvg = nothing
    try
      fvg = convert_to_finite_volume_grid(sg)
    catch convert_error
      @error "Constructing an instance of FiniteVolumeGrid type from an instance of SimplexGrid is failed."
      rethrow(convert_error)
    end
    return fvg
  end
end

function convert_to_finite_volume_grid(sg :: SimplexGrid)
  fvg = FiniteVolumeGrid()

  if sg.dimension == 0
    return fvg
  end

  fvg.dimension                  = sg.dimension
  fvg.numberOfVerticesPerSimplex = sg.dimension + 1
  fvg.numberOfEdgesPerSimplex    = div(sg.dimension * (sg.dimension + 1), 2)
  fvg.numberOfFacetsPerSimplex   = sg.dimension + 1
  fvg.numberOfVerticesPerFacet   = sg.dimension

  fvg.numberOfNodes     = sg.numberOfNodes
  fvg.nodes             = sg.nodes
  fvg.numberOfSimplices = sg.numberOfSimplices
  fvg.simplices         = sg.simplices
  fvg.numberOfFaces     = sg.numberOfFaces
  fvg.faces             = sg.faces

  if (fvg.dimension == 1)
    fvg.localEdge2Vertices = [1 2]
  elseif (fvg.dimension == 2)
    fvg.localEdge2Vertices = [2 3; 3 1; 1 2]
  end
  fvg.nodeFactors, fvg.edgeFactors = compute_internal_simplex_factors(sg)
  fvg.faceFactors = compute_face_factors(sg)

  fvg.numberOfSimplexMarkers = sg.numberOfSimplexMarkers
  fvg.simplexMarkers = sg.simplexMarkers
  permutationVector = sortperm(fvg.simplexMarkers)
  fvg.simplexMarkers = fvg.simplexMarkers[permutationVector]
  for i = 1:fvg.numberOfVerticesPerSimplex
    fvg.simplices[:, i] = fvg.simplices[permutationVector, i]
    fvg.nodeFactors[:, i] =  fvg.nodeFactors[permutationVector, i]
  end

  fvg.simplexMarkerPointers = zeros(Int64, sg.numberOfSimplexMarkers + 1)
  markerIndex = 0
  for simplexIndex = 1:fvg.numberOfSimplices
    if (markerIndex == (fvg.simplexMarkers[simplexIndex] - 1))
      markerIndex += 1
      fvg.simplexMarkerPointers[markerIndex] = simplexIndex
    end
  end
  fvg.simplexMarkerPointers[markerIndex + 1] = fvg.numberOfSimplices + 1

  fvg.numberOfFaceMarkers = sg.numberOfFaceMarkers
  fvg.faceMarkers = sg.faceMarkers
  permutationVector = sortperm(fvg.faceMarkers)
  fvg.faceMarkers = fvg.faceMarkers[permutationVector]
  for i = 1:fvg.numberOfVerticesPerFacet
    fvg.faces[:, i] = fvg.faces[permutationVector, i]
    fvg.faceFactors[:, i] = fvg.faceFactors[permutationVector, i]
  end
  fvg.faceMarkerPointers = zeros(Int64, fvg.numberOfFaceMarkers + 1)
  markerIndex = 0
  for faceIndex = 1:fvg.numberOfFaces
    if (markerIndex == (fvg.faceMarkers[faceIndex] - 1))
      markerIndex += 1
      fvg.faceMarkerPointers[markerIndex] = faceIndex
    end
  end
  fvg.faceMarkerPointers[markerIndex + 1] = fvg.numberOfFaces + 1

  return fvg
end

function compute_internal_simplex_factors(sg :: SimplexGrid)
  @info "Compute node and edge factors"
  local nodeFactors :: Array{Float64, 2}
  local edgeFactors :: Array{Float64, 2}

  if (sg.dimension == 1)
    nodeFactors, edgeFactors = compute_internal_simplex_factors_1d(sg)
  else
    @error "Unsupported dimension " * string(sg.dimension) * " in computing node and edge factors of simplex!"
    throw(ErrorException("Unsupported dimension " * string(sg.dimension) * " in computing factors of a simplex!"))
  end

  return nodeFactors, edgeFactors
end

function compute_internal_simplex_factors_1d(sg :: SimplexGrid)
  @info "Computing node and edge factors in 1D case."
  numberOfVerticesPerSimplex = 2
  numberOfEdgesPerSimplex    = 1
  nodeFactors             = Array{Float64}(undef, sg.numberOfSimplices, numberOfVerticesPerSimplex)
  edgeFactors             = Array{Float64}(undef, sg.numberOfSimplices, numberOfEdgesPerSimplex)

  for simplexIdx = 1 : sg.numberOfSimplices
    subSimplexArea = 0.5 * abs(sg.nodes[sg.simplices[simplexIdx, 1], 1] - sg.nodes[sg.simplices[simplexIdx, 2], 1])
    nodeFactors[simplexIdx, 1] = subSimplexArea
    nodeFactors[simplexIdx, 2] = subSimplexArea
    edgeFactors[simplexIdx, 1] = 0.5 / subSimplexArea
  end

  return nodeFactors, edgeFactors
end

function compute_face_factors(sg :: SimplexGrid)
  @info "Computing face factors"
  local faceFactors :: Array{Float64, 2}
  if sg.dimension == 1
    faceFactors = compute_face_factors_1d(sg)
  else 
    @error "Unsupported dimension " * string(sg.dimension) * " in computing face factors of simplex!"
    throw(ErrorException("Unsupported dimension " * string(sg.dimension) *
                         " in computing face factors of a simplex!"))
  end
  return faceFactors
end

function compute_face_factors_1d(sg :: SimplexGrid)
  numberOfVerticesPerFacet = sg.dimension
  return ones(Float64, sg.numberOfFaces, numberOfVerticesPerFacet)
end

function print(fvg :: FiniteVolumeGrid; full :: Bool = false)
  fe = Format.FormatExpr("{1:<32s}{2:8d}")
  Format.printfmtln(fe, "Dim. of grid", fvg.dimension)
  Format.printfmtln(fe, "N.o. nodes", fvg.numberOfNodes)
  if full || fvg.numberOfNodes > 0
    if fvg.dimension == 1
      Format.printfmtln("    {1:<16s}{2:>16s}", "# node", "1st coordinate")
      fe1 = FormatExpr("    {1:<16d}{2:16.8e}")
      for ii= 1 : fvg.numberOfNodes
        Format.printfmtln(fe1, ii, fvg.nodes[ii])
      end
    else
      @error "Unsupported dimension " * string(fvg.dimension) * " in printing an instance of FiniteVolumeGrid type!"
      throw(ErrorException("Unsupported dimension " * string(fvg.dimension) *
                           " in printing an instance of FiniteVolumeGrid type!"))
    end
  end
  Format.printfmtln(fe, "N.o. simplices", fvg.numberOfSimplices)
  if full || fvg.numberOfSimplices  > 0
    if fvg.dimension == 1
      Format.printfmtln("    {1:<16s}{2:>16s}{3:>16s}{4:>16s}{5:>16s}{6:>16s}{7:>16s}",
                        "# simplex", "1st vertex", "2nd vertex", "marker",
                        "1st vertex fac.", "2nd vertex fac.", "1st edge fac.")
      fe1 = "    {1:<16d}{2:16d}{3:16d}{4:16d}{5:16.4e}{6:16.4e}{7:16.4e}"
      for ii = 1 : fvg.numberOfSimplices
        Format.printfmtln(fe1, ii, fvg.simplices[ii, 1], fvg.simplices[ii, 2], fvg.simplexMarkers[ii],
                          fvg.nodeFactors[ii, 1], fvg.nodeFactors[ii, 2], fvg.edgeFactors[ii])
      end
    else
      @error "Unsupported dimension " * string(fvg.dimension) * " in printing an instance of FiniteVolumeGrid type!"
      throw(ErrorException("Unsupported dimension " * string(fvg.dimension) *
                           " in printing an instance of FiniteVolumeGrid type!"))
    end
  end
  Format.printfmtln(fe, "N.o. boundary faces", fvg.numberOfFaces)
  if full || fvg.numberOfFaces > 0
    if fvg.dimension == 1
      Format.printfmtln("    {1:<16s}{2:>16s}{3:>16s}{4:>16s}",
                        "# boundary face", "1st vertex", "marker", "1st face fac.")
      fe1 = "    {1:<16d}{2:16d}{3:16d}{4:16.8e}"
      for ii = 1 : fvg.numberOfFaces
        Format.printfmtln(fe1, ii, fvg.faces[ii], fvg.faceMarkers[ii], fvg.faceFactors[ii])
      end
    else
      @error "Unsupported dimension " * string(fvg.dimension) * " in printing an instance of FiniteVolumeGrid type!"
      throw(ErrorException("Unsupported dimension " * string(fvg.dimension) *
                           " in printing an instance of FiniteVolumeGrid type!"))
    end
  end
end

function println(fvg :: FiniteVolumeGrid; full :: Bool = false)
  print(fvg; full = full)
end
