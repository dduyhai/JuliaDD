"""
    SimplexGrid

A composite data type that holds information of a simplex grid.

It consists of following data field

* `dimension`: the dimension of the grid. It should be 1 or 2.
* `numberOfNodes`: the number of nodes of the grid.
* `numberOfSimplices`: the number of simplices of the grid.
* `numberOfFaces`: the number of boundary facets.
* `numberOfSimplexMarkers`: the number of markers that are assigned to simplices.
* `nodes`: coordinates of nodes. Its size is `(numberOfNodes, dimension)`.
* `simplices`: node-index of all simplices of the grid. Its size is `(numberOfSimplices, dimension + 1)`.
* `faces`: node-index of all boundary facets of the grid. Its size is `(numberOfFaces, dimension)`.
* `simplexMarkers`: markers of simplices. Its size is `numberOfSimplices`. The simplex marker must not be negative.
* `faceMarkers`: markers of boundary facets. It size is `numberOfFaces`. The face marker must not be negative.

It is worth to notice that the simplex markers `simplexMarkers` and face markers `faceMarkers` are completely irrelevant.
The simplex markers are used to apply different value of material parameters while the face markers
are used to impose different boundary conditions.
"""
mutable struct SimplexGrid
  dimension :: Int64
  numberOfNodes :: Int64
  numberOfSimplices :: Int64
  numberOfFaces :: Int64
  numberOfSimplexMarkers :: Int64
  numberOfFaceMarkers :: Int64
  nodes :: Array{Float64, 2}
  simplices :: Array{Int64, 2}
  faces :: Array{Int64, 2}
  simplexMarkers :: Array{Int64, 1}
  faceMarkers :: Array{Int64, 1}

  function SimplexGrid()
    this = new()
    this.dimension = 0
    this.numberOfNodes = 0
    this.numberOfSimplices = 0
    this.numberOfFaces = 0
    this.numberOfSimplexMarkers = 0
    this.numberOfFaceMarkers = 0
    return this
  end
  function SimplexGrid(sg_filepath :: String)
    grid = nothing
    try
      grid = import_simplex_grid(sg_filepath)
      println(grid.dimension)
    catch import_error
      rethrow(import_error)
    end
    return grid
  end
end

function import_simplex_grid(sgFileName::String)
  if (isfile(sgFileName) == false)
    throw(ErrorException("Unavailable the simplex grid file " * string(sgFileName)))
  end

  sgFileStream = open(sgFileName, "r");
  lineNumber = 0;
  dataBlock = 0;
  mustReadDataBlock = 5;
  dataBlockName = ["DIMENSION", "NODES", "CELLS", "FACES", "END"];
  isAvailableDataBlock = falses(mustReadDataBlock);
  grid = SimplexGrid();
  try
    while !eof(sgFileStream)
      sgFileLine = string(strip(readline(sgFileStream)));
      lineNumber += 1
      if !Utilities.is_data_line(sgFileLine)
        continue;
      end
      println(sgFileLine)
      keyword = Utilities.extract_keyword(sgFileLine);
      println(keyword)
      if (keyword == "DIMENSION")
        dataBlock += 1;
        grid.dimension = get_simplex_grid_dimension!(sgFileStream, lineNumber);
        isAvailableDataBlock[1] = true;
      elseif (keyword == "NODES")
        dataBlock += 1;
        grid.nodes = get_simplex_grid_nodes!(sgFileStream, lineNumber, grid);
        grid.numberOfNodes = size(grid.nodes)[1];
        isAvailableDataBlock[2] = true;
      elseif (keyword == "CELLS")
        dataBlock += 1;
        grid.cells, grid.cellMarkers = get_simplex_grid_cells!(sgFileStream, lineNumber, grid);
        grid.numberOfCells = size(grid.cells)[1];
        grid.numberOfCellMarkers = maximum(grid.cellMarkers);
        isAvailableDataBlock[3] = true;
      elseif (keyword == "FACES")
        dataBlock += 1;
        grid.faces, grid.faceMarkers = get_simplex_grid_faces!(sgFileStream, lineNumber, grid);
        grid.numberOfFaces = size(grid.faces)[1];
        grid.numberOfFaceMarkers = maximum(grid.faceMarkers);
        isAvailableDataBlock[4] = true;
      elseif (keyword == "END")
        dataBlock += 1;
        close(sgFileStream);
        isAvailableDataBlock[5] = true;
        break;
      else
        throw(ErrorException("Unknown keyword " * string(keyword) * "at line number " * 
                             string(lineNumber) * " in simplex grid file " * string(sgFileName) * "."));
      end
    end
  catch importSimplexGridError
    println(importSimplexGridError)
    rethrow(importSimplexGridError);
  finally
    if isopen(sgFileStream)
      close(sgFileStream);
    end
    if (dataBlock < mustReadDataBlock)
      println(isAvailableDataBlock)
      for blockIndex = 1 : mustReadDataBlock
        if !isAvailableDataBlock[blockIndex]
          throw(ErrorException("Unavailable data block " * string(dataBlockName[blockIndex]) * "."))
        end
      end
    end
  end

  return grid;
end

function get_simplex_grid_dimension!(fileStream::IOStream, lineNumber::Int64)
  fileLine = string();
  while true
    if eof(fileStream)
      throw(ErrorException("Line " * string(lineNumber) * ": there is no data in DIMENSION block!"));
    end
    fileLine = strip(readline(fileStream));
    lineNumber += 1;
    if Utilities.is_data_line(fileLine)
      break;
    end
  end
  dataVector = Utilities.extract_numeric_vector(fileLine, Int64, 1);
  return dataVector[1];
end

function get_simplex_grid_nodes!(fileStream::IOStream, lineNumber::Int64, grid::SimplexGrid)
  fileLine = string();
  while true
    if eof(fileStream)
      throw(ErrorException("Line " * string(lineNumber) * ": there is not enough data in NODES block!"));
    end
    fileLine = strip(readline(fileStream));
    lineNumber += 1;
    if Utilities.is_data_line(fileLine)
      break;
    end
  end
  dataVector = Utilities.extract_numeric_vector(fileLine, Int64, 2);
  numberOfNodes = dataVector[1];
  if (grid.dimension != dataVector[2])
    throw(ErrorException("Line " * string(lineNumber) * ": dimension of nodes and of grid are not compatible!"));
  end
  minNumberOfNodes = grid.dimension
  if (numberOfNodes < minNumberOfNodes)
    throw(ErrorException("Line " * string(lineNumber) * ": number of nodes must be greater than or equal to " * string(minNumberOfNodes) * "!"));
  end
  nodes = zeros(numberOfNodes, grid.dimension);
  for nodeIndex = 1 : numberOfNodes
    while true
      if eof(fileStream)
        throw(ErrorException("Line " * string(lineNumber) * ": there is not enough data in NODES block!"));
      end
      fileLine = strip(readline(fileStream));
      lineNumber += 1;
      if Utilities.is_data_line(fileLine)
        break;
      end
    end
    dataVector = Utilities.extract_numeric_vector(fileLine, Float64, grid.dimension);
    nodes[nodeIndex, :] = dataVector[1 : end];
  end
  return nodes;
end

function get_simplex_grid_cells!(fileStream::IOStream, lineNumber::Int64, grid::SimplexGrid)
  fileLine = string();
  while true
    if eof(fileStream)
      throw(ErrorException("Line " * string(lineNumber) * ": there is not enough data in CELLS block!"));
    end
    fileLine = strip(readline(fileStream));
    lineNumber += 1;
    if Utilities.is_data_line(fileLine)
      break;
    end
  end
  dataVector = Utilities.extract_numeric_vector(fileLine, Int64, 1);
  numberOfCells = dataVector[1];
  minNumberOfCells = 1;
  if (numberOfCells < minNumberOfCells)
    throw(ErrorException("Line " * string(lineNumber) * ": number of cells must be strictly positive!"));
  end
  cells = Array{Int64}(numberOfCells, grid.dimension + 1);
  cellMarkers = Array{Int64}(numberOfCells);
  for cellIndex = 1 : numberOfCells
    while true
      if eof(fileStream)
        throw(ErrorException("Line " * string(lineNumber) * ": there is not enough data in CELLS block!"));
      end
      fileLine = strip(readline(fileStream));
      lineNumber += 1;
      if Utilities.is_data_line(fileLine)
        break;
      end
    end
    dataVector = Utilities.extract_numeric_vector(fileLine, Int64, grid.dimension + 2);
    cells[cellIndex, :] = dataVector[1 : (end - 1)];
    cellMarkers[cellIndex] = dataVector[end];
  end
  return cells, cellMarkers;
end

function get_simplex_grid_faces!(fileStream::IOStream, lineNumber::Int64, grid::SimplexGrid)
  fileLine = string();
  while true
    if eof(fileStream)
      throw(ErrorException("Line " * string(lineNumber) * ": there is not enough data in FACES block!"));
    end
    fileLine = strip(readline(fileStream));
    lineNumber += 1;
    if Utilities.is_data_line(fileLine)
      break;
    end
  end

  dataVector = Utilities.extract_numeric_vector(fileLine, Int64, 1);
  numberOfFaces = dataVector[1];
  minNumberOfFaces = 1;
  if (numberOfFaces < minNumberOfFaces)
    throw(ErrorException("Line " * string(lineNumber) * ": number of faces must be strictly positive!"));
  end
  faces = Array{Int64}(numberOfFaces, grid.dimension);
  faceMarkers = Array{Int64}(numberOfFaces);
  for faceIndex = 1 : numberOfFaces
    while true
      if eof(fileStream)
        throw(ErrorException("Line " * string(lineNumber) * ": there is not enough data in FACES block!"));
      end
      fileLine = strip(readline(fileStream));
      lineNumber += 1;
      if Utilities.is_data_line(fileLine)
        break;
      end
    end
    dataVector = Utilities.extract_numeric_vector(fileLine, Int64, grid.dimension + 1);
    faces[faceIndex, :] = dataVector[1 : (end - 1)];
    faceMarkers[faceIndex] = dataVector[end];
  end
  return faces, faceMarkers;
end
