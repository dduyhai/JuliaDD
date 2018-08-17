module Utilities
export convert_to_ordinal_number, extract_numeric_vector, extract_keyword, extract_value, is_data_line

function convert_to_ordinal_number(cardinalNumber::Int64)
  sufficeList = ["th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th"]
  tensDigit = div(cardinalNumber, 10) % 10
  unitsDigit = cardinalNumber % 10
  local suffice::String
  if (tensDigit == 1)
    suffice = "th"
  else
    suffice = sufficeList[unitsDigit + 1]
  end
  ordinalNumber = string(cardinalNumber) * suffice
  return ordinalNumber
end

function is_data_line(fileLine::String)
  if (length(fileLine) == 0)
    return false
  elseif (fileLine[1] == '#')
    return false
  end
  return true
end

function extract_numeric_vector(fileLine::String, numericType::DataType, vectorSize::Int64)
  if (vectorSize <= 0)
    throw(ErrorException("Size of a vector must be strictly positive!"))
  end
  if ((numericType != Int64) && (numericType != Float64))
    throw(ErrorException("Data type " * string(numericType) * " has not been supported yet!"))
  end

  wordsOfLine = split(fileLine)
  if (length(wordsOfLine) < vectorSize)
    throw(ErrorException("There is not enough data to extract from line!"))
  end
  vector = Array{numericType}(vectorSize)
  for index = 1 : vectorSize
    if (isnull(tryparse(numericType, wordsOfLine[index])) == true)
      throw(ErrorException("There is an invalid data type element in line!"))
    end
    vector[index] = parse(numericType, wordsOfLine[index])
  end

  return vector
end

function extract_numeric_vector!(fileLine :: String, numericType :: DataType, vectorSize :: Int64)
  if (vectorSize <= 0)
    throw(ErrorException("Size of a vector must be strictly positive!"))
  end
  if ((numericType != Int64) && (numericType != Float64))
    throw(ErrorException("Data type " * string(numericType) * " has not been supported yet!"))
  end
  wordsOfLine = split(fileLine)
  if (length(wordsOfLine) < vectorSize)
    throw(ErrorException("There is not enough data to extract from line!"))
  end
  vector = Array{numericType}(vectorSize)
  for index = 1 : vectorSize
    if (isnull(tryparse(numericType, wordsOfLine[index])) == true)
      throw(ErrorException("There is an invalid data type element in line!"))
    end
    vector[index] = parse(numericType, wordsOfLine[index])
  end
  if (vectorSize == length(wordsOfLine))
    fileLine = "";
  else
    fileLine = join(wordsOfLine[vectorSize + 1 : end], " ");
  end
  return (vector, fileLine);
end

function extract_keyword(fileLine :: String)
  wordsOfLine = split(fileLine)
  return uppercase(wordsOfLine[1])
end

function extract_value(fileLine :: String)
  return evaluate_unit_operation(fileLine)
end
end
