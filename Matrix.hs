module Matrix (
    Matrix,
    valid,
    element,
    row,
    column,
    rowCount,
    columnCount,
    square,
    singular,
    symmetric,
    diagonal,
    adjugate,
    determinant,
    transpose,
    inverse,
    cofactor,
    complementMinor,
    complementSubmatrix,
    scalarProduct,
    matrixProduct,
    Matrix.sum,
    difference,
) where

type Matrix a = [[a]]

-- | Returns true if each row of the matrix has the same number of elements.
valid :: Matrix a -> Bool
valid [] = True
valid (xs:xss) = all (\ xs' -> length xs' == length xs) xss

-- | Returns the element at the specified row and column in the matrix.
element :: Int -> Int -> Matrix a -> a
element r c m = (m !! r) !! c

-- | Returns the row at the specified index of the matrix.
row :: Int -> Matrix a -> [a]
row i m = m !! i

-- | Returns the column at the specified index of the matrix.
column :: Int -> Matrix a -> [a]
column i m = [xs !! i | xs <- m]

-- | Returns the number of rows in the matrix.
rowCount :: Matrix a -> Int
rowCount = length

-- | Returns the number of columns in the matrix.
columnCount :: Matrix a -> Int
columnCount [] = 0
columnCount (x:_) = length x

size :: Matrix a -> (Int, Int)
size m = (rowCount m, columnCount m)

-- | Returns true if the matrix is square.
square :: Matrix a -> Bool
square m = rowCount m == columnCount m

-- | Returns true if the determinant of the matrix is 0.
singular :: (Num a, Eq a) => Matrix a -> Bool
singular m = determinant m == 0

-- | Returns true if the matrix is symmetric (i.e. equal to its transpose).
symmetric :: Eq a => Matrix a -> Bool
symmetric m = m == transpose m

-- | Returns the diagonal of the matrix (from the top left to bottom right).
diagonal :: Matrix a -> [a]
diagonal m = mapIndex (\ xs i -> element i i m) m

-- | Returns the adjugate/adjoint version of the matrix.
adjugate :: Num a => Matrix a -> Matrix a
adjugate m = transpose (mapIndex (\ xs r -> mapIndex (\ _ c -> cofactor m r c) xs) m)

-- | Returns the determinant of the matrix. Throws an error if the matrix is not square.
determinant :: Num a => Matrix a -> a
determinant [[x]] = x
determinant [[x00, x01], [x10, x11]] = (x00 * x11) - (x01 * x10)
determinant m = if square m
                then Prelude.sum (mapIndex (\ _ i -> element 0 i m * cofactor m 0 i) (row 0 m))
                else error "cannot calculate determinant of a rectangular matrix"

-- | Returns the transpose of the matrix.
transpose :: Matrix a -> Matrix a
transpose m = mapIndex (\ _ i -> column i m) (head m)

-- | Returns the inverse of the matrix. Throws an error if the matrix is singular.
inverse :: (Fractional a, Eq a) => Matrix a -> Matrix a
inverse m = if singular m
            then error "cannot invert a singular matrix"
            else scalarProduct (adjugate m) (1.0 / determinant m)

-- | Returns the cofactor of the matrix at the specified row and column.
cofactor :: Num a => Matrix a -> Int -> Int -> a
cofactor m r c = ((-1) ^ (r + c)) * complementMinor m r c

-- | Returns the complement minor of the matrix at the specified row and column.
complementMinor :: Num a => Matrix a -> Int -> Int -> a
complementMinor m r c = determinant (complementSubmatrix m r c)

-- | Returns the complement sub-matrix of the matrix at the specified row and column.
complementSubmatrix :: Num a => Matrix a -> Int -> Int -> Matrix a
complementSubmatrix m r c = exceptIndex r [exceptIndex c xs | xs <- m]

-- | Multiplies all elements of the matrix with the specified scalar.
scalarProduct :: Num a => Matrix a -> a -> Matrix a
scalarProduct m s = map (map (* s)) m

-- | Multiplies two matrices together. Throws an error if the matrices have different sizes.
matrixProduct :: Num a => Matrix a -> Matrix a -> Matrix a
matrixProduct m1 m2 = if size m1 == size (transpose m2)
                      then map (\ xs -> [scalarProductList xs ys | ys <- transpose m2]) m1
                      else error "cannot multiply matrices of different sizes"

-- | Adds two matrices. Throws an error if the matrices have different sizes.
sum :: Num a => Matrix a -> Matrix a -> Matrix a
sum m1 m2 = if size m1 == size m2
            then mapIndex (\ xs r -> mapIndex (\ _ c -> element r c m1 + element r c m2) xs) m1
            else error "cannot add matrices of different sizes"

-- | Subtracts two matrices. Throws an error if the matrices have different sizes.
difference :: Num a => Matrix a -> Matrix a -> Matrix a
difference m1 m2 = if size m1 == size m2
                   then mapIndex (\ xs r -> mapIndex (\ _ c -> element r c m1 - element r c m2) xs) m1
                   else error "cannot subtract matrices of different sizes"

{- Utility functions -}

exceptIndex :: Int -> [a] -> [a]
exceptIndex _ [] = []
exceptIndex i (x:xs) = if i == 0
                       then xs
                       else x : exceptIndex (i - 1) xs

mapIndex :: (a -> Int -> b) -> [a] -> [b]
mapIndex f xs = mapIndexInternal f xs 0

mapIndexInternal :: (a -> Int -> b) -> [a] -> Int -> [b]
mapIndexInternal _ [] _ = []
mapIndexInternal f (x:xs) i = f x i : mapIndexInternal f xs (i + 1)

scalarProductList :: Num a => [a] -> [a] -> a
scalarProductList xs ys = Prelude.sum [x * y | (x, y) <- zip xs ys]