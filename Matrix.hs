module Matrix (
    Matrix,
    valid,
    element,
    row,
    column,
    rowCount,
    columnCount,
    square,
    symmetric,
    diagonal,
    adjugate,
    determinant,
    transpose,
    cofactor,
    complementMinor,
    complementSubMatrix,
) where

type Matrix a = [[a]]

-- | Returns true if each row of the matrix has the same number of elements.
valid :: Matrix a -> Bool
valid [] = True
valid (x:xs) = all (\ r -> length r == length x) xs

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
columnCount (x:xs) = length x

-- | Returns true if the matrix is square.
square :: Matrix a -> Bool
square m = rowCount m == columnCount m

-- | Returns true if the matrix is symmetric (i.e. equal to its transpose).
symmetric :: Eq a => Matrix a -> Bool
symmetric m = m == transpose m

-- | Returns the diagonal of the matrix (from the top left to bottom right).
diagonal :: Matrix a -> [a]
diagonal m = mapIndex (\ rs r -> element r r m) m

-- | Returns the adjugate/adjoint version of the matrix.
adjugate :: Num a => Matrix a -> Matrix a
adjugate m = transpose (mapIndex (\ rs r -> mapIndex (\ cs c -> cofactor m r c) rs) m)

-- | Returns the determinant of the matrix. The matrix must be square.
determinant :: Num a => Matrix a -> a
determinant [[x]] = x
determinant [[a0, a1], [b0, b1]] = (a0 * b1) - (a1 * b0)
determinant m = sum (mapIndex (\ x i -> element 0 i m * cofactor m 0 i) (row 0 m))

-- | Returns the transpose of the matrix.
transpose :: Matrix a -> Matrix a
transpose m = mapIndex (\ xs i -> column i m) m

-- | Returns the cofactor of the matrix at the specified row and column.
cofactor :: Num a => Matrix a -> Int -> Int -> a
cofactor m r c = ((-1) ^ (r + c)) * complementMinor m r c

-- | Returns the complement minor of the matrix at the specified row and column.
complementMinor :: Num a => Matrix a -> Int -> Int -> a
complementMinor m r c = determinant (complementSubMatrix m r c)

-- | Returns the complement sub-matrix of the matrix at the specified row and column.
complementSubMatrix :: Num a => Matrix a -> Int -> Int -> Matrix a
complementSubMatrix m r c = exceptIndex r [exceptIndex c rs | rs <- m]

{- Utility functions -}

exceptIndex :: Int -> [a] -> [a]
exceptIndex _ [] = []
exceptIndex i (x:xs) = if i == 0 then xs else x : exceptIndex (i - 1) xs

mapIndex :: (a -> Int -> b) -> [a] -> [b]
mapIndex f xs = mapIndexInternal f xs 0

mapIndexInternal :: (a -> Int -> b) -> [a] -> Int -> [b]
mapIndexInternal _ [] _ = []
mapIndexInternal f (x:xs) i = f x i : mapIndexInternal f xs (i + 1)