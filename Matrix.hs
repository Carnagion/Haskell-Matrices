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
diagonal m = mapIndex (\ rs r -> element r r m) m

-- | Returns the adjugate/adjoint version of the matrix.
adjugate :: Num a => Matrix a -> Matrix a
adjugate m = transpose (mapIndex (\ rs r -> mapIndex (\ cs c -> cofactor m r c) rs) m)

-- | Returns the determinant of the matrix. Throws an error if the matrix is not square.
determinant :: Num a => Matrix a -> a
determinant [[x]] = x
determinant [[a0, a1], [b0, b1]] = (a0 * b1) - (a1 * b0)
determinant m = if square m
                then Prelude.sum (mapIndex (\ x i -> element 0 i m * cofactor m 0 i) (row 0 m))
                else error "Cannot calculate the determinant of a rectangular matrix!"

-- | Returns the transpose of the matrix.
transpose :: Matrix a -> Matrix a
transpose m = mapIndex (\ xs i -> column i m) m

-- | Returns the inverse of the matrix. Throws an error if the matrix is singular.
inverse :: (Fractional a, Eq a) => Matrix a -> Matrix a
inverse m = if singular m
            then error "Cannot invert a singular matrix!"
            else scalarProduct (adjugate m) (1.0 / determinant m)

-- | Returns the cofactor of the matrix at the specified row and column.
cofactor :: Num a => Matrix a -> Int -> Int -> a
cofactor m r c = ((-1) ^ (r + c)) * complementMinor m r c

-- | Returns the complement minor of the matrix at the specified row and column.
complementMinor :: Num a => Matrix a -> Int -> Int -> a
complementMinor m r c = determinant (complementSubmatrix m r c)

-- | Returns the complement sub-matrix of the matrix at the specified row and column.
complementSubmatrix :: Num a => Matrix a -> Int -> Int -> Matrix a
complementSubmatrix m r c = exceptIndex r [exceptIndex c rs | rs <- m]

-- | Multiplies all elements of the matrix with the specified scalar.
scalarProduct :: Num a => Matrix a -> a -> Matrix a
scalarProduct m s = map (map (* s)) m

-- | Adds two matrices. Throws an error if the matrices have different sizes.
sum :: Num a => Matrix a -> Matrix a -> Matrix a
sum m1 m2 = if size m1 == size m2
            then mapIndex (\ rs r -> mapIndex (\ cs c -> element r c m1 + element r c m2) rs) m1
            else error "Cannot add matrices of different sizes!"

-- | Subtracts two matrices. Throws an error if the matrices have different sizes.
difference :: Num a => Matrix a -> Matrix a -> Matrix a
difference m1 m2 = if size m1 == size m2
                   then mapIndex (\ rs r -> mapIndex (\ cs c -> element r c m1 - element r c m2) rs) m1
                   else error "Cannot subtract matrices of different sizes!"

{- Utility functions -}

exceptIndex :: Int -> [a] -> [a]
exceptIndex _ [] = []
exceptIndex i (x:xs) = if i == 0 then xs else x : exceptIndex (i - 1) xs

mapIndex :: (a -> Int -> b) -> [a] -> [b]
mapIndex f xs = mapIndexInternal f xs 0

mapIndexInternal :: (a -> Int -> b) -> [a] -> Int -> [b]
mapIndexInternal _ [] _ = []
mapIndexInternal f (x:xs) i = f x i : mapIndexInternal f xs (i + 1)