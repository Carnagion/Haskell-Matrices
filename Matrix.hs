module Matrix (
    Matrix,
    valid,
    Matrix.null,
    identity,
    element,
    row,
    column,
    rowCount,
    columnCount,
    size,
    square,
    singular,
    symmetric,
    triangular,
    diagonal,
    adjugate,
    determinant,
    rank,
    transpose,
    inverse,
    cofactor,
    complementMinor,
    complementSubmatrix,
    scalarProduct,
    matrixProduct,
    Matrix.sum,
    difference,
    gaussian,
) where

import Data.Ix (Ix(range))
import Data.List (sort)

-- | A matrix represented as a list of lists, where each inner list is a row.
type Matrix a = [[a]]

-- | Returns true if each row of the matrix has the same number of elements.
valid :: Matrix a -> Bool
valid [] = True
valid (xs:xss) = all (\ xs' -> length xs' == length xs) xss

-- | Constructs a null matrix of the specified size.
null :: Int -> Int -> Matrix Int
null r c = map (\ _ -> map (const 0) (range (0, c - 1))) (range (0, r - 1))

-- | Constructs an identity matrix of the specified size. Throws an error if the size is not a square.
identity :: Int -> Int -> Matrix Int
identity r c = if r == c
               then map (\ r' -> map (\ c' -> if r' == c' then 1 else 0) (range (0, c - 1))) (range (0, r - 1))
               else error "cannot construct rectangular identity matrix"

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

-- | Returns the dimensions of the matrix as (rows, columns).
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

-- | Returns true if the matrix is triangular or could be turned into a triangular square matrix by removing some columns.
triangular :: (Num a, Eq a) => [[a]] -> Bool
triangular m = sort (map (length . takeWhile (== 0)) m) == range (0, rowCount m - 1)

-- | Returns the diagonal of the matrix (from the top left to bottom right).
diagonal :: Matrix a -> [a]
diagonal m = if square m
             then mapIndex (\ _ i -> element i i m) m
             else error "cannot get the diagonal of a rectangular matrix"

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

-- | Returns the rank of the matrix.
rank :: (Num a, Eq a) => Matrix a -> Int
rank m = (if square m then rankSquare else rankRectangular) m

-- | Returns the transpose of the matrix.
transpose :: Matrix a -> Matrix a
transpose m = mapIndex (\ _ i -> column i m) (row 0 m)

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
complementSubmatrix :: Matrix a -> Int -> Int -> Matrix a
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

-- | Returns a matrix resulting from applying Gaussian elimination on the matrix.
gaussian :: (Fractional a, Eq a) => Matrix a -> Matrix a
gaussian = gaussianFrom 0

{- Utility functions -}

exceptIndex :: Int -> [a] -> [a]
exceptIndex _ [] = []
exceptIndex i (x:xs) = if i == 0
                       then xs
                       else x : exceptIndex (i - 1) xs

mapIndex :: (a -> Int -> b) -> [a] -> [b]
mapIndex f xs = mapIndexFrom f xs 0

mapIndexFrom :: (a -> Int -> b) -> [a] -> Int -> [b]
mapIndexFrom _ [] _ = []
mapIndexFrom f (x:xs) i = f x i : mapIndexFrom f xs (i + 1)

scalarProductList :: Num a => [a] -> [a] -> a
scalarProductList xs ys = Prelude.sum [x * y | (x, y) <- zip xs ys]

squareSubmatrices :: Matrix a -> [Matrix a]
squareSubmatrices m = concat (mapIndex (\ xs r -> mapIndex (\ _ c -> complementSubmatrix m r c) xs) m)

rankSquare :: (Num a, Eq a) => Matrix a -> Int
rankSquare [[0]] = 0
rankSquare [[_]] = 1
rankSquare m = if singular m
               then maximum (map rankSquare (squareSubmatrices m))
               else rowCount m

rankRectangular :: (Num a, Eq a) => Matrix a -> Int
rankRectangular [_] = 1
rankRectangular m = maximum (map rank (mapIndex (\ _ i -> exceptIndex i m') m'))
                    where m' = if rowCount m > columnCount m
                               then m
                               else transpose m

gaussianFrom :: (Fractional a, Eq a) => Int -> Matrix a -> Matrix a
gaussianFrom i m = if triangular m 
                   then m
                   else gaussianFrom (i + 1) m'
                   where m' = take (i + 1) m ++ mapIndex (\ xs r -> mapIndex (\ x c -> x + (- (element (r + i + 1) i m / element i i m)) * element i c m) xs) (drop (i + 1) m)