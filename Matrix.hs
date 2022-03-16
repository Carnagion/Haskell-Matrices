type Matrix a = [[a]]

isValid :: Matrix a -> Bool
isValid [] = True
isValid m = all (\ r -> length r == length (head m)) m

element :: Int -> Int -> Matrix a -> a
element r c m = (m !! r) !! c

row :: Int -> Matrix a -> [a]
row i m = m !! i

column :: Int -> Matrix a -> [a]
column i m = [xs !! i | xs <- m]

rowCount :: Matrix a -> Int
rowCount = length

columnCount :: Matrix a -> Int
columnCount = length . head

isSquare :: Matrix a -> Bool
isSquare m = rowCount m == columnCount m

determinant :: Num a => Matrix a -> a
determinant [[x]] = x
determinant [[a0, a1], [b0, b1]] = (a0 * b1) - (a1 * b0)
determinant m = sum (cofactorList m 0 0)

cofactor :: Num a => Matrix a -> Int -> Int -> a
cofactor m r c = ((-1) ^ (r + c)) * (complementMinor m r c)

complementMinor :: Num a => Matrix a -> Int -> Int -> a
complementMinor m r c = determinant (complementSubMatrix m r c)

complementSubMatrix :: Num a => Matrix a -> Int -> Int -> Matrix a
complementSubMatrix m r c = exceptIndex r [exceptIndex c rs | rs <- m]

transpose :: Matrix a -> Matrix a
transpose m = columnsList m 0

exceptIndex :: Int -> [a] -> [a]
exceptIndex _ [] = []
exceptIndex i (x:xs) = if i == 0 then xs else x : exceptIndex (i - 1) xs

cofactorList :: Num a => Matrix a -> Int -> Int -> [a]
cofactorList m r c = if c < columnCount m then ((element r c m) * (cofactor m r c)) : cofactorList m r (c + 1) else []

columnsList :: Matrix a -> Int -> Matrix a
columnsList m i = if i < columnCount m then (column i m) : columnsList m (i + 1) else []