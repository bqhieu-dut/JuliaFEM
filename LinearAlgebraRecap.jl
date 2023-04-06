matrixA = [1 2 3;4 5 6;7 8 9]
matrixB = [1 2 3]       # Matrix not vector
vectorC = [1,2,3]       # Column vector, Julia
vectorD = [1;2;3]
# Zeros matrix and vector
zerosmatrix = zeros(Float64,3,3)
zerosvector = zeros(Int64,3)
# Identify matrix
using LinearAlgebra
identifymatrix = Matrix{Float64}(I,3,3)
# Check equal matrix or vector
LinearAlgebra.isequal(matrixA,matrixB)
LinearAlgebra.isequal(vectorC,vectorD)
# Addition 
matrixAddition1 = matrixA + identifymatrix
matrixAddition2 = matrixA + zerosmatrix
LinearAlgebra.isequal(matrixA,matrixAddition2)
# Multiplication by a scalar
λ = 2.0
matrixA_λ = λ*matrixA
# Transpose
transposematrixA = transpose(matrixA)
LinearAlgebra.issymmetric(matrixA)
# Matrix Multiplication
A_sample = [1 2 4;-3 0 7;9 1 5]
B_sample = [-1 3;-3 1;1 0]
C_sample = A_sample*B_sample
# Determinant
det_A_sample = det(A_sample)
det_matrixA  = det(matrixA)
# Inverse
A_sample_inv = inv(A_sample)
matrixA_inv  = pinv(matrixA)

