# Determine Spring element stiffness matrix
function springstiffness(x)
    x   :: Float64
    stiffnessmatrix = zeros(2,2)
    stiffnessmatrix = [x -x;-x x]
    return stiffnessmatrix
end
# Determine Truss element stiffness matrix
function trussstiffness(E,A,L)
    E   :: Float64
    A   :: Float64
    L   :: Float64
    stiffnessmatrix = zeros(4,4)
    stiffnessmatrix[1,1] = E*A/L
    stiffnessmatrix[1,3] = -E*A/L
    stiffnessmatrix[3,1] = -E*A/L
    stiffnessmatrix[3,3] = E*A/L
    return stiffnessmatrix
end
# Determine Frame element stiffness matrix
function Frame2Dstiffness(E,A,L,momentIn)
    E           :: Float64
    A           :: Float64
    L           :: Float64
    momentIn    :: Float64
    C1          = E*A/L
    C2          = E*momentIn/L^3
    stiffnessmatrix = zeros(6,6)
    stiffnessmatrix[1,1] = C1
    stiffnessmatrix[1,4] = -C1
    stiffnessmatrix[4,1] = -C1
    stiffnessmatrix[4,4] = C1
#    
    stiffnessmatrix[2,2] = 12*C2
    stiffnessmatrix[2,3] = 6*C2*L
    stiffnessmatrix[2,5] = -12*C2
    stiffnessmatrix[2,6] = 6*C2*L
#
    stiffnessmatrix[3,2] = 6*C2*L
    stiffnessmatrix[3,3] = 4*C2*L^2
    stiffnessmatrix[3,5] = -6*C2*L
    stiffnessmatrix[3,6] = 2*C2*L^2
#
    stiffnessmatrix[5,2] = -12*C2
    stiffnessmatrix[5,3] = -6*C2*L
    stiffnessmatrix[5,5] = 12*C2
    stiffnessmatrix[5,6] = -6*C2*L
#
    stiffnessmatrix[6,2] = 6*C2*L
    stiffnessmatrix[6,3] = 2*C2*L^2
    stiffnessmatrix[6,5] = -6*C2*L
    stiffnessmatrix[6,6] = 4*C2*L^2
#
    return stiffnessmatrix
end
# Determine Geometry Property of 2D truss
function GeoPropEle2D(xcord,ycord,joint1,joint2,m)
    xcord   :: Vector{Float64}
    ycord   :: Vector{Float64}
    joint1  :: Vector{Int64}
    joint2  :: Vector{Int64}
    m       :: Int64
#
    L       = zeros(m) # Length of element member
    cs      = zeros(m) # cosin of element member
    sn      = zeros(m) # sin of element member
#
    for i in 1:m
        j1 = joint1[i]
        j2 = joint2[i]
        delx    = xcord[j2] - xcord[j1]
        dely    = ycord[j2] - ycord[j1]
        L[i]    = sqrt(delx^2+dely^2)
        cs[i]   = delx/L[i]
        sn[i]   = dely/L[i]
    end
    return L, cs, sn
end
# Transformation matrix 2D truss
function TransMatrix2DTruss(cs,sn)
    cs :: Float64
    sn :: Float64
#
    TransMatrix = zeros(4,4)
    TransMatrix[1,1] = cs
    TransMatrix[1,2] = sn
    TransMatrix[2,1] = -sn
    TransMatrix[2,2] = cs
    TransMatrix[3:4,3:4] = TransMatrix[1:2,1:2]
#
    return TransMatrix
end
# Transformation matrix 2D truss
function TransMatrix2DFrame(cs,sn)
    cs :: Float64
    sn :: Float64
#
    TransMatrix = zeros(6,6)
    TransMatrix[1,1] = cs
    TransMatrix[1,2] = sn
    TransMatrix[2,1] = -sn
    TransMatrix[2,2] = cs
    TransMatrix[3,3] = 1.0
#
    TransMatrix[4:6,4:6] = TransMatrix[1:3,1:3]
#
    return TransMatrix
end
