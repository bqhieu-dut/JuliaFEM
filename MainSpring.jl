include("Spring.jl")

# FEM for Spring systems
npoint      = 3 # Number of points
nele        = 3 # Number of spring elements
ndg         = 1 # Number of DOFs per node
nnd         = 2 # Number of nodes per elements
# Information of node and element
xsupport    = [1; 0; 0]
xload       = [0.; 0.; 1.]
joint1      = [1; 2; 1]
joint2      = [2; 3; 3]
k           = [1; 2; 3]
# 
sample1     = FemSpring(npoint,nele,ndg,nnd,xsupport,xload,joint1,joint2)
#
result_disp     = sample1[1]
result_axial    = sample1[2]