include("Frame2D.jl")

# FEM for Frame2D systems
npoint      = 9 # Number of points
nele        = 10 # Number of spring elements
ndg         = 3 # Number of DOFs per node
nnd         = 2 # Number of nodes per elements
# Information of nodes
xcord       = [0.0; 6.0; 10.0; 0.0; 6.0; 10.0; 0.0; 6.0; 10.0] 
ycord       = [0.0; 0.0; 0.0; 4.0; 4.0; 4.0; 8.0; 8.0; 8.0]
xsupport    = [1; 1; 1; 0; 0; 0; 0; 0; 0] # Support conditions: 0-->free; 1-->fix
ysupport    = [1; 1; 1; 0; 0; 0; 0; 0; 0]
rsupport    = [1; 1; 1; 0; 0; 0; 0; 0; 0]
xload       = [0.0; 0.0; 0.0; 20.0; 0.0; 0.0; 20.0; 0.0; 0.0] # Concentrated load at node in x dir
yload       = zeros(Float64,npoint)
rload       = zeros(Float64,npoint)
# Information of elements
joint1      = [1; 2; 3; 4; 5; 6; 4; 5; 7; 8] # Node 1 of element i
joint2      = [4; 5; 6; 7; 8; 9; 5; 6; 8; 9] # Node 2 of element i
A           = [0.08; 0.08; 0.08; 0.08; 0.08; 0.08; 0.09; 0.09; 0.09; 0.09] # Stiffness of element i
E           = zeros(Float64,nele)
E           .= 27.0e6
I_b         = 0.2*0.45^3/12
I_c         = 0.2*0.4^3/12
In          = [I_c; I_c; I_c; I_c; I_c; I_c; I_b; I_b; I_b; I_b]
#
SampleFrame   = FemFrame2D(npoint,nele,ndg,nnd,xcord,ycord,xsupport,ysupport,
                            rsupport,xload,yload,rload,joint1,joint2,E,A,In)
result_disp, result_force = SampleFrame
