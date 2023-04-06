include("Truss2D.jl")

# FEM for Truss2D systems
npoint      = 4 # Number of points
nele        = 5 # Number of spring elements
ndg         = 2 # Number of DOFs per node
nnd         = 2 # Number of nodes per elements
# Information of nodes
xcord       = [0.; 5000.; 5000.;10000.]
ycord       = [0.; 0.; 5000*tan(pi/3); 5000*tan(pi/3)]
xsupport    = [1; 0; 0; 0] # Support conditions: 0-->free; 1-->fix
ysupport    = [1; 1; 0; 0]
xload       = [0.; 0.; 0.; 400*cos(pi/4)] # Concentrated load at node in x dir
yload       = [0.; 0.; 0.; -400*cos(pi/4)]
# Information of elements
joint1      = [1; 1; 2; 2; 3] # Node 1 of element i
joint2      = [2; 3; 3; 4; 4] # Node 2 of element i
A           = [10000.; 15000.; 15000.; 15000.; 10000.] # Stiffness of element i
E           = zeros(Float64,nele)
E           .= 200.0
#
Bt1   = FemTruss2D(npoint,nele,ndg,nnd,xcord,ycord,xsupport,ysupport,
                    xload,yload,joint1,joint2,E,A)
result_disp     = Bt1[1]
result_axial    = Bt1[2]
