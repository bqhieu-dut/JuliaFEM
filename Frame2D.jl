using SparseArrays
include("utils.jl")
#
struct FemFrame2D
    npoint      :: Int64            # Number of points
    nele        :: Int64            # Number of elements
    ndg         :: Int64            # Number of DOFs per node
    nnd         :: Int64            # Number of nodes per element
    xcord       :: Vector{Float64}  # X coordinate of nodes
    ycord       :: Vector{Float64}  # Y coordinate of nodes
    xsupport    ::Vector{Int64}     # Support conditions in X direction
                                    # 0 is free, 1 is fix
    ysupport    ::Vector{Int64}     # Support conditions in Y direction
    rsupport    ::Vector{Int64}     # Support conditions in rotation
    xload       ::Vector{Float64}   # Point forces in X direction
    yload       ::Vector{Float64}   # Point forces in Y direction
    rload       ::Vector{Float64}   # Moment
    joint1      ::Vector{Int64}     # Node 1 of elememts
    joint2      ::Vector{Int64}     # Node 2 of elememts
    E           ::Vector{Float64}   # Young modulus of material
    A           ::Vector{Float64}   # Cross-section area of member
    In          ::Vector{Float64}   # Moment of inertial of member
    function FemFrame2D(npoint,nele,ndg,nnd,xcord,ycord,xsupport,ysupport,
                        rsupport,xload,yload,rload,joint1,joint2,E,A,In)
        #This function will return the displacment vector and member force vector
        # ****************CACULATION***********************
        skk             = spzeros(ndg*npoint,ndg*npoint)    # Stiffness matrix of the whole structures
        disp            = zeros(ndg*npoint)                 # Displacment vector
        load            = zeros(ndg*npoint)                 # Load vector
        support         = zeros(ndg*npoint)                 # Support conditions: 0-->free; 1-->fix
        nodalforces     = zeros(nele,6)                     # Nodal force vector

        # Prepare for caculation
        for m in 1:npoint
            load[3m-2:3m]     = [xload[m];yload[m];rload[m]]
            support[3m-2:3m]  = [xsupport[m];ysupport[m];rsupport[m]]
        end
        # Determine properties of 2D truss elements
        L, cs, sn = GeoPropEle2D(xcord,ycord,joint1,joint2,nele)
        # Define the stiffness matrix of the whole structure skk
        for m in 1:nele
            Ee = E[m]
            Ae = A[m]
            Ie = In[m]
            Le = L[m]
            cse = cs[m]
            sne = sn[m]
            # Define the stiffness matrix of element in local coordinate
            kee = Frame2Dstiffness(Ee,Ae,Le,Ie)
            # Define the transformation matrix of element
            Te  = TransMatrix2DFrame(cse,sne)
            # Define the stiffness matrix of element in global coordinate
            kee = Te'*kee*Te
            # Define the direction vector
            i = joint1[m]
            j = joint2[m]
            n = [3i-2;3i-1;3i;3j-2;3j-1;3j]  # It should be noted that Julia is column priority
            skk[n,n] = skk[n,n] + kee
        end
        # Define the reaction force and displacment
        indf = findall(x -> x==0,support)
        inds = findall(x -> x!=0,support)
        # Define part matrix of rearrangement structure stiffness matrix
        kff  = skk[indf,indf]
        kss  = skk[inds,inds]
        ksf  = skk[inds,indf]
        kfs  = skk[indf,inds]
        # Define unknown displacement
        disp[indf] = kff\load[indf]
        # Define reaction force
        load[inds] = ksf*disp[indf]
        # Define axial force vector
        for m in 1:nele
            Ee = E[m]
            Ae = A[m]
            Ie = In[m]
            Le = L[m]
            cse = cs[m]
            sne = sn[m]
            # Define the stiffness matrix of element in local coordinate
            kee = Frame2Dstiffness(Ee,Ae,Le,Ie)
            # Define the transformation matrix of element
            Te  = TransMatrix2DFrame(cse,sne)
            # Define the direction vector
            i = joint1[m]
            j = joint2[m]
            n = [3i-2;3i-1;3i;3j-2;3j-1;3j] 
            # Define force
            force = kee*Te*disp[n]
            # save in axial force vector
            nodalforces[m,:] = force             
        end
        return disp, nodalforces
    end
end
