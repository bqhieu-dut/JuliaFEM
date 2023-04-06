using SparseArrays
include("utils.jl")
#
struct FemSpring
    npoint  :: Int64        # Number of points
    nele    :: Int64        # Number of elements
    ndg     :: Int64        # Number of DOFs per node
    nnd     :: Int64        # Number of nodes per element
    xsupport    ::Vector{Int64} # Support conditions in X direction
                                # 0 is free, 1 is fix
    xload       ::Vector{Float64} # Point forces in X direction
    joint1      ::Vector{Int64} # Node 1 of elememts
    joint2      ::Vector{Int64} # Node 2 of elememts
    k           ::Vector{Float64}  # Stiffness of spring element
    function FemSpring(npoint,nele,ndg,nnd,xsupport,xload,joint1,joint2)
        skk         = spzeros(ndg*npoint,ndg*npoint)    # Stiffness matrix of the whole structures
        disp        = zeros(ndg*npoint)                 # Displacment vector
        load        = zeros(ndg*npoint)                 # Load vector
        support     = zeros(ndg*npoint)                 # Support conditions: 0-->free; 1-->fix
        axial       = zeros(nele)                       # Axial force vector
    # Prepare for caculation
        for m in 1:npoint
            load[m]     = xload[m]
            support[m]  = xsupport[m]
        end
    # Define the stiffness matrix of the whole structure skk
        for m in 1:nele
            ke      = k[m]
            # Define the stiffness matrix of element
            kee     = [ke -ke;-ke ke]
            # Define the direction vector
            i = joint1[m]
            j = joint2[m]
            n = [i;j]
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
            ke      = k[m]
            # Define the stiffness matrix of element
            kee     = [ke -ke;-ke ke]
            # Define the direction vector
            i = joint1[m]
            j = joint2[m]
            n = [i;j]
            # Determine axial force
            force = kee*disp[n]
            # Save in axial force vector
            axial[m] = force[2]
        end
    #
        return disp, axial
    end
end