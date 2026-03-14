#-------------------------------------------------------------------------------
# mfoil.py: class and methods for subsonic airfoil analysis (v 2023-06-28)
#
# Copyright (C) 2023 Krzysztof J. Fidkowski
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#-------------------------------------------------------------------------------

import numpy as np
from scipy.interpolate import CubicSpline, interp1d
from scipy import sparse
from scipy.sparse import linalg
import matplotlib.pyplot as plt
from matplotlib import gridspec
import copy


#-------------------------------------------------------------------------------
class Geom:   # geometry
    def __init__(S):
        S.chord = 1.                  # chord length
        S.wakelen = 1.                # wake extent length, in chords
        S.npoint = 1                  # number of geometry representation points
        S.name = 'noname'             # airfoil name, e.g. NACA XXXX
        S.xpoint = []                 # point coordinates, [2 x npoint]
        S.xref = np.array([0.25, 0])  # moment reference point

#-------------------------------------------------------------------------------
class Panel:   # paneling
    def __init__(S):
        S.N = 0            # number of nodes
        S.x = []           # node coordinates, [2 x N]
        S.s = []           # arclength values at nodes
        S.t = []           # dx/ds, dy/ds tangents at nodes

#-------------------------------------------------------------------------------
class Oper:   # operating conditions and flags
    def __init__(S):
        S.Vinf = 1.          # velocity magnitude
        S.alpha = 0.         # angle of attack, in degrees
        S.rho = 1.           # density
        S.cltgt = 0.         # lift coefficient target
        S.givencl = False    # True if cl is given instead of alpha
        S.initbl = True      # True to initialize the boundary layer
        S.viscous = False    # True to do viscous
        S.redowake = False   # True to rebuild wake after alpha changes
        S.Re = 1e5           # viscous Reynolds number
        S.Ma = 0.            # Mach number
  
#-------------------------------------------------------------------------------
class Isol:   # inviscid solution variables
    def __init__(S):
        S.AIC = []                   # aero influence coeff matrix
        S.gamref = []                # 0,90-deg alpha vortex strengths at airfoil nodes
        S.gam = []                   # vortex strengths at airfoil nodes (for current alpha)
        S.sstag = 0.                 # s location of stagnation point
        S.sstag_g = np.array([0,0])  # lin of sstag w.r.t. adjacent gammas
        S.sstag_ue = np.array([0,0]) # lin of sstag w.r.t. adjacent ue values
        S.Istag = [0,0]             # node indices before/after stagnation
        S.sgnue = []                 # +/- 1 on upper/lower surface nodes
        S.xi = []                    # distance from the stagnation at all points
        S.uewi = []                  # inviscid edge velocity in wake
        S.uewiref = []               # 0,90-deg alpha inviscid ue solutions on wake
    
#-------------------------------------------------------------------------------
class Vsol:   # viscous solution variables
    def __init__(S):
        S.th  = []              # theta = momentum thickness [Nsys]
        S.ds = []               # delta star = displacement thickness [Nsys]
        S.Is = []               # 3 arrays of surface indices
        S.wgap = []             # wake gap over wake points
        S.ue_m = []             # linearization of ue w.r.t. mass (all nodes)
        S.sigma_m = []          # d(source)/d(mass) matrix
        S.ue_sigma = []         # d(ue)/d(source) matrix
        S.turb = []             # flag over nodes indicating if turbulent (1) or lam (0) 
        S.xt = 0.               # transition location (xi) on current surface under consideration
        S.Xt = np.zeros((2,2))  # transition xi/x for lower and upper surfaces 

 
#-------------------------------------------------------------------------------
class Glob:    # global parameters
    def __init__(S):
        S.Nsys = 0      # number of equations and states
        S.U = []        # primary states (th,ds,sa,ue) [4 x Nsys]
        S.dU = []       # primary state update
        S.dalpha = 0.   # angle of attack update
        S.conv = True   # converged flag
        S.R = []        # residuals [3*Nsys x 1]
        S.R_U = []      # residual Jacobian w.r.t. primary states
        S.R_x = []      # residual Jacobian w.r.t. xi (s-values) [3*Nsys x Nsys]
        S.R_V = []      # global Jacobian [4*Nsys x 4*Nsys]
        S.realloc = False # if True, system Jacobians will be re-allocated

#-------------------------------------------------------------------------------
class Post:    # post-processing outputs, distributions
    def __init__(S):
        S.cp = []       # cp distribution
        S.cpi = []      # inviscid cp distribution
        S.cl = 0        # lift coefficient
        S.cl_ue = []    # linearization of cl w.r.t. ue [N, airfoil only]
        S.cl_alpha = 0  # linearization of cl w.r.t. alpha
        S.cm = 0        # moment coefficient
        S.cdpi = 0      # near-field pressure drag coefficient
        S.cd = 0        # total drag coefficient
        S.cd_U = []     # linearization of cd w.r.t. last wake state [4]
        S.cdf = 0       # skin friction drag coefficient
        S.cdp = 0       # pressure drag coefficient
        S.rfile = None  # results output file name

        # distributions
        S.th = []       # theta = momentum thickness distribution
        S.ds = []       # delta* = displacement thickness distribution
        S.sa = []       # amplification factor/shear lag coeff distribution
        S.ue = []       # edge velocity (compressible) distribution
        S.uei = []      # inviscid edge velocity (compressible) distribution
        S.cf = []       # skin friction distribution
        S.Ret = []      # Re_theta distribution
        S.Hk = []       # kinematic shape parameter distribution

#-------------------------------------------------------------------------------
class Param:    # parameters
    def __init__(S):
        S.verb   = 1     # printing verbosity level (higher -> more verbose)
        S.rtol   = 1e-10 # residual tolerance for Newton
        S.niglob = 50    # maximum number of global iterations
        S.doplot = True  # true to plot results after each solution
        S.axplot = []    # plotting axes (for more control of where plots go)

        # viscous parameters
        S.ncrit  = 9.0   # critical amplification factor    
        S.Cuq    = 1.0   # scales the uq term in the shear lag equation
        S.Dlr    = 0.9   # wall/wake dissipation length ratio
        S.SlagK  = 5.6   # shear lag constant

        # initial Ctau after transition
        S.CtauC  = 1.8   # Ctau constant
        S.CtauE  = 3.3   # Ctau exponent

        # G-beta locus: G = GA*sqrt(1+GB*beta) + GC/(H*Rt*sqrt(cf/2))
        S.GA     = 6.7   # G-beta A constant
        S.GB     = 0.75  # G-beta B constant
        S.GC     = 18.0  # G-beta C constant

        # operating conditions and thermodynamics
        S.Minf   = 0.    # freestream Mach number
        S.Vinf   = 0.    # freestream speed
        S.muinf  = 0.    # freestream dynamic viscosity
        S.mu0    = 0.    # stagnation dynamic viscosity
        S.rho0   = 1.    # stagnation density
        S.H0     = 0.    # stagnation enthalpy
        S.Tsrat  = 0.35  # Sutherland Ts/Tref
        S.gam    = 1.4   # ratio of specific heats
        S.KTb    = 1.    # Karman-Tsien beta = sqrt(1-Minf^2)
        S.KTl    = 0.    # Karman-Tsien lambda = Minf^2/(1+KTb)^2
        S.cps    = 0.    # sonic cp

        # station information
        S.simi   = False # true at a similarity station
        S.turb   = False # true at a turbulent station
        S.wake   = False # true at a wake station


#-------------------------------------------------------------------------------
class mfoil:
    def __init__(M, coords=None, naca='0012', npanel=199):
        M.version = '2022-02-22'   # version
        M.geom = Geom()            # geometry
        M.foil = Panel()           # airfoil panels
        M.wake = Panel()           # wake panels
        M.oper = Oper()            # operating conditions
        M.isol = Isol()            # inviscid solution variables
        M.vsol = Vsol()            # viscous solution variables
        M.glob = Glob()            # global system variables
        M.post = Post()            # post-processing quantities
        M.param = Param()          # parameters
        if coords is not None:
            set_coords(M, coords)
        else: 
            naca_points(M, naca)
        make_panels(M, npanel, None)

    # set operating conditions
    def setoper(M, alpha=None, cl=None, Re=None, Ma=None, visc=None):
        if alpha is not None: M.oper.alpha = alpha
        if cl is not None: 
            M.oper.cltgt = cl
            M.oper.givencl = True
        if Re is not None:
            M.oper.Re = Re
            M.oper.viscous = True
        if Ma is not None:
            M.oper.Ma = Ma
        if visc is not None:
            if (visc is not M.oper.viscous): clear_solution(M)
            M.oper.viscous = visc

    # solve current point
    def solve(M):
        if M.oper.viscous:
            solve_viscous(M)
        else:
            solve_inviscid(M)
        if (M.param.doplot): plot_results(M)

    # geometry functions
    def geom_flap(M, xzhinge, eta):
        mgeom_flap(M, xzhinge, eta) # add a flap
    def geom_addcamber(M, zcamb):
        mgeom_addcamber(M, zcamb) # increment camber
    def geom_derotate(M):
        mgeom_derotate(M) # derotate: make chordline horizontal
    
    # derivative pinging
    def ping(M):
        ping_test(M)



# ============ INPUT, OUTPUT, UTILITY ==============

#-------------------------------------------------------------------------------
def vprint(param, verb, *args):
    if (verb <= param.verb): print(*args)

#-------------------------------------------------------------------------------
def sind(alpha):
    return np.sin(alpha*np.pi/180.)

#-------------------------------------------------------------------------------
def cosd(alpha):
    return np.cos(alpha*np.pi/180.)

#-------------------------------------------------------------------------------
def norm2(x):
    return np.linalg.norm(x,2)

#-------------------------------------------------------------------------------
def dist(a,b):
    return np.sqrt(a**2 + b**2)

#-------------------------------------------------------------------------------
def atan2(y, x):
    return np.arctan2(y,x)


# ============ PLOTTING AND POST-PROCESSING  ==============


#-------------------------------------------------------------------------------
def plot_cpplus(ax, M):
    # makes a cp plot with outputs printed
    # INPUT
    #   M : mfoil structure
    # OUTPUT
    #   cp plot on current axes
    
    chord = M.geom.chord
    x = M.foil.x[0,:].copy()
    N = M.foil.N
    xrng = np.array([-.1,1.4])*chord;
    if (M.oper.viscous):
        x = np.concatenate((x, M.wake.x[0,:]))
        colors = ['red', 'blue', 'black']
        for si in range(3):
            Is = M.vsol.Is[si]
            ax.plot(x[Is], M.post.cp[Is], '-', color=colors[si], linewidth=2)
            ax.plot(x[Is], M.post.cpi[Is], '--', color=colors[si], linewidth=2) 
    else:
        ax.plot(x, M.post.cp, '-', color='blue', linewidth=2)

    if (M.oper.Ma > 0) and (M.param.cps > (min(M.post.cp)-.2)):
        ax.plot([xrng(1), chord], M.param.cps*[1,1], '--', color='black', linewidth=2)
        ax.text(0.8*chord, M.param.cps-0.1, r'sonic $c_p$', fontsize=18)
    
    ax.set_xlim(xrng)
    ax.invert_yaxis()
    ax.set_ylabel(r'$c_p$', fontsize=18)
    ax.tick_params(labelsize=14)
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # output text box
    textstr = '\n'.join((
        r'\underline{%s}'%(M.geom.name),
        r'$\textrm{Ma} = %.4f$'%(M.oper.Ma),
        r'$\alpha = %.2f^{\circ}$'%(M.oper.alpha),
        r'$c_{\ell} = %.4f$'%(M.post.cl),
        r'$c_{m} = %.4f$'%(M.post.cm),
        r'$c_{d} = %.6f$'%(M.post.cd)))
    ax.text(0.74, 0.97, textstr, transform=ax.transAxes, fontsize=16,
            verticalalignment='top')

    if (M.oper.viscous):
        textstr = '\n'.join((
            r'$\textrm{Re} = %.1e$'%(M.oper.Re),
            r'$c_{df} = %.5f$'%(M.post.cdf),
            r'$c_{dp} = %.5f$'%(M.post.cdp)))
        ax.text(0.74, 0.05, textstr, transform=ax.transAxes, fontsize=16,
                verticalalignment='top')


#-------------------------------------------------------------------------------
def plot_airfoil(ax, M):
    # makes an airfoil plot
    # INPUT
    #   ax : axes on which to plot
    #   M : mfoil structure
    # OUTPUT
    #   airfoil plot on given axes

    chord = M.geom.chord;
    xz = M.foil.x.copy()
    if (M.oper.viscous): xz = np.hstack((xz, M.wake.x))
    xrng = np.array([-.1,1.4])*chord;
    ax.plot(xz[0,:], xz[1,:], '-', color='black', linewidth=1)
    ax.axis('equal'); ax.set_xlim(xrng); ax.axis('off')  

#-------------------------------------------------------------------------------
def mplot_boundary_layer(ax, M):
    # makes a plot of the boundary layer
    # INPUT
    #   ax : axes on which to plot
    #   M : mfoil structure
    # OUTPUT
    #   boundary layer plot on given axes
    if (not M.oper.viscous): return
    xz = np.hstack((M.foil.x, M.wake.x))
    x = xz[0,:]; N = M.foil.N
    ds = M.post.ds  # displacement thickness
    rl = 0.5*(1+(ds[0]-ds[N-1])/ds[N]); ru = 1-rl
    t = np.hstack((M.foil.t, M.wake.t))  # tangents
    n = np.vstack((-t[1,:], t[0,:])) # outward normals
    for i in range(n.shape[1]): n[:,i] /= norm2(n[:,i])
    xzd = xz + n*ds; # airfoil + delta*
    ctype = ['red', 'blue', 'black']
    for i in range(4):
        si = i;
        if (si==2): xzd = xz + n*ds*ru
        if (si==3): xzd, si = xz - n*ds*rl, 2
        Is = M.vsol.Is[si]
        ax.plot(xzd[0,Is], xzd[1,Is], '-', color=ctype[si], linewidth=2)

#-------------------------------------------------------------------------------
def plot_results(M):
    # makes a summary results plot with cp, airfoil, BL delta, outputs
    # INPUT
    #   M : mfoil structure
    # OUTPUT
    #   summary results plot as a new figure  
    
    assert M.post.cp is not None, 'no cp for results plot'
 
    # figure parameters
    plt.rcParams["figure.figsize"] = [8, 7]
    plt.rcParams["figure.autolayout"] = True
    plt.rcParams['text.usetex'] = True

    # figure
    f = plt.figure()
    ax1 = f.add_subplot(111)
    gs = gridspec.GridSpec(4, 1)
    ax1.set_position(gs[0:3].get_position(f))
    ax1.set_subplotspec(gs[0:3])
    ax2 = f.add_subplot(gs[3])
    f.tight_layout()
    plt.show(block = M.post.rfile is None)
  
    # cp plot
    plot_cpplus(ax1, M)

    # airfoil plot
    plot_airfoil(ax2, M)
    
    # # BL thickness
    mplot_boundary_layer(ax2, M)

    if (M.post.rfile is not None): plt.savefig(M.post.rfile)


#-------------------------------------------------------------------------------
def calc_force(M):
    # calculates force and moment coefficients
    # INPUT
    #   M : mfoil structure with solution (inviscid or viscous)
    # OUTPUT
    #   M.post values are filled in
    # DETAILS
    #   lift/moment are computed from a panel pressure integration
    #   the cp distribution is stored as well
  
    chord = M.geom.chord; xref = M.geom.xref  # chord and ref moment point 
    Vinf = M.param.Vinf; rho = M.oper.rho; alpha = M.oper.alpha
    qinf = 0.5*rho*Vinf**2  # dynamic pressure
    N = M.foil.N  # number of points on the airfoil
  
    # calculate the pressure coefficient at each node
    ue = M.glob.U[3,:] if M.oper.viscous else get_ueinv(M)
    cp, cp_ue = get_cp(ue, M.param); M.post.cp = cp
    M.post.cpi, cpi_ue = get_cp(get_ueinv(M), M.param)  # inviscid cp


    # lift, moment, near-field pressure cd coefficients by cp integration  
    cl, cl_ue, cl_alpha, cm, cdpi = 0, np.zeros(N), 0, 0, 0
    for i0 in range(1,N+1):
        i,ip = (0,N-1) if (i0==N) else (i0,i0-1)
        x1, x2 = M.foil.x[:,ip], M.foil.x[:,i]  # panel points
        dxv = x2-x1; dx1 = x1-xref; dx2 = x2-xref
        dx1nds = dxv[0]*dx1[0]+dxv[1]*dx1[1]  # (x1-xref) cross n*ds
        dx2nds = dxv[0]*dx2[0]+dxv[1]*dx2[1]  # (x2-xref) cross n*ds
        dx = -dxv[0]*cosd(alpha) - dxv[1]*sind(alpha)  # minus from CW node ordering
        dz =  dxv[1]*cosd(alpha) - dxv[0]*sind(alpha)  # for drag
        cp1,cp2 = cp[ip], cp[i]; cpbar = 0.5*(cp1+cp2)  # average cp on the panel
        cl = cl + dx*cpbar
        I = [ip,i]; cl_ue[I] += dx*0.5*cp_ue[I]
        cl_alpha += cpbar*(sind(alpha)*dxv[0] - cosd(alpha)*dxv[1])*np.pi/180
        cm += cp1*dx1nds/3 + cp1*dx2nds/6 + cp2*dx1nds/6 + cp2*dx2nds/3
        cdpi = cdpi + dz*cpbar
        cl /= chord; cm /= chord**2; cdpi /= chord
        M.post.cl, M.post.cl_ue, M.post.cl_alpha = cl, cl_ue, cl_alpha
        M.post.cm, M.post.cdpi = cm, cdpi
  
    # viscous contributions
    cd, cdf = 0, 0
    if (M.oper.viscous):
        # Squire-Young relation for total drag (extrapolates theta from end of wake)
        iw = M.vsol.Is[2][-1] # station at the end of the wake
        U = M.glob.U[:,iw]
        H, H_U = get_H(U)
        uk, uk_ue = get_uk(U[3], M.param)
        cd = 2.0*U[0]*(uk/Vinf)**((5+H)/2.)
        M.post.cd_U = 2.0*U[0]*(uk/Vinf)**((5+H)/2.)*np.log(uk/Vinf)*0.5*H_U
        M.post.cd_U[0] += 2.0*(uk/Vinf)**((5+H)/2.)
        M.post.cd_U[3] += 2.0*U[0]*(5+H)/2.*(uk/Vinf)**((3+H)/2.)*(1./Vinf)*uk_ue
    
        # skin friction drag
        Df = 0.
        for si in range(2):
            Is = M.vsol.Is[si]  # surface point indices
            param = build_param(M, si) # get parameter structure
            station_param(M, param, Is[0])
            cf1 = 0 # first cf value
            ue1 = 0
            rho1 = rho
            x1 = M.isol.xstag;
            for i in range(len(Is)):
                station_param(M, param, Is[i])
                cf2, cf2_U = get_cf(M.glob.U[:,Is[i]], param) # get cf value
                ue2, ue2_ue = get_uk(M.glob.U[3,Is[i]], param)
                rho2, rho2_U = get_rho(M.glob.U[:,Is[i]], param)
                x2 = M.foil.x[:,Is[i]]; dxv = x2 - x1
                dx = dxv[0]*cosd(alpha) + dxv[1]*sind(alpha)
                Df += 0.25*(rho1*cf1*ue1**2 + rho2*cf2*ue2**2)*dx
                cf1 = cf2; ue1 = ue2; x1 = x2; rho1 = rho2
        cdf = Df/(qinf*chord)

    # store results
    M.post.cd, M.post.cdf, M.post.cdp = cd, cdf, cd-cdf

    # print out current values
    s = '  alpha=%.2fdeg, cl=%.6f, cm=%.6f, cdpi=%.6f, cd=%.6f, cdf=%.6f, cdp=%.6f'%(
        M.oper.alpha, M.post.cl, M.post.cm, M.post.cdpi, M.post.cd, M.post.cdf, M.post.cdp)
    vprint(M.param, 1, s)


#-------------------------------------------------------------------------------
def get_distributions(M):
    # computes various distributions (quantities at nodes) and stores them in M.post
    # INPUT
    #   M  : mfoil class with a valid solution in M.glob.U
    # OUTPUT
    #   M.post : distribution quantities calculated
    # DETAILS
    #   Relevant for viscous solutions
  
    assert M.glob.U is not None, 'no global solution'
  
    # quantities already in the global state
    M.post.th = M.glob.U[0,:].copy()   # theta
    M.post.ds = M.glob.U[1,:].copy()   # delta*
    M.post.sa = M.glob.U[2,:].copy()   # amp or ctau
    M.post.ue, uk_ue = get_uk(M.glob.U[3,:], M.param)  # compressible edge velocity 
    M.post.uei = get_ueinv(M)  # compressible inviscid edge velocity
  
    # derived viscous quantities
    N = M.glob.Nsys; cf = np.zeros(N); Ret = np.zeros(N); Hk = np.zeros(N)
    for si in range(3):   # loop over surfaces
        Is = M.vsol.Is[si]  # surface point indices
        param = build_param(M, si)  # get parameter structure
        for i in range(len(Is)):  # loop over points
            j = Is[i]; Uj = M.glob.U[:,j]
            station_param(M, param, j)
            uk, uk_ue = get_uk(Uj[3], param) # corrected edge speed
            cfloc, cf_u = get_cf(Uj, param)  # local skin friction coefficient
            cf[j] = cfloc * uk*uk/(param.Vinf*param.Vinf)  # free-stream-based cf
            Ret[j], Ret_U = get_Ret(Uj, param)  # Re_theta
            Hk[j], Hk_U = get_Hk(Uj, param)  # kinematic shape factor

    M.post.cf, M.post.Ret, M.post.Hk = cf, Ret, Hk


# ============ INVISCID FUNCTIONS ==============


#-------------------------------------------------------------------------------
def clear_solution(M):
    # clears inviscid/viscous solutions by re-initializing structures
    # INPUT
    #   M : mfoil structure
    # OUTPUT
    #   M : mfoil structure without inviscid or viscous solution
    # DETAILS

    M.isol = Isol()
    M.vsol = Vsol()
    M.glob = Glob()
    M.post = Post()
    M.wake.N = 0;
    M.wake.x = []
    M.wake.s = []
    M.wake.t = []


#-------------------------------------------------------------------------------
def solve_inviscid(M):
    # solves the inviscid system, rebuilds 0,90deg solutions
    # INPUT
    #   M : mfoil structure
    # OUTPUT
    #   inviscid vorticity distribution is computed
    # DETAILS
    #   Uses the angle of attack in M.oper.gamma
    #   Also initializes thermo variables for normalization

    assert M.foil.N>0, 'No panels'
    M.oper.viscous = False
    init_thermo(M)
    M.isol.sgnue = np.ones(M.foil.N)  # do not distinguish sign of ue if inviscid
    build_gamma(M, M.oper.alpha)
    #if (M.oper.givencl): cltrim_inviscid(M)
    calc_force(M)
    M.glob.conv = True; # no coupled system ... convergence is guaranteed



#-------------------------------------------------------------------------------
def get_ueinv(M):
    # computes invicid tangential velocity at every node
    # INPUT
    #   M : mfoil structure
    # OUTPUT
    #   ueinv : inviscid velocity at airfoil and wake (if exists) points
    # DETAILS
    #   The airfoil velocity is computed directly from gamma
    #   The tangential velocity is measured + in the streamwise direction

    assert len(M.isol.gam) > 0, 'No inviscid solution'
    alpha = M.oper.alpha; cs = np.array([cosd(alpha), sind(alpha)])
    uea = M.isol.sgnue * np.dot(M.isol.gamref, cs)  # airfoil
    if (M.oper.viscous) and (M.wake.N > 0):
        uew = np.dot(M.isol.uewiref, cs)  # wake
        uew[0] = uea[-1]   # ensures continuity of upper surface and wake ue
    else:
        uew = np.array([])
    ueinv = np.concatenate((uea, uew))  # airfoil/wake edge velocity
    return ueinv.transpose()

#-------------------------------------------------------------------------------
def get_ueinvref(M):
    # computes 0,90deg inviscid tangential velocities at every node
    # INPUT
    #   M : mfoil structure
    # OUTPUT
    #   ueinvref : 0,90 inviscid tangential velocity at all points (N+Nw)x2
    # DETAILS
    #   Uses gamref for the airfoil, uewiref for the wake (if exists)
  
    assert len(M.isol.gam) > 0, 'No inviscid solution'  
    uearef = M.isol.gamref * M.isol.sgnue[:, np.newaxis]
    if (M.oper.viscous) and (M.wake.N > 0):
        uewref = M.isol.uewiref; # wake
        uewref[0,:] = uearef[-1,:]; # continuity of upper surface and wake
    else:
        uewref = np.array([])
    ueinvref = np.vstack((uearef, uewref))
    return ueinvref

#-------------------------------------------------------------------------------
def build_gamma(M, alpha):
    # builds and solves the inviscid linear system for alpha=0,90,input
    # INPUT
    #   M     : mfoil structure
    #   alpha : angle of attack (degrees)
    # OUTPUT
    #   M.isol.gamref : 0,90deg vorticity distributions at each node (Nx2)
    #   M.isol.gam    : gamma for the particular input angle, alpha
    #   M.isol.AIC    : aerodynamic influence coefficient matrix, filled in
    # DETAILS
    #   Uses streamfunction approach: constant psi at each node
    #   Continuous linear vorticity distribution on the airfoil panels
    #   Enforces the Kutta condition at the TE
    #   Accounts for TE gap through const source/vorticity panels
    #   Accounts for sharp TE through gamma extrapolation

    N = M.foil.N              # number of points  
    A = np.zeros([N+1,N+1])   # influence matrix
    rhs = np.zeros([N+1,2])   # right-hand sides for 0,90
    t,hTE,dtdx,tcp,tdp = TE_info(M.foil.x) # trailing-edge info
    nogap = (abs(hTE) < 1e-10*M.geom.chord) # indicates no TE gap
  
    vprint(M.param,1, '\n <<< Solving the inviscid problem >>> \n')
  
    # Build influence matrix and rhs
    for i in range(N):            # loop over nodes
        xi = M.foil.x[:,i]        # coord of node i
        for j in range(N-1):      # loop over panels
            aij, bij = panel_linvortex_stream(M.foil.x[:,[j,j+1]], xi)
            A[i,j  ] += aij
            A[i,j+1] += bij
        A[i,N] = -1; # last unknown = streamfunction value on surf

        # right-hand sides
        rhs[i,:] = [-xi[1], xi[0]]
        # TE source influence
        a = panel_constsource_stream(M.foil.x[:,[N-1,0]], xi)
        A[i,  0] += -a*(0.5*tcp)
        A[i,N-1] +=  a*(0.5*tcp)
        # TE vortex panel
        a, b = panel_linvortex_stream(M.foil.x[:,[N-1,0]], xi)
        A[i,  0] += -(a+b)*(-0.5*tdp)
        A[i,N-1] +=  (a+b)*(-0.5*tdp)
  
    # special Nth equation (extrapolation of gamma differences) if no gap
    if (nogap):
        A[N-1,:] = 0
        A[N-1,[0,1,2,N-3,N-2,N-1]] = [1,-2,1,-1,2,-1]

    # Kutta condition
    A[N,  0] = 1
    A[N,N-1] = 1

    # Solve system for unknown vortex strengths
    M.isol.AIC = A
    g = np.linalg.solve(M.isol.AIC, rhs)

    M.isol.gamref = g[0:N,:] # last value is surf streamfunction
    M.isol.gam = M.isol.gamref[:,0]*cosd(alpha) + M.isol.gamref[:,1]*sind(alpha);



#-------------------------------------------------------------------------------
def inviscid_velocity(X, G, Vinf, alpha, x, dolin):
    # returns inviscid velocity at x due to gamma (G) on panels X, and Vinf
    # INPUT
    #   X     : coordinates of N panel nodes (N-1 panels) (Nx2)
    #   G     : vector of gamma values at each airfoil node (Nx1)
    #   Vinf  : freestream speed magnitude
    #   alpha : angle of attack (degrees)
    #   x     : location of point at which velocity vector is desired  
    #   dolin : True to also return linearization
    # OUTPUT
    #   V    : velocity at the desired point (2x1)
    #   V_G  : (optional) linearization of V w.r.t. G, (2xN)
    # DETAILS
    #   Uses linear vortex panels on the airfoil
    #   Accounts for TE const source/vortex panel
    #   Includes the freestream contribution
  
    N = X.shape[1]   # number of points  
    V = np.zeros(2)   # velocity
    if (dolin): V_G = np.zeros([2,N])
    t,hTE,dtdx,tcp,tdp = TE_info(X) # trailing-edge info
    # assume x is not a midpoint of a panel (can check for this)
    for j in range(N-1):        # loop over panels
        a, b = panel_linvortex_velocity(X[:,[j,j+1]], x, None, False)
        V += a*G[j] + b*G[j+1]
        if dolin: 
            V_G[:,j] += a
            V_G[:,j+1] += b
    # TE source influence
    a = panel_constsource_velocity(X[:,[N-1,0]], x, None)
    f1, f2 = a*(-0.5*tcp), a*0.5*tcp
    V += f1*G[0] + f2*G[N-1]
    if dolin: 
        V_G[:,0] += f1
        V_G[:,N-1] += f2
    # TE vortex influence
    a, b = panel_linvortex_velocity(X[:,[N-1,0]], x, None, False)
    f1, f2 = (a+b)*(0.5*tdp), (a+b)*(-0.5*tdp)
    V += f1*G[0] + f2*G[N-1]
    if dolin:
        V_G[:,0] += f1
        V_G[:,N-1] += f2
    # freestream influence
    V += Vinf*np.array([cosd(alpha), sind(alpha)])
    if dolin:
        return V, V_G
    else:
        return V


#-------------------------------------------------------------------------------
def build_wake(M):
    # builds wake panels from the inviscid solution
    # INPUT
    #   M     : mfoil class with a valid inviscid solution (gam)
    # OUTPUT
    #   M.wake.N  : Nw, the number of wake points
    #   M.wake.x  : coordinates of the wake points (2xNw)
    #   M.wake.s  : s-values of wake points (continuation of airfoil) (1xNw)
    #   M.wake.t  : tangent vectors at wake points (2xNw)
    # DETAILS
    #   Constructs the wake path through repeated calls to inviscid_velocity
    #   Uses a predictor-corrector method
    #   Point spacing is geometric; prescribed wake length and number of points

    assert len(M.isol.gam) > 0, 'No inviscid solution'
    N = M.foil.N  # number of points on the airfoil
    Vinf = M.oper.Vinf    # freestream speed
    Nw = int(np.ceil(N/10 + 10*M.geom.wakelen)) # number of points on wake
    S = M.foil.s   # airfoil S values
    ds1 = 0.5*(S[1]-S[0] + S[N-1]-S[N-2])  # first nominal wake panel size
    sv = space_geom(ds1, M.geom.wakelen*M.geom.chord, Nw) # geometrically-spaced points
    xyw = np.zeros([2,Nw]); tw = xyw.copy()  # arrays of x,y points and tangents on wake
    xy1, xyN = M.foil.x[:,0], M.foil.x[:,N-1] # airfoil TE points
    xyte = 0.5*(xy1 + xyN)  # TE midpoint
    n = xyN-xy1; t = np.array([n[1], -n[0]]) # normal and tangent
    assert t[0] > 0, 'Wrong wake direction; ensure airfoil points are CCW'
    xyw[:,0] = xyte + 1e-5*t*M.geom.chord  # first wake point, just behind TE
    sw = S[N-1] + sv  # s-values on wake, measured as continuation of the airfoil
  
    # loop over rest of wake
    for i in range(Nw-1):
        v1 = inviscid_velocity(M.foil.x, M.isol.gam, Vinf, M.oper.alpha, xyw[:,i], False)
        v1 = v1/norm2(v1); tw[:,i] = v1; # normalized
        xyw[:,i+1] = xyw[:,i] + (sv[i+1]-sv[i])*v1; # forward Euler (predictor) step
        v2 = inviscid_velocity(M.foil.x, M.isol.gam, Vinf, M.oper.alpha, xyw[:,i+1], False);
        v2 = v2/norm2(v2); tw[:,i+1] = v2; # normalized
        xyw[:,i+1] = xyw[:,i] + (sv[i+1]-sv[i])*0.5*(v1+v2); # corrector step
        
    # determine inviscid ue in the wake, and 0,90deg ref ue too
    uewi = np.zeros([Nw,1]); uewiref = np.zeros([Nw,2])
    for i in range(Nw):
        v = inviscid_velocity(M.foil.x, M.isol.gam, Vinf, M.oper.alpha, xyw[:,i], False)
        uewi[i] = np.dot(v,tw[:,i])
        v = inviscid_velocity(M.foil.x, M.isol.gamref[:,0], Vinf, 0, xyw[:,i], False)
        uewiref[i,0] = np.dot(v, tw[:,i])
        v = inviscid_velocity(M.foil.x, M.isol.gamref[:,1], Vinf, 90, xyw[:,i], False)
        uewiref[i,1] = np.dot(v, tw[:,i])
  
    # set values
    M.wake.N = Nw
    M.wake.x = xyw
    M.wake.s = sw
    M.wake.t = tw
    M.isol.uewi = uewi
    M.isol.uewiref = uewiref

  
#-------------------------------------------------------------------------------
def stagpoint_find(M):
    # finds the LE stagnation point on the airfoil (using inviscid solution)
    # INPUTS
    #   M  : mfoil class with inviscid solution, gam
    # OUTPUTS
    #   M.isol.sstag   : scalar containing s value of stagnation point
    #   M.isol.sstag_g : linearization of sstag w.r.t gamma (1xN)
    #   M.isol.Istag   : [i,i+1] node indices before/after stagnation (1x2)
    #   M.isol.sgnue   : sign conversion from CW to tangential velocity (1xN)
    #   M.isol.xi      : distance from stagnation point at each node (1xN)

    assert len(M.isol.gam) > 0, 'No inviscid solution'
    N = M.foil.N   # number of points on the airfoil
    for j in range(N):
        if (M.isol.gam[j] > 0): break
    assert (j<N-1), 'no stagnation point'
    I = [j-1,j]
    G = M.isol.gam[I]; S = M.foil.s[I]
    M.isol.Istag = I;  # indices of neighboring gammas
    den = (G[1]-G[0]); w1 = G[1]/den; w2 = -G[0]/den
    M.isol.sstag = w1*S[0] + w2*S[1]  # s location
    M.isol.xstag = M.foil.x[:,j-1]*w1 + M.foil.x[:,j]*w2  # x location
    st_g1 = G[1]*(S[0]-S[1])/(den*den)
    M.isol.sstag_g = np.array([st_g1, -st_g1])
    sgnue = -1*np.ones(N)  # upper/lower surface sign
    for i in range(j,N): sgnue[i] = 1
    M.isol.sgnue = sgnue
    M.isol.xi = np.concatenate((abs(M.foil.s-M.isol.sstag), M.wake.s-M.isol.sstag))


#-------------------------------------------------------------------------------
def rebuild_isol(M):
    # rebuilds inviscid solution, after an angle of attack change
    # INPUT
    #   M     : mfoil class with inviscid reference solution and angle of attack
    # OUTPUT
    #   M.isol.gam : correct combination of reference gammas
    #   New stagnation point location if inviscid
    #   New wake and source influence matrix if viscous

    assert len(M.isol.gam) > 0, 'No inviscid solution'
    vprint(M.param,2, '\n  Rebuilding the inviscid solution.')
    alpha = M.oper.alpha
    M.isol.gam = M.isol.gamref[:,0]*cosd(alpha) + M.isol.gamref[:,1]*sind(alpha)
    if (not M.oper.viscous):
        stagpoint_find(M) # viscous stag point movement is handled separately
    elif (M.oper.redowake):
        build_wake(M)
        identify_surfaces(M)
        calc_ue_m(M) # rebuild matrices due to changed wake geometry


# ============ PANELING  ==============


#-------------------------------------------------------------------------------
def make_panels(M, npanel, stgt):
    # places panels on the current airfoil, as described by M.geom.xpoint
    # INPUT
    #   M      : mfoil class
    #   npanel : number of panels
    #   stgt   : optional target s values (e.g. for adaptation), or None
    # OUTPUT
    #   M.foil.N : number of panel points
    #   M.foil.x : coordinates of panel nodes (2xN)
    #   M.foil.s : arclength values at nodes (1xN)
    #   M.foil.t : tangent vectors, not normalized, dx/ds, dy/ds (2xN)
    # DETAILS
    #   Uses curvature-based point distribution on a spline of the points

    clear_solution(M) # clear any existing solution
    Ufac = 2;  # uniformity factor (higher, more uniform paneling)
    TEfac = 0.1; # Trailing-edge factor (higher, more TE resolution)
    M.foil.x, M.foil.s, M.foil.t = spline_curvature(M.geom.xpoint, npanel+1, Ufac, TEfac, stgt)
    M.foil.N = M.foil.x.shape[1]


#-------------------------------------------------------------------------------
def TE_info(X):
    # returns trailing-edge information for an airfoil with node coords X
    # INPUT
    #   X : node coordinates, ordered clockwise (2xN)
    # OUTPUT
    #   t    : bisector vector = average of upper/lower tangents, normalized
    #   hTE  : trailing edge gap, measured as a cross-section
    #   dtdx : thickness slope = d(thickness)/d(wake x)
    #   tcp  : |t cross p|, used for setting TE source panel strength
    #   tdp  : t dot p, used for setting TE vortex panel strength
    # DETAILS
    #   p refers to the unit vector along the TE panel (from lower to upper)

    t1 = X[:,  0]-X[:, 1]; t1 = t1/norm2(t1) # lower tangent vector
    t2 = X[:, -1]-X[:,-2]; t2 = t2/norm2(t2) # upper tangent vector
    t = 0.5*(t1+t2); t = t/norm2(t)  # average tangent; gap bisector
    s = X[:,-1]-X[:,0]  # lower to upper connector vector
    hTE = -s[0]*t[1] + s[1]*t[0]  # TE gap
    dtdx = t1[0]*t2[1] - t2[0]*t1[1] # sin(theta between t1,t2) approx dt/dx
    p = s/norm2(s); # unit vector along TE panel
    tcp = abs(t[0]*p[1]-t[1]*p[0]); tdp = np.dot(t,p)

    return t, hTE, dtdx, tcp, tdp


#-------------------------------------------------------------------------------
def panel_info(Xj, xi):
    # calculates common panel properties (distance, angles)
    # INPUTS
    #   Xj    : X(:,[1,2]) = panel endpoint coordinates
    #   xi    : control point coordinates (2x1)
    #   vdir  : direction of dot product, or None
    #   onmid : true means xi is on the panel midpoint
    # OUTPUTS
    #   t, n   : panel-aligned tangent and normal vectors
    #   x, z   : control point coords in panel-aligned coord system
    #   d      : panel length
    #   r1, r2 : distances from panel left/right edges to control point
    #   theta1, theta2 : left/right angles
 
    # panel coordinates
    xj1, zj1 = Xj[0,0], Xj[1,0]
    xj2, zj2 = Xj[0,1], Xj[1,1]

    # panel-aligned tangent and normal vectors
    t = np.array([xj2-xj1, zj2-zj1]); t /= norm2(t)
    n = np.array([-t[1], t[0]])

    # control point relative to (xj1,zj1)
    xz = np.array([(xi[0]-xj1), (xi[1]-zj1)])
    x = np.dot(xz,t)  # in panel-aligned coord system
    z = np.dot(xz,n)  # in panel-aligned coord system

    # distances and angles
    d = dist(xj2-xj1, zj2-zj1)   # panel length
    r1 = dist(x, z)              # left edge to control point
    r2 = dist(x-d, z)            # right edge to control point
    theta1 = atan2(z,x)          # left angle
    theta2 = atan2(z,x-d)        # right angle

    return t, n, x, z, d, r1, r2, theta1, theta2


#-------------------------------------------------------------------------------
def panel_linvortex_velocity(Xj, xi, vdir, onmid):
    # calculates the velocity coefficients for a linear vortex panel
    # INPUTS
    #   Xj    : X(:,[1,2]) = panel endpoint coordinates
    #   xi    : control point coordinates (2x1)
    #   vdir  : direction of dot product, or None
    #   onmid : true means xi is on the panel midpoint
    # OUTPUTS
    #   a,b   : velocity influence coefficients of the panel
    # DETAILS
    #   The velocity due to the panel is then a*g1 + b*g2
    #   where g1 and g2 are the vortex strengths at the panel endpoints
    #   If vdir is None, a,b are 2x1 vectors with velocity components
    #   Otherwise, a,b are dotted with vdir
    
    # panel info
    t, n, x, z, d, r1, r2, theta1, theta2 = panel_info(Xj, xi)

    # velocity in panel-aligned coord system
    if (onmid):
        ug1, ug2 = 1/2 - 1/4,  1/4
        wg1, wg2 = -1/(2*np.pi),  1/(2*np.pi)
    else:
        temp1 = (theta2-theta1)/(2*np.pi)
        temp2 = (2*z*np.log(r1/r2) - 2*x*(theta2-theta1))/(4*np.pi*d);
        ug1 =  temp1 + temp2
        ug2 =        - temp2
        temp1 = np.log(r2/r1)/(2*np.pi);
        temp2 = (x*np.log(r1/r2) - d + z*(theta2-theta1))/(2*np.pi*d);
        wg1 =  temp1 + temp2
        wg2 =        - temp2

    # velocity influence in original coord system
    a = np.array([ug1*t[0]+wg1*n[0], ug1*t[1]+wg1*n[1]]); # point 1
    b = np.array([ug2*t[0]+wg2*n[0], ug2*t[1]+wg2*n[1]]); # point 2
    if (vdir is not None): 
        a = np.dot(a, vdir)
        b = np.dot(b, vdir)

    return a, b


#-------------------------------------------------------------------------------
def panel_linvortex_stream(Xj, xi):
    # calculates the streamfunction coefficients for a linear vortex panel
    # INPUTS
    #   Xj  : X(:,[1,2]) = panel endpoint coordinates
    #   xi  : control point coordinates (2x1)
    # OUTPUTS
    #   a,b : streamfunction influence coefficients
    # DETAILS
    #   The streamfunction due to the panel is then a*g1 + b*g2
    #   where g1 and g2 are the vortex strengths at the panel endpoints

    # panel info
    t, n, x, z, d, r1, r2, theta1, theta2 = panel_info(Xj, xi)

    # check for r1, r2 zero
    ep = 1e-9
    logr1 = 0. if (r1 < ep) else np.log(r1)
    logr2 = 0. if (r2 < ep) else np.log(r2)
  
    # streamfunction components
    P1 = (0.5/np.pi)*(z*(theta2-theta1) - d + x*logr1 - (x-d)*logr2);
    P2 = x*P1 + (0.5/np.pi)*(0.5*r2**2*logr2 - 0.5*r1**2*logr1 - r2**2/4 + r1**2/4);
  
    # influence coefficients
    a = P1-P2/d;
    b =    P2/d;
  
    return a, b


#-------------------------------------------------------------------------------
def panel_constsource_velocity(Xj, xi, vdir):
    # calculates the velocity coefficient for a constant source panel
    # INPUTS
    #   Xj    : X(:,[1,2]) = panel endpoint coordinates
    #   xi    : control point coordinates (2x1)
    #   vdir  : direction of dot product, or None
    # OUTPUTS
    #   a     : velocity influence coefficient of the panel
    # DETAILS
    #   The velocity due to the panel is then a*s
    #   where s is the panel source strength
    #   If vdir is None, a,b are 2x1 vectors with velocity components
    #   Otherwise, a,b are dotted with vdir

    # panel info
    t, n, x, z, d, r1, r2, theta1, theta2 = panel_info(Xj, xi)
  
    ep = 1e-9
    logr1, theta1, theta2 = (0, np.pi, np.pi) if (r1 < ep) else (np.log(r1), theta1, theta2)
    logr2, theta1, theta2 = (0,     0,     0) if (r2 < ep) else (np.log(r2), theta1, theta2)

    # velocity in panel-aligned coord system
    u = (0.5/np.pi)*(logr1 - logr2)
    w = (0.5/np.pi)*(theta2-theta1)
  
    # velocity in original coord system dotted with given vector
    a = np.array([u*t[0]+w*n[0], u*t[1]+w*n[1]])
    if (vdir is not None): a = np.dot(a, vdir)

    return a


#-------------------------------------------------------------------------------
def panel_constsource_stream(Xj, xi):
    # calculates the streamfunction coefficient for a constant source panel
    # INPUTS
    #   Xj    : X(:,[1,2]) = panel endpoint coordinates
    #   xi    : control point coordinates (2x1)
    # OUTPUTS
    #   a     : streamfunction influence coefficient of the panel
    # DETAILS
    #   The streamfunction due to the panel is then a*s
    #   where s is the panel source strength
  
    # panel info
    t, n, x, z, d, r1, r2, theta1, theta2 = panel_info(Xj, xi)
    
    # streamfunction
    ep = 1e-9;
    logr1, theta1, theta2 = (0, np.pi, np.pi) if (r1 < ep) else (np.log(r1), theta1, theta2)
    logr2, theta1, theta2 = (0,     0,     0) if (r2 < ep) else (np.log(r2), theta1, theta2)

    P = (x*(theta1-theta2) + d*theta2 + z*logr1 - z*logr2)/(2*np.pi)
  
    dP = d; # delta psi
    P = (P - 0.25*dP) if ((theta1+theta2) > np.pi) else (P + 0.75*dP)

    # influence coefficient
    a = P
  
    return a


#-------------------------------------------------------------------------------
def panel_linsource_velocity(Xj, xi, vdir):
    # calculates the velocity coefficients for a linear source panel
    # INPUTS
    #   Xj    : X(:,[1,2]) = panel endpoint coordinates
    #   xi    : control point coordinates (2x1)
    #   vdir  : direction of dot product
    # OUTPUTS
    #   a,b   : velocity influence coefficients of the panel
    # DETAILS
    #   The velocity due to the panel is then a*s1 + b*s2
    #   where s1 and s2 are the source strengths at the panel endpoints
    #   If vdir is None, a,b are 2x1 vectors with velocity components
    #   Otherwise, a,b are dotted with vdir

    # panel info
    t, n, x, z, d, r1, r2, theta1, theta2 = panel_info(Xj, xi)

    # velocity in panel-aligned coord system
    temp1 = np.log(r1/r2)/(2*np.pi)
    temp2 = (x*np.log(r1/r2) - d + z*(theta2-theta1))/(2*np.pi*d)
    ug1 =  temp1 - temp2
    ug2 =          temp2
    temp1 = (theta2-theta1)/(2*np.pi)
    temp2 = (-z*np.log(r1/r2) + x*(theta2-theta1))/(2*np.pi*d)
    wg1 =  temp1 - temp2
    wg2 =          temp2
  
    # velocity influence in original coord system
    a = np.array([ug1*t[0]+wg1*n[0], ug1*t[1]+wg1*n[1]]) # point 1
    b = np.array([ug2*t[0]+wg2*n[0], ug2*t[1]+wg2*n[1]]) # point 2
    if (vdir is not None): a, b = np.dot(a, vdir), np.dot(b,vdir)
  
    return a, b


#-------------------------------------------------------------------------------
def panel_linsource_stream(Xj, xi):
    # calculates the streamfunction coefficients for a linear source panel
    # INPUTS
    #   Xj  : X(:,[1,2]) = panel endpoint coordinates
    #   xi  : control point coordinates (2x1)
    # OUTPUTS
    #   a,b : streamfunction influence coefficients
    # DETAILS
    #   The streamfunction due to the panel is then a*s1 + b*s2
    #   where s1 and s2 are the source strengths at the panel endpoints

    # panel info
    t, n, x, z, d, r1, r2, theta1, theta2 = panel_info(Xj, xi)

    # make branch cut at theta = 0
    if (theta1<0): theta1 = theta1 + 2*np.pi
    if (theta2<0): theta2 = theta2 + 2*np.pi
    
    # check for r1, r2 zero
    ep = 1e-9;
    logr1, theta1, theta2 = (0, np.pi, np.pi) if (r1 < ep) else (np.log(r1), theta1, theta2)
    logr2, theta1, theta2 = (0,     0,     0) if (r2 < ep) else (np.log(r2), theta1, theta2)
    
    # streamfunction components
    P1 = (0.5/np.pi)*(x*(theta1-theta2) + theta2*d + z*logr1 - z*logr2)
    P2 = x*P1 + (0.5/np.pi)*(0.5*r2**2*theta2 - 0.5*r1**2*theta1 - 0.5*z*d)
  
    # influence coefficients
    a = P1-P2/d;
    b =    P2/d;
  
    return a, b



# ============ GEOMETRY ==============

#-------------------------------------------------------------------------------
def mgeom_flap(M, xzhinge, eta):
    # deploys a flap at hinge location xzhinge, with flap angle eta 
    # INPUTS
    #   M       : mfoil class containing an airfoil
    #   xzhinge : flap hinge location (x,z) as numpy array
    #   eta     : flap angle, positive = down, degrees
    # OUTPUTS
    #   M.foil.x : modified airfoil coordinates

    X = M.geom.xpoint; N = X.shape[1] # airfoil points
    xh = xzhinge[0]  # x hinge location

    # identify points on flap
    If = np.nonzero(X[0,:]>xh)[0]

    # rotate flap points
    R = np.array([[cosd(eta), sind(eta)], [-sind(eta), cosd(eta)]])
    for i in range(len(If)):
        X[:,If[i]] = xzhinge + R@(X[:,If[i]]-xzhinge)
  
    # remove flap points to left of hinge
    I = If[X[0,If]<xh]; I = np.setdiff1d(np.arange(N),I)
  
    # re-assemble the airfoil; note, chord length is *not* redefined
    M.geom.xpoint = X[:,I]; M.geom.npoint = M.geom.xpoint.shape[1] 
    
    # repanel
    if (M.foil.N > 0): make_panels(M, M.foil.N-1,None)
  
    # clear the solution
    clear_solution(M)


#-------------------------------------------------------------------------------
def mgeom_addcamber(M, xzcamb):
    # adds camber to airfoil from given coordinates
    # INPUTS
    #   M       : mfoil class containing an airfoil
    #   xzcamb  : (x,z) points on camberline increment, 2 x Nc
    # OUTPUTS
    #   M.foil.x : modified airfoil coordinates

    if (xzcamb.shape[0] > xzcamb.shape[1]): xzcamb = np.transpose(xzcamb)  

    X = M.geom.xpoint # airfoil points

    # interpolate camber delta, add to X
    dz = interp1d(xzcamb[0,:], xzcamb[1,:], 'cubic', )(X[0,:])
    X[1,:] += dz;
    
    # store back in M.geom
    M.geom.xpoint = X; M.geom.npoint = M.geom.xpoint.shape[1] 

    # repanel
    if (M.foil.N > 0): make_panels(M, M.foil.N-1,None)
  
    # clear the solution
    clear_solution(M)


#-------------------------------------------------------------------------------
def mgeom_derotate(M):
    # derotates airfoil about leading edge to make chordline horizontal
    # INPUTS
    #   M       : mfoil class containing an airfoil
    # OUTPUTS
    #   M.foil.x : modified airfoil coordinates

    X = M.geom.xpoint; N = X.shape[1] # airfoil points
  
    xLE = X[:,np.argmin(X[0,:])] # LE point
    xTE = 0.5*(X[:,0] + X[:,N-1])  # TE "point"
  
    theta = atan2(xTE[1]-xLE[1], xTE[0]-xLE[0]) # rotation angle
    R = np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
    for i in range(N):
        X[:,i] = xLE + R@(X[:,i]-xLE)

    # store back in M.geom
    M.geom.xpoint = X; M.geom.npoint = M.geom.xpoint.shape[1] 

    # repanel
    if (M.foil.N > 0): make_panels(M, M.foil.N-1,None)
  
    # clear the solution
    clear_solution(M)


#-------------------------------------------------------------------------------
def space_geom(dx0, L, Np):
    # spaces Np points geometrically from [0,L], with dx0 as first interval
    # INPUTS
    #   dx0 : first interval length
    #   L   : total domain length
    #   Np  : number of points, including endpoints at 0,L
    # OUTPUTS
    #   x   : point locations (1xN)
    
    assert Np>1, 'Need at least two points for spacing.'
    N = Np - 1 # number of intervals
    # L = dx0 * (1 + r + r^2 + ... r^{N-1}) = dx0*(r^N-1)/(r-1)
    # let d = L/dx0, and for a guess, consider r = 1 + s
    # The equation to solve becomes d*s  = (1+s)^N - 1
    # Initial guess: (1+s)^N ~ 1 + N*s + N*(N-1)*s^2/2 + N*(N-1)*(N-2)*s^3/3
    d = L/dx0; a = N*(N-1.)*(N-2.)/6.; b = N*(N-1.)/2.; c = N-d
    disc = max(b*b-4.*a*c, 0.); r = 1 + (-b+np.sqrt(disc))/(2*a)
    for k in range(10):
        R = r**N -1-d*(r-1); R_r = N*r**(N-1)-d; dr = -R/R_r
        if (abs(dr)<1e-6): break
        r -= R/R_r
    return np.r_[0,np.cumsum(dx0*r**(np.array(range(N))))]
  

#-------------------------------------------------------------------------------
def set_coords(M, X):
    # sets geometry from coordinate matrix
    # INPUTS
    #   M : mfoil class
    #   X : matrix whose rows or columns are (x,z) points, CW or CCW
    # OUTPUTS
    #   M.geom.npoint : number of points
    #   M.geom.xpoint : point coordinates (2 x npoint)
    #   M.geom.chord  : chord length
    # DETAILS
    #   Coordinates should start and end at the trailing edge
    #   Trailing-edge point must be repeated if sharp
    #   Points can be clockwise or counter-clockwise (will detect and make CW)

    if (X.shape[0] > X.shape[1]): X = X.transpose()

    # ensure CCW
    A = 0.;
    for i in range(X.shape(1)-1): A += (X[0,i+1]-X[0,i])*(X[1,i+1]+X[1,i])
    if (A<0): X = np.fliplr(X)
  
    # store points in M
    M.geom.npoint = X.shape(1)
    M.geom.xpoint = X;
    M.geom.chord = max(X[0,:]) - min(X[0,:])


#-------------------------------------------------------------------------------
def naca_points(M, digits):
    # calculates coordinates of a NACA 4-digit airfoil, stores in M.geom
    # INPUTS
    #   M      : mfoil class
    #   digits : character array containing NACA digits
    # OUTPUTS
    #   M.geom.npoint : number of points
    #   M.geom.xpoint : point coordinates (2 x npoint)
    #   M.geom.chord  : chord length
    # DETAILS
    #   Uses analytical camber/thickness formulas

    M.geom.name = 'NACA ' + digits
    N, te = 100, 1.5  # points per side and trailing-edge bunching factor
    f = np.linspace(0,1,N+1)  # linearly-spaced points between 0 and 1
    x = 1 - (te+1)*f*(1-f)**te - (1-f)**(te+1)  # bunched points, x, 0 to 1

    # normalized thickness, gap at trailing edge (use -.1035*x**4 for no gap)
    t = .2969*np.sqrt(x) - .126*x - .3516*x**2 + .2843*x**3 - .1015*x**4;
    tmax = float(digits[-2:])*.01 # max thickness
    t = t*tmax/.2;
  
    if (len(digits)==4):
        # 4-digit series
        m, p = float(digits[0])*.01, float(digits[1])*.1
        c = m/(1-p)**2 * ((1-2.*p)+2.*p*x-x**2)
        for i in range(len(x)): 
            if x[i] < p: c[i] = m/p**2*(2*p*x[i]-x[i]**2)
    elif (len(digits)==5):
        # 5-digit series
        n = float(digits[1])
        valid = digits[0]=='2' and digits[2]=='0' and n>0 and n<6
        assert valid, '5-digit NACA must begin with 2X0, X in 1-5'
        mv = [.058, .126, .2025, .29, .391]; m = mv(n);
        cv = [361.4, 51.64, 15.957, 6.643, 3.23]; cc = cv(n);
        c = (cc/6.)*(x**3 - 3*m*x**2 + m**2*(3-m)*x);
        for i in range(len(x)): 
            if x[i] > m: c[i] = (cc/6.)*m**3*(1-x(i));
    else:
        raise ValueError('Provide 4 or 5 NACA digits')
    
    zu = c + t; zl = c - t  # upper and lower surfaces
    xs = np.concatenate((np.flip(x), x[1:])) # x points
    zs = np.concatenate((np.flip(zl), zu[1:])) # z points

    # store points in M
    M.geom.npoint = len(xs)
    M.geom.xpoint = np.vstack((xs,zs))
    M.geom.chord = max(xs) - min(xs)


#-------------------------------------------------------------------------------
def spline_curvature(Xin, N, Ufac, TEfac, stgt):
    # Splines 2D points in Xin and samples using curvature-based spacing 
    # INPUT
    #   Xin   : points to spline
    #   N     : number of points = one more than the number of panels
    #   Ufac  : uniformity factor (1 = normal; > 1 means more uniform distribution)
    #   TEfac : trailing-edge resolution factor (1 = normal; > 1 = high; < 1 = low)
    #   stgt  : optional target s values
    # OUTPUT
    #   X  : new points (2xN)
    #   S  : spline s values (N)
    #   XS : spline tangents (2xN)
  
    # min/max of given points (x-coordinate)
    xmin, xmax = min(Xin[0,:]), max(Xin[0,:])
  
    # spline given points
    PP = spline2d(Xin)

    # curvature-based spacing on geom
    nfine = 501
    s = np.linspace(0,PP['X'].x[-1],nfine)
    xyfine = splineval(PP, s)
    PPfine = spline2d(xyfine)
  
    if (stgt == None):
        s = PPfine['X'].x
        sk = np.zeros(nfine)
        xq, wq = quadseg()
        for i in range(nfine-1):
            ds = s[i+1]-s[i]
            st = xq*ds
            px = PPfine['X'].c[:,i]
            xss = 6.0*px[0]*st + 2.0*px[1]
            py = PPfine['Y'].c[:,i]
            yss = 6.0*py[0]*st + 2.0*py[1]
            skint = 0.01*Ufac+0.5*np.dot(wq, np.sqrt(xss*xss + yss*yss))*ds
      
            # force TE resolution
            xx = (0.5*(xyfine[0,i]+xyfine[0,i+1])-xmin)/(xmax-xmin) # close to 1 means at TE
            skint = skint + TEfac*0.5*np.exp(-100*(1.0-xx))

            # increment sk
            sk[i+1] = sk[i] + skint
    
        # offset by fraction of average to avoid problems with zero curvature
        sk = sk + 2.0*sum(sk)/nfine
    
        # arclength values at points
        skl = np.linspace(min(sk), max(sk), N)
        s = interp1d(sk, s, 'cubic')(skl)
    else:
        s = stgt
    
    # new points
    X, S, XS  = splineval(PPfine, s), s, splinetan(PPfine, s)

    return X, S, XS


#-------------------------------------------------------------------------------
def spline2d(X):
    # splines 2d points
    # INPUT
    #   X : points to spline (2xN)
    # OUTPUT
    #   PP : two-dimensional spline structure 
  
    N = X.shape[1]; S, Snew = np.zeros(N), np.zeros(N)

    # estimate the arclength and spline x, y separately
    for i in range(1,N): S[i] = S[i-1] + norm2(X[:,i]-X[:,i-1])
    PPX = CubicSpline(S,X[0,:])
    PPY = CubicSpline(S,X[1,:])
  
    # re-integrate to true arclength via several passes
    xq, wq = quadseg()
    for ipass in range(10):
        serr = 0
        Snew[0] = S[0]
        for i in range(N-1):
            ds = S[i+1]-S[i]
            st = xq*ds
            px = PPX.c[:,i]
            xs = 3.0*px[0]*st*st + 2.0*px[1]*st + px[2]
            py = PPY.c[:,i]
            ys = 3.0*py[0]*st*st + 2.0*py[1]*st + py[2]
            sint = np.dot(wq, np.sqrt(xs*xs + ys*ys))*ds
            serr = max(serr, abs(sint-ds))
            Snew[i+1] = Snew[i] + sint
        S[:] = Snew
        PPX = CubicSpline(S,X[0,:])
        PPY = CubicSpline(S,X[1,:])

    return {'X':PPX, 'Y':PPY}
  
  
#-------------------------------------------------------------------------------
def splineval(PP, S):
    # evaluates 2d spline at given S values
    # INPUT
    #   PP : two-dimensional spline structure 
    #   S  : arclength values at which to evaluate the spline
    # OUTPUT
    #   XY : coordinates on spline at the requested s values (2xN)

    return np.vstack((PP['X'](S), PP['Y'](S)))


#-------------------------------------------------------------------------------
def splinetan(PP, S):
    # evaluates 2d spline tangent (not normalized) at given S values
    # INPUT
    #   PP  : two-dimensional spline structure 
    #   S   : arclength values at which to evaluate the spline tangent
    # OUTPUT
    #   XYS : dX/dS and dY/dS values at each point (2xN)
    
    DPX = PP['X'].derivative()
    DPY = PP['Y'].derivative()
    return np.vstack((DPX(S), DPY(S)))
  

#-------------------------------------------------------------------------------
def quadseg():
    # Returns quadrature points and weights for a [0,1] line segment
    # INPUT
    # OUTPUT
    #   x : quadrature point coordinates (1d)
    #   w : quadrature weights
    
    x = np.array([ 0.046910077030668, 0.230765344947158, 0.500000000000000,
                   0.769234655052842, 0.953089922969332])
    w = np.array([ 0.118463442528095, 0.239314335249683, 0.284444444444444,
                   0.239314335249683, 0.118463442528095])

    return x, w



# ============ VISCOUS FUNCTIONS ==============

#-------------------------------------------------------------------------------
def calc_ue_m(M):
    # calculates sensitivity matrix of ue w.r.t. transpiration BC mass sources
    # INPUT
    #   M : mfoil class with wake already built
    # OUTPUT
    #   M.vsol.sigma_m : d(source)/d(mass) matrix, for computing source strengths
    #   M.vsol.ue_m    : d(ue)/d(mass) matrix, for computing tangential velocity
    # DETAILS
    #   "mass" flow refers to area flow (we exclude density)
    #   sigma_m and ue_m return values at each node (airfoil and wake)
    #   airfoil panel sources are constant strength
    #   wake panel sources are two-piece linear
  
    assert len(M.isol.gam) > 0, 'No inviscid solution'
    N, Nw = M.foil.N, M.wake.N  # number of points on the airfoil/wake
    assert Nw>0, 'No wake'
  
    # Cgam = d(wake uei)/d(gamma)   [Nw x N]   (not sparse)
    Cgam = np.zeros([Nw,N])
    for i in range(Nw):
        [v, v_G] = inviscid_velocity(M.foil.x, M.isol.gam, 0, 0, M.wake.x[:,i], True)
        Cgam[i,:] = v_G[0,:]*M.wake.t[0,i] + v_G[1,:]*M.wake.t[1,i]
  
    # B = d(airfoil surf streamfunction)/d(source)  [(N+1) x (N+Nw-2)]  (not sparse)
    B = np.zeros([N+1,N+Nw-2])  # note, N+Nw-2 = # of panels

    for i in range(N):  # loop over points on the airfoil
        xi = M.foil.x[:,i] # coord of point i
        for j in range(N-1): # loop over airfoil panels
            B[i,j] = panel_constsource_stream(M.foil.x[:,[j,j+1]], xi)
        for j in range(Nw-1): # loop over wake panels
            Xj = M.wake.x[:,[j,j+1]] # panel endpoint coordinates
            Xm = 0.5*(Xj[:,0] + Xj[:,1]) # panel midpoint
            Xj = np.transpose(np.vstack((Xj[:,0], Xm, Xj[:,1]))) # left, mid, right coords on panel
            if (j==(Nw-2)): Xj[:,2] = 2*Xj[:,2] - Xj[:,1]  # ghost extension at last point
            a, b = panel_linsource_stream(Xj[:,[0,1]], xi)  # left half panel
            if (j > 0):
                B[i,N-1+j] += 0.5*a + b
                B[i,N-1+j-1] += 0.5*a
            else:
                B[i,N-1+j] += b
            a, b = panel_linsource_stream(Xj[:,[1,2]], xi) # right half panel
            B[i,N-1+j] += a + 0.5*b
            if (j<Nw-2):
                B[i,N-1+j+1] += 0.5*b
            else:
                B[i,N-1+j] += 0.5*b

    # Bp = - inv(AIC) * B   [N x (N+Nw-2)]  (not sparse)
    # Note, Bp is d(airfoil gamma)/d(source)
    Bp = -np.linalg.solve(M.isol.AIC, B)  # this has N+1 rows, but the last one is zero
    Bp = Bp[:-1,:]  # trim the last row

    # Csig = d(wake uei)/d(source) [Nw x (N+Nw-2)]  (not sparse)
    Csig = np.zeros([Nw, N+Nw-2])
    for i in range(Nw):
        xi, ti = M.wake.x[:,i], M.wake.t[:,i]  # point, tangent on wake
    
        # first/last airfoil panel effects on i=0 wake point handled separately
        jstart, jend = 0 + (i==0), N-1 - (i==0)
        for j in range(jstart, jend):   # constant sources on airfoil panels
            Csig[i,j] = panel_constsource_velocity(M.foil.x[:,[j,j+1]], xi, ti)
    
        # piecewise linear sources across wake panel halves (else singular)
        for j in range(Nw):  # loop over wake points
            I = [max(j-1,0), j, min(j+1,Nw-1)]  # left, self, right

            Xj = M.wake.x[:,I]  # point coordinates
            Xj[:,0] = 0.5*(Xj[:,0] + Xj[:,1])  # left midpoint
            Xj[:,2] = 0.5*(Xj[:,1] + Xj[:,2])  # right midpoint

            if (j==Nw-1): Xj[:,2] = 2*Xj[:,1] - Xj[:,0]  # ghost extension at last point
            d1 = norm2(Xj[:,1]-Xj[:,0]) # left half-panel length
            d2 = norm2(Xj[:,2]-Xj[:,1]) # right half-panel length
            if (i==j):
                if (j==0): # first point: special TE system (three panels meet)
                    dl = norm2(M.foil.x[:,  1]-M.foil.x[:,  0]) # lower surface panel length
                    du = norm2(M.foil.x[:,N-1]-M.foil.x[:,N-2]) # upper surface panel length
                    Csig[i,  0] += (0.5/np.pi)*(np.log(dl/d2) + 1) # lower panel effect
                    Csig[i,N-2] += (0.5/np.pi)*(np.log(du/d2) + 1) # upper panel effect
                    Csig[i,N-1] += - 0.5/np.pi # self effect
                elif (j==Nw-1): # last point: no self effect of last pan (ghost extension)
                    Csig[i,N-1+j-1] += 0 # hence the 0
                else: # all other points
                    aa = (0.25/np.pi)*np.log(d1/d2)
                    Csig[i,N-1+j-1] += aa + 0.5/np.pi
                    Csig[i,N-1+j  ] += aa - 0.5/np.pi
            else:
                if (j==0): # first point only has a half panel on the right
                    a, b = panel_linsource_velocity(Xj[:,[1,2]], xi, ti)
                    Csig[i,N-1] += b  # right half panel effect
                    Csig[i,  0] += a  # lower airfoil panel effect
                    Csig[i,N-2] += a  # upper airfoil panel effect
                elif (j==Nw-1): # last point has a constant source ghost extension
                    a = panel_constsource_velocity(Xj[:,[0,2]], xi, ti)
                    Csig[i,N+Nw-3] += a  # full const source panel effect
                else: # all other points have a half panel on left and right
                    [a1,b1] = panel_linsource_velocity(Xj[:,[0,1]], xi, ti) # left half-panel ue contrib
                    [a2,b2] = panel_linsource_velocity(Xj[:,[1,2]], xi, ti) # right half-panel ue contrib
                    Csig[i,N-1+j-1] += a1 + 0.5*b1
                    Csig[i,N-1+j  ] += 0.5*a2 + b2
  
    # compute ue_sigma = d(unsigned ue)/d(source) [(N+Nw) x (N+Nw-2)] (not sparse)
    # Df = +/- Bp = d(foil uei)/d(source)  [N x (N+Nw-2)]  (not sparse)
    # Dw = (Cgam*Bp + Csig) = d(wake uei)/d(source)  [Nw x (N+Nw-2)]  (not sparse)
    Dw = np.dot(Cgam,Bp) + Csig
    Dw[0,:] = Bp[-1,:]  # ensure first wake point has same ue as TE
    M.vsol.ue_sigma = np.vstack((Bp, Dw))  # store combined matrix
  
    # build ue_m from ue_sigma, using sgnue
    rebuild_ue_m(M)

#-------------------------------------------------------------------------------
def rebuild_ue_m(M):
    # rebuilds ue_m matrix after stagnation panel change (new sgnue)
    # INPUT
    #   M : mfoil class with calc_ue_m already called once
    # OUTPUT
    #   M.vsol.sigma_m : d(source)/d(mass) matrix, for computing source strengths
    #   M.vsol.ue_m    : d(ue)/d(mass) matrix, for computing tangential velocity
    # DETAILS
    #   "mass" flow refers to area flow (we exclude density)
    #   sigma_m and ue_m return values at each node (airfoil and wake)
    #   airfoil panel sources are constant strength
    #   wake panel sources are two-piece linear

    assert len(M.vsol.ue_sigma) > 0, 'Need ue_sigma to build ue_m'
  
    # Dp = d(source)/d(mass)  [(N+Nw-2) x (N+Nw)]  (sparse)
    N, Nw = M.foil.N, M.wake.N  # number of points on the airfoil/wake
    if (type(M.vsol.sigma_m) == list) or not (M.vsol.sigma_m.shape == (N+Nw-2, N+Nw)):
        alloc_Dp = True
        M.vsol.sigma_m = sparse.lil_matrix((N+Nw-2,N+Nw))  # empty matrix
    else:
        alloc_Dp = False
        M.vsol.sigma_m *= 0.
    for i in range(N-1):
        ds = M.foil.s[i+1]-M.foil.s[i]
        # Note, at stagnation: ue = K*s, dstar = const, m = K*s*dstar
        # sigma = dm/ds = K*dstar = m/s (separate for each side, +/-)
        M.vsol.sigma_m[i,[i,i+1]] = M.isol.sgnue[[i,i+1]]*np.array([-1.,1.])/ds
    for i in range(Nw-1):
        ds = M.wake.s[i+1]-M.wake.s[i]
        M.vsol.sigma_m[N-1+i,[N+i,N+i+1]] = np.array([-1.,1.])/ds
    if (alloc_Dp):
        M.vsol.sigma_m = M.vsol.sigma_m.tocsr()
  
    # sign of ue at all points (wake too)
    sgue = np.concatenate((M.isol.sgnue, np.ones(Nw)))
  
    # ue_m = ue_sigma * sigma_m [(N+Nw) x (N+Nw)] (not sparse)
    M.vsol.ue_m = sparse.spdiags(sgue,0,N+Nw,N+Nw,'csr') @ M.vsol.ue_sigma @ M.vsol.sigma_m
  
#-------------------------------------------------------------------------------
def init_thermo(M):
    # initializes thermodynamics variables in param structure
    # INPUT
    #   M  : mfoil class with oper structure set
    # OUTPUT
    #   M.param fields filled in based on M.oper
    #   Gets ready for compressibilty corrections if M.oper.Ma > 0

    g = M.param.gam; gmi = g-1
    rhoinf = M.oper.rho  # freestream density
    Vinf = M.oper.Vinf; M.param.Vinf = Vinf # freestream speed
    M.param.muinf = rhoinf*Vinf*M.geom.chord/M.oper.Re # freestream dyn viscosity 
    Minf = M.oper.Ma; M.param.Minf = Minf # freestream Mach
    if (Minf > 0):
        M.param.KTb = np.sqrt(1-Minf**2)  # Karman-Tsien beta
        M.param.KTl = Minf**2/(1+M.param.KTb)**2 # Karman-Tsien lambda
        M.param.H0 = (1+0.5*gmi*Minf**2)*Vinf**2/(gmi*Minf**2) # stagnation enthalpy
        Tr = 1-0.5*Vinf**2/M.param.H0  # freestream/stagnation temperature ratio
        finf = Tr**1.5*(1+M.param.Tsrat)/(Tr + M.param.Tsrat)  # Sutherland's ratio
        M.param.cps = 2/(g*Minf**2)*(((1+0.5*gmi*Minf**2)/(1+0.5*gmi))**(g/gmi) - 1)
    else:
        finf = 1  # incompressible case

    M.param.mu0 = M.param.muinf/finf   # stag visc (Sutherland ref temp is stag)
    M.param.rho0 = rhoinf*(1+0.5*gmi*Minf**2)**(1/gmi)  # stag density

#-------------------------------------------------------------------------------
def identify_surfaces(M):
    # identifies lower/upper/wake surfaces
    # INPUT
    #   M  : mfoil class with stagnation point found
    # OUTPUT
    #   M.vsol.Is : list of of node indices for lower(1), upper(2), wake(3)
    
    M.vsol.Is = [range(M.isol.Istag[0],-1,-1), 
                 range(M.isol.Istag[1],M.foil.N),
                 range(M.foil.N, M.foil.N+M.wake.N)]


#-------------------------------------------------------------------------------
def set_wake_gap(M):
    # sets height (delta*) of dead air in wake
    # INPUT
    #   M  : mfoil class with wake built and stagnation point found
    # OUTPUT
    #   M.vsol.wgap : wake gap at each wake point
    # DETAILS
    #   Uses cubic function to extrapolate the TE gap into the wake
    #   See Drela, IBL for Blunt Trailing Edges, 1989, 89-2166-CP
  
    t,hTE,dtdx,tcp,tdp = TE_info(M.foil.x) # trailing-edge info
    flen = 2.5 # length-scale factor
    dtdx = min(max(dtdx,-3./flen), 3./flen)  # clip TE thickness slope
    Lw = flen*hTE
    wgap = np.zeros(M.wake.N)
    for i in range(M.wake.N):
        xib = (M.isol.xi[M.foil.N+i] - M.isol.xi[M.foil.N])/Lw
        if (xib <= 1): wgap[i] = hTE*(1+(2+flen*dtdx)*xib)*(1-xib)**2
    M.vsol.wgap = wgap


#-------------------------------------------------------------------------------
def stagpoint_move(M):
    # moves the LE stagnation point on the airfoil using the global solution ue
    # INPUT
    #   M  : mfoil class with a valid solution in M.glob.U
    # OUTPUT
    #   New sstag, sstag_ue, xi in M.isol
    #   Possibly new stagnation panel, Istag, and hence new surfaces and matrices
  
    N = M.foil.N   # number of points on the airfoil
    I = M.isol.Istag  # current adjacent node indices
    ue = M.glob.U[3,:]  # edge velocity
    sstag0 = M.isol.sstag  # original stag point location
  
    newpanel = True;  # are we moving to a new panel?
    if (ue[I[1]] < 0):
        # move stagnation point up (larger s, new panel)
        vprint(M.param,2, '  Moving stagnation point up')
        for j in range(I[1], N):
            if (ue[j] > 0): break
        assert (j<N), 'no stagnation point'
        I1 = j
        for j in range(I[1], I1): ue[j] *= -1.
        I[0], I[1] = I1-1, I1  # new panel
    elif (ue[I[0]] < 0):
        # move stagnation point down (smaller s, new panel)
        vprint(M.param,2, '  Moving stagnation point down')
        for j in range(I[0], -1, -1):
            if (ue[j] > 0): break
        assert j>0, 'no stagnation point'
        I0 = j
        for j in range(I0+1,I[0]+1): ue[j] *= -1.
        I[0], I[1] = I0, I0+1  # new panel
    else: newpanel = False; # staying on the current panel

    # move point along panel
    ues, S = ue[I], M.foil.s[I]
    assert (ues[0] > 0) and (ues[1] > 0), 'stagpoint_move: velocity error'
    den = ues[0] + ues[1]; w1 = ues[1]/den; w2 = ues[0]/den
    M.isol.sstag = w1*S[0] + w2*S[1]  # s location
    M.isol.xstag = np.dot(M.foil.x[:,I], np.r_[w1,w2]) # x location
    M.isol.sstag_ue = np.r_[ues[1], -ues[0]]*(S[1]-S[0])/(den*den)
    vprint(M.param,2, '  Moving stagnation point: s=%.15e -> s=%.15e'%(sstag0, M.isol.sstag))
  
    # set new xi coordinates for every point
    M.isol.xi = np.concatenate((abs(M.foil.s-M.isol.sstag), M.wake.s-M.isol.sstag))
  
    # matrices need to be recalculated if on a new panel
    if (newpanel):
        vprint(M.param,2, '  New stagnation panel = %d %d'%(I[0], I[1]))
        M.isol.Istag = I  # new panel indices
        for i in range(I[0]+1): M.isol.sgnue[i] = -1
        for i in range(I[0]+1,N): M.isol.sgnue[i] = 1
        identify_surfaces(M)  # re-identify surfaces
        M.glob.U[3,:] = ue  # sign of ue changed on some points near stag
        M.glob.realloc = True
        rebuild_ue_m(M)


#-------------------------------------------------------------------------------
def solve_viscous(M):
    # solves the viscous system (BL + outer flow concurrently)
    # INPUT
    #   M  : mfoil class with an airfoil
    # OUTPUT
    #   M.glob.U : global solution
    #   M.post   : post-processed quantities

    solve_inviscid(M)
    M.oper.viscous = True
    init_thermo(M)
    build_wake(M)
    stagpoint_find(M)  # from the inviscid solution
    identify_surfaces(M)
    set_wake_gap(M)  # blunt TE dead air extent in wake
    calc_ue_m(M)
    init_boundary_layer(M) # initialize boundary layer from ue
    stagpoint_move(M) # move stag point, using viscous solution  
    solve_coupled(M) # solve coupled system
    calc_force(M)
    get_distributions(M)



#-------------------------------------------------------------------------------
def solve_coupled(M):
    # Solves the coupled inviscid and viscous system
    # INPUT
    #   M  : mfoil class with an inviscid solution
    # OUTPUT
    #   M.glob.U : global coupled solution
    # DETAILS
    #   Inviscid solution should exist, and BL variables should be initialized
    #   The global variables are [th, ds, sa, ue] at every node
    #   th = momentum thickness; ds = displacement thickness
    #   sa = amplification factor or sqrt(ctau); ue = edge velocity
    #   Nsys = N + Nw = total number of unknowns
    #   ue is treated as a separate variable for improved solver robustness
    #   The alternative is to eliminate ue, ds and use mass flow (not done here):
    #     Starting point: ue = uinv + D*m -> ue_m = D
    #     Since m = ue*ds, we have ds = m/ue = m/(uinv + D*m)
    #     So, ds_m = diag(1/ue) - diag(ds/ue)*D
    #     The residual linearization is then: R_m = R_ue*ue_m + R_ds*ds_m
  
    # Newton loop
    nNewton = M.param.niglob  # number of iterations
    M.glob.conv = False
    M.glob.realloc = True  # reallocate Jacobian on first iter
    vprint(M.param,1, '\n <<< Beginning coupled solver iterations >>>')

    for iNewton in range(nNewton):

        # set up the global system
        vprint(M.param, 2, 'Building global system')
        build_glob_sys(M)
    
        # compute forces
        vprint(M.param, 2, 'Calculating force')
        calc_force(M)
    
        # convergence check
        Rnorm = norm2(M.glob.R)
        vprint(M.param,1, '\nNewton iteration %d, Rnorm = %.10e'%(iNewton, Rnorm))
        if (Rnorm < M.param.rtol): 
            M.glob.conv = True
            break
    
        # solve global system
        vprint(M.param, 2, 'Solving global system')
        solve_glob(M)
    
        # update the state
        vprint(M.param, 2, 'Updating the state')
        update_state(M)

        M.glob.realloc = False  # assume Jacobian will not get reallocated
        
        # update stagnation point; Newton still OK; had R_x effects in R_U
        vprint(M.param, 2, 'Moving stagnation point')
        stagpoint_move(M)
    
        # update transition
        vprint(M.param, 2, 'Updating transition')
        update_transition(M)
                
  
    if (not M.glob.conv): vprint(M.param,1, '\n** Global Newton NOT CONVERGED **\n')


#-------------------------------------------------------------------------------
def update_state(M):
    # updates state, taking into account physical constraints
    # INPUT
    #   M  : mfoil class with a valid solution (U) and proposed update (dU)
    # OUTPUT
    #   M.glob.U : updated solution, possibly with a fraction of dU added
    # DETAILS
    #   U = U + omega * dU; omega = under-relaxation factor
    #   Calculates omega to prevent big changes in the state or negative values
    
    if (any(np.iscomplex(M.glob.U[2,:]))): raise ValueError('imaginary amp in U')
    if (any(np.iscomplex(M.glob.dU[2,:]))): raise ValueError('imaginary amp in dU')

    # max ctau
    It = np.nonzero(M.vsol.turb==True)[0]
    ctmax = max(M.glob.U[2,It])
  
    # starting under-relaxation factor
    omega = 1.0
      
    # first limit theta and delta*
    for k in range(2):
        Uk = M.glob.U[k,:]; dUk = M.glob.dU[k,:]
        # prevent big decreases in th, ds
        fmin = min(dUk/Uk)  # find most negative ratio
        om = abs(0.5/fmin) if (fmin < -0.5) else 1.
        if (om<omega): 
            omega = om
            vprint(M.param,3, '  th/ds decrease: omega = %.5f'%(omega))
  
    # limit negative amp/ctau
    Uk = M.glob.U[2,:]; dUk = M.glob.dU[2,:]
    for i in range(len(Uk)):
        if (not M.vsol.turb[i]) and (Uk[i]<.2): continue  # do not limit very small amp (too restrictive)
        if (M.vsol.turb[i]) and (Uk[i]<0.1*ctmax): continue # do not limit small ctau
        if (Uk[i]==0.) or (dUk[i]==0.): continue
        if (Uk[i]+dUk[i] < 0):
            om = 0.8*abs(Uk[i]/dUk[i]); 
            if (om<omega): 
                omega = om
                vprint(M.param,3, '  neg sa: omega = %.5f'%(omega))
  
    # prevent big changes in amp
    I = np.nonzero(M.vsol.turb==False)[0]
    if (any(np.iscomplex(Uk[I]))): raise ValueError('imaginary amplification')
    dumax = max(abs(dUk[I]))
    om = abs(2./dumax) if (dumax>0) else 1.
    if (om<omega):
        omega = om
        vprint(M.param,3, '  amp: omega = %.5f'%(omega))
    
    # prevent big changes in ctau
    I = np.nonzero(M.vsol.turb==True)[0]
    dumax = max(abs(dUk[I]))
    om = abs(.05/dumax) if (dumax>0) else 1.
    if (om<omega):
        omega = om
        vprint(M.param,3, '  ctau: omega = %.5f'%(omega))
  
    # prevent large ue changes
    dUk = M.glob.dU[3,:]
    fmax = max(abs(dUk)/M.oper.Vinf)
    om = 0.2/fmax if (fmax > 0) else 1.
    if (om<omega):
        omega = om
        vprint(M.param,3, '  ue: omega = %.5f'%(omega))
  
    # prevent large alpha changes
    if (abs(M.glob.dalpha) > 2):
        omega = min(omega, abs(2/M.glob.dalpha))

    # take the update
    vprint(M.param,2, '  state update: under-relaxation = %.5f'%(omega))
    M.glob.U += omega*M.glob.dU
    M.oper.alpha += omega*M.glob.dalpha
  
    # fix bad Hk after the update
    for si in range(3): # loop over surfaces
        Hkmin = 1.00005 if (si==2) else 1.02
        Is = M.vsol.Is[si]  # surface point indices
        param = build_param(M, si); # get parameter structure
        for i in range(len(Is)): # loop over points
            j = Is[i]
            Uj = M.glob.U[:,j]
            station_param(M, param, j)
            Hk, Hk_U = get_Hk(Uj, param)
            if (Hk < Hkmin):
                M.glob.U[1,j] += 2*(Hkmin-Hk)*M.glob.U[1,j]
  
    # fix negative ctau after the update
    for ii in range(len(It)):
        i = It[ii]
        if (M.glob.U[2,i] < 0): M.glob.U[2,i] = 0.1*ctmax
  
    # rebuild inviscid solution (gam, wake) if angle of attack changed
    if (abs(omega*M.glob.dalpha) > 1e-10): rebuild_isol(M)
    
#-------------------------------------------------------------------------------
def solve_glob(M):
    # solves global system for the primary variable update dU
    # INPUT
    #   M  : mfoil class with residual and Jacobian calculated
    # OUTPUT
    #   M.glob.dU : proposed solution update
    # DETAILS
    #   Uses the augmented system: fourth residual = ue equation
    #   Supports lift-constrained mode, with an extra equation: cl - cltgt = 0
    #   Extra variable in cl-constrained mode is angle of attack
    #   Solves sparse matrix system for for state/alpha update  
  
    Nsys = M.glob.Nsys  # number of dofs
    docl = M.oper.givencl # 1 if in cl-constrained mode
  
    # get edge velocity and displacement thickness
    ue = M.glob.U[3,:]; ds = M.glob.U[1,:]
    uemax = max(abs(ue))
    for i in range(len(ue)): ue[i] = max(ue[i], 1e-10*uemax) # avoid 0/negative ue
    
    # use augmented system: variables = th, ds, sa, ue
    
    # inviscid edge velocity on the airfoil and wake
    ueinv = get_ueinv(M)
  
    # initialize the global variable Jacobian
    NN = 4*Nsys + docl
    if (M.glob.realloc) or (type(M.glob.R_V) == list) or not (M.glob.R_V.shape == (NN, NN)):
        alloc_R_V = True
        M.glob.R_V = sparse.lil_matrix((NN,NN)) # +1 for cl-alpha constraint
    else:
        alloc_R_V = False
        M.glob.R_V *= 0. # matrix already allocated, just zero it out
  
    # state indices in the global system
    Ids = slice(1,4*Nsys,4) # delta star indices
    Iue = slice(3,4*Nsys,4) # ue indices
  
    # assemble the residual
    R = np.concatenate((M.glob.R, ue - (ueinv + M.vsol.ue_m @ (ds*ue))))
    #print('first Norm(R) = %.10e'%(norm2(R)))
  
    # assemble the Jacobian
    M.glob.R_V[0:3*Nsys,0:4*Nsys] = M.glob.R_U
    I = slice(3*Nsys,4*Nsys,1)
    M.glob.R_V[I,Iue] = sparse.identity(Nsys) - M.vsol.ue_m @ np.diag(ds)
    M.glob.R_V[I,Ids] = -M.vsol.ue_m @ np.diag(ue)
  
    if (docl):
        # include cl-alpha residual and Jacobian
        Rcla, Ru_alpha, Rcla_U = clalpha_residual(M)
        R = np.append(R, Rcla) #np.concatenate((R, Rcla))
        M.glob.R_V[I,4*Nsys] = Ru_alpha
        M.glob.R_V[4*Nsys,:] = Rcla_U
  
    # solve system for dU, dalpha
    if (alloc_R_V): M.glob.R_V = M.glob.R_V.tocsr()
    dV = -sparse.linalg.spsolve(M.glob.R_V, R)
  
    # store dU, reshaped, in M
    M.glob.dU = np.reshape(dV[0:4*Nsys],(4,Nsys),order='F')
    if (docl): M.glob.dalpha = dV[-1]

   
#-------------------------------------------------------------------------------
def clalpha_residual(M):
    # computes cl constraint (or just prescribed alpha) residual and Jacobian
    # INPUT
    #   M  : mfoil class with inviscid solution and post-processed cl_alpha, cl_ue
    # OUTPUT
    #   Rcla     : cl constraint residual = cl - cltgt (scalar)
    #   Ru_alpha : lin of ue residual w.r.t. alpha (Nsys x 1)
    #   Rcla_U   : lin of cl residual w.r.t state (1 x 4*Nsys)
    # DETAILS
    #   Used for cl-constrained mode, with alpha as the extra variable
    #   Should be called with up-to-date cl and cl linearizations
  
    Nsys = M.glob.Nsys    # number of dofs
    N = M.foil.N          # number of points (dofs) on airfoil
    alpha = M.oper.alpha  # angle of attack (deg)
  
    if (M.oper.givencl):  # cl is prescribed, need to trim alpha    
        Rcla = M.post.cl - M.oper.cltgt  # cl constraint residual
        Rcla_U = np.zeros(4*Nsys+1); Rcla_U[-1] = M.post.cl_alpha
        Rcla_U[3:4*N:4] = M.post.cl_ue  # only airfoil nodes affected 
        # Ru = ue - [uinv + ue_m*(ds.*ue)], uinv = uinvref*[cos(alpha);sin(alpha)]
        Ru_alpha = -get_ueinvref(M) @ np.r_[-sind(alpha), cosd(alpha)]*np.pi/180
    else:     # alpha is prescribed, easy
        Rcla = 0  # no residual
        Ru_alpha = np.zeros(Nsys,1)  # not really, but alpha is not changing
        Rcla_U = np.zeros(4*Nsys+1); Rcla_U[-1] = 1
        
    return Rcla, Ru_alpha, Rcla_U
  

#-------------------------------------------------------------------------------
def build_glob_sys(M):
    # builds the primary variable global residual system for the coupled problem
    # INPUT
    #   M  : mfoil class with a valid solution in M.glob.U
    # OUTPUT
    #   M.glob.R   : global residual vector (3*Nsys x 1)
    #   M.glob.R_U : residual Jacobian matrix (3*Nsys x 4*Nsys, sparse)
    #   M.glob.R_x : residual linearization w.r.t. x (3*Nsys x Nsys, sparse)
    # DETAILS
    #   Loops over nodes/stations to assemble residual and Jacobian
    #   Transition dicated by M.vsol.turb, which should be consistent with the state
    #   Accounts for wake initialization and first-point similarity solutions
    #   Also handles stagnation point on node via simple extrapolation
  
    Nsys = M.glob.Nsys

    # allocate matrices if [], if size changed, or if global realloc flag is true
    if (M.glob.realloc) or (type(M.glob.R) == list) or not (M.glob.R.shape[0] == 3*Nsys):
        M.glob.R = np.zeros(3*Nsys)
    else:
        M.glob.R *= 0.
    if (M.glob.realloc) or (type(M.glob.R_U) == list) or not(M.glob.R_U.shape == (3*Nsys,4*Nsys)):
        alloc_R_U = True
        M.glob.R_U = sparse.lil_matrix((3*Nsys,4*Nsys))
    else:
        alloc_R_U = False
        M.glob.R_U *= 0.

    if (M.glob.realloc) or (type(M.glob.R_x) == list) or not (M.glob.R_x == (3*Nsys, Nsys)):
        alloc_R_x = True
        M.glob.R_x = sparse.lil_matrix((3*Nsys,Nsys))
    else:
        alloc_R_x = False
        M.glob.R_x *= 0.
  
    for si in range(3):     # loop over surfaces
        Is = M.vsol.Is[si]  # surface point indices
        xi = M.isol.xi[Is]  # distance from LE stag point
        N = len(Is)         # number of points on this surface
        U = M.glob.U[:,Is]  # [th, ds, sa, ue] states at all points on this surface
        Aux = np.zeros(N)   # auxiliary data at all points: [wgap]
    
        # get parameter structure
        param = build_param(M, si);
    
        # set auxiliary data
        if (si == 2): Aux[:] = M.vsol.wgap
    
        # special case of tiny first xi -- will set to stagnation state later
        i0 = 1 if (si < 2) and (xi[0] < 1e-8*xi[-1]) else 0
    
        # first point system 
        if (si < 2): 
      
            # calculate the stagnation state, a function of U1 and U2
            Ip = [i0,i0+1]
            Ust, Ust_U, Ust_x, xst = stagnation_state(U[:,Ip], xi[Ip]) # stag state
            param.turb, param.simi = False, True  # similarity station flag  
            R1, R1_Ut, R1_x = residual_station(param, np.r_[xst,xst], np.stack((Ust, Ust),axis=-1), Aux[[i0,i0]])
            param.simi = False
            R1_Ust = R1_Ut[:,0:4] + R1_Ut[:,4:8]
            R1_U = np.dot(R1_Ust, Ust_U)
            R1_x = np.dot(R1_Ust, Ust_x)
            J = [Is[i0], Is[i0+1]]
            
            if (i0 == 1): 
                # i0=0 point landed right on stagnation: set value to Ust
                vprint(param, 2, 'hit stagnation!')
                Ig = slice(3*Is[0], 3*Is[0]+3)
                M.glob.R[Ig] = U[0:3,0] - Ust[0:3]
                M.glob.R_U[Ig, 4*Is[0]:(4*Is[0]+4)] += np.eye(3,4)
                M.glob.R_U[Ig, 4*J[0]:(4*J[0]+4)] -= Ust_U[0:3,0:4]
                M.glob.R_U[Ig, 4*J[1]:(4*J[1]+4)] -= Ust_U[0:3,4:8]
                M.glob.R_x[Ig, J] = -Ust_x[0:3,:]
      
        else:
            # wake initialization
            R1, R1_U, J = wake_sys(M, param)
            R1_x = []  # no xi dependence of first wake residual
            param.turb, param.wake = True, True  # force turbulent in wake if still laminar

        # store first point system in global residual, Jacobian
        Ig = slice(3*Is[i0], 3*Is[i0]+3)
        M.glob.R[Ig] = R1
        if (alloc_R_U): R1_U += 1e-15  # hack: force lil sparse format to allocate
        if (alloc_R_x) and (len(R1_x)>0): R1_x += 1e-15  # hack: force lil sparse format to allocate
        for j in range(len(J)):
            M.glob.R_U[Ig, 4*J[j]:(4*J[j]+4)] += R1_U[:,4*j:(4*j+4)]
            if (len(R1_x)>0): M.glob.R_x[Ig,J[j]] += R1_x[:,j:(j+1)]

        # march over rest of points
        for i in range(i0+1,N):
            Ip = [i-1,i]  # two points involved in the calculation
      
            tran = M.vsol.turb[Is[i-1]] ^ M.vsol.turb[Is[i]]  # transition flag
      
            # residual, Jacobian for point i
            if (tran):
                Ri, Ri_U, Ri_x = residual_transition(M, param, xi[Ip], U[:,Ip], Aux[Ip])
                store_transition(M, si, i)
            else:
                Ri, Ri_U, Ri_x = residual_station(param, xi[Ip], U[:,Ip], Aux[Ip])
        
            # store point i contribution in global residual, Jacobian
            Ig = slice(3*Is[i], 3*Is[i]+3)
            if (alloc_R_U): Ri_U += 1e-15  # hack: force lil sparse format to allocate
            if (alloc_R_x): Ri_x += 1e-15  # hack: force lil sparse format to allocate
            M.glob.R[Ig] += Ri
            M.glob.R_U[Ig, 4*Is[i-1]:(4*Is[i-1]+4)] += Ri_U[:,0:4]
            M.glob.R_U[Ig, 4*Is[i  ]:(4*Is[i  ]+4)] += Ri_U[:,4:8]
            M.glob.R_x[Ig, [Is[i-1],Is[i]]] += Ri_x
      
            # following transition, all stations will be turbulent
            if (tran): param.turb = True
      
    # include effects of R_x into R_U: R_ue += R_x*x_st*st_ue
    #   The global residual Jacobian has a column for ue sensitivity
    #   ue, the edge velocity, also affects the location of the stagnation point
    #   The location of the stagnation point (st) dictates the x value at each node
    #   The residual also depends on the x value at each node (R_x)
    #   We use the chain rule (formula above) to account for this
    Nsys = M.glob.Nsys  # number of dofs
    Iue = range(3,4*Nsys,4)  # ue indices in U
    x_st = -M.isol.sgnue  # st = stag point [Nsys x 1]
    x_st = np.concatenate((x_st, -np.ones(M.wake.N)))  # wake same sens as upper surface
    R_st = M.glob.R_x @ x_st[:,np.newaxis]  # [3*Nsys x 1]
    Ist, st_ue = M.isol.Istag, M.isol.sstag_ue # stag points, sens
    if (alloc_R_x) or (alloc_R_U): R_st += 1e-15; # hack to avoid sparse matrix warning
    M.glob.R_U[:,Iue[Ist[0]]] += R_st*st_ue[0]
    M.glob.R_U[:,Iue[Ist[1]]] += R_st*st_ue[1]

    if (alloc_R_U): M.glob.R_U = M.glob.R_U.tocsr()
    if (alloc_R_x): M.glob.R_x = M.glob.R_x.tocsr()


#-------------------------------------------------------------------------------
def stagnation_state(U, x):
    # extrapolates two states in U, first ones in BL, to stagnation
    # INPUT
    #   U  : [U1,U2] = states at first two nodes (4x2)
    #   x  : [x1,x2] = x-locations of first two nodes (2x1)
    # OUTPUT
    #   Ust    : stagnation state (4x1)
    #   Ust_U  : linearization of Ust w.r.t. U1 and U2 (4x8)
    #   Ust_x  : linearization of Ust w.r.t. x1 and x2 (4x2)
    #   xst    : stagnation point location ... close to 0
    # DETAILS
    #   fits a quadratic to the edge velocity: 0 at x=0, then through two states
    #   linearly extrapolates other states in U to x=0, from U1 and U2

    # pull off states
    U1, U2, x1, x2 = U[:,0], U[:,1], x[0], x[1]
    dx = x2-x1; dx_x = np.r_[-1, 1]
    rx = x2/x1; rx_x = np.r_[-rx,1]/x1
  
    # linear extrapolation weights and stagnation state
    w1 =  x2/dx; w1_x = -w1/dx*dx_x + np.r_[ 0,1]/dx
    w2 = -x1/dx; w2_x = -w2/dx*dx_x + np.r_[-1,0]/dx
    Ust = U1*w1 + U2*w2
  
    # quadratic extrapolation of the edge velocity for better slope, ue=K*x
    wk1 = rx/dx; wk1_x = rx_x/dx - wk1/dx*dx_x
    wk2 = -1/(rx*dx); wk2_x = -wk2*(rx_x/rx + dx_x/dx) 
    K = wk1*U1[3] + wk2*U2[3]
    K_U = np.r_[0,0,0,wk1, 0,0,0,wk2]
    K_x = U1[3]*wk1_x + U2[3]*wk2_x
  
    # stagnation coord cannot be zero, but must be small
    xst = 1e-6
    Ust[3] = K*xst  # linear dep of ue on x near stagnation
    Ust_U = np.block([[w1*np.eye(3,4), w2*np.eye(3,4)], [K_U*xst]])
    Ust_x = np.vstack((np.outer(U1[0:3],w1_x) + np.outer(U2[0:3],w2_x), K_x*xst))
 
    return Ust, Ust_U, Ust_x, xst


#-------------------------------------------------------------------------------
def thwaites_init(K, nu):
    # uses Thwaites correlation to initialize first node in stag point flow
    # INPUT
    #   K  : stagnation point constant
    #   nu : kinematic viscosity
    # OUTPUT
    #   th : momentum thickness
    #   ds : displacement thickness
    # DETAILS
    #   ue = K*x -> K = ue/x = stag point flow constant
    #   th^2 = ue^(-6) * 0.45 * nu * int_0^x ue^5 dx = 0.45*nu/(6*K)
    #   ds = Hstag*th = 2.2*th

    th = np.sqrt(0.45*nu/(6.*K)) # momentum thickness
    ds = 2.2*th  # displacement thickness

    return th, ds

#-------------------------------------------------------------------------------
def wake_sys(M, param):
    # constructs residual system corresponding to wake initialization
    # INPUT
    #   param  : parameters
    # OUTPUT
    #   R   : 3x1 residual vector for th, ds, sa
    #   R_U : 3x12 residual linearization, as three 3x4 blocks
    #   J   : indices of the blocks of U in R_U (lower, upper, wake)

    il = M.vsol.Is[0][-1]; Ul = M.glob.U[:,il] # lower surface TE index, state
    iu = M.vsol.Is[1][-1]; Uu = M.glob.U[:,iu] # upper surface TE index, state
    iw = M.vsol.Is[2][ 0]; Uw = M.glob.U[:,iw] # first wake index, state
    t,hTE,dtdx,tcp,tdp = TE_info(M.foil.x) # trailing-edge gap is hTE

    # Obtain wake shear stress from upper/lower; transition if not turb
    param.turb = True; param.wake = False # calculating turbulent quantities right before wake
    if (M.vsol.turb[il]): ctl = Ul[2]; ctl_Ul = np.r_[0,0,1,0] # already turb; use state
    else: ctl, ctl_Ul = get_cttr(Ul, param) # transition shear stress, lower
    if (M.vsol.turb[iu]): ctu = Uu[2]; ctu_Uu = np.r_[0,0,1,0] # already turb; use state
    else: ctu, ctu_Uu = get_cttr(Uu, param) # transition shear stress, upper
    thsum = Ul[0] + Uu[0]  # sum of thetas
    ctw = (ctl*Ul[0] + ctu*Uu[0])/thsum  # theta-average
    ctw_Ul = (ctl_Ul*Ul[0] + (ctl - ctw)*np.r_[1,0,0,0])/thsum
    ctw_Uu = (ctu_Uu*Uu[0] + (ctu - ctw)*np.r_[1,0,0,0])/thsum

    # residual; note, delta star in wake includes the TE gap, hTE
    R = np.r_[Uw[0]-(Ul[0]+Uu[0]), Uw[1]-(Ul[1]+Uu[1]+hTE), Uw[2]-ctw]
    J = [il, iu, iw]  # R depends on states at these nodes
    R_Ul = np.vstack((-np.eye(2,4), -ctw_Ul))
    R_Uu = np.vstack((-np.eye(2,4), -ctw_Uu)) 
    R_Uw = np.eye(3,4)
    R_U = np.hstack((R_Ul, R_Uu, R_Uw))
    return R, R_U, J
  

#-------------------------------------------------------------------------------
def wake_init(M, ue):
    # initializes the first point of the wake, using data in M.glob.U
    # INPUT
    #   ue  : edge velocity at the wake point
    # OUTPUT
    #   Uw  : 4x1 state vector at the wake point
  
    iw = M.vsol.Is[2][0]; Uw = M.glob.U[:,iw]  # first wake index, state
    [R, R_U, J] = wake_sys(M, M.param) # construct the wake system
    Uw[0:3] -= R; Uw[3] = ue  # solve the wake system, use ue
    return Uw


#-------------------------------------------------------------------------------
def build_param(M, si):
    # builds a parameter structure for side is
    # INPUT
    #   si  : side number, 0 = lower, 1 = upper, 2 = wake
    # OUTPUT
    #   param : M.param structure with side information

    param = copy.deepcopy(M.param)  
    param.wake = (si == 2)
    param.turb = param.wake # the wake is fully turbulent
    param.simi = False # true for similarity station  
    return param


#-------------------------------------------------------------------------------
def station_param(M, param, i):
    # modifies parameter structure to be specific for station i
    # INPUT
    #   i  : station number (node index along the surface)
    # OUTPUT
    #   param : modified parameter structure
    param.turb = M.vsol.turb[i] # turbulent
    param.simi = i in M.isol.Istag # similarity


#-------------------------------------------------------------------------------
def init_boundary_layer(M):
    # initializes BL solution on foil and wake by marching with given edge vel, ue
    # INPUT
    #   The edge velocity field ue must be filled in on the airfoil and wake
    # OUTPUT
    #   The state in M.glob.U is filled in for each point

    Hmaxl = 3.8 # above this shape param value, laminar separation occurs
    Hmaxt = 2.5 # above this shape param value, turbulent separation occurs 
  
    ueinv = get_ueinv(M) # get inviscid velocity

    M.glob.Nsys = M.foil.N + M.wake.N; # number of global variables (nodes)
  
    # do we need to initialize?
    if (not M.oper.initbl) and (M.glob.U.shape[1]==M.glob.Nsys):
        vprint(M.param,1, '\n <<< Starting with current boundary layer >>> \n')
        M.glob.U[3,:] = ueinv  # do set a new edge velocity
        return
    
    vprint(M.param,1, '\n <<< Initializing the boundary layer >>> \n')
  
    M.glob.U = np.zeros((4,M.glob.Nsys)) # global solution matrix
    M.vsol.turb = np.zeros(M.glob.Nsys,dtype=int) # node flag: 0 = lam, 1 = turb
  
    for si in range(3):  # loop over surfaces
    
        vprint(M.param, 3, '\nSide is = %d:\n'%(si))
    
        Is = M.vsol.Is[si] # surface point indices
        xi = M.isol.xi[Is] # distance from LE stag point
        ue = ueinv[Is] # edge velocities
        N = len(Is) # number of points
        U = np.zeros([4,N]) # states at all points: [th, ds, sa, ue]
        Aux = np.zeros(N) # auxiliary data at all points: [wgap]
    
        # ensure edge velocities are not tiny
        uemax = max(abs(ue))
        for i in range(N): ue[i] = max(ue[i], 1e-8*uemax)
     
        # get parameter structure
        param = build_param(M, si)
    
        # set auxiliary data
        if (si == 2): Aux[:]= M.vsol.wgap
    
        # initialize state at first point
        i0 = 0
        if (si < 2): 

            # Solve for the stagnation state (Thwaites initialization + Newton)
            if (xi[0]<1e-8*xi[-1]): K, hitstag = ue[1]/xi[1], True 
            else: K, hitstag = ue[0]/xi[0], False
            th, ds = thwaites_init(K, param.mu0/param.rho0)
            xst = 1.e-6; # small but nonzero
            Ust = np.array([th, ds, 0, K*xst])
            nNewton = 20
            for iNewton in range(nNewton):
                # call residual at stagnation
                param.turb, param.simi = False, True  # similarity station flag 
                R, R_U, R_x = residual_station(param, np.r_[xst,xst], np.stack((Ust,Ust),axis=-1), np.zeros(2))
                param.simi = False
                if (norm2(R) < 1e-10): break
                A = R_U[:,4:7] + R_U[:,0:3]; b = -R; dU = np.append(np.linalg.solve(A,b), 0)
                # under-relaxation
                dm = max(abs(dU[0]/Ust[0]), abs(dU[1]/Ust[1]))
                omega = 1 if (dm < 0.2) else 0.2/dm
                dU = dU*omega
                Ust = Ust + dU
      
            # store stagnation state in first one (rarely two) points
            if (hitstag):
                U[:,0] = Ust; U[3,0] = ue[0]; i0=1
            U[:,i0] = Ust; U[3,i0] = ue[i0]
      
        else: # wake
            U[:,0] = wake_init(M, ue[0])  # initialize wake state properly
            param.turb = True # force turbulent in wake if still laminar
            M.vsol.turb[Is[0]] = True # wake starts turbulent

        # march over rest of points
        tran = False # flag indicating that we are at transition
        i = i0+1
        while (i<N):
            Ip = [i-1,i]; # two points involved in the calculation
            U[:,i] = U[:,i-1]; U[3,i] = ue[i]  # guess = same state, new ue
            if (tran): # set shear stress at transition interval
                ct, ct_U = get_cttr(U[:,i], param); U[2,i] = ct
            M.vsol.turb[Is[i]] = (tran or param.turb)  # flag node i as turbulent
            direct = True  # default is direct mode
            nNewton, iNswitch = 30, 12
            for iNewton in range(nNewton):
        
                # call residual at this station
                if (tran): # we are at transition
                    vprint(param, 4, 'i=%d, residual_transition (iNewton = %d) \n'%(i, iNewton))
                    try:
                        R, R_U, R_x = residual_transition(M, param, xi[Ip], U[:,Ip], Aux[Ip])
                    except:
                        vprint(param, 1, 'Transition calculation failed in BL init. Continuing.')
                        M.vsol.xt = 0.5*sum(xi[Ip])
                        U[:,i] = U[:,i-1]; U[3,i] = ue[i]; U[2,i] = ct
                        R = 0  # so we move on
                else:
                    vprint(param, 4, 'i=%d, residual_station (iNewton = %d)'%(i, iNewton))
                    R, R_U, R_x = residual_station(param, xi[Ip], U[:,Ip], Aux[Ip])
                if (norm2(R) < 1e-10): break
        
                if (direct): # direct mode => ue is prescribed => solve for th, ds, sa
                    A = R_U[:,4:7]; b = -R; dU = np.append(np.linalg.solve(A,b), 0)
                else: # inverse mode => Hk is prescribed 
                    Hk, Hk_U = get_Hk(U[:,i], param);
                    A = np.vstack((R_U[:,4:8], Hk_U)); b = np.r_[-R, Hktgt-Hk]
                    dU = np.linalg.solve(A,b)
          
                # under-relaxation
                dm = max(abs(dU[0]/U[0,i-1]), abs(dU[1]/U[1,i-1]))
                if (not direct): dm = max(dm, abs(dU[3]/U[3,i-1]))
                if (param.turb): dm = max(dm, abs(dU[2]/U[2,i-1]))
                elif (direct): dm = max(dm, abs(dU[2]/10)) 

                omega = 0.3/dm if (dm > 0.3) else 1
                dU = dU*omega

                # trial update
                Ui = U[:,i] + dU
        
                # clip extreme values
                if (param.turb): Ui[2] = max(min(Ui[2], .3), 1e-7)
                #Hklim = 1.02; if (param.wake), Hklim = 1.00005; end
                #[Hk,Hk_U] = get_Hk(Ui, param);
                #dH = max(0,Hklim-Hk); Ui(2) = Ui(2) + dH*Ui(1);
        
                # check if about to separate
                Hmax = Hmaxt if (param.turb) else Hmaxl
                Hk, Hk_U = get_Hk(Ui, param)

                if (direct) and ((Hk>Hmax) or (iNewton > iNswitch)):
                    # no update; need to switch to inverse mode: prescribe Hk
                    direct = False
                    vprint(param, 2, '** switching to inverse: i=%d, iNewton=%d'%(i, iNewton))
                    [Hk,Hk_U] = get_Hk(U[:,i-1], param); Hkr = (xi[i]-xi[i-1])/U[0,i-1]
                    if (param.wake):
                        H2 = Hk 
                        for k in range(6): H2 -= (H2+.03*Hkr*(H2-1)**3-Hk)/(1+.09*Hkr*(H2-1)**2)
                        Hktgt = max(H2, 1.01)
                    elif (param.turb): Hktgt = Hk - .15*Hkr  # turb: decrease in Hk
                    else: Hktgt = Hk + .03*Hkr # lam: increase in Hk 
                    if (not param.wake): Hktgt = max(Hktgt, Hmax)
                    if (iNewton > iNswitch): # reinit 
                        U[:,i] = U[:,i-1]
                        U[3,i] = ue[i]
                else: U[:,i] = Ui  # take the update


            if (iNewton >= nNewton-1):
                vprint(param, 1, '** BL init not converged: si=%d, i=%d **\n'%(si, i))
                # extrapolate values
                U[:,i] = U[:,i-1]; U[3,i] = ue[i]
                if (si<3):
                    U[0,i] = U[0,i-1]*(xi[i]/xi[i-1])**.5
                    U[1,i] = U[1,i-1]*(xi[i]/xi[i-1])**.5
                else:
                    rlen = (xi[i]-xi[i-1])/(10.*U[1,i-1])
                    U[1,i] = (U[1,i-1] + U[0,i-1]*rlen)/(1.+rlen)  # TODO check on this extrap

      
            # check for transition
            if (not param.turb) and (not tran) and (U[2,i]>param.ncrit):
                vprint(param,2, 'Identified transition at (si=%d, i=%d): n=%.5f, ncrit=%.5f\n'%(si, i, U[2,i], param.ncrit))
                tran = True 
                continue # redo station with transition
      
            if (tran):
                store_transition(M, si, i)  # store transition location
                param.turb = True; tran = False  # turbulent after transition

            i += 1  # next point

        # store states
        M.glob.U[:,Is] = U

#-------------------------------------------------------------------------------
def store_transition(M, si, i):
    # stores xi and x transition locations using current M.vsol.xt 
    # INPUT
    #   si,i : side,station number
    # OUTPUT
    #   M.vsol.Xt stores the transition location s and x values
  
    xt = M.vsol.xt;
    i0, i1 = M.vsol.Is[si][i-1], M.vsol.Is[si][i] # pre/post transition nodes
    xi0, xi1 = M.isol.xi[i0], M.isol.xi[i1]  # xi (s) locations at nodes
    assert ((i0<M.foil.N) and (i1<M.foil.N)), 'Can only store transition on airfoil'
    x0, x1 = M.foil.x[0,i0], M.foil.x[0,i1]  # x locations at nodes
    if ((xt<xi0) or (xt>xi1)):
        vprint(M.param,1, 'Warning: transition (%.3f) off interval (%.3f,%.3f)!'%(xt, xi0, xi1))
    M.vsol.Xt[si,0] = xt  # xi location
    M.vsol.Xt[si,1] = x0 + (xt-xi0)/(xi1-xi0)*(x1-x0) # x location
    slu = ['lower', 'upper']
    vprint(M.param,1, '  transition on %s side at x=%.5f'%(slu[si], M.vsol.Xt[si,1]))


#-------------------------------------------------------------------------------
def update_transition(M):
    # updates transition location using current state
    # INPUT
    #   a valid state in M.glob.U
    # OUTPUT
    #   M.vsol.turb : updated with latest lam/turb flags for each node
    #   M.glob.U    : updated with amp factor or shear stress as needed at each node
  
    for si in range(2):  # loop over lower/upper surfaces
    
        Is = M.vsol.Is[si]  # surface point indices
        N = len(Is) # number of points

        # get parameter structure
        param = build_param(M, si)
    
        # current last laminar station
        for ilam0 in range(N):
            if (M.vsol.turb[Is[ilam0]]): 
                ilam0 -=1
                break
        
        # current amp/ctau solution (so we do not change it unnecessarily)
        sa = M.glob.U[2,Is].copy()
    
        # march amplification equation to get new last laminar station
        ilam = march_amplification(M, si)
    
        if (ilam == ilam0):
            M.glob.U[2,Is] = sa[:]  # no change
            continue
    
        vprint(param, 2, '  Update transition: last lam [%d]->[%d]'%(ilam0, ilam))
    
        if (ilam < ilam0):
            # transition is now earlier: fill in turb between [ilam+1, ilam0]
            param.turb = True
            sa0, temp = get_cttr(M.glob.U[:,Is[ilam+1]], param)
            sa1 = M.glob.U[2,Is[ilam0+1]] if (ilam0<N-1) else sa0
            xi = M.isol.xi[Is]
            dx = xi[min(ilam0+1,N-1)]-xi[ilam+1]
            for i in range(ilam+1, ilam0+1):
                f = 0 if (dx==0) or (i==ilam+1) else (xi[i]-xi[ilam+1])/dx
                if ((ilam+1) == ilam0): f = 1
                M.glob.U[2,Is[i]] = sa0 + f*(sa1-sa0)
                assert M.glob.U[2,Is[i]] > 0, 'negative ctau in update_transition'
                M.vsol.turb[Is[i]] = True

        elif (ilam > ilam0):
            # transition is now later: lam already filled in; leave turb alone
            for i in range(ilam0,ilam+1): M.vsol.turb[Is[i]] = False


#-------------------------------------------------------------------------------
def march_amplification(M, si):
    # marches amplification equation on surface si
    # INPUT
    #   si : surface number index
    # OUTPUT
    #   ilam : index of last laminar station before transition
    #   M.glob.U : updated with amp factor at each (new) laminar station
  
    Is = M.vsol.Is[si]  # surface point indices
    N = len(Is)  # number of points
    param = build_param(M, si)  # get parameter structure
    U = M.glob.U[:,Is]  # states
    turb = M.vsol.turb[Is]  # turbulent station flag
  
    # loop over stations, calculate amplification
    U[2,0] = 0.; # no amplification at first station
    param.turb, param.wake = False, False
    i = 1
    while (i < N):
        U1, U2 = U[:,i-1], U[:,i].copy()  # states
        if turb[i]: U2[2] = U1[2]*1.01 # initialize amp if turb
        dx = M.isol.xi[Is[i]] - M.isol.xi[Is[i-1]] # interval length
        
        # Newton iterations, only needed if adding extra amplification in damp
        nNewton = 20;
        for iNewton in range(nNewton):
            # amplification rate, averaged
            damp1, damp1_U1 = get_damp(U1, param)
            damp2, damp2_U2 = get_damp(U2, param)
            damp, damp_U = upwind(0.5, 0, damp1, damp1_U1, damp2, damp2_U2)
          
            Ramp = U2[2] - U1[2] - damp*dx

            if (iNewton > 11):
                vprint(param,3,'i=%d, iNewton=%d, sa = [%.5e, %.5e], damp = %.5e, Ramp = %.5e'%(
                    i, iNewton, U1[2], U2[2], damp, Ramp))
      
            if (abs(Ramp)<1e-12): break  # converged
            Ramp_U = np.r_[0,0,-1,0, 0,0,1,0] - damp_U*dx
            dU = -Ramp/Ramp_U[6]
            omega = 1; dmax = 0.5*(1.01-iNewton/nNewton)
            if (abs(dU) > dmax): omega = dmax/abs(dU)
            U2[2] += omega*dU

        if (iNewton >= nNewton-1): vprint(param,1, 'march amp Newton unconverged!')

        # check for transition
        if (U2[2]>param.ncrit):
            vprint(param, 2,'  march_amplification (si,i=%d,%d): %.5e is above critical.'%(si, i, U2[2]))
            break
        else:
            M.glob.U[2,Is[i]] = U2[2]  # store amplification in M.glob.U (also seen in view U)
            U[2,i] = U2[2]
            if (np.iscomplex(U[2,i])): raise ValueError('imaginary amp during march')

    
        i += 1  # next station
  
    return (i-1) # return last laminar station
  

#-------------------------------------------------------------------------------
def residual_transition(M, param, x, U, Aux):
    # calculates the combined lam + turb residual for a transition station
    # INPUT
    #   param : parameter structure
    #   x     : 2x1 vector, [x1, x2], containing xi values at the points
    #   U     : 4x2 matrix, [U1, U2], containing the states at the points
    #   Aux   : ()x2 matrix, [Aux1, Aux2] of auxiliary data at the points
    # OUTPUT
    #   R     : 3x1 transition residual vector
    #   R_U   : 3x8 residual Jacobian, [R_U1, R_U2]
    #   R_x   : 3x2 residual linearization w.r.t. x, [R_x1, R_x2]
    # DETAILS
    #   The state U1 should be laminar; U2 should be turbulent
    #   Calculates and linearizes the transition location in the process
    #   Assumes linear variation of th and ds from U1 to U2
  
    # states
    U1 = U[:,0]; U2 = U[:,1]; sa = U[2,:]
    I1 = range(4); I2 = range(4,8); Z = np.zeros(4)
  
    # interval
    x1 = x[0]; x2 = x[1]; dx = x2-x1
  
    # determine transition location (xt) using amplification equation
    xt = x1 + 0.5*dx  # guess
    ncrit = param.ncrit  # critical amp factor
    nNewton = 20
    vprint(param, 3, '  Transition interval = [%.5e, %.5e]'%(x1, x2))
    #  U1, U2
    for iNewton in range(nNewton):
        w2 = (xt-x1)/dx; w1 = 1.-w2  # weights
        Ut = w1*U1 + w2*U2; Ut_xt = (U2-U1)/dx  # state at xt
        Ut[2] = ncrit; Ut_xt[2] = 0.; # amplification at transition
        damp1, damp1_U1 = get_damp(U1, param)
        dampt, dampt_Ut = get_damp(Ut, param); dampt_Ut[2] = 0.
        Rxt = ncrit - sa[0] - 0.5*(xt-x1)*(damp1 + dampt)
        Rxt_xt = -0.5*(damp1+dampt) - 0.5*(xt-x1)*np.dot(dampt_Ut,Ut_xt)
        dxt = -Rxt/Rxt_xt
        vprint(param, 4, '   Transition: iNewton,Rxt,xt = %d,%.5e,%.5e'%(iNewton,Rxt,xt))
        dmax = 0.2*dx*(1.1-iNewton/nNewton)
        if (abs(dxt)>dmax): dxt = dxt*dmax/abs(dxt)
        if (abs(Rxt) < 1e-10): break
        if (iNewton<nNewton): xt += dxt

    if (iNewton >= nNewton): vprint(param, 1, 'Transition location calculation failed.')
    M.vsol.xt = xt  # save transition location
  
    # prepare for xt linearizations
    Rxt_U = -0.5*(xt-x1)*np.concatenate((damp1_U1+dampt_Ut*w1, dampt_Ut*w2)); Rxt_U[2] -= 1.
    Ut_x1 = (U2-U1)*(w2-1)/dx; Ut_x2 = (U2-U1)*(-w2)/dx  # at fixed xt
    Ut_x1[2] = 0; Ut_x2[2] = 0  # amp at xt is always ncrit
    Rxt_x1 = 0.5*(damp1+dampt) - 0.5*(xt-x1)*np.dot(dampt_Ut,Ut_x1)
    Rxt_x2 =                   - 0.5*(xt-x1)*np.dot(dampt_Ut,Ut_x2)
  
    # sensitivity of xt w.r.t. U,x from Rxt(xt,U,x) = 0 constraint
    xt_U = -Rxt_U/Rxt_xt;  xt_U1 = xt_U[I1]; xt_U2 = xt_U[I2]
    xt_x1 = -Rxt_x1/Rxt_xt; xt_x2 = -Rxt_x2/Rxt_xt

    # include derivatives w.r.t. xt in Ut_x1 and Ut_x2
    Ut_x1 += Ut_xt*xt_x1
    Ut_x2 += Ut_xt*xt_x2
    
    # sensitivity of Ut w.r.t. U1 and U2
    Ut_U1 = w1*np.eye(4) + np.outer((U2-U1),xt_U1)/dx # w1*I + U1*w1_xt*xt_U1 + U2*w2_xt*xt_U1;
    Ut_U2 = w2*np.eye(4) + np.outer((U2-U1),xt_U2)/dx # w2*I + U1*w1_xt*xt_U2 + U2*w2_xt*xt_U2;
  
    # laminar and turbulent states at transition
    Utl = Ut.copy(); Utl_U1 = Ut_U1.copy(); Utl_U2 = Ut_U2.copy(); Utl_x1 = Ut_x1.copy(); Utl_x2 = Ut_x2.copy();
    Utl[2] = ncrit; Utl_U1[2,:] = Z; Utl_U2[2,:] = Z; Utl_x1[2] = 0; Utl_x2[2] = 0;
    Utt = Ut.copy(); Utt_U1 = Ut_U1.copy(); Utt_U2 = Ut_U2.copy(); Utt_x1 = Ut_x1.copy(); Utt_x2 = Ut_x2.copy();
  
    # parameter structure
    param = build_param(M, 0)
  
    # set turbulent shear coefficient, sa, in Utt
    param.turb = True
    cttr, cttr_Ut = get_cttr(Ut, param)
    Utt[2] = cttr; Utt_U1[2,:] = np.dot(cttr_Ut, Ut_U1); Utt_U2[2,:] = np.dot(cttr_Ut, Ut_U2)
    Utt_x1[2] = np.dot(cttr_Ut,Ut_x1); Utt_x2[2] = np.dot(cttr_Ut,Ut_x2)
  
    # laminar/turbulent residuals and linearizations
    param.turb = False
    Rl, Rl_U, Rl_x = residual_station(param, np.r_[x1,xt], np.stack((U1,Utl),axis=-1), Aux)
    Rl_U1 = Rl_U[:,I1]; Rl_Utl = Rl_U[:,I2]
    param.turb = True
    Rt, Rt_U, Rt_x = residual_station(param, np.r_[xt,x2], np.stack((Utt,U2),axis=-1), Aux)
    Rt_Utt = Rt_U[:,I1]; Rt_U2 = Rt_U[:,I2]

    # combined residual and linearization
    R = Rl + Rt
    if (any(np.imag(R) != 0)): raise ValueError('imaginary transition residual')
    R_U1 = Rl_U1 + np.dot(Rl_Utl,Utl_U1) + np.outer(Rl_x[:,1],xt_U1) + np.dot(Rt_Utt,Utt_U1) + np.outer(Rt_x[:,0], xt_U1)
    R_U2 = np.dot(Rl_Utl,Utl_U2) + np.outer(Rl_x[:,1], xt_U2) + np.dot(Rt_Utt,Utt_U2) + Rt_U2 + np.outer(Rt_x[:,0], xt_U2)
    R_U = np.hstack((R_U1, R_U2))
    R_x = np.stack((Rl_x[:,0] + Rl_x[:,1]*xt_x1 + Rt_x[:,0]*xt_x1 + np.dot(Rl_Utl,Utl_x1) + np.dot(Rt_Utt,Utt_x1),
                    Rt_x[:,1] + Rl_x[:,1]*xt_x2 + Rt_x[:,0]*xt_x2 + np.dot(Rl_Utl,Utl_x2) + np.dot(Rt_Utt,Utt_x2)), axis=-1)
                   
    return R, R_U, R_x
  

#-------------------------------------------------------------------------------
def residual_station(param, x, Uin, Aux):
    # calculates the viscous residual at one non-transition station
    # INPUT
    #   param : parameter structure
    #   x     : 2x1 vector, [x1, x2], containing xi values at the points
    #   U     : 4x2 matrix, [U1, U2], containing the states at the points
    #   Aux   : ()x2 matrix, [Aux1, Aux2] of auxiliary data at the points
    # OUTPUT
    #   R     : 3x1 residual vector (mom, shape-param, amp/lag)
    #   R_U   : 3x8 residual Jacobian, [R_U1, R_U2]
    #   R_x   : 3x2 residual linearization w.r.t. x, [R_x1, R_x2]
    # DETAILS
    #   The input states are U = [U1, U2], each with th,ds,sa,ue

    # so that we do not overwrite Uin
    U = Uin.copy()

    # modify ds to take out wake gap (in Aux) for all calculations below
    U[1,:] -= Aux
  
    # states
    U1 = U[:,0]; U2 = U[:,1]; Um = 0.5*(U1+U2)
    th = U[0,:]; ds = U[1,:]; sa = U[2,:] 
  
    # speed needs compressibility correction
    uk1, uk1_u = get_uk(U1[3],param)
    uk2, uk2_u = get_uk(U2[3],param)
    
    # log changes
    thlog = np.log(th[1]/th[0])
    thlog_U = np.r_[-1./th[0],0,0,0, 1./th[1],0,0,0]
    uelog = np.log(uk2/uk1)
    uelog_U = np.r_[0,0,0,-uk1_u/uk1, 0,0,0,uk2_u/uk2]
    xlog = np.log(x[1]/x[0]); xlog_x = np.r_[-1./x[0], 1./x[1]]
    dx = x[1]-x[0]; dx_x = np.r_[-1., 1.]
  
    # upwinding factor
    upw, upw_U = get_upw(U1, U2, param)
  
    # shape parameter
    H1, H1_U1 = get_H(U1)
    H2, H2_U2 = get_H(U2)
    H = 0.5*(H1+H2)
    H_U = 0.5*np.r_[H1_U1, H2_U2]
  
    # Hstar = KE shape parameter, averaged
    Hs1, Hs1_U1 = get_Hs(U1, param)
    Hs2, Hs2_U2 = get_Hs(U2, param)
    Hs, Hs_U = upwind(0.5, 0, Hs1, Hs1_U1, Hs2, Hs2_U2)  
  
    # log change in Hstar
    Hslog = np.log(Hs2/Hs1)
    Hslog_U = np.r_[-1./Hs1*Hs1_U1, 1./Hs2*Hs2_U2]
  
    # similarity station is special: U1 = U2, x1 = x2
    if (param.simi):
        thlog = 0.; thlog_U *= 0.
        Hslog = 0.; Hslog_U *= 0.
        uelog = 1.; uelog_U *= 0.
        xlog = 1.; xlog_x = np.r_[0., 0.]
        dx = 0.5*(x[0]+x[1]); dx_x = np.r_[0.5,0.5]
    
    # Hw = wake shape parameter
    Hw1, Hw1_U1 = get_Hw(U1, Aux[0])
    Hw2, Hw2_U2 = get_Hw(U2, Aux[1])
    Hw = 0.5*(Hw1 + Hw2)
    Hw_U = 0.5*np.r_[Hw1_U1, Hw2_U2]
  
    # set up shear lag or amplification factor equation
    if (param.turb):

        # log change of root shear stress coeff
        salog = np.log(sa[1]/sa[0])
        salog_U = np.r_[0,0,-1./sa[0],0, 0,0,1./sa[1],0]
    
        # BL thickness measure, averaged
        de1, de1_U1 = get_de(U1, param)
        de2, de2_U2 = get_de(U2, param)
        de, de_U = upwind(0.5, 0, de1, de1_U1, de2, de2_U2)
    
        # normalized slip velocity, averaged
        Us1, Us1_U1 = get_Us(U1, param)
        Us2, Us2_U2 = get_Us(U2, param)
        Us, Us_U = upwind(0.5, 0, Us1, Us1_U1, Us2, Us2_U2)
    
        # Hk, upwinded
        Hk1, Hk1_U1 = get_Hk(U1, param)
        Hk2, Hk2_U2 = get_Hk(U2, param)
        Hk, Hk_U = upwind(upw, upw_U, Hk1, Hk1_U1, Hk2, Hk2_U2)
        
        # Re_theta, averaged
        Ret1, Ret1_U1 = get_Ret(U1, param)
        Ret2, Ret2_U2 = get_Ret(U2, param)
        Ret, Ret_U = upwind(0.5, 0, Ret1, Ret1_U1, Ret2, Ret2_U2)
    
        # skin friction, upwinded
        cf1, cf1_U1 = get_cf(U1, param)
        cf2, cf2_U2 = get_cf(U2, param)
        cf, cf_U = upwind(upw, upw_U, cf1, cf1_U1, cf2, cf2_U2)
    
        # displacement thickness, averaged
        dsa = 0.5*(ds[0] + ds[1])
        dsa_U = 0.5*np.r_[0,1,0,0, 0,1,0,0]
    
        # uq = equilibrium 1/ue * due/dx
        uq, uq_U = get_uq(dsa, dsa_U, cf, cf_U, Hk, Hk_U, Ret, Ret_U, param)
    
        # cteq = root equilibrium wake layer shear coeficient: (ctau eq)^.5
        cteq1, cteq1_U1 = get_cteq(U1, param)
        cteq2, cteq2_U2 = get_cteq(U2, param)
        cteq, cteq_U = upwind(upw, upw_U, cteq1, cteq1_U1, cteq2, cteq2_U2)

        # root of shear coefficient (a state), upwinded
        saa, saa_U = upwind(upw, upw_U, sa[0], np.r_[0,0,1,0], sa[1], np.r_[0,0,1,0])
    
        # lag coefficient
        Klag = param.SlagK
        beta = param.GB
        Clag = Klag/beta*1./(1.+Us)
        Clag_U = -Clag/(1.+Us)*Us_U
    
        # extra dissipation in wake
        ald = 1.0;
        if (param.wake): ald = param.Dlr
    
        # shear lag equation
        Rlag = Clag*(cteq-ald*saa)*dx - 2*de*salog + 2*de*(uq*dx-uelog)*param.Cuq
        Rlag_U = Clag_U*(cteq-ald*saa)*dx + Clag*(cteq_U-ald*saa_U)*dx \
            - 2*de_U*salog - 2*de*salog_U \
            + 2*de_U*(uq*dx-uelog)*param.Cuq + 2*de*(uq_U*dx-uelog_U)*param.Cuq
        Rlag_x = Clag*(cteq-ald*saa)*dx_x + 2*de*uq*dx_x
        
    else:
        # laminar, amplification factor equation
    
        if (param.simi):
            # similarity station
            Rlag = sa[0] + sa[1]  # no amplification
            Rlag_U = np.array([0,0,1,0, 0,0,1,0])
            Rlag_x = np.array([0,0])
        else:
            # amplification factor equation in Rlag

            # amplification rate, averaged
            damp1, damp1_U1 = get_damp(U1, param)
            damp2, damp2_U2 = get_damp(U2, param)
            damp, damp_U = upwind(0.5, 0, damp1, damp1_U1, damp2, damp2_U2)
      
            Rlag = sa[1] - sa[0] - damp*dx
            Rlag_U = np.array([0,0,-1,0, 0,0,1,0]) - damp_U*dx
            Rlag_x = -damp*dx_x
  
    # squared mach number, symmetrical average
    Ms1, Ms1_U1 = get_Mach2(U1, param)
    Ms2, Ms2_U2 = get_Mach2(U2, param)
    Ms, Ms_U = upwind(0.5, 0, Ms1, Ms1_U1, Ms2, Ms2_U2)
  
    # skin friction * x/theta, symmetrical average
    cfxt1, cfxt1_U1, cfxt1_x1 = get_cfxt(U1, x[0], param)
    cfxt2, cfxt2_U2, cfxt2_x2 = get_cfxt(U2, x[1], param)
    cfxtm, cfxtm_Um, cfxtm_xm = get_cfxt(Um, 0.5*(x[0]+x[1]), param)
    cfxt = 0.25*cfxt1 + 0.5*cfxtm + 0.25*cfxt2
    cfxt_U = 0.25*np.concatenate((cfxt1_U1+cfxtm_Um, cfxtm_Um+cfxt2_U2))
    cfxt_x = 0.25*np.array([cfxt1_x1+cfxtm_xm, cfxtm_xm+cfxt2_x2])
  
    # momentum equation
    Rmom = thlog + (2+H+Hw-Ms)*uelog - 0.5*xlog*cfxt
    Rmom_U = thlog_U + (H_U+Hw_U-Ms_U)*uelog + (2+H+Hw-Ms)*uelog_U - 0.5*xlog*cfxt_U
    Rmom_x = -0.5*xlog_x*cfxt - 0.5*xlog*cfxt_x

    # dissipation function times x/theta: cDi = (2*cD/H*)*x/theta, upwinded
    cDixt1, cDixt1_U1, cDixt1_x1 = get_cDixt(U1, x[0], param)
    cDixt2, cDixt2_U2, cDixt2_x2 = get_cDixt(U2, x[1], param)
    cDixt, cDixt_U = upwind(upw, upw_U, cDixt1, cDixt1_U1, cDixt2, cDixt2_U2)
    cDixt_x = np.array([(1.-upw)*cDixt1_x1, upw*cDixt2_x2])
  
    # cf*x/theta, upwinded
    cfxtu, cfxtu_U = upwind(upw, upw_U, cfxt1, cfxt1_U1, cfxt2, cfxt2_U2)
    cfxtu_x = np.array([(1.-upw)*cfxt1_x1, upw*cfxt2_x2])
  
    # Hss = density shape parameter, averaged
    [Hss1, Hss1_U1] = get_Hss(U1, param)
    [Hss2, Hss2_U2] = get_Hss(U2, param)
    [Hss, Hss_U] = upwind(0.5, 0, Hss1, Hss1_U1, Hss2, Hss2_U2)
  
    Rshape = Hslog + (2*Hss/Hs + 1-H-Hw)*uelog + xlog*(0.5*cfxtu - cDixt)
    Rshape_U = Hslog_U + (2*Hss_U/Hs - 2*Hss/Hs**2*Hs_U -H_U - Hw_U)*uelog \
        + (2*Hss/Hs + 1-H-Hw)*uelog_U + xlog*(0.5*cfxtu_U - cDixt_U)
    Rshape_x = xlog_x*(0.5*cfxtu - cDixt) + xlog*(0.5*cfxtu_x - cDixt_x)
  
    # put everything together
    R = np.array([Rmom, Rshape, Rlag])
    R_U = np.vstack((Rmom_U, Rshape_U, Rlag_U))
    R_x = np.vstack((Rmom_x, Rshape_x, Rlag_x))

    return R, R_U, R_x



# ============ GET FUNCTIONS ==============


#-------------------------------------------------------------------------------
def get_upw(U1, U2, param):
    # calculates a local upwind factor (0.5 = trap; 1 = BE) based on two states
    # INPUT
    #   U1,U2 : first/upwind and second/downwind states (4x1 each)
    #   param : parameter structure
    # OUTPUT
    #   upw   : scalar upwind factor
    #   upw_U : 1x8 linearization vector, [upw_U1, upw_U2]
    # DETAILS
    #   Used to ensure a stable viscous discretization
    #   Decision to upwind is made based on the shape factor change
  
    Hk1, Hk1_U1 = get_Hk(U1, param)
    Hk2, Hk2_U2 = get_Hk(U2, param)
    Z = np.zeros(Hk1_U1.shape)
    Hut = 1.  # triggering constant for upwinding
    C = 1. if (param.wake) else 5.
    Huc = C*Hut/Hk2**2; # only depends on U2
    Huc_U = np.concatenate((Z, -2*Huc/Hk2*Hk2_U2))
    aa = (Hk2-1.)/(Hk1-1.); sga = np.sign(aa)
    la = np.log(sga*aa)
    la_U = np.concatenate((-1./(Hk1-1.)*Hk1_U1, 1./(Hk2-1.)*Hk2_U2))
    Hls = la**2; Hls_U = 2*la*la_U
    if (Hls > 15): Hls, Hls_U = 15, Hls_U*0.
    upw = 1. - 0.5*np.exp(-Hls*Huc)
    upw_U = -0.5*np.exp(-Hls*Huc)*(-Hls_U*Huc-Hls*Huc_U)

    return upw, upw_U
  

#-------------------------------------------------------------------------------
def upwind(upw, upw_U, f1, f1_U1, f2, f2_U2):
    # calculates an upwind average (and derivatives) of two scalars
    # INPUT
    #   upw, upw_U : upwind scalar and its linearization w.r.t. U1,U2
    #   f1, f1_U   : first scalar and its linearization w.r.t. U1
    #   f2, f2_U   : second scalar and its linearization w.r.t. U2
    # OUTPUT
    #   f    : averaged scalar
    #   f_U  : linearization of f w.r.t. both states, [f_U1, f_U2]
  
    f = (1-upw)*f1 + upw*f2
    f_U = (-upw_U)*f1 + upw_U*f2 + np.concatenate(((1-upw)*f1_U1, upw*f2_U2))
  
    return f, f_U


#-------------------------------------------------------------------------------
def get_uq(ds, ds_U, cf, cf_U, Hk, Hk_U, Ret, Ret_U, param):
    # calculates the equilibrium 1/ue*due/dx
    # INPUT
    #   ds, ds_U   : delta star and linearization (1x4)
    #   cf, cf_U   : skin friction and linearization (1x4)
    #   Hk, Hk_U   : kinematic shape parameter and linearization (1x4)
    #   Ret, Ret_U : theta Reynolds number and linearization (1x4)
    #   param      : parameter structure
    # OUTPUT
    #   uq, uq_U   : equilibrium 1/ue*due/dx and linearization w.r.t. state (1x4)

    beta, A, C = param.GB, param.GA, param.GC
    if (param.wake): A, C = A*param.Dlr, 0.
    # limit Hk (TODO smooth/eliminate)
    if (param.wake) and (Hk < 1.00005): Hk, Hk_U = 1.00005, Hk_U*0.
    if (not param.wake) and (Hk < 1.05): Hk, Hk_U = 1.05, Hk_U*0.
    Hkc = Hk - 1. - C/Ret
    Hkc_U = Hk_U + C/Ret**2*Ret_U
  
    if (Hkc < .01): Hkc, Hkc_U = .01, Hkc_U*0.
    ut = 0.5*cf - (Hkc/(A*Hk))**2
    ut_U = 0.5*cf_U - 2*(Hkc/(A*Hk))*(Hkc_U/(A*Hk) - Hkc/(A*Hk**2)*Hk_U)
    uq = ut/(beta*ds);
    uq_U = ut_U/(beta*ds) - uq/ds * ds_U
  
    return uq, uq_U


#-------------------------------------------------------------------------------
def get_cttr(U, param):
    # calculates root of the shear stress coefficient at transition
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   param : parameter structure
    # OUTPUT
    #   cttr, cttr_U : sqrt(shear stress coeff) and its lin w.r.t. U (1x4)
    # DETAILS
    #   used to initialize the first turb station after transition
  
    param.wake = False  # transition happens just before the wake starts
    cteq, cteq_U = get_cteq(U, param)
    Hk, Hk_U = get_Hk(U, param)
    if (Hk < 1.05): Hk, Hk_U = 1.05, Hk_U*0.
    C, E = param.CtauC, param.CtauE
    c = C*np.exp(-E/(Hk-1.)); c_U = c*E/(Hk-1)**2*Hk_U
    cttr = c*cteq; cttr_U = c_U*cteq + c*cteq_U

    return cttr, cttr_U


#-------------------------------------------------------------------------------
def get_cteq(U, param):
    # calculates root of the equilibrium shear stress coefficient: sqrt(ctau_eq)
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   param : parameter structure
    # OUTPUT
    #   cteq, cteq_U : sqrt(equilibrium shear stress) and its lin w.r.t. U (1x4)
    # DETAILS
    #   uses equilibrium shear stress correlations
    CC, C = 0.5/(param.GA**2*param.GB), param.GC
    Hk, Hk_U = get_Hk(U, param)
    Hs, Hs_U = get_Hs(U, param)
    H, H_U = get_H(U)
    Ret, Ret_U = get_Ret(U, param)
    Us, Us_U = get_Us(U, param)
    if (param.wake):
        if (Hk < 1.00005): Hk, Hk_U = 1.00005, Hk_U*0.
        Hkc = Hk - 1.; Hkc_U = Hk_U
    else:
        if (Hk < 1.05): Hk, HK_U = 1.05, Hk_U*0.
        Hkc = Hk - 1. - C/Ret
        Hkc_U = Hk_U + C/Ret**2*Ret_U
        if (Hkc < 0.01): Hkc, Hkc_U = 0.01, Hkc_U*0.

    num = CC*Hs*(Hk-1)*Hkc**2
    num_U = CC*(Hs_U*(Hk-1)*Hkc**2 + Hs*Hk_U*Hkc**2 + Hs*(Hk-1)*2*Hkc*Hkc_U)
    den = (1-Us)*H*Hk**2
    den_U = (-Us_U)*H*Hk**2 + (1-Us)*H_U*Hk**2 + (1-Us)*H*2*Hk*Hk_U
    cteq = np.sqrt(num/den)
    cteq_U = 0.5/cteq*(num_U/den - num/den**2*den_U)

    return cteq, cteq_U



#-------------------------------------------------------------------------------
def get_Hs(U, param):
    # calculates Hs = Hstar = K.E. shape parameter, from U
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   param : parameter structure
    # OUTPUT
    #   Hs, Hs_U : Hstar and its lin w.r.t. U (1x4)
    # DETAILS
    #   Hstar is the ratio theta*/theta, where theta* is the KE thicknes
    Hk, Hk_U = get_Hk(U, param)
  
    # limit Hk (TODO smooth/eliminate)
    if (param.wake) and (Hk < 1.00005): Hk, Hk_U = 1.00005, Hk_U*0.
    if (not param.wake) and (Hk < 1.05): Hk, Hk_U = 1.05, Hk_U*0.
  
    if (param.turb): # turbulent
        Hsmin, dHsinf = 1.5, .015
        Ret, Ret_U = get_Ret(U, param)
        # limit Re_theta and dependence
        Ho = 4.; Ho_U = 0.
        if (Ret > 400): Ho, Ho_U = 3 + 400./Ret, -400./Ret**2*Ret_U
        Reb, Reb_U = Ret, Ret_U
        if (Ret < 200): Reb, Reb_U = 200, Reb_U*0.
        if (Hk < Ho):  # attached branch
            Hr = (Ho-Hk)/(Ho-1)
            Hr_U = (Ho_U - Hk_U)/(Ho-1) - (Ho-Hk)/(Ho-1)**2*Ho_U
            aa = (2-Hsmin-4/Reb)*Hr**2
            aa_U = (4/Reb**2*Reb_U)*Hr**2 + (2-Hsmin-4/Reb)*2*Hr*Hr_U
            Hs = Hsmin + 4/Reb + aa * 1.5/(Hk+.5)
            Hs_U = -4/Reb**2*Reb_U + aa_U*1.5/(Hk+.5) - aa*1.5/(Hk+.5)**2*Hk_U
        else:  # separated branch
            lrb = np.log(Reb); lrb_U = 1/Reb*Reb_U
            aa = Hk - Ho + 4/lrb
            aa_U = Hk_U - Ho_U - 4/lrb**2*lrb_U
            bb = .007*lrb/aa**2 + dHsinf/Hk
            bb_U = .007*(lrb_U/aa**2 - 2*lrb/aa**3*aa_U) - dHsinf/Hk**2*Hk_U
            Hs = Hsmin + 4/Reb + (Hk-Ho)**2*bb
            Hs_U = -4/Reb**2*Reb_U + 2*(Hk-Ho)*(Hk_U-Ho_U)*bb + (Hk-Ho)**2*bb_U
        # slight Mach number correction
        M2, M2_U = get_Mach2(U, param) # squared edge Mach number
        den = 1+.014*M2; den_M2 = .014
        Hs = (Hs+.028*M2)/den
        Hs_U = (Hs_U+.028*M2_U)/den - Hs/den*den_M2*M2_U
    else: # laminar
        a = Hk-4.35
        if (Hk < 4.35):
            num = .0111*a**2 - .0278*a**3
            Hs = num/(Hk+1) + 1.528 - .0002*(a*Hk)**2
            Hs_Hk = (.0111*2*a - .0278*3*a**2)/(Hk+1) - num/(Hk+1)**2 - .0002*2*a*Hk*(Hk+a)
        else:
            Hs = .015*a**2/Hk + 1.528
            Hs_Hk = .015*2*a/Hk - .015*a**2/Hk**2
        Hs_U = Hs_Hk*Hk_U

    return Hs, Hs_U



#-------------------------------------------------------------------------------
def get_cp(u, param):
    # calculates pressure coefficient from speed, with compressibility correction
    # INPUT
    #   u     : speed
    #   param : parameter structure
    # OUTPUT
    #   cp, cp_U : pressure coefficient and its linearization w.r.t. u
    # DETAILS
    #   Karman-Tsien correction is included

    Vinf = param.Vinf
    cp = 1-(u/Vinf)**2; cp_u = -2*u/Vinf**2
    if (param.Minf > 0):
        l, b = param.KTl, param.KTb
        den = b+0.5*l*(1+b)*cp; den_cp = 0.5*l*(1+b)
        cp /= den; cp_u *= (1-cp*den_cp)/den

    return cp, cp_u



#-------------------------------------------------------------------------------
def get_uk(u, param):
    # calculates Karman-Tsien corrected speed
    # INPUT
    #   u     : incompressible speed
    #   param : parameter structure
    # OUTPUT
    #   uk, uk_u : compressible speed and its linearization w.r.t. u
    # DETAILS
    #   Uses the Karman-Tsien correction, Minf from param
  
    if (param.Minf > 0):
        l, Vinf = param.KTl, param.Vinf
        den = 1-l*(u/Vinf)**2; den_u = -2*l*u/Vinf**2
        uk = u*(1-l)/den; uk_u = (1-l)/den - (uk/den)*den_u
    else:
        uk, uk_u = u, 1.

    return uk, uk_u
  

#-------------------------------------------------------------------------------
def get_Mach2(U, param):
    # calculates squared Mach number
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   param : parameter structure
    # OUTPUT
    #   M2, M2_U : squared Mach number and its linearization w.r.t. U (1x4)
    # DETAILS
    #   Uses constant total enthalpy from param.H0
    #   The speed of sound varies; depends on enthalpy, which depends on speed
    #   The compressible edge speed must be used
  
    if (param.Minf > 0):
        H0, g = param.H0, param.gam
        uk, uk_u = get_uk(U[3], param)
        c2 = (g-1)*(H0-0.5*uk**2); c2_uk = (g-1)*(-uk)  # squared speed of sound
        M2 = uk**2/c2; M2_uk = 2*uk/c2 - M2/c2*c2_uk; M2_U = np.array([0,0,0,M2_uk*uk_u])
    else:
        M2 = 0.; M2_U = np.zeros(4)
  
    return M2, M2_U


#-------------------------------------------------------------------------------
def get_H(U):
    # calculates H = shape parameter = delta*/theta, from U
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   param : parameter structure
    # OUTPUT
    #   H, H_U : shape parameter and its linearization w.r.t. U (1x4)
    # DETAILS
    #   H is the ratio of the displacement thickness to the momentum thickness
    #   In U, the ds entry should be (delta*-wgap) ... i.e wake gap taken out
    #   When the real H is needed with wake gap, Hw is calculated and added
  
    H = U[1]/U[0]
    H_U = np.array([-H/U[0], 1/U[0], 0, 0])
  
    return H, H_U


#-------------------------------------------------------------------------------
def get_Hw(U, wgap):
    # calculates Hw = wake gap shape parameter = wgap/theta
    # INPUT
    #   U    : state vector [th; ds; sa; ue]
    #   wgap : wake gap
    # OUTPUT
    #   Hw, Hw_U : wake gap shape parameter and its linearization w.r.t. U (1x4)
    # DETAILS
    #   Hw is the ratio of the wake gap to the momentum thickness
    #   The wake gap is the TE gap extrapolated into the wake (dead air region)
  
    Hw = wgap/U[0] # wgap/th
    Hw_U = np.array([-Hw/U[0],0,0,0])
  
    return Hw, Hw_U


#-------------------------------------------------------------------------------
def get_Hk(U, param):
    # calculates Hk = kinematic shape parameter, from U
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   param : parameter structure
    # OUTPUT
    #   Hk, Hk_U : kinematic shape parameter and its linearization w.r.t. U (1x4)
    # DETAILS
    #   Hk is like H but with no density in the integrals defining th and ds
    #   So it is exactly the same when density is constant (= freestream)
    #   Here, it is computed from H with a correlation using the Mach number
  
    H, H_U = get_H(U)
  
    if (param.Minf > 0):
        M2, M2_U = get_Mach2(U, param)  # squared edge Mach number
        den = (1+0.113*M2); den_M2 = 0.113
        Hk = (H-0.29*M2)/den
        Hk_U = (H_U-0.29*M2_U)/den - Hk/den*den_M2*M2_U
    else:
        Hk, Hk_U = H, H_U

    return Hk, Hk_U


#-------------------------------------------------------------------------------
def get_Hss(U, param):
    # calculates Hss = density shape parameter, from U
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   param : parameter structure
    # OUTPUT
    #   Hss, Hss_U : density shape parameter and its linearization w.r.t. U (1x4)
    # DETAILS
  
    M2, M2_U = get_Mach2(U, param)  # squared edge Mach number
    Hk, Hk_U = get_Hk(U,param) 
    num = 0.064/(Hk-0.8) + 0.251; num_U = -.064/(Hk-0.8)**2*Hk_U
    Hss = M2*num; Hss_U = M2_U*num + M2*num_U
    
    return Hss, Hss_U


#-------------------------------------------------------------------------------
def get_de(U, param):
    # calculates simplified BL thickness measure
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   param : parameter structure
    # OUTPUT
    #   de, de_U : BL thickness "delta" and its linearization w.r.t. U (1x4)
    # DETAILS
    #   delta is delta* incremented with a weighted momentum thickness, theta
    #   The weight on theta depends on Hk, and there is an overall cap
  
    Hk, Hk_U = get_Hk(U, param)
    aa = 3.15 + 1.72/(Hk-1); aa_U = -1.72/(Hk-1)**2*Hk_U
    de = U[0]*aa + U[1]; de_U = np.array([aa,1,0,0]) + U[0]*aa_U
    dmx = 12.0
    if (de > dmx*U[0]): de, de_U = dmx*U[0], np.array([dmx,0,0,0])
    
    return de, de_U


#-------------------------------------------------------------------------------
def get_rho(U, param):
    # calculates the density (useful if compressible)
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   param : parameter structure
    # OUTPUT
    #   rho, rho_U : density and linearization
    # DETAILS
    #   If compressible, rho is calculated from stag rho + isentropic relations
  
    if (param.Minf > 0):
        M2, M2_U = get_Mach2(U, param) # squared edge Mach number
        uk, uk_u = get_uk(U[3], param) # corrected speed
        H0, gmi = param.H0, param.gam-1
        den = 1+0.5*gmi*M2; den_M2 = 0.5*gmi
        rho = param.rho0/den**(1/gmi); rho_U = (-1/gmi)*rho/den*den_M2*M2_U
    else:
        rho = param.rho0 
        rho_U = np.zeros(4)
  
    return rho, rho_U


#-------------------------------------------------------------------------------
def get_Ret(U, param):
    # calculates theta Reynolds number, Re_theta, from U
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   param : parameter structure
    # OUTPUT
    #   Ret, Ret_U : Reynolds number based on the momentum thickness, linearization
    # DETAILS
    #   Re_theta = rho*ue*theta/mu
    #   If compressible, rho is calculated from stag rho + isentropic relations
    #   ue is the edge speed and must be comressibility corrected
    #   mu is the dynamic viscosity, from Sutherland's law if compressible
  
    if (param.Minf > 0):
        M2, M2_U = get_Mach2(U, param) # squared edge Mach number
        uk, uk_u = get_uk(U[3], param) # corrected speed
        H0, gmi, Ts = param.H0, param.gam-1, param.Tsrat
        Tr = 1-0.5*uk**2/H0; Tr_uk = -uk/H0  # edge/stagnation temperature ratio
        f = Tr**1.5*(1+Ts)/(Tr+Ts); f_Tr = 1.5*f/Tr-f/(Tr+Ts) # Sutherland's ratio
        mu = param.mu0*f; mu_uk = param.mu0*f_Tr*Tr_uk # local dynamic viscosity
        den = 1+0.5*gmi*M2; den_M2 = 0.5*gmi
        rho = param.rho0/den**(1/gmi); rho_U = (-1/gmi)*rho/den*den_M2*M2_U # density
        Ret = rho*uk*U[0]/mu
        Ret_U = rho_U*uk*U[0]/mu + (rho*U[0]/mu-Ret/mu*mu_uk)*np.array([0,0,0,uk_u]) + rho*uk/mu*np.array([1,0,0,0])
    else:
        Ret = param.rho0*U[0]*U[3]/param.mu0
        Ret_U = np.array([U[3], 0, 0, U[0]])/param.mu0
  
    return Ret, Ret_U


#-------------------------------------------------------------------------------
def get_cf(U, param):
    # calculates cf = skin friction coefficient, from U
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   param : parameter structure
    # OUTPUT
    #   cf, cf_U : skin friction coefficient and its linearization w.r.t. U (1x4)
    # DETAILS
    #   cf is the local skin friction coefficient = tau/(0.5*rho*ue^2)
    #   Correlations are used based on Hk and Re_theta

    if (param.wake): return 0, np.zeros(4) # zero cf in wake
    Hk, Hk_U = get_Hk(U, param)
    Ret, Ret_U = get_Ret(U, param)
  
    # TODO: limit Hk

    if (param.turb): # turbulent cf
        M2, M2_U = get_Mach2(U, param)  # squared edge Mach number
        Fc = np.sqrt(1+0.5*(param.gam-1)*M2)
        Fc_U = 0.5/Fc*0.5*(param.gam-1)*M2_U
        aa = -1.33*Hk; aa_U = -1.33*Hk_U
        #if (aa < -20), aa = -20; aa_U = aa_U*0; warning('aa in cfturb'); end
        # smooth limiting of aa
        if (aa < -17): 
            aa = -20+3*np.exp((aa+17)/3)
            aa_U = (aa+20)/3*aa_U  # TODO: ping me  
        bb = np.log(Ret/Fc); bb_U = Ret_U/Ret - Fc_U/Fc
        if (bb < 3): bb, bb_U = 3, bb_U*0
        bb /= np.log(10); bb_U /= np.log(10)
        cc = -1.74 - 0.31*Hk; cc_U = -0.31*Hk_U
        dd = np.tanh(4.0-Hk/0.875); dd_U = (1-dd**2)*(-Hk_U/0.875)
        cf0 = 0.3*np.exp(aa)*bb**cc
        cf0_U = cf0*aa_U + 0.3*np.exp(aa)*cc*bb**(cc-1)*bb_U + cf0*np.log(bb)*cc_U
        cf = (cf0 + 1.1e-4*(dd-1))/Fc
        cf_U = (cf0_U + 1.1e-4*dd_U)/Fc - cf/Fc*Fc_U
    else:  # laminar cf
        if (Hk < 5.5):
            num = .0727*(5.5-Hk)**3/(Hk+1) - .07
            num_Hk = .0727*(3*(5.5-Hk)**2/(Hk+1)*(-1) - (5.5-Hk)**3/(Hk+1)**2)
        else:
            num = .015*(1-1./(Hk-4.5))**2 - .07
            num_Hk = .015*2*(1-1./(Hk-4.5))/(Hk-4.5)**2
        cf = num/Ret
        cf_U = num_Hk/Ret*Hk_U - num/Ret**2*Ret_U

    return cf, cf_U


#-------------------------------------------------------------------------------
def get_cfxt(U, x, param):
    # calculates cf*x/theta from the state
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   x     : distance along wall (xi)
    #   param : parameter structure
    # OUTPUT
    #   cfxt,  : the combination cf*x/theta (calls cf function)
    #   cfxt_U : linearization w.r.t. U (1x4)
    #   cfxt_x : linearization w.r.t x (scalar)  
    # DETAILS
    #   This combination appears in the momentum and shape parameter equations
  
    cf, cf_U = get_cf(U, param)
    cfxt = cf*x/U[0]
    cfxt_U = cf_U*x/U[0]; cfxt_U[0] = cfxt_U[0] - cfxt/U[0]
    cfxt_x = cf/U[0]
  
    return cfxt, cfxt_U, cfxt_x


#-------------------------------------------------------------------------------
def get_cfutstag(Uin, param):
    # calculates cf*ue*theta, used in stagnation station calculations
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   param : parameter structure
    # OUTPUT
    #   F, F_U : value and linearization of cf*ue*theta
    # DETAILS
    #   Only for stagnation and laminar

    U = Uin.copy()
    U[3] = 0
    Hk, Hk_U = get_Hk(U, param)

    if (Hk < 5.5):
        num = .0727*(5.5-Hk)**3/(Hk+1) - .07
        num_Hk = .0727*(3*(5.5-Hk)**2/(Hk+1)*(-1) - (5.5-Hk)**3/(Hk+1)**2)
    else:
        num = .015*(1-1./(Hk-4.5))**2 - .07
        num_Hk = .015*2*(1-1./(Hk-4.5))/(Hk-4.5)**2
    nu = param.mu0/param.rho0
    F = nu*num
    F_U = nu*num_Hk*Hk_U

    return F, F_U


#-------------------------------------------------------------------------------
def get_cdutstag(Uin, param):
    # calculates cDi*ue*theta, used in stagnation station calculations
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   param : parameter structure
    # OUTPUT
    #   D, D_U : value and linearization of cDi*ue*theta
    # DETAILS
    #   Only for stagnation and laminar

    U = Uin.copy()
    U[3] = 0.
    Hk, Hk_U = get_Hk(U, param)
  
    if (Hk<4):
        num = .00205*(4-Hk)**5.5 + .207
        num_Hk = .00205*5.5*(4-Hk)**4.5*(-1)
    else:
        Hk1 = Hk-4;
        num = -.0016*Hk1**2/(1+.02*Hk1**2) + .207
        num_Hk = -.0016*(2*Hk1/(1+.02*Hk1**2) - Hk1**2/(1+.02*Hk1**2)**2*.02*2*Hk1)
  
    nu = param.mu0/param.rho0
    D = nu*num
    D_U = nu*num_Hk*Hk_U

    return D, D_U


#-------------------------------------------------------------------------------
def get_cDixt(U, x, param):
    # calculates cDi*x/theta from the state
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   x     : distance along wall (xi)
    #   param : parameter structure
    # OUTPUT
    #   cDixt,  : the combination cDi*x/theta (calls cDi function)
    #   cDixt_U : linearization w.r.t. U (1x4)
    #   cDixt_x : linearization w.r.t x (scalar)  
    # DETAILS
    #   cDi is the dissipation function
  
    cDi, cDi_U = get_cDi(U, param)
    cDixt = cDi*x/U[0] 
    cDixt_U = cDi_U*x/U[0]; cDixt_U[0] = cDixt_U[0] - cDixt/U[0]
    cDixt_x = cDi/U[0]
 
    return cDixt, cDixt_U, cDixt_x
 

#-------------------------------------------------------------------------------
def get_cDi(U, param):
    # calculates cDi = dissipation function = 2*cD/H*, from the state
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   param : parameter structure
    # OUTPUT
    #   cDi, cDi_U : dissipation function and its linearization w.r.t. U (1x4)
    # DETAILS
    #   cD is the dissipation coefficient, int(tau*du/dn*dn)/(rho*ue^3)
    #   The combination with H* appears in the shape parameter equation
  
    if (param.turb):  # turbulent includes wake
    
        # initialize to 0; will add components that are needed
        cDi, cDi_U = 0, np.zeros(4)
    
        if (not param.wake):
            # turbulent wall contribution (0 in the wake) 
            cDi0, cDi0_U = get_cDi_turbwall(U, param)
            cDi = cDi + cDi0; cDi_U = cDi_U + cDi0_U
            cDil, cDil_U = get_cDi_lam(U, param) # for max check
        else:
            cDil, cDil_U = get_cDi_lamwake(U, param) # for max check
    
        # outer layer contribution
        cDi0, cDi0_U = get_cDi_outer(U, param)
        cDi = cDi + cDi0; cDi_U = cDi_U + cDi0_U
    
        # laminar stress contribution
        cDi0, cDi0_U = get_cDi_lamstress(U, param)
        cDi = cDi + cDi0; cDi_U = cDi_U + cDi0_U
        
        # maximum check
        if (cDil > cDi): cDi, cDi_U = cDil, cDil_U
    
        # double dissipation in the wake
        if (param.wake): cDi, cDi_U = 2*cDi, 2*cDi_U
    else:
        # just laminar dissipation
        [cDi, cDi_U] = get_cDi_lam(U, param)
  
    return cDi, cDi_U



#-------------------------------------------------------------------------------
def get_cDi_turbwall(U, param):
    # calculates the turbulent wall contribution to cDi
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   param : parameter structure
    # OUTPUT
    #   cDi, cDi_U : dissipation function and its linearization w.r.t. U (1x4)
    # DETAILS
    #   This is one contribution to the dissipation function cDi = 2*cD/H*
  
    if (param.wake): return 0, np.zeros(4)
  
    # get cf, Hk, Hs, Us
    cf, cf_U = get_cf(U, param)
    Hk, Hk_U = get_Hk(U, param)
    Hs, Hs_U = get_Hs(U, param)
    Us, Us_U = get_Us(U, param)
    Ret, Ret_U = get_Ret(U, param)
  
    lr = np.log(Ret); lr_U = Ret_U/Ret
    Hmin = 1 + 2.1/lr; Hmin_U = -2.1/lr**2*lr_U
    aa = np.tanh((Hk-1)/(Hmin-1)); fac = 0.5 + 0.5*aa
    fac_U = 0.5*(1-aa**2)*(Hk_U/(Hmin-1)-(Hk-1)/(Hmin-1)**2*Hmin_U)

    cDi = 0.5*cf*Us*(2/Hs)*fac
    cDi_U = cf_U*Us/Hs*fac + cf*Us_U/Hs*fac - cDi/Hs*Hs_U + cf*Us/Hs*fac_U

    return cDi, cDi_U


#-------------------------------------------------------------------------------
def get_cDi_lam(U, param):
    # calculates the laminar dissipation function cDi
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   param : parameter structure
    # OUTPUT
    #   cDi, cDi_U : dissipation function and its linearization w.r.t. U (1x4)
    # DETAILS
    #   This is one contribution to the dissipation function cDi = 2*cD/H*
  
    # first get Hk and Ret
    Hk, Hk_U = get_Hk(U, param)
    Ret, Ret_U = get_Ret(U, param)
  
    if (Hk<4):
        num = .00205*(4-Hk)**5.5 + .207
        num_Hk = .00205*5.5*(4-Hk)**4.5*(-1)
    else:
        Hk1 = Hk-4
        num = -.0016*Hk1**2/(1+.02*Hk1**2) + .207
        num_Hk = -.0016*(2*Hk1/(1+.02*Hk1**2) - Hk1**2/(1+.02*Hk1**2)**2*.02*2*Hk1)
    
    cDi = num/Ret;
    cDi_U = num_Hk/Ret*Hk_U - num/Ret**2*Ret_U
  
    return cDi, cDi_U


#-------------------------------------------------------------------------------
def get_cDi_lamwake(U, paramin):
    # laminar wake dissipation function cDi
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   param : parameter structure
    # OUTPUT
    #   cDi, cDi_U : dissipation function and its linearization w.r.t. U (1x4)
    # DETAILS
    #   This is one contribution to the dissipation function cDi = 2*cD/H*
    
    param = copy.deepcopy(paramin)
    param.turb = False  # force laminar
  
    # dependencies
    Hk, Hk_U = get_Hk(U, param)
    Hs, Hs_U = get_Hs(U, param)
    Ret, Ret_U = get_Ret(U, param)
    HsRet = Hs*Ret
    HsRet_U = Hs_U*Ret + Hs*Ret_U
  
    num = 2*1.1*(1-1/Hk)**2*(1/Hk)
    num_Hk = 2*1.1*(2*(1-1/Hk)*(1/Hk**2)*(1/Hk)+(1-1/Hk)**2*(-1/Hk**2))
    cDi = num/HsRet
    cDi_U = num_Hk*Hk_U/HsRet - num/HsRet**2*HsRet_U
      
    return cDi, cDi_U


#-------------------------------------------------------------------------------
def get_cDi_outer(U, param):
    # turbulent outer layer contribution to dissipation function cDi
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   param : parameter structure
    # OUTPUT
    #   cDi, cDi_U : dissipation function and its linearization w.r.t. U (1x4)
    # DETAILS
    #   This is one contribution to the dissipation function cDi = 2*cD/H*
  
    if (not param.turb): return 0, np.zeros(4) # for pinging
  
    # first get Hs, Us
    [Hs, Hs_U] = get_Hs(U, param)
    [Us, Us_U] = get_Us(U, param)

    # shear stress: note, state stores ct^.5
    ct = U[2]**2; ct_U = np.array([0,0,2*U[2],0])
  
    cDi = ct*(0.995-Us)*2/Hs
    cDi_U = ct_U*(0.995-Us)*2/Hs + ct*(-Us_U)*2/Hs - ct*(0.995-Us)*2/Hs**2*Hs_U

    return cDi, cDi_U


#-------------------------------------------------------------------------------
def get_cDi_lamstress(U, param):
    # laminar stress contribution to dissipation function cDi
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   param : parameter structure
    # OUTPUT
    #   cDi, cDi_U : dissipation function and its linearization w.r.t. U (1x4)
    # DETAILS
    #   This is one contribution to the dissipation function cDi = 2*cD/H*
  
    # first get Hs, Us, and Ret
    Hs, Hs_U = get_Hs(U, param)
    Us, Us_U = get_Us(U, param)
    Ret, Ret_U = get_Ret(U, param)
    HsRet = Hs*Ret
    HsRet_U = Hs_U*Ret + Hs*Ret_U
  
    num = 0.15*(0.995-Us)**2*2
    num_Us = 0.15*2*(0.995-Us)*(-1)*2
    cDi = num/HsRet
    cDi_U = num_Us*Us_U/HsRet - num/HsRet**2*HsRet_U
  
    return cDi, cDi_U


#-------------------------------------------------------------------------------
def get_Us(U, param):
    # calculates the normalized wall slip velocity Us
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   param : parameter structure
    # OUTPUT
    #   Us, Us_U : normalized wall slip velocity and its linearization w.r.t. U (1x4)
  
    [Hs, Hs_U] = get_Hs(U, param)
    [Hk, Hk_U] = get_Hk(U, param)
    [H, H_U] = get_H(U)
  
    # limit Hk (TODO smooth/eliminate)
    if (param.wake) and (Hk < 1.00005): Hk, Hk_U = 1.00005, Hk_U*0
    if (not param.wake) and (Hk < 1.05): Hk, Hk_U = 1.05, Hk_U*0
  
    beta = param.GB; bi = 1./beta
    Us = 0.5*Hs*(1-bi*(Hk-1)/H)
    Us_U = 0.5*Hs_U*(1-bi*(Hk-1)/H) + 0.5*Hs*(-bi*(Hk_U)/H +bi*(Hk-1)/H**2*H_U)
    # limits
    if ((not param.wake) and (Us>0.95   )): Us, Us_U = 0.98, Us_U*0
    if ((not param.wake) and (Us>0.99995)): Us, Us_U = 0.99995, Us_U*0
  
    return Us, Us_U



#-------------------------------------------------------------------------------
def get_damp(U, param):
    # calculates the amplification rate, dn/dx, used in predicting transition
    # INPUT
    #   U     : state vector [th; ds; sa; ue]
    #   param : parameter structure
    # OUTPUT
    #   damp, damp_U : amplification rate and its linearization w.r.t. U (1x4)
    # DETAILS
    #   damp = dn/dx is used in the amplification equation, prior to transition
  
    [Hk, Hk_U] = get_Hk(U, param)
    [Ret, Ret_U] = get_Ret(U, param)
    th = U[0]
    
    # limit Hk (TODO smooth/eliminate)
    if (Hk < 1.05): Hk, Hk_U = 1.05, Hk_U*0
  
    Hmi = 1./(Hk-1); Hmi_U = -Hmi**2*Hk_U
    aa = 2.492*Hmi**0.43; aa_U = 0.43*aa/Hmi*Hmi_U
    bb = np.tanh(14*Hmi-9.24); bb_U = (1-bb**2)*14*Hmi_U
    lrc = aa + 0.7*(bb+1); lrc_U = aa_U + 0.7*bb_U
    lten = np.log(10); lr = np.log(Ret)/lten; lr_U = (1/Ret)*Ret_U/lten
    dl = .1;  # changed from .08 to make smoother
    damp = 0; damp_U = np.zeros(len(U))  # default no amplification
    if (lr >= lrc-dl):
        rn = (lr-(lrc-dl))/(2*dl); rn_U = (lr_U - lrc_U)/(2*dl)
        if (rn >= 1):
            rf = 1; rf_U = np.zeros(len(U))
        else:
            rf = 3*rn**2-2*rn**3; rf_U = (6*rn-6*rn**2)*rn_U
        ar = 3.87*Hmi-2.52; ar_U = 3.87*Hmi_U
        ex = np.exp(-ar**2); ex_U = ex*(-2*ar*ar_U)
        da = 0.028*(Hk-1)-0.0345*ex; da_U = 0.028*Hk_U-0.0345*ex_U
        af = -0.05+2.7*Hmi-5.5*Hmi**2+3*Hmi**3+0.1*np.exp(-20*Hmi)
        af_U = (2.7-11*Hmi+9*Hmi**2-1*np.exp(-20*Hmi))*Hmi_U
        damp = rf*af*da/th
        damp_U = (rf_U*af*da + rf*af_U*da + rf*af*da_U)/th - damp/th*np.array([1,0,0,0])
    
    # extra amplification to ensure dn/dx > 0 near ncrit
    ncrit = param.ncrit
  
    Cea = 5; nx = Cea*(U[2]-ncrit); nx_U = Cea*np.array([0,0,1,0])
    eex = 1+np.tanh(nx); eex_U = (1-np.tanh(nx)**2)*nx_U
  
    ed = eex*.001/th
    ed_U = eex_U*.001/th - ed/th*np.array([1,0,0,0])
    damp = damp + ed
    damp_U = damp_U + ed_U
    
    return damp, damp_U


#-------------------------------------------------------------------------------
def check_ping(ep, v, v_u, sname):
# checks convergence of 3 values/derivatives
# INPUT
#   v     : list of three function evaluations at 0,+ep,+2*ep
#   v_u   : list of three derivative evaluations at 0,+ep,+2*ep
#   sname : descriptive name of where values came from for printing
# OUTPUT
#   E     : error values for two finite-difference comparisons
#   rate  : convergence rate, also printed
  
    E = np.zeros(2)
    for i in range(2): E[i] = norm2((v[1+i]-v[0])/(ep*(i+1.)) - 0.5*(v_u[0] + v_u[1+i]))
    rate = np.log2(E[1]/E[0])
    print('%s ping error convergence rate = %.4f'%(sname, rate))
    return E, rate

  
#-------------------------------------------------------------------------------
def ping_test(M):
# checks derivatives of various functions through finite-difference pinging
# INPUT
#   M : mfoil class
# OUTPUT
#   printouts of rates (2 = second order expected).
  
    M.oper.alpha = 3 # angle, in degrees
    M.oper.Ma = 0.4 # Mach number
    M.oper.viscous = True  # tests are viscous
    np.random.seed(17) # for consistent pseudo random numbers
    M.param.verb = 2  # to minimize prints to screen
  
    # freestream Reynolds numbers
    Rev = np.r_[2e3, 1e5]
  
    # laminar/turbulent test states: th, ds, sa, ue
    Uv = [np.r_[0.01, 0.02, 8.4, 0.9], np.r_[0.023, 0.05, .031, 1.1]]
  
    # functions to test
    fv = [get_Hk, get_Ret, get_cf, get_cDi, get_Hs, get_Us, \
          get_cDi_turbwall, get_cDi_lam, get_cDi_lamwake, get_cDi_outer, \
          get_cDi_lamstress, get_cteq, get_cttr, get_de, get_damp, \
          get_Mach2, get_Hss, residual_station]
  
    # ping tests
    sturb = ['lam', 'turb', 'wake']
    for iRe in range(len(Rev)):  # loop over Reynolds numbers
        M.oper.Re = Rev[iRe]
        init_thermo(M)
        param = build_param(M, 1)
        for it in range(3):  # loop over lam, turb, wake
            param.turb, param.wake = (it>0), (it==2)
            for ih in range(len(fv)): # loop over functions
                U, srates, smark, serr, f = Uv[min(it,1)], '', '', '', fv[ih]
                if (f==residual_station): U = np.concatenate((U,U*np.r_[1.1,.8,.9,1.2]))
                for k in range(len(U)): # test all state component derivatives
                    ep, E = 1e-2*U[k], np.zeros(2)
                    if (f==residual_station):
                        xi, Aux, dx = np.r_[0.7,0.8], np.r_[.002, .0018], np.r_[-.2,.3]
                        v0, v_U0, v_x0 = f(param,xi,np.stack((U[0:4],U[4:8]),axis=-1),Aux)
                        for iep in range(2): # test with two epsilons
                            U[k] += ep; xi += ep*dx
                            v1, v_U1, v_x1 = f(param,xi,np.stack((U[0:4],U[4:8]),axis=-1),Aux)
                            U[k] -= ep; xi -= ep*dx
                            E[iep] = norm2((v1-v0)/ep - 0.5*(v_U1[:,k] + v_U0[:,k] + np.dot(v_x0+v_x1, dx)))
                            ep /= 2
                    else:
                        [v0, v_U0] = f(U, param)
                        for iep in range(2): # test with two epsilons
                            U[k] += ep; v1, v_U1 = f(U, param); U[k] -= ep
                            E[iep] = abs((v1-v0)/ep - 0.5*(v_U1[k] + v_U0[k]))
                            ep /= 2
                    srate = ' N/A'; skip = False
                    if (not skip) and (E[0]>5e-11) and (E[1]>5e-11):
                        m = np.log2(E[0]/E[1]); srate = '%4.1f'%(m)
                        if (m<1.5): smark = '<==='
                    srates += ' ' + srate
                    serr += ' %.2e->%.2e'%(E[0],E[1])
                vprint(param, 0, '%-18s %-5s err=[%s]  rates=[%s] %s'%(f.__name__, sturb[it], serr, srates, smark))
  
    # transition residual ping
    M.oper.Re = 2e6; init_thermo(M); param = build_param(M,1)
    U, x, Aux = np.transpose(np.array([[0.01, 0.02, 8.95, 0.9], [0.013, 0.023, .028, 0.85]])), np.r_[0.7,0.8], np.r_[0,0]
    dU, dx, ep, v, v_u = np.random.rand(4,2), np.random.rand(2), 1e-4, [], []
    for ie in range(3):
        R, R_U, R_x = residual_transition(M, param, x, U, Aux)
        v.append(R); v_u.append(np.dot(R_U, np.reshape(dU,8,order='F')) + np.dot(R_x,dx))
        U += ep*dU; x += ep*dx
    check_ping(ep, v, v_u, 'transition residual');  

    # stagnation residual ping
    M.oper.Re = 1e6; M.oper.alpha = 1; init_thermo(M); param = build_param(M,1)
    U, x, Aux = np.array([0.00004673616, 0.000104289, 0, 0.11977917547]), 4.590816441485401e-05, [0,0]
    dU, dx, ep, v, v_u = np.random.rand(4), np.random.rand(1), 1e-6, [], []
    for ie in range(3):
        param.simi = True
        R, R_U, R_x = residual_station(param, np.r_[x,x], np.stack((U,U),axis=-1), Aux)
        param.simi = False
        v.append(R); v_u.append(np.dot(R_U[:,range(0,4)] + R_U[:,range(4,8)], dU) + (R_x[:,0] + R_x[:,1])*dx[0])
        U += ep*dU; x += ep*dx[0]
    check_ping(ep, v, v_u, 'stagnation residual')
    
    # need a viscous solution for the next tests
    solve_viscous(M)
    # entire system ping
    #M.param.niglob = 10
    Nsys = M.glob.U.shape[1]
    dU, dx, ep = np.random.rand(4,Nsys), 0.1*np.random.rand(Nsys), 1e-6
    for ix in range(2): # ping with explicit and implicit (baked-in) R_x effects
        if (ix==1): 
            dx *= 0
            stagpoint_move(M) # baked-in check
        v, v_u = [], []
        for ie in range(3):
            build_glob_sys(M)
            v.append(M.glob.R.copy()); v_u.append(M.glob.R_U @ np.reshape(dU,4*Nsys,order='F') + M.glob.R_x @ dx)
            M.glob.U += ep*dU; M.isol.xi += ep*dx 
            if (ix==1): stagpoint_move(M) # baked-in check: stagnation point moves
        M.glob.U -= 3*ep*dU; M.isol.xi -= 3*ep*dx
        check_ping(ep, v, v_u, 'global system, ix=%d'%(ix))

    # wake system ping  
    dU, ep, v, v_u = np.random.rand(M.glob.U.shape[0], M.glob.U.shape[1]), 1e-5, [], []
    for ie in range(3):
        R, R_U, J = wake_sys(M, param)
        v.append(R); v_u.append(np.dot(R_U, np.reshape(dU[:,J],4*len(J),order='F')))
        M.glob.U += ep*dU
    M.glob.U -= 2*ep*dU
    check_ping(ep, v, v_u, 'wake system')
  
    # stagnation state ping
    M.oper.Re = 5e5; init_thermo(M); param = build_param(M,1)
    U, x, = np.transpose(np.array([[5e-5, 1.1e-4, 0, .0348], [4.9e-5, 1.09e-4, 0, .07397]])), np.r_[5.18e-4, 1.1e-3]
    dU, dx, ep, v, v_u = np.random.rand(4, 2), np.random.rand(2), 1e-6, [], []
    for ie in range(3):
        Ust, Ust_U, Ust_x, xst = stagnation_state(U, x)
        v.append(Ust); v_u.append(np.dot(Ust_U, np.reshape(dU,8,order='F')) + np.dot(Ust_x,dx))
        U += ep*dU; x += ep*dx
    check_ping(ep, v, v_u, 'stagnation state')
  
    # force calculation ping
    Nsys, N, v, v_u = M.glob.U.shape[1], M.foil.N, [], []
    due = np.random.rand(N); dU = np.zeros((4,Nsys)); dU[3,0:N] = due; da = 10; ep = 1e-2
    for ie in range(3):
        calc_force(M)
        v.append(np.array([M.post.cl])); v_u.append(np.array([np.dot(M.post.cl_ue, due) + M.post.cl_alpha*da]))
        M.glob.U += ep*dU; M.oper.alpha += ep*da
    M.glob.U -= 3*ep*dU; M.oper.alpha -= 3*ep*da
    check_ping(ep, v, v_u, 'lift calculation');


#-------------------------------------------------------------------------------
def main():
    # make a NACA 2412 airfoil
    m = mfoil(naca='0012', npanel=199)
    print('NACA geom name =', m.geom.name, '  num. panels =', m.foil.N)
    # add a flap
    #m.geom_flap(np.r_[0.8,0],5)
    # derotate the geometry
    #m.geom_derotate()
    # add camber
    #m.geom_addcamber(np.array([[0,0.3,0.7,1],[0,-.03,.01,0]]))
    # set up a compressible viscous run (note, alpha is in degrees)
    m.setoper(alpha=2, Re=5*10**6, Ma=0.0)
    # request plotting, specify the output file for the plot
    m.param.doplot, m.post.rfile = True, 'results.pdf'
    # set the verbosity (higher -> more output to stdout)
    m.param.verb = 1
    # run the solver
    print('Running the solver.')
    m.solve()
    # run the derivative ping check
    #print('Derivative ping check.')
    #m.ping()


if __name__ == "__main__":
    main()
