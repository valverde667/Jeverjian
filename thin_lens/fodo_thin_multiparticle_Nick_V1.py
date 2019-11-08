import numpy as np
import matplotlib.pyplot as plt
# Program for thin lens FODO quad transport

# Unit conversions
mm = 1.e-3 # unit mm
mr = 1.e-3 # unit mrad

# Lattice parameters
d = 1.  # half drift length [m]
f = 4.      # focal length [m]
#f = 0.9999  # focal length [m] unstable value ... just beyond stability

Lp = 4*d  # Lattice period [m]


#n_advance = int(input('How many times through the FODO Cell: ')) # lattice periods to advance
#n_particles = int(input('How many particles: ')) # number of particles to send through
n_advance = 100
n_particles = 35
# Setup lattice elements
M0 = np.array( [[1.,d],[0,1.]] )         # drift matrix
MF = np.array( [[1.,0],[-1./f,1.]] )     # focus matrix
MD = np.array( [[1.,0],[ 1./f,1.]] )     # defocus matrix

# Advance xv through element
def xv_adv(M,xv):
    return( np.dot(M,xv) )


fig, axes = plt.subplots(nrows = 2, ncols = 1, figsize = (8,8), sharex = True) #Initiate Plot
ax_pos = axes[0] #Name x_history plot object
ax_phase_space = axes[1] #Name phase-space plot object
for particles in range(n_particles):

    # Particle initial conditions
    xi  = np.random.uniform(-1,1)*mm  # inital x-coordinate
    xpi = np.random.uniform(0.1,0.6)*mr  # inital xp_history
    s  = 0.     # inital axial coordinate


    #Initial array for particle
    xv   = np.array( [xi,xpi] ) # coordinate vector set by initial condition
    x_history = [xv[0]]  # List of x_history coordinates to be updated
    xp_history = [xv[1]]     # List of xp_history coordinates to be updated
    s_history = [s]        # List of lognitudinal coordinates to be updated

    #Loop for going through FODO cell
    for turn in range(n_advance):

        #First half of drift
        xv = xv_adv(M0, xv)
        s += d
        x_history.append(xv[0])
        xp_history.append(xv[1])
        s_history.append(s)

        #Focus quadrupole
        xv = xv_adv(MF, xv)
        x_history.append(xv[0])
        xp_history.append(xv[1])
        s_history.append(s) #keep sizes the same

        # two half drifts between quadrupoles
        xv = xv_adv(M0,xv)
        xv = xv_adv(M0,xv)
        s += 2.*d
        x_history.append(xv[0])
        xp_history.append(xv[1])
        s_history.append(s)

        # Defocus quadrupole
        xv = xv_adv(MD,xv)
        x_history.append(xv[0])
        xp_history.append(xv[1])
        s_history.append(s)

        #Final Half drift
        xv = xv_adv(M0, xv)
        s += d
        x_history.append(xv[0])
        xp_history.append(xv[1])
        s_history.append(s)

    #Scale lists and turn them into arrays
    x_history = np.array(x_history)/mm
    xp_history = np.array(xp_history)/mr
    s_history = np.array(s_history)/Lp

    #Plot particle's transverse x_history
    ax_pos.plot(s_history, x_history, color = 'k')

    #Plot particle's phase space
    ax_phase_space.plot(s_history, xp_history, color = 'k')


#Label x_history plot
ax_pos.set_title(r"$x$-$x'$ Phase-Space Evolution, %g FODO Cells and %g Particles" %(n_advance, n_particles), fontsize = 16)
ax_pos.set_ylabel(r"Position, $x$ [mm]", fontsize = 14)

#Label Phase-space plot
ax_phase_space.set_ylabel(r"Angle, $x'$[mr]", fontsize = 14)

ax_phase_space.set_xlabel(r"Axial Coordinate, $s/L_p$ [Lattice Periods]", fontsize = 14)


plt.tight_layout()
plt.savefig('Fodo_example_multiparticles', dpi = 300)
