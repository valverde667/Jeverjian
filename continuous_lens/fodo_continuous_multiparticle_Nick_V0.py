import numpy as np
import matplotlib.pyplot as plt

# Unit conversions
mm = 1.e-3 # unit mm
mr = 1.e-3 # unit mrad

# Lattice parameters
d = 1.  # half drift length [m]
f = 4.      # focal length [m]
#f = 0.9999  # focal length [m] unstable value ... just beyond stability

Lp = 4*d  # Lattice period [m]
sigma_0 = np.pi/3 #phase advance [rad]
kbeta_0 = sigma_0/Lp #[rad/Lp]

# Setup Transfer Matrix
def transfer_matrix(kbeta_0, s, si):
    """Create Transfer matrix"""
    #create individual elements for matrix
    m11 = np.cos(np.sqrt(kbeta_0)*(s - si))
    m12 = np.sin(np.sqrt(kbeta_0)*(s-si))/np.sqrt(kbeta_0)
    m21 = -np.sin(np.sqrt(kbeta_0)*(s-si))*np.sqrt(kbeta_0)
    m22 = np.cos(np.sqrt(kbeta_0)*(s - si))
    matrix = np.array([[m11, m12], [m21, m22]])

    return matrix

# Advance xv through element
def xv_adv(M,xv):
    return( np.dot(M,xv) )

#Create Plots
fig, axes = plt.subplots(nrows = 2, ncols = 1, figsize=(8,8), sharex = True)
ax_pos = axes[0]
ax_phase_space = axes[1]

#n_particles = int(input("Enter number of particles: ")) #User input number of particles
n_particles = 25
#Main loop
for particles in range(n_particles):
    # Particle initial conditions
    xi  = np.random.uniform(-1,1)*mm  # inital x-coordinate
    xpi = np.random.uniform(0.1,0.6)*mr  # inital xp_history

    #Initialize arrays lists
    xv = np.array([xi,xpi]) #x-vector containing position and phase
    x_history = [xv[0]]  # List of x_history coordinates to be updated
    xp_history = [xv[1]]     # List of xp_history coordinates to be updated
    s_pos = [i/Lp for i in range(0, 90)] #create array of s position

    #Advance particles through lattice
    for s in range(len(s_pos)-1):
        M = transfer_matrix(kbeta_0,s_pos[s+1], s_pos[s]) #Create transfer matrix
        xv = xv_adv(M, xv) #advance particle
        x_history.append(xv[0]) #append new particle position to history
        xp_history.append(xv[1]) #append new particle phase to history

    x_history = np.array(x_history) #vectorize x_history
    xp_history = np.array(xp_history) #vecotrize xp_history

    #Plotting
    ax_pos.plot(s_pos, x_history/mm) #position plot
    ax_phase_space.plot(s_pos, xp_history/mr) #phase-space plot



#Label Plots
ax_pos.set_title(r"$x$-$x'$ Phase-Space Evolution, Continuous FODO Cell For %g Particles" %(n_particles), fontsize = 14)
ax_pos.axhline(y = 0, alpha = 0.5, color = 'k', linewidth = 0.5)
ax_pos.set_ylabel(r"Transverse Position $x$ [mm]", fontsize = 14)

ax_phase_space.axhline(y = 0, alpha = 0.5, color = 'k', linewidth = 0.5)
ax_phase_space.set_xlabel(r"Longitudinal Position $s/L_p$", fontsize = 14)
ax_phase_space.set_ylabel(r"Phase Coordinate $x'$ [mrad]", fontsize = 14)

plt.tight_layout()
plt.savefig('Fodo_continous_multiparticle', dpi = 300)
