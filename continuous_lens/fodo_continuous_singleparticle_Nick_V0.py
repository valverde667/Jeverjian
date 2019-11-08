import numpy as np
import matplotlib.pyplot as plt


# Program for thin lens FODO quad transport
#--Add Mod1:  Add code to make continuous focusing version of program.
#--In this case the lattice period is just a length scale and set k_beta0 = sigma_0/Lp where
#--sigma_0 is an input “phase-advance” in radians per lattice period.

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



# Particle initial conditions
xi  = 1.0  # inital x-coordinate
xpi = 0.  # inital xp_history
s_pos = np.array([i for i in range(0, 11)])/Lp #create array of s positions


#Initial array for particle
xv   = np.array( [xi,xpi] ) # coordinate vector set by initial condition
x_history = [xv[0]]  # List of x_history coordinates to be updated
xp_history = [xv[1]]     # List of xp_history coordinates to be updated
s_pos = [i/Lp for i in range(0, 90)] #create array of s position

#Create Plots
fig, axes = plt.subplots(nrows = 2, ncols = 1, figsize=(8,8), sharex = True)
ax_pos = axes[0]
ax_phase_space = axes[1]

#Main for-loop
for s in range(len(s_pos)-1):
    M = transfer_matrix(kbeta_0,s_pos[s+1], s_pos[s]) #Create transfer matrix
    xv = xv_adv(M, xv) #advance particle
    x_history.append(xv[0]) #append new particle position to history
    xp_history.append(xv[1]) #append new particle phase to history


#Add to plots
ax_pos.plot(s_pos, x_history) #position plot
ax_pos.axhline(y = 0, alpha = 0.5, color = 'k', linewidth = 0.5)
ax_pos.set_ylabel(r"Transverse Position $x$ [mm]", fontsize = 14)

ax_phase_space.plot(s_pos, xp_history) #phase-space plot
ax_phase_space.axhline(y = 0, alpha = 0.5, color = 'k', linewidth = 0.5)
ax_phase_space.set_xlabel(r"Longitudinal Position $s/L_p$", fontsize = 14)
ax_phase_space.set_ylabel(r"Phase Coordinate $x'$ [mrad]", fontsize = 14)
