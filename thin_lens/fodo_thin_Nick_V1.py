# Program for thin lens FODO quad transport

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

# Particle initial conditions and advance control
xi  = 1.*mm  # inital x-coordinate
xpi = 0.0*mr  # inital angle
si  = 0.     # inital axial coordinate

n_advance = 75 # lattice periods to advance

# Setup coordinate array
#   xv = [x,x'] = [x-position,x-angle] == [x,xp]
#    s = Axial coordinate
xv  = np.array( [xi,xpi] ) # coordinate vector set by initial condition
s   = si                    # axial coordinate vector

# Setup lattice elements
M0 = np.array( [[1.,d],[0,1.]] )         # drift matrix
MF = np.array( [[1.,0],[-1./f,1.]] )     # focus matrix
MD = np.array( [[1.,0],[ 1./f,1.]] )     # defocus matrix

# Advance xv through element
def xv_adv(M,xv):
    return( np.dot(M,xv) )

# Initialize history "_h" arrays that will be used to plot x and x' evolution
x_history  = [xv[0]]
xp_history = [xv[1]]
s_history  = [s]

for turn in range(n_advance):

    # First half-drift
    # --- advance coordinates
    xv = xv_adv(M0, xv)
    s += d
    # -- save history appending arrays
    x_history.append(xv[0])
    xp_history.append(xv[1])
    s_history.append(s)

    # Focus quadrupole
    # --- advance coordinates
    xv = xv_adv(MF, xv)
    # --- save history
    x_history.append(xv[0])
    xp_history.append(xv[1])
    s_history.append(s) #keep sizes the same, thin kick does not advance s

    # two half drifts between quadrupoles
    # --- advance coordinates
    xv = xv_adv(M0,xv)
    xv = xv_adv(M0,xv)
    s += 2.*d
    # --- save history
    x_history.append(xv[0])
    xp_history.append(xv[1])
    s_history.append(s)

    # Defocus quadrupole
    # --- advance coordinates
    xv = xv_adv(MD,xv)
    # --- save history
    x_history.append(xv[0])
    xp_history.append(xv[1])
    s_history.append(s)

    # Final Half drift
    # --- advance coordinates
    xv = xv_adv(M0, xv)
    s += d
    # --- save history
    x_history.append(xv[0])
    xp_history.append(xv[1])
    s_history.append(s)

# Convert to arrays
x_history  = np.array(x_history)
xp_history = np.array(xp_history)
s_history  = np.array(s_history)

fig, axes = plt.subplots(nrows = 2, ncols = 1, figsize = (8,8), sharex = True)
ax_pos = axes[0]
ax_phase_space = axes[1]

ax_pos.plot(s_history/Lp, x_history/mm)
ax_pos.set_title(r"$x$-$x'$ Phase-Space Evolution, %g FODO CElls" %n_advance, fontsize = 16)
ax_pos.set_ylabel(r'Position, $x$ [mm]', fontsize = 14)

ax_phase_space.plot(s_history/Lp, xp_history/mr)
#ax_phase_space.set_title('Phase Space Evolution through %g FODO Cells' %n_advance)
ax_phase_space.set_ylabel(r"Angle, $x'$ [mr]", fontsize = 14)
ax_phase_space.set_xlabel(r"Axial Coordinate, $s/L_p$ [Lattice Periods]", fontsize = 14)


plt.tight_layout()
plt.savefig('fodo_thin_xxpps', dpi = 300)
