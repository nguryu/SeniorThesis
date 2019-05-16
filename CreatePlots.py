from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
import time
start = time.time()

# Initialize lists to store data in.
times = []
x = []
y = []
z = []
r = []
x_hold = []
y_hold = []
z_hold = []
timestamp = []
mass = []
x_pos = []
y_pos = []
z_pos = []
vx = []
vy = []
vz = []
vx_gr = []
vy_gr = []
vz_gr = []
epsilon_gr = []
minus_ut = []
minus_uth = []
x_inter = []
y_inter = []
z_inter = []
t_mid = []
x_actual = []
y_actual = []
z_actual = []
x_coord = []
y_coord = []
z_coord = []
relCorrection = []
x_axis = []

# Read in data from files.
with open('ParticleTrajectories.dat', 'r') as data_file:
    for line in data_file:
        line = line.strip()  # Remove whitespace.
        if not line:  # Skip empty lines.
            continue
        if not line.startswith('#'):
            col = line.split(' ')
            times.append(float(col[0]))
            x.append(float(col[1]))
            y.append(float(col[2]))
            z.append(float(col[3]))

with open('AllParticlesFromOutflowLev0.dat', 'r') as data_file2:
    # data_slice = [next(data_file2) for x in xrange(6000)]  # Read in file up to line N.
    for line in data_file2:
        line = line.strip()  # Remove whitespace.
        if not line:  # Skip empty lines.
            continue
        if not line.startswith('#'):
            col = line.split(' ')
            timestamp.append(float(col[0]))
            mass.append(float(col[1]))
            x_pos.append(float(col[2]))
            y_pos.append(float(col[3]))
            z_pos.append(float(col[4]))
            minus_ut.append(float(col[8]))
            minus_uth.append(float(col[9]))
            vx.append(float(col[10]))
            vy.append(float(col[11]))
            vz.append(float(col[12]))

def boundCond(x_pos, y_pos, z_pos, vx, vy, vz, minus_ut):
    for i in range(len(timestamp)):
        r = float(sqrt(x_pos[i] ** 2 + y_pos[i] ** 2 + z_pos[i] ** 2))
        v = sqrt(vx[i] ** 2 + vy[i] ** 2 + vz[i] ** 2)
        A = sqrt(2 * (minus_ut[i] + 2.6/r - 1) / v**2)  # Transformation term
        if (2 * (minus_ut[i] + 2.6/r - 1) / v**2) < 1:
            A = 0
        vx_gr.append(A * vx[i])
        vy_gr.append(A * vy[i])
        vz_gr.append(A * vz[i])
        gr_velocity = sqrt(vx_gr[i]**2 + vy_gr[i]**2 + vz_gr[i]**2)
        epsilon_gr.append(1 + (gr_velocity ** 2) / 2 - 2.6 / r)
    return epsilon_gr

epsilon_gr = boundCond(x_pos, y_pos, z_pos, vx, vy, vz, minus_ut)

# ======================================================================= #
# Plot Particles
# ======================================================================= #

# Figure parameters
SMALL_SIZE = 8
MEDIUM_SIZE = 12
LARGE_SIZE = 16
plt.rc('font', size = SMALL_SIZE)          # Controls default text sizes
plt.rc('axes', titlesize = LARGE_SIZE)     # Font size of the axes title
plt.rc('axes', labelsize = SMALL_SIZE)     # Font size of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)    # Font size of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)    # Font size of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)    # Legend font size
# Axes ticks
majorFormatter = FormatStrFormatter('%1.1f')
majorLocator = MultipleLocator(1000)  # Major tick intervals
minorLocator = MultipleLocator(250)  # Minor tick intervals

# Figure set-up
fig1 = plt.figure(1)
# ax1 = fig1.add_subplot(111, projection = '3d')
ax1 = fig1.add_subplot(111)
ax1.set_xlabel('X axis')
ax1.set_ylabel('Y axis')
# ax1.set_zlabel('Z axis')

# Create new list of coordinates for simple plotting purposes.
n = 1  # Choose the first (2800) in the list of solutions to plot [2800, 4500, 8000, 10000]
x_coord = x[(n-1)::n]
y_coord = y[(n-1)::n]
z_coord = z[(n-1)::n]

for i in range(len(timestamp)):
    if epsilon_gr[i] > 0:
        # ax1.scatter(x_coord[i], y_coord[i], z_coord[i], marker = 'o', color = 'red', s = 1)
        ax1.scatter(x_coord[i], y_coord[i], marker = 'o', color = 'red', s = 1)
    if epsilon_gr[i] < 0:
        # ax1.scatter(x_coord[i], y_coord[i], z_coord[i], marker = 'o', color = 'blue', s = 1)
        ax1.scatter(x_coord[i], y_coord[i], marker = 'o', color = 'blue', s = 1)

    # Legend
    custom_lines = [Line2D([0], [0], color='red', lw=1),
                    Line2D([0], [0], color='blue', lw=1)]
    ax1.legend(custom_lines, ['Unbound', 'Bound'])

# Output
fig1.savefig('outflow2800-2d.eps', format='eps', bbox_inches='tight', pad_inches=0.1, dpi=1000)

end = time.time()
print "\nRuntime:", end - start, "s"