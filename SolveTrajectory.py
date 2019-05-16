from numpy import *
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import time
start = time.time()

# Initialize lists to store data in.
timestamp = []
mass = []
x_pos = []
y_pos = []
z_pos = []
minus_ut = []
minus_uth = []
vx = []
vy = []
vz = []
c_epsilon = []
gr_epsilon = []
vx_gr = []
vy_gr = []
vz_gr = []
solx = []
soly = []
solz = []

# Read in data from file.
with open('AllParticlesFromOutflowLev0.dat', 'r') as data_file:
    # data_slice = [next(data_file2) for x in xrange(6000)]  # Read in file up to line N.
    for line in data_file:
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

# Define matrix to store (x, y, z) spatial data in.
matrix = zeros((len(timestamp),3,4))  # Row, inner row, column.

# Differential equations of motion.
def diffEq(t, x):
    G = 1
    M = 1.2 + 1.4
    tau = e**log(t)
    # Vectorized RHS of differential equation (vx, vy, vz, a_x, a_y, a_z).
    f = [x[3],
             x[4],
             x[5],
            -G*M*x[0] / (x[0]**2 + x[1]**2 + x[2]**2)**(3./2),
            -G*M*x[1] / (x[0]**2 + x[1]**2 + x[2]**2)**(3./2),
            -G*M*x[2] / (x[0]**2 + x[1]**2 + x[2]**2)**(3./2)]
    return f

# Transform classical velocities to relativistic velocities.
def veloTransform(x_pos, y_pos, z_pos, vx, vy, vz, minus_ut):
    for i in range(len(timestamp)):
        r = float(sqrt(x_pos[i] ** 2 + y_pos[i] ** 2 + z_pos[i] ** 2))
        v = sqrt(vx[i] ** 2 + vy[i] ** 2 + vz[i] ** 2)
        A = sqrt(2 * (minus_ut[i] + 2.6/r - 1) / v**2)  # Transformation term
        vx_gr.append(A * vx[i])
        vy_gr.append(A * vy[i])
        vz_gr.append(A * vz[i])
    return vx_gr, vy_gr, vz_gr

# Relativistic Trajectories.
def relativisticCalc(x_pos, y_pos, z_pos, vx, vy, vz, minus_ut, timestamp, t_eval):
    v_gr = veloTransform(x_pos, y_pos, z_pos, vx, vy, vz, minus_ut)
    for i in range(len(timestamp)):
        # Relativistic correction to velocity
        gr_velocity = sqrt(v_gr[0][i]**2 + v_gr[1][i]**2 + v_gr[2][i]**2)
        # Initial conditions
        initial_cond = [x_pos[i], y_pos[i], z_pos[i], v_gr[0][i], v_gr[1][i], v_gr[2][i]]
        # Solve differential equation
        sol = solve_ivp(diffEq, (timestamp[i], 10000), initial_cond, method = 'RK45', t_eval = t_eval)
        # Check for bound conditions
        r = float(sqrt(x_pos[i] ** 2 + y_pos[i] ** 2 + z_pos[i] ** 2))
        gr_epsilon.append(1 + (gr_velocity**2)/2 - 2.6/r)
        # Temporarily store position solutions
        solx = ndarray.tolist(sol.y[0])
        soly = ndarray.tolist(sol.y[1])
        solz = ndarray.tolist(sol.y[2])
        # Store data into NxMx3 array: {([x,y,z], [x,y,z], [x,y,z]), ([...], [...], [...])}.
        for j in range(len(solx)):
            matrix[i][0][j] = solx[j]
            matrix[i][1][j] = soly[j]
            matrix[i][2][j] = solz[j]
        # Clear lists for next iteration.
        solx[:] = []
        soly[:] = []
        solz[:] = []

    # Write relevant data into new file
    with open('ParticleTrajectories.dat', 'w') as outfile:
        outfile.write("# Evaluated Time\n# X Position\n# Y Position\n# Z Position\n")
        for i in range(len(timestamp)):
            data_lists = [sol.t, matrix[i][0], matrix[i][1], matrix[i][2]]  # Each entry is a list itself.
            for x in zip(*data_lists):
                outfile.write('{0} {1} {2} {3}\n'.format(*x))

t_eval = [2800, 4500, 8000, 10000]
relativisticCalc(x_pos, y_pos, z_pos, vx, vy, vz, minus_ut, timestamp, t_eval)

end = time.time()
print "\nRuntime:", end - start, "s"