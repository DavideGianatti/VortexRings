#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import pandas as pd

# Define simulation parameters
Nt = 10                # Number of time evolution steps
dt = 5 / 10            # Time discretization
Nv = 2                 # Number of tubes
M = 100                # Number of discretization points per tube
strength = np.full(Nv, 5)       # Vorticity flux constant for each tube
length = np.ones(Nv, dtype=float)  # Length of each tube
length_0 = np.copy(length)      # Initial length of each tube (fixed)
sigma = np.full(Nv, 0.5)        # Width of each tube
sigma_0 = np.copy(sigma)        # Initial width of each tube (fixed)
position = np.zeros((3 * Nv, M), dtype=float)  # Positions of discretization points for each tube
velocity = np.zeros((3 * Nv, M), dtype=float)  # Velocities of discretization points
df = pd.DataFrame({"time": np.empty(Nv * M * Nt), "x": np.empty(Nv * M * Nt), "y": np.empty(Nv * M * Nt), "z": np.empty(Nv * M * Nt)})

# Smoothing function parameters
K = 10**4              # Number of discretization points for smoothing function
dx = 3 / K             # Smoothing function integration step size
smoothing_values = np.zeros(K, dtype=float)  # Smoothing function values

# Calculate smoothing function values
integral = 0
for i in range(K):
    integral += (i * dx) ** 2 * np.exp(-(i * dx) ** 2) / np.pi ** (3 / 2) * dx
    smoothing_values[i] = integral
smoothing_values *= 4 * np.pi

def smoothing(x):
    """Smoothing function to avoid divergences in the Biot-Savart integral.

    Args:
        x (float): The distance for which the smoothing function value is needed.

    Returns:
        float: The smoothing function value for the given distance.
    """
    if x >= 3:
        return 1
    else:
        return smoothing_values[round(x / dx) - 1]

def generate_circle(r, cx, cy, cz, phi, n):
    """Generates a matrix with the coordinates of the points of a circle.

    Args:
        r (float): Radius of the circle.
        cx (float): X-coordinate of the center of the circle.
        cy (float): Y-coordinate of the center of the circle.
        cz (float): Z-coordinate of the center of the circle.
        phi (float): Angle of inclination of the circle.
        n (int): Number of discretization points along the circle.

    Returns:
        np.ndarray: A 3 x n array with the coordinates of the points on the circle.
    """
    dtheta = 2 * np.pi / n
    coordinates = np.zeros((3, n), dtype=float)
    for i in range(n):
        coordinates[0, i] = cx + r * np.cos(i * dtheta)
        coordinates[1, i] = cy + r * np.sin(i * dtheta) * np.cos(phi)
        coordinates[2, i] = cz - r * np.sin(i * dtheta) * np.sin(phi)
    return coordinates

def initialize_positions(r, cx, cy, cz, phi, n, deltax, deltay, deltaz):
    """Initialize the positions of the tubes.

    Args:
        r (float): Radius of the circle used to initialize the positions.
        cx (float): X-coordinate of the center of the initial circle.
        cy (float): Y-coordinate of the center of the initial circle.
        cz (float): Z-coordinate of the center of the initial circle.
        phi (float): Angle of inclination of the initial circle.
        n (int): Number of discretization points along the initial circle.
        deltax (float): Increment in the x-direction for each tube.
        deltay (float): Increment in the y-direction for each tube.
        deltaz (float): Increment in the z-direction for each tube.
    """
    for i in range(Nv):
        coordinates = generate_circle(r, cx, cy, cz, phi, n)
        length_0[i] = 2 * np.pi * r  # Fix the initial lengths
        for j in range(M):
            position[3 * i, j] = coordinates[0, j] + i * deltax
            position[3 * i + 1, j] = coordinates[1, j] + i * deltay
            position[3 * i + 2, j] = coordinates[2, j] + i * deltaz
            df.loc[i * M + j, ['x', 'y', 'z']] = position[3 * i:3 * i + 3, j]
        df.loc[i * M:(i + 1) * M, 'time'] = 0

def calculate_velocity():
    """Calculate the velocity at each discretization point using the Biot-Savart law."""
    u = np.zeros(3, dtype=float)  # Velocity at a point
    dr = np.zeros(3, dtype=float)  # Displacement vector between two adjacent points on the tube
    xr = np.zeros(3, dtype=float)  # Displacement vector from the point of interest to a point on the tube
    for i in range(Nv):
        for j in range(M):
            for k in range(Nv):
                for l in range(M):
                    if k == i and l == j:
                        continue  # Avoid division by zero
                    dr[:] = position[3 * k:3 * k + 3, l] - position[3 * k:3 * k + 3, l - 1]
                    xr[:] = position[3 * i:3 * i + 3, j] - position[3 * k:3 * k + 3, l]
                    dxr = np.linalg.norm(xr)
                    const = strength[k] * smoothing(dxr / sigma[k]) / dxr ** 3
                    u += np.cross(xr, dr) * const
            u = -u / (4 * np.pi)
            velocity[3 * i:3 * i + 3, j] = u
            u[:] = 0

def update_positions():
    """Update the positions of the tubes based on their velocities."""
    for i in range(Nv):
        for j in range(M):
            position[3 * i:3 * i + 3, j] += velocity[3 * i:3 * i + 3, j] * dt
            df.loc[i * M + j + M * Nv * m, ['x', 'y', 'z']] = position[3 * i:3 * i + 3, j]
    df.loc[M * Nv * m:M * Nv * (m + 1), 'time'] = m

def update_width():
    """Update the width of each tube to conserve volume."""
    length[:] = 0
    dr = np.zeros(3, dtype=float)  # Displacement vector between two adjacent points on the tube
    for k in range(Nv):
        for l in range(M):
            dr[:] = position[3 * k:3 * k + 3, l] - position[3 * k:3 * k + 3, l - 1]
            length[k] += np.linalg.norm(dr)
        sigma[k] = sigma_0[k] * np.sqrt(length_0[k] / length[k])

def animate(frame):
    """Animation function to update the plot.

    Args:
        frame (int): The current frame number to display.
    """
    data = df[df['time'] == frame]
    graph._offsets3d = (data.x, data.y, data.z)

# Initialize positions
initialize_positions(5, 0, 0, 0, 0, M, 0, 0, 0.2)

# Time evolution loop
for m in range(1, Nt):
    print(f"t = {m}")
    calculate_velocity()
    update_positions()
    update_width()

# Plot and animate the results
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_xlim3d(-5, 5)
ax.set_ylim3d(-5, 5)
ax.set_zlim3d(-4, 8)
data = df[df['time'] == 0]
graph = ax.scatter(data.x, data.y, data.z, color='red')

anim = animation.FuncAnimation(fig, animate, Nt, interval=10, blit=False)
ax.view_init(30, 30)

writervideo = animation.FFMpegWriter(fps=20)

# Uncomment the line below to save the animation as a video file
# anim.save("formazione.mp4", writer=writervideo)

plt.show()
