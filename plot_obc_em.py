# plot the results using matplotlib
import matplotlib.pyplot as plt
import numpy as np
import argparse

# Set up argument parsing, there is one input for the output file, but there will be three output files, output_file_q, output_file_ohm, output_file_com
parser = argparse.ArgumentParser()
parser.add_argument(
    'output_file',
    nargs='?',          # makes the argument optional
    default='parsed',   # value to use if not provided
    help='Path to the output file (default: parsed)'
)
args = parser.parse_args()

# Load the data from the output file
meas_q = np.loadtxt(f"{args.output_file}_q.txt")
meas_ohm = np.loadtxt(f"{args.output_file}_ohm.txt")
meas_com = np.loadtxt(f"{args.output_file}_com.txt")
eclipse = np.loadtxt(f"{args.output_file}_eclipse.txt")

mask = eclipse[:,1].astype(bool)

# labels for plots
label_q = ["i", "j", "k", "w"]
label_ohm_com = ["x", "y", "z"]

# estimate quaternions plot
fig1, ax1 = plt.subplots(figsize=(8, 4))

# Grey background where in eclipse (full height of axes)
ax1.fill_between(
    meas_q[:,0], 0, 1,
    where=mask,
    transform=ax1.get_xaxis_transform(),  # span full y axis in axes coords
    color='0.9', edgecolor='0.9'
)

# Plot each of the 4 columns
for i in range(4):
    ax1.plot(meas_q[:,0], meas_q[:, i+1]/1000, label=label_q[i])

ax1.set_xlabel("Time")
ax1.set_ylabel("Quaternion estimates")
ax1.set_title("Estimated Quaternions from OBC sim emulator")
ax1.legend()
ax1.grid(True)
plt.tight_layout()

# estimate ohm plot
fig2, ax2 = plt.subplots(figsize=(8, 4))

# Grey background where in eclipse (full height of axes)
ax2.fill_between(
    meas_ohm[:,0], 0, 1,
    where=mask,
    transform=ax2.get_xaxis_transform(),  # span full y axis in axes coords
    color='0.9', edgecolor='0.9'
)

# Plot each of the 3 columns
for i in range(3):
    ax2.plot(meas_ohm[:,0], meas_ohm[:, i+1]/1000, label=label_ohm_com[i])
ax2.set_xlabel("Time")
ax2.set_ylabel("Ohm estimates")
ax2.set_title("Estimated Ohm from OBC sim emulator")
ax2.legend()
ax2.grid(True)

plt.tight_layout()

# estimate commanded dipole plot
fig3, ax3 = plt.subplots(figsize=(8, 4))

# Grey background where in eclipse (full height of axes)
ax3.fill_between(
    meas_com[:,0], 0, 1,
    where=mask,
    transform=ax3.get_xaxis_transform(),  # span full y axis in axes coords
    color='0.9', edgecolor='0.9'
)

# Plot each of the 3 columns
for i in range(3):
    ax3.plot(meas_com[:,0], meas_com[:, i+1]/1000, label=label_ohm_com[i])
ax3.set_xlabel("Time")
ax3.set_ylabel("Commanded dipole estimates")
ax3.set_title("Estimated Commanded Dipole from OBC sim emulator")
ax3.legend()
ax3.grid(True)
plt.tight_layout()
plt.show()