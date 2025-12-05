import numpy as np
from scipy.spatial.transform import Rotation as R
import os
import subprocess
import math

# Function to convert a number to short scientific notation
def short_sci(x, decimals=2):
    if x == 0:
        return f"0E00"
    exp = int(math.floor(math.log10(abs(float(x)))))
    mantissa = float(x) / (10 ** exp)
    mantissa_str = f"{mantissa:.{decimals}f}"
    exp_str = f"{exp:02d}" if exp >= 0 else f"{exp:03d}"  # Handles negative exponents if needed
    return f"{mantissa_str}E{exp_str}"

# Generate kp values from 10^4 to 9*10^40
kp_vals = []
for exp in range(4, 5):  # 10^4 to 10^40
    for mul in range(4, 6):  # 1*10^exp, 1.25*10^exp, ..., 2*10^exp, ... 9*10^exp
        kp_vals.append(0.25*mul * 10**exp)
kp_vals = np.array(kp_vals)

# Write header to results file
header = f'{"kp gain":>9} {"kd gain":>9} {"mean":>12} {"median":>12} {"min":>9} {"max":>12} {"% below 20":>12}\n'
with open('gain_search_results.txt', 'w') as fout:
    fout.write(header)

# Loop over each kp value, run simulation, and record results
for i in range(len(kp_vals)):
    print(f'Running simulation {i+1} of {len(kp_vals)} with kp={kp_vals[i]}')
    kp = kp_vals[i]
    kd = 10 * kp

    with open('../../adcs-simulator/AOCS_oDyn_sim/data/ADCS_param.txt', 'r') as f:
        lines = f.readlines()

    targ_line = 5
    fields = lines[targ_line].strip().split(',')

    # Update the parameters in the target line
    fields[1] = str(kp)
    fields[2] = str(kd)
    lines[targ_line] = ','.join(fields) + '\n'

    # Write the updated lines back to the file
    with open('../../adcs-simulator/AOCS_oDyn_sim/data/ADCS_param.txt', 'w') as f:
        f.writelines(lines)

    # Run the simulation
    exe_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'adcs-simulator', 'AOCS_oDyn_sim', 'build', 'AOCS_oDyn_sim.exe'))
    exe_dir = os.path.dirname(exe_path)  # This will be your working dir for .exe

    subprocess.run([exe_path], cwd=exe_dir, stdout=subprocess.DEVNULL,) # devnull to suppress output

    # Read the output quaternion data
    with open('../../adcs-simulator/AOCS_oDyn_sim/data/output/quatR_B.txt', 'r') as fh:
        quatR_B = [list(map(float, line.strip().split(','))) for line in fh]

    quatR_B = np.array(quatR_B)  # list-of-lists, shape (N, 4)
    DCMO2B = np.zeros((3,3,len(quatR_B)))
    angleIWant = np.zeros(len(quatR_B))

    # Process each quaternion to compute the poionting error angle (I kept the terrible variable name)
    for i in range(len(quatR_B)):
        # MATLAB: quat2dcm([quatR_B(i,4),quatR_B(i,1),quatR_B(i,2),quatR_B(i,3)]);
        # scipy expects quaternion as [x, y, z, w]
        q = [quatR_B[i][0], quatR_B[i][1], quatR_B[i][2], quatR_B[i][3]]
        rot = R.from_quat(q)
        DCMO2B[:,:,i] = rot.as_matrix()
        zBody = DCMO2B[:,:,i].T @ np.array([0,0,1])
        angleIWant[i] = np.degrees(np.arccos(np.dot(zBody, np.array([0,0,1]))))

    # Compute statistics
    mean_val = np.mean(angleIWant)
    median_val = np.median(angleIWant)
    min_val = np.min(angleIWant)
    max_val = np.max(angleIWant)

    # Find the percentage of errors that are below accepted 20 degrees in the non-eclipse phase
    with open('../../adcs-simulator/AOCS_oDyn_sim/data/output/Eclipse_Matlab.txt', 'r') as fe:
        eclipse_data = [float(line.strip()) for line in fe if line.strip()] # 1 is non-eclipse, 0 is eclipse

    eclipse_data = np.array(eclipse_data)
    non_eclipse_indices = np.where(eclipse_data == 1)[0]
    non_eclipse_angles = angleIWant[non_eclipse_indices]
    acc_err = np.sum(non_eclipse_angles <= 20)
    per_acc_err = (acc_err / len(non_eclipse_angles)) * 100 if len(non_eclipse_angles) > 0 else 0

    kp_sci = short_sci(kp)
    kd_sci = short_sci(kd)
    line = f"{kp_sci:>9} {kd_sci:>9} {mean_val:12.4f} {median_val:12.4f} {min_val:9.4f} {max_val:12.4f} {per_acc_err:12.4f}\n"

    # Append results to the txt file
    with open('gain_search_results.txt', 'a') as fout:
        fout.write(line)

# plot the statistics
import matplotlib.pyplot as plt

data = np.loadtxt('gain_search_results.txt', skiprows=1)
rng = range(0, len(data))  # change the range as needed
kp_vals_plot = data[rng,0]
mean_vals = data[rng,2]
median_vals = data[rng,3]
min_vals = data[rng,4]
max_vals = data[rng,5]
perc_vals = data[rng,6]

# figure
fig = plt.figure(figsize=(20, 8))
plt.plot(rng, mean_vals, label='Mean Pointing Error')
plt.plot(rng, median_vals, label='Median Pointing Error')
plt.plot(rng, min_vals, label='Min Pointing Error')
plt.plot(rng, max_vals, label='Max Pointing Error')
plt.plot(rng, perc_vals, label='% Below 20 Degrees')
plt.xlabel('Kp Gain')
plt.ylabel('Pointing Error Statistics')
plt.title('Pointing Error Statistics vs Kp Gain')
plt.grid(True)

# put legend outside, on the right
ax = plt.gca()
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1.02, 0.6), borderaxespad=0.)

# give extra room on the right for the legend + text
plt.tight_layout(rect=[0, 0, 0.75, 1])

# ticks and labels - use the short_sci function, python cannot handle printing these kinds of numbers in a readible and sensible way
am_xticks = 21
#xticks_rep = np.linspace(np.min(kp_vals_plot), np.max(kp_vals_plot), am_xticks)
xticks = np.linspace(0, len(data)-1, am_xticks)
xlabels = [short_sci(v) for v in kp_vals_plot[xticks.astype(int)]]
yticks = np.linspace(0, 180, 19)
ax.set_xticks(xticks)
ax.set_xticklabels(xlabels, rotation=45, ha='right')
ax.set_yticks(yticks)

# side text: top Kp values
top_n = min(10, len(data))
lines = ["Top Kp (%<=20°):"]
for i in range(top_n):
    idx = np.argsort(-perc_vals)[i]
    lines.append(f"{i+1}. Kp {short_sci(kp_vals_plot[idx])} | {perc_vals[idx]:5.2f}% | max {max_vals[idx]:.2f}")
fig.text(0.78, 0.05, "\n".join(lines), ha='right', va='bottom', fontsize=9, family='monospace')

plt.show()