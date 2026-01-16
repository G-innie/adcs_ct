
import numpy as np
from scipy.spatial.transform import Rotation as R
import os
import subprocess
import matplotlib.pyplot as plt


# Generate kp values from 10^4 to 9*10^40
satdip_vals =  np.arange(0.01, 0.31, 0.01)


# Write header to results file
header = f'{"satdip gain":>9} {"mean":>12} {"median":>12} {"min":>9} {"max":>12} {"% below 20 in sun":>12} {"% below 40 in ecl":>12} \n'
with open('sat_dip_search_results.txt', 'w') as fout:
    fout.write(header)

# Loop over each kp value, run simulation, and record results
for i in range(len(satdip_vals)):
    print(f'Running simulation {i+1} of {len(satdip_vals)} with satdip={satdip_vals[i]}')
    satdip = satdip_vals[i]

    with open('../../../adcs-simulator/AOCS_oDyn_sim/data/Param.txt', 'r') as f:
        lines = f.readlines()

    # If sat_dip is on line 42 (use 0-indexed, so line 42 = index 41)
    lines[10] = str(satdip)

    with open('../../../adcs-simulator/AOCS_oDyn_sim/data/Param.txt', 'w') as f:
        f.writelines(lines)
        f.flush()


    # Run the simulation
    exe_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', 'adcs-simulator', 'AOCS_oDyn_sim', 'build', 'AOCS_oDyn_sim.exe'))
    exe_dir = os.path.dirname(exe_path)  # This will be your working dir for .exe

    # Run the simulation
    subprocess.run([exe_path], cwd=exe_dir, stdout=subprocess.DEVNULL) # devnull to suppress output

    # Read the output quaternion data
    with open('../../../adcs-simulator/AOCS_oDyn_sim/data/output/quatR_B.txt', 'r') as fh:
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
    with open('../../../adcs-simulator/AOCS_oDyn_sim/data/output/Eclipse_Matlab.txt', 'r') as fe:
        eclipse_data = [float(line.strip()) for line in fe if line.strip()] # 1 is non-eclipse, 0 is eclipse

    eclipse_data = np.array(eclipse_data)
    sun_indices = np.where(eclipse_data == 1)[0]
    eclipse_indices = np.where(eclipse_data == 0)[0]
    sun_angles = angleIWant[sun_indices]
    eclipse_angles = angleIWant[eclipse_indices]
    sun_err = np.sum(sun_angles <= 20)
    ecl_err = np.sum(eclipse_angles <= 40)
    per_sun_err = (sun_err / len(sun_angles)) * 100 if len(sun_angles) > 0 else 0
    per_ecl_err = (ecl_err / len(eclipse_angles)) * 100 if len(eclipse_angles) > 0 else 0

    line = f"{satdip:9.2f} {mean_val:12.4f} {median_val:12.4f} {min_val:9.4f} {max_val:12.4f} {per_sun_err:12.4f} {per_ecl_err:12.4f}\n"

    # Append results to the txt file
    with open('sat_dip_search_results.txt', 'a') as fout:
        fout.write(line)

    # plot the pointing error for this satdip value for future reference if needed
    plt.figure()
    plt.plot(angleIWant)
    plt.xlabel('Time')
    plt.ylabel('Pointing Error (degrees)')
    plt.title(f'Pointing Error for satdip={satdip}')
    plt.grid(True)
    plt.savefig(f'pointing_error_satdip_{satdip:.2f}.png')
    plt.close()

# plot the statistics

data = np.loadtxt('sat_dip_search_results.txt', skiprows=1)
rng = range(0, len(data))  # change the range as needed
satdip_plot = data[rng,0]
mean_vals = data[rng,1]
median_vals = data[rng,2]
min_vals = data[rng,3]
max_vals = data[rng,4]
perc_sun_vals = data[rng,5]
perc_ecl_vals = data[rng,6]

# figure
fig = plt.figure(figsize=(20, 8))
plt.plot(rng, mean_vals, label='Mean Pointing Error')
plt.plot(rng, median_vals, label='Median Pointing Error')
plt.plot(rng, min_vals, label='Min Pointing Error')
plt.plot(rng, max_vals, label='Max Pointing Error')
plt.plot(rng, perc_sun_vals, label='% Below 20 Degrees (Sun)')
plt.plot(rng, perc_ecl_vals, label='% Below 40 Degrees (Eclipse)')
plt.xlabel('Sat Dip')
plt.ylabel('Pointing Error Statistics')
plt.title('Pointing Error Statistics vs Sat Dip')
plt.grid(True)

# put legend outside, on the right
ax = plt.gca()
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1.02, 0.6), borderaxespad=0.)

plt.show()