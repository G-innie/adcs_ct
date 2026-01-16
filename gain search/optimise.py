# This function is meant to find as good kp, kd, and sat_dip values as possible within defined ranges. The functionality is similar to the other files in this folder, but the new part comes from combnining the kp, and kd search with the sat_dip search. On top of that, bayesian optimisation is used to find the best values rather than a brute-force grid search.
# last edidted: 16th January 2026 by Ginnie

import os
import subprocess
import numpy as np
from scipy.spatial.transform import Rotation as R
from skopt import gp_minimize
from skopt.space import Real
from skopt.callbacks import VerboseCallback
import csv

exe_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', 'adcs-simulator', 'AOCS_oDyn_sim', 'build', 'AOCS_oDyn_sim.exe'))
exe_dir = os.path.dirname(exe_path)  # The folder where the .exe is located

def write_params(kp, kd, satdip):
    # read existing ADCS_param.txt file, modify the kp and kd values, and write back
    with open('../../../adcs-simulator/AOCS_oDyn_sim/data/ADCS_param.txt', 'r') as f:
        lines = f.readlines()
    
    targ_line = 5
    fields = lines[targ_line].strip().split(',')

    # Update the parameters in the target line
    fields[1] = str(kp)
    fields[2] = str(kd)
    lines[targ_line] = ','.join(fields) + '\n'

    # Write the updated lines back to the file
    with open('../../../adcs-simulator/AOCS_oDyn_sim/data/ADCS_param.txt', 'w') as f:
        f.writelines(lines)

    # read existing Param.txt file, modify the sat_dip value, and write back
    with open('../../../adcs-simulator/AOCS_oDyn_sim/data/Param.txt', 'r') as f:
        lines = f.readlines()
    
    lines[10] = str(satdip)

    with open('../../../adcs-simulator/AOCS_oDyn_sim/data/Param.txt', 'w') as f:
        f.writelines(lines)

def run_simulation():
    subprocess.run([exe_path], cwd=exe_dir, stdout=subprocess.DEVNULL) # devnull to suppress output

def evaluate():
    # Read the output quaternion data
    with open('../../../adcs-simulator/AOCS_oDyn_sim/data/output/quatR_B.txt', 'r') as fh:
        quatR_B = [list(map(float, line.strip().split(','))) for line in fh]

    quatR_B = np.array(quatR_B)  # list-of-lists, shape (N, 4)
    DCMO2B = np.zeros((3,3,len(quatR_B)))
    angleIWant = np.zeros(len(quatR_B))

    # Process each quaternion to compute the pointing error angle
    for i in range(len(quatR_B)):
        q = [quatR_B[i][0], quatR_B[i][1], quatR_B[i][2], quatR_B[i][3]]
        rot = R.from_quat(q)
        DCMO2B[:,:,i] = rot.as_matrix()
        zBody = DCMO2B[:,:,i].T @ np.array([0,0,1])
        angleIWant[i] = np.degrees(np.arccos(np.dot(zBody, np.array([0,0,1]))))

    # Only consider data after index 5000 (assumes data after 5000 seconds)
    start_idx = 5000
    angleIWant_filtered = angleIWant[start_idx:]
    
    # Compute statistics on filtered data
    median_val = np.median(angleIWant_filtered)
    max_val = np.max(angleIWant_filtered)

    # Find the percentage of errors that are below accepted 20 degrees in the non-eclipse phase
    with open('../../../adcs-simulator/AOCS_oDyn_sim/data/output/Eclipse_Matlab.txt', 'r') as fe:
        eclipse_data = [float(line.strip()) for line in fe if line.strip()] # 1 is non-eclipse, 0 is eclipse

    eclipse_data = np.array(eclipse_data)
    eclipse_data_filtered = eclipse_data[start_idx:]
    
    sun_indices = np.where(eclipse_data_filtered == 1)[0]
    ecl_indices = np.where(eclipse_data_filtered == 0)[0]
    sun_point_err = angleIWant_filtered[sun_indices]
    ecl_point_err = angleIWant_filtered[ecl_indices]
    sun_err = np.sum(sun_point_err <= 20)
    ecl_err = np.sum(ecl_point_err <= 40)
    per_sun_err = (sun_err / len(sun_point_err)) if len(sun_point_err) > 0 else 0
    per_ecl_err = (ecl_err / len(ecl_point_err)) if len(ecl_point_err) > 0 else 0
    c_per_sun = 1 - per_sun_err
    c_per_ecl = 1 - per_ecl_err
    c_median = min(median_val / 15, 1)  # cap at 1 for median errors above 20 degrees
    c_max = max_val / 125 if max_val < 125 else 1 + (max_val - 125) / 125  # linear penalty for max errors above 90 degrees

    cost = 4 * c_per_sun + 2 * c_median + 1 * c_max + 0.5 * c_per_ecl

    return cost, per_sun_err, median_val, max_val, per_ecl_err

def objective(params):
    kp, kd, satdip = params

    write_params(kp, kd, satdip)
    run_simulation()
    cost, per_sun_err, median_val, max_val, per_ecl_err = evaluate()

    # log everything for later analysis
    results_log.append({
        'kp': kp,
        'kd': kd,
        'satdip': satdip,
        'cost': cost,
        'per_sun_err': per_sun_err,
        'median_val': median_val,
        'max_val': max_val,
        'per_ecl_err': per_ecl_err,
    })

    return cost

results_log = []

space = [
    Real(1e5, 1e9, prior='log-uniform', name='kp'),
    Real(1e6, 1e10, prior='log-uniform', name='kd'),
    Real(0.01, 0.3, prior='uniform', name='satdip')
]

result = gp_minimize(
    func = objective,
    dimensions = space,
    n_calls=20, 
    n_initial_points=10,
    acq_func='EI',
    random_state=0,
    callback=[VerboseCallback(n_total=20)]
    )

best_kp, best_kd, best_satdip = result.x
print(f'Best kp: {best_kp}, Best kd: {best_kd}, Best sat_dip: {best_satdip}')

# Save results log to a CSV for further analysis
with open('opt_results_log.csv', 'w', newline='') as csvfile:
    fieldnames = ['kp', 'kd', 'satdip', 'cost', 'per_sun_err', 'median_val', 'max_val', 'per_ecl_err']

    writer = csv.DictWriter(
        csvfile, 
        fieldnames=fieldnames)

    writer.writeheader()
    for entry in results_log:
        writer.writerow(entry)