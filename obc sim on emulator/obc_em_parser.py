import re
import argparse

# Regular expressions to match the desired lines
re_day = re.compile(r"sim_time:\s*year\s*(\d+)\s*,\s*month\s*(\d+)\s*,\s*day\s*(\d+)")
re_hour = re.compile(r"hours\s*(\d+)\s*,\s*minutes\s*(\d+)\s*,\s*seconds\s*(\d+)")
re_q = re.compile(r"q_out:\s*(-?\d+)\s*(-?\d+)\s*(-?\d+)")
re_q3 = re.compile(r"q4\s*(-?\d+)")
re_ohm = re.compile(r"ohm_out:\s*(-?\d+)\s*(-?\d+)\s*(-?\d+)")
re_com = re.compile(r"com_dip_out:\s*(-?\d+)\s*(-?\d+)\s*(-?\d+)")
re_eclipse = re.compile(r"eclipseflag:\s*(\d+)")


# Set up argument parsing, there is one input for the output file, but there will be three output files, output_file_q, output_file_ohm, output_file_com
parser = argparse.ArgumentParser()
parser.add_argument(
    'input_file',
    nargs='?',          # makes the argument optional
    default='out.txt',   # value to use if not provided
    help='Path to the input file (default: out.txt)'
)
parser.add_argument(
    'output_file',
    nargs='?',          # makes the argument optional
    default='parsed',   # value to use if not provided
    help='Path to the output file (default: parsed)'
)
args = parser.parse_args()

# states are set such that the first line we are looking for is q_out, then we are looking for the next line for q3, then ohm_out, then com_dip_out, then back to q_out
state = 'default'
pending_q = None

# variables for storing current readings
curr_q = []
curr_ohm = []
curr_com = []
curr_day = []
curr_eclipse = 0

# lists for storing all measurements
meas_q = []
meas_ohm = []
meas_com = []
meas_ecl = []

# t is telling us which reading we are encountring (probably will need to be adjusted to collect time readings from mock time) - TO DO
t = -1

print(f"Parsing input file: {args.input_file}")
print(f"Output will be written to: {args.output_file}_q.txt, {args.output_file}_ohm.txt, {args.output_file}_com.txt, {args.output_file}_eclipse.txt")

with open(args.input_file, 'r') as fin:
    for line in fin:
        # getting the day information
        if state == 'default':
            match_day = re_day.search(line)
            if match_day:
                curr_day = list(map(int, match_day.groups()))
                state = 'got day'

        # getting the time information
        elif state == 'got day':
            match_hour = re_hour.search(line)
            if match_hour:
                curr_day += list(map(int, match_hour.groups()))
                
                if t < 0:
                    t = curr_day[5]
                    prev_day = curr_day
                else:
                    # time difference calculation, could be improved for month lengths and leap years
                    t += 365*86400*(curr_day[0] - prev_day[0])  # years to seconds
                    t += 30*86400*(curr_day[1] - prev_day[1])   # months to seconds
                    t += 86400*(curr_day[2] - prev_day[2])    # days to seconds
                    t += 3600*(curr_day[3] - prev_day[3])    # hours to seconds
                    t += 60*(curr_day[4] - prev_day[4])      # minutes to seconds
                    t += (curr_day[5] - prev_day[5])       # seconds to seconds
                    prev_day = curr_day
                state = 'got time'

        # collecting q_out measurements
        elif state == 'got time': 
            if pending_q is None:
                match_q = re_q.search(line)
                if match_q:
                    pending_q = list(map(float, match_q.groups()))
            # collecting the last component of q_out
            else:
                match_q3 = re_q3.search(line)
                if match_q3:
                    pending_q.append(float(match_q3.group(1)))
                    curr_q = pending_q
                    meas_q.append((t, *curr_q))
                    pending_q = None
                    state = 'got quats'

        # collecting ohm_out measurements
        elif state == 'got quats':
            match_ohm = re_ohm.search(line)
            if match_ohm:
                curr_ohm = list(map(float, match_ohm.groups()))
                meas_ohm.append((t, *curr_ohm))
                state = 'got ohm'

        # collecting com_dip_out measurements
        elif state == 'got ohm':
            match_com = re_com.search(line)
            if match_com:
                curr_com = list(map(float, match_com.groups()))

                # Store the complete measurement
                meas_com.append((t, *curr_com))
                state = 'got com'

        # collecting eclipseflag measurements
        elif state == 'got com':
            match_eclipse = re_eclipse.search(line)
            if match_eclipse:
                curr_eclipse = int(match_eclipse.group(1))
                meas_ecl.append((t, curr_eclipse))
                state = 'default'  # ready for the next measurement set


# Write the collected measurements to output files
with open(f"{args.output_file}_q.txt", 'w') as fq, \
     open(f"{args.output_file}_ohm.txt", 'w') as fohm, \
     open(f"{args.output_file}_com.txt", 'w') as fcom, \
     open(f"{args.output_file}_eclipse.txt", 'w') as fecl:
    
    for entry in meas_q:
        fq.write(" ".join(map(str, entry)) + "\n")
    
    for entry in meas_ohm:
        fohm.write(" ".join(map(str, entry)) + "\n")
    
    for entry in meas_com:
        fcom.write(" ".join(map(str, entry)) + "\n")

    for entry in meas_ecl:
        fecl.write(" ".join(map(str, entry)) + "\n")
        

print(f"Parsing complete. Output files generated. Found {t} measurements.")

