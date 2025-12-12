import random
import logging

def log_out(message):
    line = random.randint(1, 100)
    logger.debug(f"{line} {message}")

# Configure logging
LOG_FMT = "[%(levelname)s] [%(asctime)s] subsystems/adcs/adcs.tim.c:%(message)s"
DATE_FMT = "%Y-%m-%d %H:%M:%S"

logging.basicConfig(
    level=logging.DEBUG,
    format=LOG_FMT,
    datefmt=DATE_FMT,
)

logger = logging.getLogger(__name__)

# loop over 1000 iterations to simulate data generation
eclipse = 0

for i in range(10000):
    # Generate dummy data for simulation output
    q_out = [random.random()*1000 for _ in range(4)]
    ohm_out = [random.random()*1000 for _ in range(3)]
    com_dip_out = [random.random()*1000 for _ in range(3)]

    if i % 100 == 0:
        eclipse = 1 if eclipse == 0 else 0

    if i % 10 == 0:
        log_out("No ADCS parameters have changed, skipping initialization.")

    if i % 5 == 0:
        log_out("internal delay for 2 ms")

    if i % 7 == 0:
        log_out("main delay for 1ms, will accept at most 0 tc")

    log_out(f"q_out: {q_out[0]} {q_out[1]} {q_out[2]}")
    log_out(f"{q_out[3]}")
    log_out(f"ohm_out: {ohm_out[0]} {ohm_out[1]} {ohm_out[2]}")
    log_out(f"com_dip_out: {com_dip_out[0]} {com_dip_out[1]} {com_dip_out[2]}")
    log_out(f"eclipseflag: {eclipse}")
    