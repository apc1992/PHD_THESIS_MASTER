import os
from datetime import datetime

# Define the custom log directory

LOG_DIR = os.path.join(os.getcwd(), "logs", "json")

# Ensure the log directory exists
os.makedirs(LOG_DIR, exist_ok=True)

# Generate a timestamped log file name
timestamp = datetime.now().strftime("%Y-%m-%dT%H%M%S.%f")
log_file_path = f"{LOG_DIR}/{timestamp}.snakemake.log"

# Function to handle log records
def log_handler(msg):
    """
    Snakemake will call this function with log messages.
    This function preserves the original format of the messages.
    """
    # Convert the message to a string (if necessary)
    msg_str = str(msg)

    # Append the log message to the log file
    with open(log_file_path, "a") as log_file:
        log_file.write(msg_str + "\n")