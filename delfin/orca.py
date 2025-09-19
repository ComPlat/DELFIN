import platform
import subprocess
import logging
import sys

def find_orca_executable():
    command = "where" if platform.system() == "Windows" else "which"
    result = subprocess.run([command, "orca"], capture_output=True, text=True)
    orca_path = result.stdout.strip()
    if not orca_path:
        logging.error("ORCA executable not found. Please ensure ORCA is installed and in your PATH.")
        return None
    return orca_path

def run_orca(input_file_path, output_log):
    orca_path = find_orca_executable()
    if not orca_path:
        return
    with open(output_log, "w") as output_file:
        try:
            subprocess.run([orca_path, input_file_path], check=True, stdout=output_file, stderr=output_file)
            logging.info(f"ORCA run successful for '{input_file_path}'")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running ORCA: {e}")

def run_orca_IMAG(input_file_path, iteration):
    orca_path = find_orca_executable()
    output_log = f"output_{iteration}.out"
    with open(output_log, "w") as output_file:
        try:
            subprocess.run([orca_path, input_file_path], check=True, stdout=output_file, stderr=output_file)
            logging.info(f"ORCA run successful for '{input_file_path}', output saved to '{output_log}'")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running ORCA: {e}")
            sys.exit(1)

def run_orca_plot(homo_index):
    for index in range(homo_index - 10, homo_index + 11):
        process = subprocess.Popen(['orca_plot', 'input.gbw', '-i'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        inputs = f'1\n1\n4\n100\n5\n7\n2\n{index}\n10\n11\n'.encode()
        stdout, stderr = process.communicate(input=inputs)
        if process.returncode == 0:
            logging.info(f"orca_plot ran successfully for index {index}")
        else:
            logging.error(f"orca_plot encountered an error for index {index}: {stderr.decode()}")
