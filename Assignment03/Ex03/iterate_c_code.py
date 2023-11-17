import subprocess
import numpy as np
import os

ROOT = os.getcwd()


def rename_file(folder_path, old_filename, new_filename):
    # Construct the full paths for the old and new files
    old_file_path = os.path.join(folder_path, old_filename)
    new_file_path = os.path.join(folder_path, new_filename)

    try:
        # Rename the file
        os.rename(old_file_path, new_file_path)
    except FileNotFoundError:
        print(f"Error: File '{old_filename}' not found.")
    except FileExistsError:
        print(f"Error: File '{new_filename}' already exists.")


def execute_c_program(program: str, input_values: list) -> str:
    c_program_path = os.path.join(ROOT, program)

    # Construct the command to execute the C program
    command = [c_program_path]

    # Run the C program using Popen to interact with stdin and stdout
    with subprocess.Popen(
        command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True
    ) as process:
        # Provide input to the subprocess
        process.stdin.write("\n".join(input_values))
        process.stdin.close()  # Close stdin to signal the end of input

        # Capture and print the output
        output = process.stdout.read()

    return output


def parse_diff_output(output: str, wanted_fields: list, datatypes: list) -> str:
    lines = output.splitlines()
    lines = lines[-8:]
    wanted_output = []
    for line in lines:
        for i, wanted_field in enumerate(wanted_fields):
            if wanted_field in line:
                wanted_output.append(
                    tuple(
                        [
                            datatypes[i](x.strip()) if j == 1 else x.strip()
                            for j, x in enumerate(line.split(":"))
                        ]
                    )
                )
    return wanted_output


field_to_optimize = "MSE"

noisy_image_path = "noisy.pgm"
original_image_path = "original.pgm"
output_image_path = "output_image.pgm"
num_iterations = 10

# Start values for each variable
ts_start = 0.23
lambda_start = 12.7
sigma_start = 1.7

# Max value for each variable
time_step_max = 0.23
lambda_max = 12.7
sigma_max = 1.7

# Step for each variable
ts_step = 0.01
lambda_step = 0.1
sigma_step = 0.01

current_state = []
min_value = 1000000
optimized_params = []
counter = 1
for ts in np.arange(ts_start, time_step_max + ts_step, ts_step):
    for lambda_ in np.arange(lambda_start, lambda_max + lambda_step, lambda_step):
        for sigma in np.arange(sigma_start, sigma_max + sigma_step, sigma_step):
            input_gauss_conv = [
                noisy_image_path,
                original_image_path,
                str(lambda_),
                str(sigma),
                str(ts),
                str(num_iterations),
                output_image_path,
            ]

            output = execute_c_program("iso_non_diff", input_gauss_conv)

            parsed_output = parse_diff_output(
                output, [field_to_optimize], [float, float]
            )

            curr_value = parsed_output[0][1]
            if curr_value < min_value:
                optimized_params = [ts, lambda_, sigma]
                min_value = curr_value
                print(
                    f"Found better params in try #{counter} [ts = {ts}, la = {lambda_}, sig = {sigma}] with {field_to_optimize} = {curr_value}"
                )
                rename_file(ROOT, "output_image.pgm", "best_output_image.pgm")

            counter += 1
