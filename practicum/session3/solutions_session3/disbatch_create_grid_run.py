import os
import shutil

# --- User settings ---
base_dir = "/mnt/home/lvanson/ASBE-course/practicum/session3/solutions_session3/data"
mesa_base = os.path.join(base_dir, "HW3_base")
output_dir = os.path.join(base_dir, "HW3_output")
inlist_template = os.path.join(mesa_base, "inlist_gridrun")

# Define your mass grid and Dutch wind scaling
masses = [50, 55, 60]          # np.arange(15, 105,5)
dutch_scaling = 1.0                # Could later loop over different values if needed

# --- Prepare Taskfile ---
taskfile = os.path.join(base_dir, "Taskfile")
with open(taskfile, "w") as tf:
    # following line is effectively executed before each disBatch command
    # It sets up the environment for MESA runs
    tf.write(f"#DISBATCH PREFIX mesa-12778; export MESA_BASE=\"{mesa_base}\"; export MESA_INLIST=\"$MESA_BASE/inlist\";\n\n")

    # for every mass: set up a run directory and write disBatch command
    for m in masses:
        folder_name = f"M{m}_DW{int(dutch_scaling)}"
        run_dir = os.path.join(output_dir, folder_name)
        os.makedirs(run_dir, exist_ok=True)

        # Copy inlist_gridrun template
        inlist_new = os.path.join(run_dir, "inlist_gridrun")
        shutil.copy(inlist_template, inlist_new)

        # Edit log_directory and parameters
        with open(inlist_new, "r") as f:
            lines = f.readlines()

        new_lines = []
        for line in lines:
            if "initial_mass" in line and "=" in line:
                new_lines.append(f"      initial_mass = {m}\n")
            elif "log_directory" in line and "=" in line:
                new_lines.append(f"      log_directory = 'LOGS/M{m}_DW{int(dutch_scaling)}'\n")
            elif "Dutch_scaling_factor" in line and "=" in line:
                new_lines.append(f"      Dutch_scaling_factor = {dutch_scaling}\n")
            else:
                new_lines.append(line)

        with open(inlist_new, "w") as f:
            f.writelines(new_lines)

        # Write the corresponding disBatch line
        tf.write(
            "ml openblas; export OMP_NUM_THREADS=32; "
            f"cd {run_dir}; "
            f"/usr/bin/time -p {mesa_base}/star &> MESArun.out \n"
        )

print(f"Taskfile created at: {taskfile}")
print(f"Set up {len(masses)} runs in {output_dir}")
print("You can now submit the jobs using: sbatch -p cca -n 4 -c 32 -t 0-04:00 disBatch Taskfile")
print("!! Keep in mind inlists need to point to absolute paths for MESA to find the necessary files !!")
# -n N = Number of parallel disBatch workers (each worker handles one task at a time).
# -c C = Number of CPU cores allocated to each worker/task.
# 