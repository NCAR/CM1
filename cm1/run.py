import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Dict
from dataclasses import dataclass, field

# You may need to install f90nml: pip install f90nml
import f90nml


# --- Best Practice Improvements ---
@dataclass
class PBS:
    """
    Configuration for a PBS batch job.
    REFAC-NOTE: This class now only contains PBS-specific scheduler information.
    """

    name: str
    account: str
    walltime: str
    nodes: int
    queue: str = "main@desched1"
    cpus_per_node: int = 128
    mpi_procs_per_node: int = 128
    omp_threads: int = 1
    env_vars: Dict[str, str] = field(
        default_factory=lambda: {
            "PALS_PPN": "128",
            "PALS_DEPTH": "1",
            "PALS_CPU_BIND": "depth",
        }
    )


class CM1Run:
    """
    A class to manage the setup, configuration, and execution of a CM1 model run.
    """

    def __init__(
        self,
        cm1_path: Path,
        run_dir: Path,
        executable_path: Path,
        # REFAC-NOTE: pbs_config is now optional.
        pbs_config: Optional[PBS] = None,
        printout: str = "cm1.print.out",
        sounding: Optional[object] = None,
        serial: bool = False,
        background: bool = False,
    ):
        """
        Initialize a CM1 model run.

        :param cm1_path: Path to the root of the CM1 repository.
        :param run_dir: The directory where the model will be run. Its name determines the config to use.
        :param executable_path: Path to the CM1 executable.
        :param pbs_config: An instance of the PBS configuration dataclass. If provided, PBS mode is enabled.
        :param printout: Filename for the standard output log.
        :param sounding: An optional sounding object with a `to_txt()` method.
        :param serial: If True, the model will be run serially on the command line.
        :param background: If True (and serial is True), run the serial job in the background.
        """
        self.cm1_path = Path(cm1_path)
        self.pbs_config = pbs_config
        self.printout = printout
        self.sounding = sounding
        self.serial = serial
        self.background = background
        self.pbs = pbs_config is not None

        if self.serial and self.pbs:
            raise ValueError("Cannot select both serial and PBS execution modes.")

        self.run_dir = Path(run_dir)
        self.run_dir.mkdir(parents=True, exist_ok=True)
        self.executable_path = Path(executable_path)

        if not self.executable_path.is_file():
            raise FileNotFoundError(f"Executable not found at: {self.executable_path}")
        if not os.access(self.executable_path, os.X_OK):
            raise PermissionError(
                f"Executable is not executable: {self.executable_path}"
            )

        # REFAC-NOTE: Derive test case and paths from the run_dir name.
        if not self.run_dir.name.startswith("run_"):
            raise ValueError(
                "run_dir name must start with 'run_' to derive the test case."
            )
        testcase = self.run_dir.name.split("run_")[-1]
        self.config_source_dir = self.cm1_path / "run" / "config_files" / testcase
        self.original_namelist_path = self.config_source_dir / "namelist.input"

        self.readme = self._load_readme()
        self.namelist = self._load_namelist()

    def _load_readme(self) -> str:
        """
        Loads the README file from the derived config source directory.
        """
        readme_path = self.config_source_dir / "README"
        if readme_path.is_file():
            return readme_path.read_text()
        logging.warning(f"No README file found at {readme_path}")
        return ""

    def _load_namelist(self) -> f90nml.Namelist:
        """
        Loads the namelist from the derived original_namelist_path.
        """
        if not self.original_namelist_path.is_file():
            raise FileNotFoundError(
                f"Derived namelist file not found: {self.original_namelist_path}"
            )

        logging.info(f"Loading namelist from: {self.original_namelist_path}")
        return f90nml.read(self.original_namelist_path)

    def generate_pbs_script(self, script_path: str = "job.pbs") -> Path:
        """
        Generate a PBS job script for the CM1 model run.

        :param script_path: Name of the PBS script file to generate.
        :return: Path object of the generated PBS script.
        """
        if not self.pbs_config:
            raise ValueError("Cannot generate PBS script without pbs_config.")

        run_dir_abs = self.run_dir.resolve()
        executable_abs = self.executable_path.resolve()

        env_setup = "\n".join(
            [f"export {key}={value}" for key, value in self.pbs_config.env_vars.items()]
        )

        script_content = f"""#!/bin/bash
# Job Name:
#PBS -N {self.pbs_config.name}
# Account:
#PBS -A {self.pbs_config.account}
# Number of nodes:
#PBS -l select={self.pbs_config.nodes}:ncpus={self.pbs_config.cpus_per_node}:mpiprocs={self.pbs_config.mpi_procs_per_node}:ompthreads={self.pbs_config.omp_threads}
# Wall-clock time:
#PBS -l walltime={self.pbs_config.walltime}
# Queue:
#PBS -q {self.pbs_config.queue}
# Redirect output and error streams:
#PBS -o {run_dir_abs / self.printout}
#PBS -e {run_dir_abs / (self.printout + '.err')}

# Environment setup
module purge
module load ncarenv/24.12
module reset

# Create a temporary directory if needed
mkdir -p $TMPDIR

# MPI environment variables
{env_setup}

# Change to the run directory before execution
cd {run_dir_abs}

# Execute the model
mpiexec --cpu-bind depth {executable_abs}
"""
        script_full_path = self.run_dir / script_path
        script_full_path.write_text(script_content)
        script_full_path.chmod(0o755)
        logging.info(f"Generated PBS script: {script_full_path}")
        return script_full_path

    def prepare_run_dir(self) -> None:
        """
        Prepare the run directory by creating it and copying necessary files.
        Errors out if expected output files already exist.
        """
        run_dir = self.run_dir

        # Check for existing output files before creating/populating the directory.
        expected_outputs = ["cm1out.nc", "cm1out_stats.nc"]
        for filename in expected_outputs:
            output_file = run_dir / filename
            if output_file.exists():
                raise FileExistsError(
                    f"Output file {output_file} already exists. Aborting to prevent overwrite."
                )

        logging.info(f"Preparing run directory: {run_dir}")
        run_dir.mkdir(parents=True, exist_ok=True)

        # REFAC-NOTE: Copy from the derived config source directory.
        if self.config_source_dir.is_dir():
            shutil.copytree(self.config_source_dir, run_dir, dirs_exist_ok=True)
        else:
            logging.warning(
                f"Config source directory not found, skipping copy: {self.config_source_dir}"
            )

        common_files_path = self.cm1_path / "run"
        common_files = ["RRTMG_LW_DATA", "RRTMG_SW_DATA"]
        for file_name in common_files:
            source_file = common_files_path / file_name
            if source_file.is_file():
                shutil.copy(source_file, run_dir / file_name)
            else:
                logging.warning(f"Common file not found, skipping: {source_file}")

        shutil.copy(self.executable_path, run_dir)

        namelist_out = run_dir / "namelist.input"
        # The namelist object loaded from the original path is written to the run directory.
        self.namelist.write(str(namelist_out), force=True)
        logging.info(f"Wrote namelist to {namelist_out}")

        if self.sounding:
            sounding_path = run_dir / "input_sounding"
            sounding_path.write_text(self.sounding.to_txt())
            logging.info(f"Wrote input sounding to {sounding_path}")

    def run(self) -> None:
        """
        Prepares the run directory and executes the model based on instance attributes.
        """
        self.prepare_run_dir()

        if self.pbs:
            script_path = self.generate_pbs_script()
            logging.info(f"Submitting job script: {script_path}")
            subprocess.run(["qsub", str(script_path)], check=True, cwd=self.run_dir)
        elif self.serial:
            executable_local_name = self.executable_path.name
            run_dir = self.run_dir
            printout_path = run_dir / self.printout

            logging.info(
                f"Running serially in {run_dir}. Output will be in {printout_path}"
            )

            with open(printout_path, "w") as f_out:
                command = ["./" + executable_local_name]
                if self.background:
                    subprocess.Popen(
                        command, stdout=f_out, stderr=subprocess.STDOUT, cwd=run_dir
                    )
                    logging.info(f"Process started in background in {run_dir}.")
                else:
                    result = subprocess.run(
                        command,
                        stdout=f_out,
                        stderr=subprocess.STDOUT,
                        check=False,
                        cwd=run_dir,
                    )
                    if result.returncode != 0:
                        logging.error(
                            f"Serial run failed with exit code {result.returncode}. Check {printout_path}."
                        )
                    else:
                        logging.info("Serial run completed successfully.")
