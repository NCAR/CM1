import logging
import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Dict
import f90nml
from cm1.input.sounding import Sounding


@dataclass
class PBS:
    """
    Configuration for a PBS batch job.
    This class only contains PBS-specific scheduler information.
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
        pbs_config: Optional[PBS] = None,
        sounding: Optional[Sounding] = None,
    ):
        """
        Initialize a CM1 model run.

        :param cm1_path: Path to the root of the CM1 repository.
        :param run_dir: The directory where the model will be run. Its name determines the config to use.
        :param executable_path: Path to the CM1 executable.
        :param pbs_config: An instance of the PBS configuration dataclass. If provided, PBS mode is enabled.
        :param sounding: An optional sounding object with a `to_txt()` method.
        """
        self.cm1_path = Path(cm1_path)
        self.pbs_config = pbs_config
        self.sounding = sounding

        self.run_dir = Path(run_dir)
        self.run_dir.mkdir(parents=True, exist_ok=True)
        self.executable_path = Path(executable_path)

        if not self.executable_path.is_file():
            raise FileNotFoundError(f"Executable not found at: {self.executable_path}")
        if not os.access(self.executable_path, os.X_OK):
            raise PermissionError(
                f"Executable is not executable: {self.executable_path}"
            )

        # Derive test case and paths from the run_dir name.
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
#PBS -o {run_dir_abs / 'cm1.print.out'}
#PBS -e {run_dir_abs / 'cm1.print.err'}

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
        logging.info(f"Preparing {self.run_dir}")

        # Copy from the derived config source directory.
        if self.config_source_dir.is_dir():
            shutil.copytree(self.config_source_dir, self.run_dir, dirs_exist_ok=True)
        else:
            logging.warning(
                f"Config source directory not found, skipping copy: {self.config_source_dir}"
            )

        # Check for existing output files before creating/populating the directory.
        expected_outputs = ["cm1out.nc", "cm1out_stats.nc"]
        for filename in expected_outputs:
            output_file = self.run_dir / filename
            if output_file.exists():
                raise FileExistsError(
                    f"Output file {output_file} already exists. Aborting to prevent overwrite."
                )

        self.run_dir.mkdir(parents=True, exist_ok=True)

        for common_file in ["RRTMG_LW_DATA", "RRTMG_SW_DATA"]:
            src = self.cm1_path / "run" / common_file
            if src.is_file():
                shutil.copy(src, self.run_dir)
            else:
                logging.warning(f"Common file {common_file} not found.")

        shutil.copy(self.executable_path, self.run_dir)

        if self.namelist:
            # The namelist object loaded from the original path is written to the run directory.
            namelist_out = self.run_dir / "namelist.input"
            self.namelist.write(namelist_out, force=True)
            logging.info(f"Wrote namelist to {namelist_out}")

        if self.sounding:
            sounding_path = self.run_dir / "input_sounding"
            sounding_path.write_text(self.sounding.to_txt())
            logging.info(f"Wrote input sounding to {sounding_path}")

    def run(self) -> None:
        """
        Prepares the run directory and executes the model.
        If pbs_config is present, it submits via qsub.
        Otherwise, it runs the executable directly.
        """
        self.prepare_run_dir()

        if self.pbs_config:
            script_path = self.generate_pbs_script()
            logging.info(f"Submitting job script: {script_path}")
            subprocess.run(["qsub", str(script_path)], check=True, cwd=self.run_dir)
        else:
            printout_path = self.run_dir / "cm1.print.out"

            logging.info(
                f"Running locally in {self.run_dir}. Output will be in {printout_path}"
            )

            with open(printout_path, "w") as f_out:
                result = subprocess.run(
                    self.executable_path,
                    stdout=f_out,
                    stderr=subprocess.STDOUT,
                    cwd=self.run_dir,
                    check=False,
                )
                if result.returncode != 0:
                    logging.error(
                        f"Serial run failed with exit code {result.returncode}. Check {printout_path}."
                    )
                else:
                    logging.info("Serial run completed successfully.")
