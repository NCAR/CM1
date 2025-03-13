import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional

import f90nml


class PBS:
    def __init__(
        self,
        name: str,
        account: str,
        walltime: str,
        nodes: int,
        run_dir: Path,
        executable_path: Path,
    ):
        """
        Initialize PBS job configuration.

        :param name: Job name.
        :param account: Account number for the PBS job.
        :param walltime: Wall clock time for the job (e.g., '02:00:00' for 2 hours).
        :param nodes: Number of nodes to request.
        :param run_dir: Directory where the model run will be executed.
        :param executable_path: Executable path.
        """
        self.name: str = name
        self.account: str = account
        self.walltime: str = walltime
        self.nodes: int = nodes
        self.run_dir: Path = run_dir if isinstance(run_dir, Path) else Path(run_dir)
        self.executable_path: Path = (
            executable_path if isinstance(executable_path, Path) else Path(executable_path)
        )


class CM1Run:
    def __init__(
        self,
        cm1_path: Path,
        pbs_config: PBS,
        namelist_path: Optional[Path] = None,
        printout: Optional[str] = "cm1.print.out",
    ):
        """
        Initialize a CM1 model run.
        Assign a default namelist and readme.

        :param cm1_path: Path to CM1 repository.
        :param pbs_config: Instance of PBS.
        :param namelist_path: Path to namelist.
        :param printout: Filename for the standard output log.
        """
        self.cm1_path: Path = cm1_path
        self.pbs: PBS = pbs_config
        self.printout: str = printout
        self.readme: str = None

        # Path to default namelist and readme
        defaults_path = cm1_path / "run"
        # Change defaults_path if pbs.name matches an existing config_files directory.
        if os.path.exists(defaults_path / "config_files" / self.pbs.name):
            defaults_path = defaults_path / "config_files" / self.pbs.name
        if namelist_path is None:
            logging.warning(f"assign default {self.pbs.name} namelist")
            namelist_path = defaults_path / "namelist.input"
        self.namelist: f90nml.Namelist = f90nml.read(namelist_path)
        with open(defaults_path / "README", "r") as file:
            self.readme = file.read()

    def generate_pbs_script(self, script_path: str = "pbs.job") -> str:
        """
        Generate a PBS job script for the CM1 model run.

        :param script_path: Name of the PBS script file.
        :return: Full path of the generated PBS script.
        """
        script_content = f"""#!/bin/bash
# job name:
#PBS -N {self.pbs.name}
#PBS -A {self.pbs.account}
# Number of nodes:
#PBS -l select={self.pbs.nodes}:ncpus=128:mpiprocs=128:ompthreads=1
# Maximum wall-clock time:
#PBS -l walltime={self.pbs.walltime}
# Queue:
#PBS -q main@desched1
# Redirect output:
#PBS -o {self.pbs.run_dir}
#PBS -e {self.pbs.run_dir}

# Use ncarenv version before it is loaded by default. Executable must have
# been compiled with same modules loaded.
module purge
module load ncarenv/24.12
module reset
module load intel/2025.0.3

export TMPDIR={os.getenv("TMPDIR", Path(os.getenv("SCRATCH")) / "tmp")}
mkdir -p $TMPDIR

export PALS_PPN=128
export PALS_DEPTH=1
export PALS_CPU_BIND=depth

cd {self.pbs.run_dir}

mpiexec --cpu-bind depth {self.pbs.executable_path} >& {self.printout}
"""
        script_full_path = self.pbs.run_dir / script_path
        with open(script_full_path, "w") as script_file:
            script_file.write(script_content)
        return str(script_full_path)

    def prepare_run_dir(self) -> None:
        """
        Prepare the run directory by copying necessary files.
        """
        if self.namelist is None:
            raise ValueError("Namelist must be defined before preparing the run directory.")

        run_dir = self.pbs.run_dir
        run_dir.mkdir(parents=True, exist_ok=True)

        # Copy testcase configuration files to run_dir
        #  input_sounding
        #  LANDUSE.TBL
        #  namelist.input (overwritten below)
        #  README
        shutil.copytree(
            self.cm1_path / "run/config_files" / self.pbs.name, run_dir, dirs_exist_ok=True
        )

        # Copy RRTMG radiation data from cm1_path/run
        for file in ["RRTMG_LW_DATA", "RRTMG_SW_DATA"]:
            shutil.copy(self.cm1_path / "run" / file, run_dir / file)
        # Copy executable to run_dir
        shutil.copy(self.pbs.executable_path, run_dir)

        # Write namelist to run directory
        namelist_out = run_dir / "namelist.input"
        self.namelist.write(namelist_out, force=True)
        logging.info(f"Wrote namelist to {namelist_out}")

        # TODO: function to get expected output files from namelist.
        # Check for existing output file(s). Right now assume one netCDF.
        output_file = self.pbs.run_dir / "cm1out.nc"
        if output_file.exists():
            raise FileExistsError(f"Output file {output_file} already exists.")

    def submit_job(self) -> None:
        """Submit the PBS job."""
        self.prepare_run_dir()
        script_path = self.generate_pbs_script()
        subprocess.run(["qsub", script_path], check=True)

    def run_serial(self, background: bool = False) -> None:
        """
        Run executable serially as from the command line.

        :param background: Whether to run the process in the background.
        """
        self.prepare_run_dir()
        os.chdir(self.pbs.run_dir)
        logging.info(f"stdout and stderr directed to {self.printout}")
        with open(self.printout, "w") as f:
            if background:
                subprocess.Popen([self.pbs.executable_path], stdout=f, stderr=f)
            else:
                subprocess.run([self.pbs.executable_path], stdout=f, stderr=f, check=True)
