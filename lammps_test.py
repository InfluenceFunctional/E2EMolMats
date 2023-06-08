import lammps
import os
from pathlib import Path

# path = Path(r"C:\Users\mikem\OneDrive\NYU\CSD\lammps_test\deposit")
# infile = r"in.deposit.molecule.rigid-nvt-small.lmp"
path = Path(r"C:\Users\mikem\crystals\clusters\For_Michael")
infile = r"run_NNout.lmp"

os.chdir(path)
lmp = lammps.lammps()
lmp.file(infile)