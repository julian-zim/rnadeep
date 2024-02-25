import os
import subprocess

for file in os.listdir('.'):
    if file.startswith('M-'):
        subprocess.run(['sbatch', file])
