import os
import subprocess

for file in os.listdir('.'):
    if file.startswith('AT'):
        subprocess.run(['sbatch', file])
