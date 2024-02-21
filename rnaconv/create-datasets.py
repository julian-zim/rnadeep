import os
import sys
import numpy as np
import shutil
import alignment_generator
import subprocess


def script(seqlength, alicount):
	tag = 'l' + str(seqlength) + '-c' + str(alicount)
	filename = 'AT-' + tag

	subprocess.run(['rm', '-rf', str(os.path.join('../examples/slurm', filename + '.slrm'))])
	subprocess.run(['rm', '-rf', str(os.path.join('../examples/slurm', filename + '.slrm'))])

	with open(os.path.join('../examples/slurm', filename + '.slrm'), 'w') as scriptfile:
		scriptfile.write('#!/bin/sh\n'
						 '#SBATCH -J ' + tag + '\n'
						 '#SBATCH --partition=zen2_0256_a40x2\n'
						 '#SBATCH --qos zen2_0256_a40x2\n'
						 '#SBATCH --gres=gpu:2\n'
						 '#SBATCH --mail-user=a12144285@unet.univie.ac.at\n'
						 '#SBATCH --mail-type=BEGIN,END,FAIL\n'
						 '#SBATCH --output=/home/fs71475/julianz123/workspace/rnadeep/examples/slurm/out/' + filename + '-bs5.%A.out\n'
						 '#SBATCH --error=/home/fs71475/julianz123/workspace/rnadeep/examples/slurm/out/' + filename + '-bs5.%A.err\n'
						 '\n'
						 'module load miniconda3\n'
						 'eval \"$(conda shell.bash hook)\"\n'
						 'conda activate rnadeep\n'
						 '\n'
						 'python ../train_ali.py --ali-dir ../../rnaconv/data/sissi/' + filename + '/alignments/ \\\n'
						 '--dbn-dir ../../rnaconv/data/rfam/' + filename + '/seed_neighbourhoods/dbn/ \\\n'
						 '--model-log-dir ../models \\\n'
						 '--data-tag sm3-' + tag + '-sissi --smodel 3 --batch-size 5 --epochs 4')


def copy(datasetdir, filenames):
	to_path = os.path.join('data/rfam/', str(datasetdir))
	from_path = 'data/rfam/full'

	tree_orig_path = 'seed_trees/original'
	tree_fixed_path = 'seed_trees/fixed'

	ali_path = 'seed_alignments'

	neigh_wuss_path = 'seed_neighbourhoods/wuss'
	neigh_dbn_path = 'seed_neighbourhoods/dbn'
	neigh_ct_path = 'seed_neighbourhoods/ct'
	neigh_nei_path = 'seed_neighbourhoods/nei'

	single_freq_path = 'seed_frequencies/single'
	doublet_freq_path = 'seed_frequencies/doublet'

	os.makedirs(os.path.join(to_path, tree_orig_path), exist_ok=True)
	os.makedirs(os.path.join(to_path, tree_fixed_path), exist_ok=True)
	os.makedirs(os.path.join(to_path, ali_path), exist_ok=True)
	os.makedirs(os.path.join(to_path, neigh_wuss_path), exist_ok=True)
	os.makedirs(os.path.join(to_path, neigh_dbn_path), exist_ok=True)
	os.makedirs(os.path.join(to_path, neigh_ct_path), exist_ok=True)
	os.makedirs(os.path.join(to_path, neigh_nei_path), exist_ok=True)
	os.makedirs(os.path.join(to_path, single_freq_path), exist_ok=True)
	os.makedirs(os.path.join(to_path, doublet_freq_path), exist_ok=True)

	for filename in filenames:
		condition = [os.path.exists(os.path.join(from_path, tree_orig_path, filename + '.seed_tree')),
				os.path.exists(os.path.join(from_path, tree_fixed_path, filename + '.seed_tree')),
				os.path.exists(os.path.join(from_path, ali_path, filename + '.aln')),
				os.path.exists(os.path.join(from_path, neigh_wuss_path, filename + '.wuss')),
				os.path.exists(os.path.join(from_path, neigh_dbn_path, filename + '.dbn')),
				os.path.exists(os.path.join(from_path, neigh_ct_path, filename + '.ct')),
				os.path.exists(os.path.join(from_path, neigh_nei_path, filename + '.nei')),
				os.path.exists(os.path.join(from_path, single_freq_path, filename + '.freq')),
				os.path.exists(os.path.join(from_path, doublet_freq_path, filename + '.freq'))]
		if all(condition):
			shutil.copy(os.path.join(from_path, tree_orig_path, filename + '.seed_tree'),
						os.path.join(to_path, tree_orig_path, filename + '.seed_tree'))
			shutil.copy(os.path.join(from_path, tree_fixed_path, filename + '.seed_tree'),
						os.path.join(to_path, tree_fixed_path, filename + '.seed_tree'))
			shutil.copy(os.path.join(from_path, ali_path, filename + '.aln'),
						os.path.join(to_path, ali_path, filename + '.aln'))
			shutil.copy(os.path.join(from_path, neigh_wuss_path, filename + '.wuss'),
						os.path.join(to_path, neigh_wuss_path, filename + '.wuss'))
			shutil.copy(os.path.join(from_path, neigh_dbn_path, filename + '.dbn'),
						os.path.join(to_path, neigh_dbn_path, filename + '.dbn'))
			shutil.copy(os.path.join(from_path, neigh_ct_path, filename + '.ct'),
						os.path.join(to_path, neigh_ct_path, filename + '.ct'))
			shutil.copy(os.path.join(from_path, neigh_nei_path, filename + '.nei'),
						os.path.join(to_path, neigh_nei_path, filename + '.nei'))
			shutil.copy(os.path.join(from_path, single_freq_path, filename + '.freq'),
						os.path.join(to_path, single_freq_path, filename + '.freq'))
			shutil.copy(os.path.join(from_path, doublet_freq_path, filename + '.freq'),
						os.path.join(to_path, doublet_freq_path, filename + '.freq'))


def create(seqlength, alicount):
	foldername = 'AT-l' + str(seqlength) + '-c' + str(alicount)
	subprocess.run(['rm', '-rf', str(os.path.join('data/rfam', foldername))])
	os.makedirs(os.path.join('data/rfam', foldername), exist_ok=False)
	subprocess.run(['rm', '-rf', str(os.path.join('data/sissi', foldername))])
	os.makedirs(os.path.join('data/sissi', foldername), exist_ok=False)

	filenames = os.listdir('data/rfam/full/seed_alignments')
	alilist = list()
	while len(alilist) < 1:
		id = np.random.choice(range(len(filenames)), 1)[0]
		with open(os.path.join('data/rfam/full/seed_alignments', filenames[id]), 'r') as candidate:
			candidate.readline()
			line = candidate.readline().split()[1]
		if len(line) == seqlength:
			alilist.append(filenames[id].split('.')[0])
		del filenames[id]

	copy(foldername, alilist)
	alignment_generator.generate_alignments(alicount, os.path.join('data/rfam', foldername), os.path.join('data/sissi', foldername))
	script(seqlength, alicount)


def main():
	seqlengths = [100, 150, 200, 250]
	alicounts = [5, 10, 20, 40]

	for seqlength in seqlengths:
		for alicount in alicounts:
			create(seqlength, alicount)


if __name__ == '__main__':
	main()
