import os
import subprocess
import shutil
import numpy as np
import family_generator


def get_paths(rfam_path):
	tree_dirpath = os.path.join(rfam_path, 'seed_trees', 'rescaled')
	sfreq_dirpath = os.path.join(rfam_path, 'seed_frequencies', 'single')
	dfreq_dirpath = os.path.join(rfam_path, 'seed_frequencies', 'doublet')
	return [tree_dirpath, sfreq_dirpath, dfreq_dirpath]


def script(batchsize, epoch, seqlength, alicount):
	tag = 'L' + str(seqlength) + '-b' + str(batchsize) + '-ac' + str(alicount) + '-e' + str(epoch)
	filename = 'AT-' + tag
	filename_lcc = 'AT-' + 'L' + str(seqlength) + '-ac' + str(alicount)
	jobname = 'e' + str(epoch) + '-l' + str(seqlength)

	subprocess.run(['rm', '-rf', str(os.path.join('../examples/slurm', filename + '.slrm'))])
	subprocess.run(['rm', '-rf', str(os.path.join('../examples/script', filename + '.sh'))])

	commands = ('module load miniconda3\n'
				'eval \"$(conda shell.bash hook)\"\n'
				'conda activate rnadeep\n'
				'\n'
				'python ../train_ali.py --ali-dir ../../rnaconv/data/generated/family/' + filename_lcc + '/alignments/ \\\n'
				'\t   --dbn-dir ../../rnaconv/data/generated/family/' + filename_lcc + '/neighbourhoods/dbn/ \\\n'
				'\t   --model-log-dir ../models \\\n'
				'\t   --data-tag sm3-' + tag + '-sissi --smodel 3 --batch-size ' + str(batchsize) + ' --epochs ' + str(epoch) + '\n')

	with open(os.path.join('../examples/slurm', filename + '.slrm'), 'w') as scriptfile:
		scriptfile.write('#!/bin/sh\n'
						 '#SBATCH -J ' + jobname + '\n'
						 '#SBATCH --partition=zen3_0512_a100x2\n'
						 '#SBATCH --qos zen3_0512_a100x2\n'
						 '#SBATCH --gres=gpu:2\n'
						 '#SBATCH --mail-user=a12144285@unet.univie.ac.at\n'
						 '#SBATCH --mail-type=BEGIN,END,FAIL\n'
						 '#SBATCH --output=/home/fs71475/julianz123/workspace/rnadeep/examples/slurm/out/' + filename + '.%A.out\n'
						 '#SBATCH --error=/home/fs71475/julianz123/workspace/rnadeep/examples/slurm/out/' + filename + '.%A.err\n'
						 '\n' +
						 commands)

	with open(os.path.join('../examples/script', filename + '.sh'), 'w') as scriptfile:
		scriptfile.write(commands)


def copy(from_path, to_path, filenames):
	tree_orig_path = 'seed_trees/original'
	tree_fixed_path = 'seed_trees/fixed'
	tree_rescaled_path = 'seed_trees/rescaled'

	ali_path = 'seed_alignments'

	neigh_wuss_path = 'seed_neighbourhoods/wuss'
	neigh_dbn_path = 'seed_neighbourhoods/dbn'
	neigh_ct_path = 'seed_neighbourhoods/ct'
	neigh_nei_path = 'seed_neighbourhoods/nei'

	single_freq_path = 'seed_frequencies/single'
	doublet_freq_path = 'seed_frequencies/doublet'

	os.makedirs(os.path.join(to_path, tree_orig_path), exist_ok=True)
	os.makedirs(os.path.join(to_path, tree_fixed_path), exist_ok=True)
	os.makedirs(os.path.join(to_path, tree_rescaled_path), exist_ok=True)
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
				os.path.exists(os.path.join(from_path, tree_rescaled_path, filename + '.seed_tree')),
				os.path.exists(os.path.join(from_path, ali_path, filename + '.aln')),
				os.path.exists(os.path.join(from_path, neigh_wuss_path, filename + '.wuss')),
				os.path.exists(os.path.join(from_path, neigh_dbn_path, filename + '.dbn')),
				os.path.exists(os.path.join(from_path, neigh_ct_path, filename + '.ct')),
				os.path.exists(os.path.join(from_path, neigh_nei_path, filename + '.nei')),
				os.path.exists(os.path.join(from_path, single_freq_path, filename + '.sfreq')),
				os.path.exists(os.path.join(from_path, doublet_freq_path, filename + '.dfreq'))]
		if all(condition):
			shutil.copy(os.path.join(from_path, tree_orig_path, filename + '.seed_tree'),
						os.path.join(to_path, tree_orig_path, filename + '.seed_tree'))
			shutil.copy(os.path.join(from_path, tree_fixed_path, filename + '.seed_tree'),
						os.path.join(to_path, tree_fixed_path, filename + '.seed_tree'))
			shutil.copy(os.path.join(from_path, tree_rescaled_path, filename + '.seed_tree'),
						os.path.join(to_path, tree_rescaled_path, filename + '.seed_tree'))
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
			shutil.copy(os.path.join(from_path, single_freq_path, filename + '.sfreq'),
						os.path.join(to_path, single_freq_path, filename + '.sfreq'))
			shutil.copy(os.path.join(from_path, doublet_freq_path, filename + '.dfreq'),
						os.path.join(to_path, doublet_freq_path, filename + '.dfreq'))
		else:
			print('Skipping \'' + filename + '\' as it is missing data. (Raw output: ' + str(condition) + ')')


def create(alicount, seqlength, batch_sizes, epochs):
	rfam_path = 'data/rfam'
	generated_path = 'data/generated/family'
	foldername = 'AT-' + 'L' + str(seqlength) + '-ac' + str(alicount)
	subprocess.run(['rm', '-rf', str(os.path.join(rfam_path, foldername))])
	os.makedirs(os.path.join(rfam_path, foldername), exist_ok=False)
	subprocess.run(['rm', '-rf', str(os.path.join(generated_path, foldername))])
	os.makedirs(os.path.join(generated_path, foldername), exist_ok=False)

	filenames = os.listdir(os.path.join(rfam_path, 'full/seed_trees/rescaled'))
	filename = [filenames[np.random.randint(len(filenames))].split('.')[0]]

	copy(os.path.join(rfam_path, 'full'), os.path.join(rfam_path, foldername), filename)
	family_generator.generate_family_set('./sissi099',
										 alicount,
										 seqlength,
										 *get_paths(os.path.join(rfam_path, foldername)),
										 os.path.join(generated_path, foldername))
	for batch_size in batch_sizes:
		for epoch in epochs:
			script(batch_size, epoch, seqlength, alicount)


def main():
	alicounts = [5]
	seqlengths = [285, 290, 295, 300]  # [100, 200, 301, 401, 503, 604, 673]
	batch_sizes = [4]  # [1, 2, 3, 4, 5]
	epochs = [1]

	for alicount in alicounts:
		for seqlength in seqlengths:
			create(alicount, seqlength, batch_sizes, epochs)


if __name__ == '__main__':
	main()
