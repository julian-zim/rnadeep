import os
import subprocess
import shutil
import numpy as np
import alignment_generator


def get_paths(rfam_path):
	tree_dirpath = os.path.join(rfam_path, 'seed_trees', 'rescaled')
	neigh_dirpath = os.path.join(rfam_path, 'seed_neighbourhoods', 'nei')
	sfreq_dirpath = os.path.join(rfam_path, 'seed_frequencies', 'single')
	dfreq_dirpath = os.path.join(rfam_path, 'seed_frequencies', 'doublet')
	ali_dirpath = os.path.join(rfam_path, 'seed_alignments')
	return [tree_dirpath, neigh_dirpath, sfreq_dirpath, dfreq_dirpath, ali_dirpath]


def script(batchsize, epoch, seqlength, alicount):
	tag_nc = 'L' + str(seqlength) + '-b' + str(batchsize)
	tag = 'L' + str(seqlength) + '-b' + str(batchsize) + '-c' + str(alicount) + '-e' + str(epoch)
	filename = 'AT-' + 'L' + str(seqlength) + '-b' + str(batchsize) + '-c' + str(alicount) + '-e' + str(epoch)
	filename_nb = 'AT-' + 'L' + str(seqlength) + '-c' + str(alicount)

	subprocess.run(['rm', '-rf', str(os.path.join('../examples/slurm', filename + '.slrm'))])
	subprocess.run(['rm', '-rf', str(os.path.join('../examples/script', filename + '.sh'))])

	commands = ('module load miniconda3\n'
				'eval \"$(conda shell.bash hook)\"\n'
				'conda activate rnadeep\n'
				'\n'
				'python ../train_ali.py --ali-dir ../../rnaconv/data/generated/alignment/' + filename_nb + '/alignments/ \\\n'
				'\t   --dbn-dir ../../rnaconv/data/rfam/' + filename_nb + '/seed_neighbourhoods/dbn/ \\\n'
				'\t   --model-log-dir ../models \\\n'
				'\t   --data-tag sm3-' + tag + '-sissi --smodel 3 --batch-size ' + str(batchsize) + ' --epochs ' + str(epoch) + '\n')

	with open(os.path.join('../examples/slurm', filename + '.slrm'), 'w') as scriptfile:
		scriptfile.write('#!/bin/sh\n'
						 '#SBATCH -J ' + tag_nc + '\n'
						 '#SBATCH --partition=zen2_0256_a40x2\n'
						 '#SBATCH --qos zen2_0256_a40x2\n'
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
				os.path.exists(os.path.join(from_path, single_freq_path, filename + '.freq')),
				os.path.exists(os.path.join(from_path, doublet_freq_path, filename + '.freq'))]
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
			shutil.copy(os.path.join(from_path, single_freq_path, filename + '.freq'),
						os.path.join(to_path, single_freq_path, filename + '.freq'))
			shutil.copy(os.path.join(from_path, doublet_freq_path, filename + '.freq'),
						os.path.join(to_path, doublet_freq_path, filename + '.freq'))
		else:
			print('Skipping \'' + filename + '\' as it is missing data. (Raw output: ' + str(condition) + ')')


def create(seqlength, alicount, batch_sizes, epochs):
	rfam_path = 'data/rfam'
	generated_path = 'data/generated/alignment'
	foldername = 'AT-L' + str(seqlength) + '-c' + str(alicount)
	subprocess.run(['rm', '-rf', str(os.path.join(rfam_path, foldername))])
	os.makedirs(os.path.join(rfam_path, foldername), exist_ok=False)
	subprocess.run(['rm', '-rf', str(os.path.join(generated_path, foldername))])
	os.makedirs(os.path.join(generated_path, foldername), exist_ok=False)

	filenames = os.listdir(os.path.join(rfam_path, 'full/seed_alignments'))
	alilist = list()
	while len(alilist) < 1:
		id = np.random.choice(range(len(filenames)), 1)[0]
		with open(os.path.join(os.path.join(rfam_path, 'full/seed_alignments'), filenames[id]), 'r') as candidate:
			candidate.readline()
			line = candidate.readline().split()[1]
		if len(line) == seqlength:
			alilist.append(filenames[id].split('.')[0])
		del filenames[id]

	copy(os.path.join(rfam_path, 'full'), os.path.join(rfam_path, foldername), alilist)
	alignment_generator.generate_alignment_set('./sissi099',
											   alicount,
											   *get_paths(os.path.join(rfam_path, foldername)),
											   os.path.join(generated_path, foldername))
	for batch_size in batch_sizes:
		for epoch in epochs:
			script(batch_size, epoch, seqlength, alicount)


def main():
	seqlengths = [100, 200, 301, 401, 503, 604, 673]
	alicounts = [101, 200, 309, 407, 493]
	batch_sizes = [1] # [1, 2, 3, 4, 5]
	epochs = [2] # [1, 2, 3, 4, 5]

	for seqlength in seqlengths:
		for alicount in alicounts:
			create(seqlength, alicount, batch_sizes, epochs)


if __name__ == '__main__':
	main()
