import sys
import os
import subprocess


def generate_alignment(n, directory, filename, outpath):
	#ct_filepath = os.path.join(directory, 'seed_neighbourhoods', 'ct', filename + '.ct')
	neigh_filepath = os.path.join(directory, 'seed_neighbourhoods', 'nei', filename + '.nei')
	single_freq_filepath = os.path.join(directory, 'seed_frequencies', 'single', filename + '.freq')
	doublet_freq_filepath = os.path.join(directory, 'seed_frequencies', 'doublet', filename + '.freq')
	tree_filepath = os.path.join(directory, 'seed_trees', 'fixed', filename + '.seed_tree')
	out_filepath = os.path.join(outpath, filename + '.ali')

	if (not os.path.exists(neigh_filepath)
			or not os.path.exists(single_freq_filepath)
			or not os.path.exists(doublet_freq_filepath)
			or not os.path.exists(tree_filepath)):
		print('Warning: Skipping \'' + filename + '\' as it is missing a neighbourhood, frequency or tree file.')
		return

	with open(neigh_filepath) as nei_file:
		seq_length = sum(1 for _ in nei_file)
	with open(single_freq_filepath) as single_freq_file:
		single_frequencies = single_freq_file.readline()[:-1]
	with open(doublet_freq_filepath) as doublet_freq_file:
		doublet_frequencies = doublet_freq_file.readline()[:-1]
	command = ('./sissi099' +
			  #' -nr' + ct_filepath +
			   ' -nn ' + str(neigh_filepath) +
			   ' -fs ' + single_frequencies +
			   ' -fd ' + doublet_frequencies +
			   ' -l' + str(seq_length) +
			   ' -a' + str(n) +
			   ' ' + str(tree_filepath) +
			   ' > ' + str(out_filepath))
	result = subprocess.run(command, text=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if os.path.getsize(out_filepath) == 0:
		#print('Running \'' + command + '\'')
		print('Error in \'' + filename + '\':\n' + result.stderr)
		os.remove(out_filepath)
	else:
		print('Successfully generated alignment' + ('s' if n > 1 else '') + ' for \'' + filename + '\'.')


def generate_alignments(n, directory, outpath, alignments=None):
	try:
		os.makedirs(os.path.join(outpath, 'alignments'))
	except FileExistsError:
		pass

	if alignments is None:
		filenames = os.listdir(os.path.join(directory, 'seed_trees', 'fixed'))
		for filename in filenames:
			generate_alignment(n, directory, filename.split('.')[0], os.path.join(outpath, 'alignments'))
	else:
		for filename in alignments:
			generate_alignment(n, directory, filename, os.path.join(outpath, 'alignments'))

	print('Done.')


def main():
	if len(sys.argv) < 4:
		print('Usage: ./alignment_generator.py <amount> <rfam path> <outpath> [alignment list]')
		return -1

	n = int(sys.argv[1])
	rfam_path = sys.argv[2]
	outpath = sys.argv[3]
	alignments = None
	if len(sys.argv) > 4:
		alignments = sys.argv[4:]

	print('========== GENERATING SISSI ALIGNMENTS ==========')
	generate_alignments(n, rfam_path, outpath, alignments)


if __name__ == '__main__':
	main()
