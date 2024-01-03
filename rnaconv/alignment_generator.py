import sys
import os
import subprocess
from rfam_converter import convert_rfam_data


def generate_alignment(n, directory, filename, outpath):
	#ct_filepath = base_path + 'Rfam.seed_neighbourhood/ct/' + filename + 'fa.ct'
	neigh_filepath = directory + 'seed_neighbourhood/nei/' + filename + '.nei'
	freq_filepath = directory + 'seed_frequency/freq/' + filename + '.freq'
	tree_filepath = directory + 'seed_tree/fixed/' + filename + '.seed_tree'
	out_filepath = outpath + filename + '.ali'

	with open(neigh_filepath) as nei_file:
		seq_length = sum(1 for _ in nei_file)
	with open(freq_filepath) as freq_file:
		frequencies = freq_file.readline()[:-1]
	command = ('./sissi099' +
			  #' -nr' + ct_filepath +
			   ' -nn ' + neigh_filepath +
			   ' -fs ' + frequencies +
			   ' -l' + str(seq_length) +
			   ' -a' + str(n) +
			   ' ' + tree_filepath +
			   ' > ' + out_filepath)
	result = subprocess.run(command, text=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if os.path.getsize(out_filepath) == 0:
		#print('Running \'' + command + '\'')
		print('Error in \'' + filename + '\':\n' + result.stderr)
		os.remove(out_filepath)
	else:
		print('Successfully generated alignment(s) for \'' + filename + '\'.')


def generate_alignments(n, directory, outpath):
	try:
		os.makedirs(outpath)
	except FileExistsError:
		pass

	filenames = os.listdir(directory + 'seed_tree/fixed/')
	for filename in filenames:
		filename_base = filename.split('.')[0]
		if os.path.exists(directory + 'seed_neighbourhood/nei/' + filename_base + '.nei') \
			and os.path.exists(directory + 'seed_frequency/freq/' + filename_base + '.freq'):
			generate_alignment(n, directory, filename_base, outpath)
		else:
			print('Warning: Skipping \'' + filename_base + '\' as it is missing a neighbourhood or frequency file.')

	print('Done.')


def main():
	if len(sys.argv) != 4:
		print('Usage: ./alignment_generator.py <amount> <rfam path> <outpath>')
		return -1

	n = sys.argv[1]
	rfam_path = sys.argv[2]
	outpath = sys.argv[3]

	if not os.path.exists(rfam_path + 'Rfam.seed') or not os.path.exists(rfam_path + 'seed_tree/original/'):
		print('Required rfam folder structure:\n- rfam\n\tRfam.seed\n\t- seed_tree\n\t\t- original\n\t\t\t<tree files>')
		return -1

	print('========== CONVERTING RFAM DATA POINTS ==========')
	convert_rfam_data(rfam_path + 'Rfam.seed',
					  rfam_path + 'seed_tree/original/',
					  rfam_path + 'seed_tree/fixed/',
					  rfam_path + 'seed_frequency/',
					  rfam_path + 'seed_neighbourhood/')

	print('\n========== GENERATING SISSI ALIGNMENTS ==========')
	generate_alignments(n, rfam_path, outpath)


if __name__ == '__main__':
	main()
