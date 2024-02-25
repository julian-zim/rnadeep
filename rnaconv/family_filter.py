import sys
import os
import numpy as np
import RNA


default_max_dbrs_deviation = 20  # in percent


def filter_alignments(path, max_dbrs_deviation):

	for ali_filename in os.listdir(os.path.join(path, 'alignments')):
		filename = ali_filename.split('.')[0]
		with open(os.path.join(path, 'neighbourhoods/dbn', filename + '.dbn'), 'r') as dbnfile:
			cons_dbrs = dbnfile.readlines()[-1].split()[0]

		with open(os.path.join(path, 'alignments', ali_filename), 'r') as alifile:
			ali = alifile.readlines()

		blacklist = list()
		for i, seq in enumerate([''.join(line.split()[1:]) for line in ali[1:]]):
			dbrs, _ = RNA.fold(seq)

			bp_diff = RNA.bp_distance(dbrs, cons_dbrs)
			if bp_diff / len(dbrs) > max_dbrs_deviation / 100:
				blacklist.append(i)
				print('Deleting sequence ' + str(i) + ' of alignment ' + filename + ' due to its secondary structure '
									'predicted by RNAfold deviating '
									'by over ' + str(max_dbrs_deviation) + '%.')
		for i in reversed(blacklist):
			del ali[i+1]

		if len(ali) == 1:
			os.remove(os.path.join(path, 'alignments', ali_filename))
			os.remove(os.path.join(path, 'sequences', filename + '.seq'))
			os.remove(os.path.join(path, 'neighbourhoods/dbn', filename + '.dbn'))
			os.remove(os.path.join(path, 'neighbourhoods/ct', filename + '.ct'))
		else:
			with open(os.path.join(path, 'alignments', ali_filename), 'w') as alifile:
				alifile.write(''.join(ali))

		'''command = 'RNAalifold --noPS < ' + str(os.path.join(sissi_path, 'alignments', ali_filename))
		result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		if len(result.stdout) == 0:
			raise RuntimeError('RNAalifold failed to predict the secondary structure for \'' + filename + '\':\n' + result.stderr)
		mfe_dbrs = result.stdout.split('\n')[1].split()[0]

		dbrs_diff = RNA.bp_distance(dbrs, mfe_dbrs)
		if dbrs_diff / len(dbrs) > max_dbrs_deviation / 100:
			os.remove(os.path.join(sissi_path, 'alignments', ali_filename))
			print('Removed alignment \'' + filename + '\' due to the secondary structure predicted by RNAalifold deviating '
									'by over ' + str(max_dbrs_deviation) + '%.')'''

	# TODO: filter alignments who show a frequence distance aboe a certain threshold
	# or, you know, just dont allow the difference values to be higher or equal to 0.1

	#diffs = (list(), list())
	#diffs[0].append(read_frequencies(os.path.join(sissi_path, 'frequencies/differences/single', filename + '.freq')))
	#diffs[1].append(read_frequencies(os.path.join(sissi_path, 'frequencies/differences/doublet', filename + '.freq')))
	#diffs_flattened = ([vector for diffset in diffs[0] for vector in diffset],
	#					[vector for diffset in diffs[1] for vector in diffset])
	pass


def main():
	if not 2 <= len(sys.argv) <= 3:
		print('Usage: ./generator_filter.py <path> [<max_ss_deviation>]')
		return -1
	path = sys.argv[1]
	max_dbrs_deviation = None
	if len(sys.argv) == 3:
		max_dbrs_deviation = sys.argv[2]

	print('========== FILTERING SISSI ALIGNMENTS ==========')
	#obtain_sissi_frequencies(rfam_path, sissi_path)
	filter_alignments(path, default_max_dbrs_deviation if max_dbrs_deviation is None else max_dbrs_deviation)
	print('Done.')


if __name__ == '__main__':
	main()
