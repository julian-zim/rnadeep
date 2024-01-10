import random

import numpy as np


def profile_vec_encode(char1, char2):
	char1 = char1.upper()
	char2 = char2.upper()

	base_to_ids = {
		'A': [0],
		'C': [1],
		'G': [2],
		'U': [3],
		'-': [4],
		'M': [0, 1],
		'R': [0, 2],
		'W': [0, 3],
		'S': [1, 2],
		'Y': [1, 3],
		'K': [2, 3],
		'V': [0, 1, 2],
		'H': [0, 1, 3],
		'D': [0, 2, 3],
		'B': [1, 2, 3],
		'N': [0, 1, 2, 3],
	}

	try:
		ids1 = base_to_ids[char1]
	except KeyError:
		raise ValueError('Char ' + char1 + ' is not allowed in an alignment sequence!')
	try:
		ids2 = base_to_ids[char2]
	except KeyError:
		raise ValueError('Char ' + char2 + ' is not allowed in an alignment sequence!')

	idxs = list()
	for id1 in ids1:
		for id2 in ids2:
			idxs.append(id1 * 5 + id2)

	vector = np.zeros(25)
	for idx in idxs:
		vector[idx] = 1 / len(idxs)

	return vector.flatten()


def alt_met(char1, char2):
	char1 = char1.upper()
	char2 = char2.upper()

	base_to_ids = {
		'A': np.array([1, 0, 0, 0, 0]),
		'C': np.array([0, 1, 0, 0, 0]),
		'G': np.array([0, 0, 1, 0, 0]),
		'U': np.array([0, 0, 0, 1, 0]),
		'-': np.array([0, 0, 0, 0, 1]),
		'M': np.array([0.5, 0.5, 0, 0, 0]),
		'R': np.array([0.5, 0, 0.5, 0, 0]),
		'W': np.array([0.5, 0, 0, 0.5, 0]),
		'S': np.array([0, 0.5, 0.5, 0, 0]),
		'Y': np.array([0, 0.5, 0, 0.5, 0]),
		'K': np.array([0, 0, 0.5, 0.5, 0]),
		'V': np.array([1./3., 1./3., 1./3., 0, 0]),
		'H': np.array([1./3., 1./3., 0, 1./3., 0]),
		'D': np.array([1./3., 0, 1./3., 1./3., 0]),
		'B': np.array([0, 1./3., 1./3., 1./3., 0]),
		'N': np.array([0.25, 0.25, 0.25, 0.25, 0]),
	}

	return np.outer(base_to_ids[char1], base_to_ids[char2]).flatten()


def main():
	codes = ['C', 'G', 'U', '-', 'M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'N']

	for i in range(50):
		char1 = codes[random.randint(0, len(codes)-1)]
		char2 = codes[random.randint(0, len(codes)-1)]

		res_old = profile_vec_encode(char1, char2)
		res_new = alt_met(char1, char2)

		print(res_old)
		print(res_new)
		print('\n')
		assert all(res_old) == all(res_new)


if __name__ == '__main__':
	main()
