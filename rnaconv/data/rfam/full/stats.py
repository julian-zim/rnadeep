import os
from matplotlib import pyplot as plt


def main():
	path = 'seed_alignments'
	filenames = os.listdir(path)
	lengths = list()
	for filename in filenames:
		with open(os.path.join(path, filename), 'r') as file:
			file.readline()
			length = len(file.readline().split()[1])
			lengths.append((length, filename.split('.')[0]))

	lengths.sort(key=lambda x: x[0])

	print(lengths[0])
	print(lengths[-1])

	plt.boxplot([e[0] for e in lengths])
	plt.violinplot([e[0] for e in lengths])
	plt.show()


if __name__ == '__main__':
	main()
