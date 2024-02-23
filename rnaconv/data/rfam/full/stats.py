import os
from matplotlib import pyplot as plt


def main():
	path = 'seed_alignments'
	filenames = os.listdir(path)
	lengths = list()
	seqcounts = list()
	for filename in filenames:
		with open(os.path.join(path, filename), 'r') as file:
			line = file.readline()
			length = len(file.readline().split()[1])
			lengths.append((length, filename.split('.')[0]))
			seqcount = 1
			while line != '':
				line = file.readline()
				seqcount += 1
			seqcounts.append((seqcount, filename.split('.')[0]))

	lengths.sort(key=lambda x: x[0])
	seqcounts.sort(key=lambda x: x[0])

	print(lengths[0])
	print(lengths[-1])

	print(seqcounts[0])
	print(seqcounts[-1])

	plt.boxplot([e[0] for e in lengths])
	plt.violinplot([e[0] for e in lengths])
	plt.savefig('lenplot')
	plt.clf()
	plt.boxplot([e[0] for e in seqcounts])
	plt.violinplot([e[0] for e in seqcounts])
	plt.savefig('countplot')
	plt.clf()

	cutcounts = [e for e in seqcounts if e[0] <= 500]
	print(cutcounts[-1])
	plt.boxplot([e[0] for e in cutcounts])
	plt.violinplot([e[0] for e in cutcounts])
	plt.show()


if __name__ == '__main__':
	main()
