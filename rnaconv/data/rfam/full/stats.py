import os
from matplotlib import pyplot as plt


def main():
	path = 'seed_alignments'
	filenames = os.listdir(path)
	lengths = list()
	seqcounts = list()
	for filename in filenames:
		with open(os.path.join(path, filename), 'r') as file:
			file.readline()
			line = file.readline().split()[1]
			seqcount = 0
			length = len(line)
			lengths.append((length, filename.split('.')[0]))
			while line != '':
				seqcount += 1
				line = file.readline()
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

	cutlengths = [e for e in lengths if e[0] <= 300]
	plt.boxplot([e[0] for e in cutlengths])
	plt.violinplot([e[0] for e in cutlengths])
	plt.savefig('len300plot')
	plt.clf()
	print(len(cutlengths))

	"""
	len85 = [e for e in lengths if abs(e[0] - 85) < 1]
	print(len85)
	len85counts = list()
	for name in [e[1] for e in len85]:
		new = [esq for esq in seqcounts if esq[1] == name]
		if len(new) != 1:
			raise ValueError()
		len85counts.append(*new)

	len85counts.sort(key=lambda x: x[0])
	print(len85counts)

	names11 = [e[1] for e in len85counts if e[0] == 11]
	print(names11)
	names20 = [e[1] for e in len85counts if e[0] == 20]
	print(names20)
	names30 = [e[1] for e in len85counts if e[0] == 30]
	print(names30)
	names40 = [e[1] for e in len85counts if e[0] == 40]
	print(names40)
	names62 = [e[1] for e in len85counts if e[0] == 62]
	print(names62)
	names72 = [e[1] for e in len85counts if e[0] == 72]
	print(names72)

	counts = dict()
	for e in [v[0] for v in len85counts]:
		if e not in counts.keys():
			counts[e] = 1
		else:
			counts[e] += 1

	print(counts)

	# format: seqcount, amount
	res = list(counts.items())
	res.sort(key=lambda x: x[0])
	print(res)
	pass
	"""

	"""cutcounts = [e for e in seqcounts if e[0] <= 500]
	print(cutcounts[-1])
	plt.boxplot([e[0] for e in cutcounts])
	plt.violinplot([e[0] for e in cutcounts])
	plt.show()"""


if __name__ == '__main__':
	main()
