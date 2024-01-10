#!/usr/bin/env python

import os
import numpy as np
from tensorflow.keras.models import load_model

from rnadeep.metrics import mcc, f1, sensitivity
from rnadeep.data_generators import PaddedAlignmentMatrixEncoding
from rnadeep.sampling_ali import draw_ali_sets

import absl.logging

absl.logging.set_verbosity(absl.logging.ERROR)


def main():
	ali_dir = '../rnaconv/data/small/rfam/seed_frequency/ali/'
	dbn_dir = '../rnaconv/data/small/rfam/seed_neighbourhood/dbn/'
	model = './models/sm0_rfam-sissi_004'  # choose a trained model here
	outdir = f"predictions/{model}-rfamseed"

	evaluate = True
	predict = False
	batch_size = 1

	# Choose a testset.
	[data] = list(draw_ali_sets(ali_dir, dbn_dir))
	[alis, dbrs] = zip(*data)
	tgen = PaddedAlignmentMatrixEncoding(batch_size, alis, dbrs)

	m = load_model(model,
				   custom_objects={"mcc": mcc, "f1": f1,
								   "sensitivity": sensitivity})
	if evaluate:
		print(m.evaluate(tgen))

	if predict:
		mrxs = m.predict(tgen)
		# Convert tf.RaggedTensor to numpy arrays.
		if not isinstance(mrxs, np.ndarray):
			mrxs = mrxs.numpy()
			for i in range(len(mrxs)):
				mrxs[i] = np.asarray([x for x in mrxs[i]], dtype=np.float32)

		if not os.path.exists(outdir):
			os.makedirs(outdir)
		np.save(os.path.join(outdir, 'alignments'), alis)
		np.save(os.path.join(outdir, 'structures'), dbrs)
		np.save(os.path.join(outdir, 'matrices'), mrxs)


if __name__ == '__main__':
	main()
