#!/usr/bin/env python

import os
import argparse
from tensorflow.keras.models import load_model
from tensorflow.keras.callbacks import CSVLogger, ModelCheckpoint

from rnadeep import __version__
from rnadeep.models import spotrna_alignment_models
from rnadeep.metrics import mcc, f1, sensitivity
from rnadeep.data_generators import PaddedAlignmentMatrixEncoding
from rnadeep.sampling_ali import draw_ali_sets

import absl.logging

absl.logging.set_verbosity(absl.logging.ERROR)


def training(datatag, dbn_dir, ali_dir,
			 spotmodel=None,
			 basemodel=None,
			 savedir='.',
			 epochs=50,
			 epoch0=0,
			 batch_size=4):
	assert '_' not in datatag, f"forbidden character '_' in {datatag = }"
	assert (spotmodel is not None) or (basemodel is not None)
	#
	# Model Setup
	#
	if basemodel is None:
		assert epoch0 == 0
		logname = os.path.join(savedir, f"sm{spotmodel}_{datatag}")
		model = spotrna_alignment_models(spotmodel, True)
		model.compile(optimizer="adam",
					  loss="binary_crossentropy",
					  metrics=["acc", mcc, f1, sensitivity],
					  run_eagerly=True)
	else:
		bm, dt, ep = os.path.basename(basemodel).split('_')
		assert (spotmodel is None) or (bm == f'sm{spotmodel}')
		if dt == datatag:
			assert int(ep) == epoch0
			logname = os.path.join(savedir, f"{bm}_{dt}")
		else:
			assert epoch0 == 0
			logname = os.path.join(savedir, f"{bm}_{dt}_{ep}_{datatag}")
		model = load_model(basemodel,
						   custom_objects={"mcc": mcc, "f1": f1,
										   "sensitivity": sensitivity})
	#
	# Model training
	#
	# Callback functions for training.
	csv_logger = CSVLogger(f"{logname}.csv", separator=';', append=True)
	model_checkpoint = ModelCheckpoint(filepath=logname + '_{epoch:03d}',
									   save_weights_only=False)
	# Get the data for analysis
	[train, valid] = list(draw_ali_sets(ali_dir, dbn_dir, [0.8, 0.2]))
	[train_alis, train_dbrs] = zip(*train)
	train_generator = PaddedAlignmentMatrixEncoding(batch_size, train_alis, train_dbrs)
	[valid_alis, valid_dbrs] = zip(*valid)
	valid_generator = PaddedAlignmentMatrixEncoding(batch_size, valid_alis, valid_dbrs)

	model.fit(
		x=train_generator,
		validation_data=valid_generator,
		initial_epoch=epoch0,
		epochs=epochs,
		shuffle=True,
		verbose=1,
		callbacks=[csv_logger, model_checkpoint]
	)
	#
	# Save final model
	#
	model.save(f"RNAdeep-{logname}-ep{epochs}")
	return model


def parse_rnadeep_args(p):
	"""Arguments that are used by RNAdeep."""
	p.add_argument('--version', action='version',
				   version='%(prog)s ' + __version__)
	p.add_argument("-d", "--data-tag", action='store', required=True,
				   metavar='<str>', default=None,
				   help="Provide a tag to describe training data.")
	p.add_argument("-dd", "--dbn-dir", action='store', required=True,
				   metavar='<str>', default=None,
				   help="Provide a directory with dot bracket notation data for training and validation.")
	p.add_argument("-sd", "--ali-dir", action='store', required=True,
				   metavar='<str>', default=None,
				   help="Provide a directory with alignment data for training and validation.")
	p.add_argument("-m", "--smodel", type=int,
				   metavar='<int>', default=None,
				   choices=(0, 1, 2, 3, 4),
				   help="Select a model for training.")
	p.add_argument("-l", "--load-model", action='store',
				   metavar='<str>', default=None,
				   help="Provide a pretrained model for continued training.")
	p.add_argument("-s", "--model-log-dir", action='store',
				   metavar='<str>', default='.',
				   help="Store intermediate models and csv output in this directory.")
	p.add_argument("-e", "--epochs", type=int,
				   metavar='<int>', default=50,
				   help="Train for --epochs.")
	p.add_argument("-e0", "--epoch0", type=int,
				   metavar='<int>', default=0,
				   help="Start at this epoch.")
	p.add_argument("-b", "--batch-size", type=int,
				   metavar='<int>', default=4,
				   help="Batch size.")


def main():
	""" RNAdeep training interface.

	"""
	parser = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description='RNAdeep: RNA folding using deep learning.')
	parse_rnadeep_args(parser)
	args = parser.parse_args()
	#
	# Argument Parsing
	#
	if '_' in args.data_tag:
		raise SystemExit('--data-tag must not contain the following character: "_"')
	if args.smodel is args.load_model is None:
		raise SystemExit(('Load exisiting model with --load-model ',
						  'or specify which model to train with --smodel'))
	m = training(args.data_tag,
				 args.dbn_dir,
				 args.ali_dir,
				 spotmodel=args.smodel,
				 basemodel=args.load_model,
				 savedir=args.model_log_dir,
				 epochs=args.epochs,  # select the epoch where you want to stop training
				 epoch0=args.epoch0,  # select the epoch where you want to pick up training
				 batch_size=args.batch_size)


if __name__ == '__main__':
	main()
