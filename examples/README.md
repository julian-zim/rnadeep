### A typical single sequences example workflow would consist of:

 - *generate_data.py*: generate some random sequence/structure training pairs using RNAfold.

 - *train.py*: select one of the models from rnadeep/models.py (currently those are all SPOTRNA reimplementations.) and train with your generated test/validation data, using pairs of sequences and secondary strucutres.

 - *predict.py*: either test the performance of the model on different data or generate the output matrices from your trained model.

 - *mlforensics.py*: postprocess output matrices to generate secondary structures.


### A typical sequence alignments example workflow would consist of:

 - Generate alignment/consensus structure training pairs using rnaconv or use the converted rfam database directly.

 - *train_ali.py*: select one of the models from rnadeep/models.py (currently those are all SPOTRNA reimplementations.) and train with your generated test/validation data, using pairs of sequence alignments and consensus structures.

 - *predict_ali.py*: either test the performance of the model on different data or generate the output matrices from your trained model.
