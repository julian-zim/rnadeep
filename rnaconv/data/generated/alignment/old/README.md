Contains old datasets that were used for training on the zen2 VSC cluster (which ran out of time) together with the respective job outputs.<br>
<br>
These datasets were generated with 10 alignments per datapoint and using an old implementation of the alignment generator that didn't make use of RNAinverse yet and used random ancestral sequences, resulting the data filter to delete much more data, only leaving about 6500 alignments for the 39520 generated alignments left in the filtered dataset. Also, in the old implementation, filtered alignments with only one sequences were still left, so the usable data is even smaller.