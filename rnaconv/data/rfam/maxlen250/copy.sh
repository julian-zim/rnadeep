#!/usr/bin/bash

# shellcheck disable=SC2046
python ../copy.py ./ ../full $(python get_filenames.py)
