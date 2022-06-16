#!/bin/bash

# git clone https://github.com/Genomicsplc/variantkey.git
cd variantkey/python
make wheel
wheel=$(find . -type f -name "*.whl")
pip install $wheel
