# turtles
TURTLES: Time-sensitive Untemplated Recording using TdT for Local Environmental Signals

### Getting Started

Clone this repo using

`git clone https://github.com/tyo-nu/turtles`

Then, make a virtual environment. If using Conda, make sure to first install pip (`conda install pip`). Then, do

`pip install -r requirements.txt`

Then, make sure you are in the turtles/ directory and run

`pytest`

to run all tests.

Note: To import the turtles package, make sure the directory containing the turtles folder is in your PATH.

### Layout

* filter_and_trim_TdT.sh - shell script used to cut and trim all NGS sequences. It requires fastqc and cutadapt (v1.16) to be installed.

* Notebooks/ - Jupyter notebooks which were used to make most of the Figures. Each notebook's filename lists the 
figures in that notebook.

* turtles_utils.py - contains the source code for the custom methods used in the notebooks.

* tests/ - directory for tests. Tests are written in test_turtles_utils.py using pytest. Test datafiles in the tests/ directory were taken from the 0 control, 1 control, and 20 min switch (0 to 1) fastq files for Cobalt (paper Figure 2). Only the first 1000 sequences were copied to each test datafile so that tests run more quickly.

### Paper

This repository provides all the code and notebooks used in [this paper](https://www.biorxiv.org/content/10.1101/634790v3), 
currently on bioRxiv.

**Title**: Recording temporal data onto DNA with minutes resolution

**Abstract**: Recording biological signals can be difficult in three-dimensional matrices, such as
tissue. We present a DNA polymerase-based strategy that records temporal biosignals locally onto DNA to be read out
later, which could obviate the need to extract information from tissue on the fly. We use a template-independent DNA polymerase, 
terminal deoxynucleotidyl transferase (TdT) that probabilistically adds dNTPs to single- stranded DNA (ssDNA) substrates without a 
template. We show that in vitro, the dNTP- incorporation preference of TdT changes with the presence of Co2+, Ca2+, Zn2+ and 
temperature. Extracting the signal profile over time is possible by examining the dNTP incorporation preference along the length of 
synthesized ssDNA strands like a molecular ticker tape. We call this TdT-based untemplated recording of temporal local
environmental signals (TURTLES). We show that we can determine the time of Co2+ addition to within two minutes over a 
60-minute period. Further, TURTLES has the capability to record multiple fluctuations. We can estimate the rise and fall of an 
input Co2+ pulse to within three minutes. TURTLES has at least 200-fold better temporal resolution than all previous DNA-based 
recording techniques.

### Cite
Bhan, N. J. et al. Recording temporal data onto DNA with minutes resolution. bioRxiv 634790 (2019). doi:10.1101/634790
