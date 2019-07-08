# turtles
TURTLES: Time-sensitive Untemplated Recording using TdT for Local Environmental Signals

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

### Layout
For the shell script used to cut and trim all NGS sequences, see filter_and_trim_TdT.sh. It requires fastqc and cutadapt (v1.16) to be installed.

For Jupyter notebooks, which were used to make most of the Figures, see the Notebooks directory. Each notebook's filename lists the 
figures in that notebook.

Finally, see turtles_utils.py, which contains the source code for the custom methods used in the notebooks.

Unit tests in progress.

### Cite
Bhan, N. J. et al. Recording temporal data onto DNA with minutes resolution. bioRxiv 634790 (2019). doi:10.1101/634790
