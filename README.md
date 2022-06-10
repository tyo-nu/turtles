# turtles
TURTLES: Time-sensitive Untemplated Recording using TdT for Local Environmental Signals

### Paper

This repository provides all the code and notebooks used in [this paper](https://pubs.acs.org/doi/full/10.1021/jacs.1c07331).

**Title**: Recording temporal data onto DNA with minutes resolution

**Abstract**: Employing DNA as a high-density data storage medium has paved the way for next-generation digital storage and biosensing technologies. However, the multipart architecture of current DNA-based recording techniques renders them inherently slow and incapable of recording fluctuating signals with subhour frequencies. To address this limitation, we developed a simplified system employing a single enzyme, terminal deoxynucleotidyl transferase (TdT), to transduce environmental signals into DNA. TdT adds nucleotides to the 3′-ends of single-stranded DNA (ssDNA) in a template-independent manner, selecting bases according to inherent preferences and environmental conditions. By characterizing TdT nucleotide selectivity under different conditions, we show that TdT can encode various physiologically relevant signals such as Co2+, Ca2+, and Zn2+ concentrations and temperature changes in vitro. Further, by considering the average rate of nucleotide incorporation, we show that the resulting ssDNA functions as a molecular ticker tape. With this method we accurately encode a temporal record of fluctuations in Co2+ concentration to within 1 min over a 60 min period. Finally, we engineer TdT to allosterically turn off in the presence of a physiologically relevant concentration of calcium. We use this engineered TdT in concert with a reference TdT to develop a two-polymerase system capable of recording a single-step change in the Ca2+ signal to within 1 min over a 60 min period. This work expands the repertoire of DNA-based recording techniques by developing a novel DNA synthesis-based system that can record temporal environmental signals into DNA with a resolution of minutes.

### Cite
Recording Temporal Signals with Minutes Resolution Using Enzymatic DNA Synthesis.
Namita Bhan, Alec Callisto, Jonathan Strutz, Joshua Glaser, Reza Kalhor, Edward S. Boyden, George Church, Konrad Kording, and Keith E. J. Tyo
Journal of the American Chemical Society 2021 143 (40), 16630-16640
DOI: 10.1021/jacs.1c07331

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
