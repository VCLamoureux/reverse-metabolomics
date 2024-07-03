# Reverse-metabolomics

This is the code repository to perform reverse metabolomics.

----------------------------------------
### System requirements
This code has been tested on `Apple MacBook Pro` Sonoma version 14.5 (specifications: Apple M2 max, 64 GB RAM, 38 cores GPU, 12 cores CPU) and `Windows 11` version 23H2 (specifications: 13th Gen Intel(R) core (TM) i7-13850HX, 2100 Mhz, 32 GB RAM, 20 cores, 28 logical processors). 

----------------------------------------
### Installation guide
#### R version
The following packages should be installed: 
- data.table (>= 1.15.4)
- tidyverse (>= 2.0.0)
- pheatmap (>= 1.0.12)

The installation time should be about 30 seconds

#### Python version
A Python version >=3.8 is required.
Additionally, the following packages should be installed: 
- pandas >= 1.2
- numpy >= 1.20, !=1.24.0
- seaborn >= 0.13.2
- matplotlib >= 3.4, !=3.6.1

Installation can be performed using the `pip` installer with the command `pip install -r requirements.txt`. Make sure that the `requirements.txt` file is within the working directory when this is performed, or edit the command to include the path to this file. This command should either be entered within the terminal or a Jupyter notebook code cell. If entered within a code cell, the command should be prefaced with `!`.

The installation time should be less than 20 seconds

For users new to Python, we recommend that the code be run on a cloud-based platform such as Google Colab, as this tends to make plot rendering easier and more standardized, in addition to providing the user with more memory, allowing for a greater portion of the metadata table to be read in at once



It should be noted that package installation times may vary depending on the user's network latency and processor
