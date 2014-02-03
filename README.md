IDentify EPItopes README
=====================

IDEPI is a domain-specific and extensible software library for supervised
learning of models that relate genotype to phenotype for HIV-1 and other
organisms. 

IDEPI provides a programming interface to allow users to engineer sequence
features and select machine algorithms appropriate for their application.

 
For the impatient
=================

A virtual machine has been created with Ubuntu 13 and python environment for
use with idepi located at http://hyphy.org/pubs/idepi-vm.tar.gz

Virtual Machine Credentials : 
- user: user
- password: reverse

To use the virtual machine:
- Decompress the downloaded file and take note of the file path.
- Open VirtualBox <https://www.virtualbox.org/>
- Click on the 'New' button 
- Type the desired name, select 'Linux' as your operating system type, and 'Ubuntu(64 bit)' as your version.
- Select your desired memory size(at least 1024MB)
- Under Hard Drive, please select 'Use an existing virtual hard drive file' and navigate to the idepi-vm.vdi file that you have just downloaded

Dependencies
============

- Python 3.3
- BioPython, http://biopython.org/wiki/Main_Page
- BioExt, http://github.com/veg/BioExt 
- hppy, http://github.com/veg/hppy
- sklrmmr, http://github.com/nlhepler/sklmrmr
- NumPy, http://www.numpy.org 
- SciPy, http://www.scipy.org 
- six, http://pypi.python.org/pypi/six
- hmmer, http://hmmer.janelia.org/ 


Using IDEPI
=================
Please look inside the <code>examples</code> directory for examples of how
to use IDEPI to learn and apply phenotype prediction models based on viral
sequences and to identify most predictive features of viral genomes.
