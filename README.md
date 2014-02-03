IDentify EPItopes README
=====================

IDEPI is a domain-specific and extensible software library for supervised
learning of models that relate genotype to phenotype for HIV-1 and other
organisms. 

IDEPI provides a programming interface to allow users to engineer sequence
features and select machine algorithms appropriate for their application.

Using IDEPI
=================
Please look inside the [examples](/examples) directory for examples on how
to use IDEPI to learn and apply phenotype prediction models based on viral
sequences and to identify most predictive features of viral genomes.

 
Virtual Machine
=================

A virtual machine has been created with a Python environment preinstalled with packages that are depended on by IDEPI. The virtual disk image can be found [here](http://hyphy.org/pubs/idepi-vm.tar.gz).

The default user and password for the virtual machine is: 
- user: user
- password: reverse

##Using the virtual machine disk image
- Decompress the downloaded file and take note of the file path
- Open VirtualBox (https://www.virtualbox.org/)
- Create a new virtual machine
- Enter a desired name, select 'Linux' as your operating system type, and 'Ubuntu(64 bit)' as your version.
- Select the desired memory size
- When prompted for hard drive options, please select *Use an existing virtual hard drive file* and navigate to the idepi-vm.vdi file that you have just downloaded and decompressed.
- When the machine has started, double-click the idepi launcher on the desktop

###Note
A Python virtual environment was used to create the IDEPI environment on the virtual machine. The root directory of the Python virtual environment is <code>/home/user/Programming/env/</code>

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



