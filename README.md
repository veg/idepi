# IDentify EPItopes (IDEPI)

IDEPI is a domain-specific and extensible software library for supervised
learning of models that relate genotype to phenotype for HIV-1 and other
organisms.  IDEPI makes use of open source libraries for machine learning
(scikit- learn, scikit-learn.org/), sequence alignment (HMMER,
hmmer.janelia.org/), sequence manipulation (BioPython, biopython.org), and
parallelization (joblib, pythonhosted.org/joblib), and provides a programming
interface to allow the users to engineer sequence features and select machine
learning algorithms appropriate for their application.

For more information, a corresponding paper is located at
http://hyphy.org/pubs/idepi/paper.pdf. 


##Using IDEPI

Please look inside the [examples](/examples) directory for examples on how
to use IDEPI to learn and apply phenotype prediction models based on viral
sequences and to identify most predictive features of viral genomes.

 
##Quickstart

Because some of the dependencies of IDEPI (see below) can be somewhat tricky to
install, we've created a virtual machine based on Ubuntu 13. It includes a
complete preconfigured python environment for use with IDEPI and can be
downloaded from <http://hyphy.org/pubs/idepi/idepi-vm.tar.gz>

> Please note that this VM file is ~2GB in size and may take a bit of 
> time to download

Virtual Machine Credentials : 
- user: user
- password: reverse

To use the virtual machine:

- Decompress the downloaded file and take note of the file path on the host system
- Open VirtualBox (which can be downloaded from <https://www.virtualbox.org/>)
- Create a new virtual machine
- Enter a desired name, select 'Linux' as your operating system type, and 'Ubuntu(64 bit)' as your version.
- Select the desired memory size
- When prompted for hard drive options, please select *Use an existing virtual hard drive file* and navigate to the idepi-vm.vdi file that you have just downloaded and decompressed.
- When the machine has started, double-click the _idepi_ launcher on the desktop

> A Python virtual environment was used to create the IDEPI environment on the
> virtual machine. The root directory of the Python virtual environment is
> <code>/home/user/Programming/env/</code>


##Dependencies

- Python 3.3
- BioPython (>=1.58), <http://biopython.org/wiki/Main_Page>
- BioExt (>=0.14.0), <http://github.com/veg/BioExt>
- hppy (>=0.9.5), <http://github.com/veg/hppy>
- sklrmmr (>=0.2.0), <http://github.com/nlhepler/sklmrmr>
- NumPy, <http://www.numpy.org>
- SciPy, <http://www.scipy.org>
- six, <http://pypi.python.org/pypi/six>
- HMMer, <http://hmmer.janelia.org/>
- sklearn (>=0.14.0), <http://scikit-learn.org/> 

