# IDentify EPItopes (IDEPI)


IDEPI is a domain-specific and extensible software library for supervised
learning of models that relate genotype to phenotype for HIV-1 and other
organisms. 

IDEPI provides a programming interface to allow users to engineer sequence
features and select machine algorithms appropriate for their application.

 
##Quickstart

Because some of the dependancies of IDEPI (see below) can be somewhat tricky
to install, we created a virtual machine based on Ubuntu 13. It includes a complete
preconfigured python environment for
use with IDEPI and can be downloaded from [http://hyphy.org/pubs/idepi-vm.tar.gz](http://hyphy.org/pubs/idepi-vm.tar.gz)

> Please note that this VM file is ~2GB in size and make take a bit of 
> time to download

  User: user  
  password: reverse   

To use the virtual machine:

- Decompress the downloaded file and take note of the file path on the host system
- Open VirtualBox (which can be downloaded from <https://www.virtualbox.org/>)
- Click on the 'New' button 
- Type the desired name, select 'Linux' as your operating system type, and 'Ubuntu(64 bit)' as your version.
- Select your desired memory size(at least 1024MB)
- Under Hard Drive, please select 'Use an existing virtual hard drive file' and navigate to the idepi-vm.vdi file that you have just downloaded

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


##Using IDEPI

Please look inside the **examples** directory for examples of how
to use IDEPI to learn and apply phenotype prediction models based on viral
sequences and to identify most predictive features of viral genomes.
