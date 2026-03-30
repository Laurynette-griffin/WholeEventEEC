# PYTHIA configuration file.
# Generated on Tue May 20 10:22:38 AM CDT 2025 with the following command:
# ./configure --prefix=/data/rke_group/Software/new_PYTHIA/pythia_install --with-root=/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcc12/root/6.28.06 --with-fastjet3=/data/rke_group/Software/new_RIVET/local/ --with-hepmc2=/data/rke_group/Software/new_RIVET/local/ --with-hepmc3=/data/rke_group/Software/new_RIVET/local/ --with-python --with-hdf5=/data/rke_group/Software/new_RIVET/local/ --with-lhapdf6=/data/rke_group/Software/new_RIVET/local/

# Install directory prefixes.
PREFIX_BIN=/data/rke_group/Software/new_PYTHIA/pythia_install/bin
PREFIX_INCLUDE=/data/rke_group/Software/new_PYTHIA/pythia_install/include
PREFIX_LIB=/data/rke_group/Software/new_PYTHIA/pythia_install/lib
PREFIX_SHARE=/data/rke_group/Software/new_PYTHIA/pythia_install/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=g++
CXX_DTAGS=-Wl,--disable-new-dtags
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=

EVTGEN_USE=false
EVTGEN_CONFIG=
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

FASTJET3_USE=true
FASTJET3_CONFIG=fastjet-config
FASTJET3_BIN=/panfs/accrepfs.vampire/data/rke_group/Software/new_RIVET/local/bin/
FASTJET3_INCLUDE=-I/panfs/accrepfs.vampire/data/rke_group/Software/new_RIVET/local/include
FASTJET3_LIB=-L/panfs/accrepfs.vampire/data/rke_group/Software/new_RIVET/local/lib -Wl,-rpath,/panfs/accrepfs.vampire/data/rke_group/Software/new_RIVET/local/lib -lfastjet

HEPMC2_USE=false
HEPMC2_CONFIG=
HEPMC2_BIN=
HEPMC2_INCLUDE=
HEPMC2_LIB=

HEPMC3_USE=true
HEPMC3_CONFIG=HepMC3-config
HEPMC3_BIN=/panfs/accrepfs.vampire/data/rke_group/Software/new_RIVET/local/bin/
HEPMC3_INCLUDE=-I/panfs/accrepfs.vampire/data/rke_group/Software/new_RIVET/local/include
HEPMC3_LIB=-L/panfs/accrepfs.vampire/data/rke_group/Software/new_RIVET/local/lib -Wl,-rpath,/panfs/accrepfs.vampire/data/rke_group/Software/new_RIVET/local/lib -lHepMC3

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=true
LHAPDF6_CONFIG=lhapdf-config
LHAPDF6_BIN=/panfs/accrepfs.vampire/data/rke_group/Software/new_RIVET/local/bin/
LHAPDF6_INCLUDE=-I/panfs/accrepfs.vampire/data/rke_group/Software/new_RIVET/local/include
LHAPDF6_LIB=-L/panfs/accrepfs.vampire/data/rke_group/Software/new_RIVET/local/lib -Wl,-rpath,/panfs/accrepfs.vampire/data/rke_group/Software/new_RIVET/local/lib -lLHAPDF

POWHEG_USE=false
POWHEG_CONFIG=
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=false
RIVET_CONFIG=
RIVET_BIN=
RIVET_INCLUDE=
RIVET_LIB=

YODA_USE=false
YODA_CONFIG=
YODA_BIN=
YODA_INCLUDE=
YODA_LIB=

ROOT_USE=true
ROOT_CONFIG=root-config
ROOT_BIN=/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcc12/root/6.28.06/bin/
ROOT_INCLUDE=-I/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcc12/root/6.28.06/include
ROOT_LIB=-L/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcc12/root/6.28.06/lib -Wl,-rpath,/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcc12/root/6.28.06/lib -lCore

GZIP_USE=false
GZIP_CONFIG=
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

PYTHON_USE=true
PYTHON_CONFIG=python-config
PYTHON_BIN=
PYTHON_INCLUDE=-I/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/python/3.11.5/include/python3.11 -I/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/python/3.11.5/include/python3.11
PYTHON_LIB=

MG5MES_USE=false
MG5MES_CONFIG=
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=

MPICH_USE=false
MPICH_CONFIG=
MPICH_BIN=
MPICH_INCLUDE=
MPICH_LIB=

HDF5_USE=true
HDF5_CONFIG=
HDF5_BIN=/panfs/accrepfs.vampire/data/rke_group/Software/new_RIVET/local/bin/
HDF5_INCLUDE=-I/panfs/accrepfs.vampire/data/rke_group/Software/new_RIVET/local/include
HDF5_LIB=-L/panfs/accrepfs.vampire/data/rke_group/Software/new_RIVET/local/lib -Wl,-rpath,/panfs/accrepfs.vampire/data/rke_group/Software/new_RIVET/local/lib -lhdf5

HIGHFIVE_USE=false
HIGHFIVE_CONFIG=
HIGHFIVE_BIN=
HIGHFIVE_INCLUDE=
HIGHFIVE_LIB=
