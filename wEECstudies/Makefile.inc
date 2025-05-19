# PYTHIA configuration file.
# Generated on Mon Apr 28 08:36:31 AM CDT 2025 with the following command:
# ./configure --with-fastjet3=/data/rke_group/Software/RIVET410/local --with-hepmc3=/data/rke_group/Software/RIVET410/local --with-rivet=/data/rke_group/Software/RIVET410/local --with-yoda=/data/rke_group/Software/RIVET410/local --with-hdf5=/data/rke_group/Software/RIVET410/local --with-highfive=/data/rke_group/Software/RIVET410/local --with-root=/cvmfs/sft.cern.ch/lcg/releases/LCG_107/ROOT/6.34.02/x86_64-el9-gcc13-opt/

# Install directory prefixes.
PREFIX_BIN=/data/rke_group/Software/PYTHIA/pythia8314/bin
PREFIX_INCLUDE=/data/rke_group/Software/PYTHIA/pythia8314/include
PREFIX_LIB=/data/rke_group/Software/PYTHIA/pythia8314/lib
PREFIX_SHARE=/data/rke_group/Software/PYTHIA/pythia8314/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=/cvmfs/sft.cern.ch/lcg/releases/gcc/13.1.0-b3d18/x86_64-el9/bin/g++
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
FASTJET3_BIN=/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/bin/
FASTJET3_INCLUDE=-I/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/include
FASTJET3_LIB=-L/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/lib -Wl,-rpath,/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/lib -lfastjet

HEPMC2_USE=false
HEPMC2_CONFIG=
HEPMC2_BIN=
HEPMC2_INCLUDE=
HEPMC2_LIB=

HEPMC3_USE=true
HEPMC3_CONFIG=HepMC3-config
HEPMC3_BIN=/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/bin/
HEPMC3_INCLUDE=-I/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/include
HEPMC3_LIB=-L/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/lib -Wl,-rpath,/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/lib -lHepMC3

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=false
LHAPDF6_CONFIG=
LHAPDF6_BIN=
LHAPDF6_INCLUDE=
LHAPDF6_LIB=

POWHEG_USE=false
POWHEG_CONFIG=
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=true
RIVET_CONFIG=rivet-config
RIVET_BIN=/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/bin/
RIVET_INCLUDE=-I/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/include
RIVET_LIB=-L/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/lib -Wl,-rpath,/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/lib -lRivet

YODA_USE=true
YODA_CONFIG=yoda-config
YODA_BIN=/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/bin/
YODA_INCLUDE=-I/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/include
YODA_LIB=-L/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/lib -Wl,-rpath,/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/lib -lYODA

ROOT_USE=true
ROOT_CONFIG=root-config
ROOT_BIN=/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.34.02-18eb6/x86_64-el9-gcc13-opt/bin/
ROOT_INCLUDE=-I/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.34.02-18eb6/x86_64-el9-gcc13-opt/include
ROOT_LIB=-L/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.34.02-18eb6/x86_64-el9-gcc13-opt/lib -Wl,-rpath,/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.34.02-18eb6/x86_64-el9-gcc13-opt/lib -lCore

GZIP_USE=false
GZIP_CONFIG=
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

PYTHON_USE=false
PYTHON_CONFIG=
PYTHON_BIN=
PYTHON_INCLUDE=
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
HDF5_BIN=/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/bin/
HDF5_INCLUDE=-I/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/include
HDF5_LIB=-L/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/lib -Wl,-rpath,/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/lib -lhdf5

HIGHFIVE_USE=true
HIGHFIVE_CONFIG=
HIGHFIVE_BIN=/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/bin/
HIGHFIVE_INCLUDE=-I/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/include
HIGHFIVE_LIB=-L/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/lib -Wl,-rpath,/panfs/accrepfs.vampire/data/rke_group/Software/RIVET410/local/lib
