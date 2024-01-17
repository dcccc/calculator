# 1. download orca 
# download orca from https://orcaforum.kofo.mpg.de/app.php/dlext/, before that, 
# you need to register an account. orca is free for academic use.
# here orca_5_0_4_linux_x86-64_shared_openmpi411.tar.xz is downloaded.


# 2.  install openmpi
# openmpi 4.1.1 is required by orca_5_0_4_linux_x86-64_shared_openmpi411.
# before installing openmpi, you need to make sure gcc, g++ and gfortran are
# installed.

wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.1.tar.gz
tar -xf openmpi-4.1.1.tar.gz
cd openmpi-4.1.1
./configure --prefix=/home/username/openmpi411  CC=gcc CXX=g++ FC=gfortran
make -j 4
make install

# the path of openmpi execution file mpirun is in  /home/username/openmpi411/openmpi411/bin/
# and the path of openmpi library is in /home/username/openmpi411/openmpi411/lib/
# add the path of openmpi execution file mpirun to the PATH environment 
# if you want to use mpirun directly, you need to add the paths of openmpi 
# execution file and openmpi library to the PATH environment variable
# for example:
# pwd=`pwd`
# export PATH=/home/username/openmpi411/bin:$PATH
# export LD_LIBRARY_PATH=/home/username/openmpi411/lib:$LD_LIBRARY_PATH


# 3. install orca
# unzip the downloaded file
tar -xf orca_5_0_4_linux_x86-64_shared_openmpi411.tar.xz 

# add the path of orca execution files and library files to the PATH environment variable
pwd=`pwd`
export PATH=$pwd/orca_5_0_4_linux_x86-64_shared_openmpi411:$PATH
export LD_LIBRARY_PATH=$pwd/orca_5_0_4_linux_x86-64_shared_openmpi411:$LD_LIBRARY_PATH

