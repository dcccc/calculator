

# currenty directory

pwd=`pwd`

# install openmpi

wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.1.tar.gz
tar -zxf openmpi-4.1.1.tar.gz
cd openmpi-4.1.1
./configure --prefix=$pwd/gromacs_2022_install/openmpi_1.1.1 CC=gcc CXX=g++ 
make -j 16
make install 
cd ../
# set enviroment variables for openmpi

export PATH=$PATH:$pwd/gromacs_2022_install/openmpi_1.1.1/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$pwd/gromacs_2022_install/openmpi_1.1.1/lib


# install cmake

wget https://github.com/Kitware/CMake/releases/download/v3.24.2/cmake-3.24.2-linux-x86_64.tar.gz
tar -zxf cmake-3.24.2-linux-x86_64.tar.gz

cmake=$pwd/cmake-3.24.2-linux-x86_64/bin/cmake


# compile gromacs 

wget ftp://ftp.gromacs.org/gromacs/gromacs-2022.tar.gz
tar -zxf gromacs-2022.tar.gz

cd gromacs-2022

# ntmpi mixed precision verison
mkdir build_sp 
cd build_sp
$cmake .. -DCMAKE_INSTALL_PREFIX=$pwd/gromacs_2022_install/sp -DCMAKE_C_COMPILER=gcc  -DCMAKE_CXX_COMPILER=g++  -DGMX_MPI=off  -DGMX_GPU=off -DGMX_DOUBLE=off  -DGMX_BUILD_OWN_FFTW=ON 
make -j 16
make install 
cd ../


# ntmpi double precision verison
mkdir build_dp 
cd build_dp
$cmake .. -DCMAKE_INSTALL_PREFIX=$pwd/gromacs_2022_install/dp -DCMAKE_C_COMPILER=gcc  -DCMAKE_CXX_COMPILER=g++  -DGMX_MPI=off  -DGMX_GPU=off -DGMX_DOUBLE=on  -DGMX_BUILD_OWN_FFTW=ON 
make -j 16
make install 
cd ../

# openmpi mixed precision verison
mkdir build_sp_mpi 
cd build_sp_mpi
$cmake .. -DCMAKE_INSTALL_PREFIX=$pwd/gromacs_2022_install/sp_mpi -DCMAKE_C_COMPILER=mpicc  -DCMAKE_CXX_COMPILER=mpicxx  -DGMX_MPI=on  -DGMX_GPU=off -DGMX_DOUBLE=on  -DGMX_BUILD_OWN_FFTW=ON 
make -j 16
make install 
cd ../


# ntmpi double precision verison
mkdir build_dp_mpi 
cd build_dp_mpi 
$cmake .. -DCMAKE_INSTALL_PREFIX=$pwd/gromacs_2022_install/dp_mpi -DCMAKE_C_COMPILER=mpicc  -DCMAKE_CXX_COMPILER=mpicxx  -DGMX_MPI=on  -DGMX_GPU=off -DGMX_DOUBLE=on  -DGMX_BUILD_OWN_FFTW=ON 
make -j 16
make install 
cd ../
cd ../

# install nvcc

wget https://developer.download.nvidia.com/compute/cuda/11.3.0/local_installers/cuda_11.3.0_465.19.01_linux.run
sh cuda_11.3.0_465.19.01_linux.run  --silent --toolkit --toolkitpath=$pwd/cuda_11.3/


# ntmmpi  gpu  verison
cd gromacs_2022_install/
mkdir build_gpu 
cd build_gpu 
$cmake .. -DCMAKE_INSTALL_PREFIX=$pwd/gromacs_2022_install/gpu -DCMAKE_C_COMPILER=gcc  -DCMAKE_CXX_COMPILER=g++  -DGMX_MPI=off  -DGMX_GPU=CUDA  -DCUDA_TOOLKIT_ROOT_DIR=$pwd/cuda_11.3/ -DGMX_DOUBLE=off  -DGMX_BUILD_OWN_FFTW=ON 
make -j 16
make install 
cd ../

# openmpi  gpu  verison

mkdir build_mpi_gpu 
cd build_mpi_gpu 
$cmake .. -DCMAKE_INSTALL_PREFIX=$pwd/gromacs_2022_install/gpu_mpi -DCMAKE_C_COMPILER=mpicc  -DCMAKE_CXX_COMPILER=mpicxx   -DGMX_GPU=CUDA  -DCUDA_TOOLKIT_ROOT_DIR=$pwd/cuda_11.3/  -DGMX_GPU=off -DGMX_DOUBLE=on  -DGMX_BUILD_OWN_FFTW=ON 
make -j 16
make install 
cd ../
cd ../






