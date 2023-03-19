# Kokkos Project
Kokkos is a computation model for high performance computing (HPC). These abstraction layers are designed specifically to isolate software developers from fluctuation and diversity in hardware details yet provide portability and high levels of performance across many architectures. We use its implementation as an embedded C++ library.

## Important Notes About Kokkos
Remember there are just two important functions (classes) in Kokkos, `Kokkos::parallel_for`(or corresponding `Kokkos::parallel_reduce`) and `Kokkos::View<YourType*>`. The first is to replace for loop, the second is to replace vector or array of c++. See documentation [Kokkos documentation](https://kokkos.org/documentation/) for more details.

## Download and Install on Linux

### Platform
You need to download OpenMP with `sudo apt install libomp-dev`. If you want to use CUDA, it's more recommended to download from `developer.nvidia.com` because CUDA version of apt is 11.5 and you will see some bugs when compile (we will talk about it later). You need to make sure gpu driver has been downloaded before or may you will fail to install CUDA. I used `.run` file to install NVIDIA driver and CUDA.

### NVIDIA Driver and CUDA
You need to find the corresponding file of your gpu on the official website. If you install NVIDIA driver from `.run` file, you may meet some problems. First you need to disable nouveau (the default and open source driver of NVIDIA gpu on Ubuntu). It's recommended to generate your private key of secure boot before you install. There maybe another warning tell that you are using X window. Use:
```sh
# disable nouveau
sudo echo "blacklist nouveau" >> /etc/modprobe.d/blacklist.conf
sudo echo "options nouveau modeset=0" >> /etc/modprobe.d/blacklist.conf
# generate your key, input following the hint
openssl req -new -x509 -newkey rsa:2048 -keyout MOK.priv -outform DER -out MOK.der -nodes -days 36500
# if you haven't install openssl, use sudo apt install openssl
sudo mokutil --import MOK.der
# the password can be any, remember it because it will be used when reboot
sudo bash ./NVIDIA-Linux-x86_64-525.89.02.run -no-x-check
```
After that, choose `sign` and select the generated `MOK.priv`. Reboot and you will see an interface, choose `enroll` and input your password. See more [details](https://wiki.debian.org/SecureBoot). The driver will be load by itself.

To install CUDA, just run `sudo bash ./cuda_12.1.0_530.30.02_linux.run`. Remeber don't choose any driver to install again.

### Kokkos
You can use `git clone https://github.com/kokkos/kokkos` to download and comfiguration with such command:

```sh
cmake .. -DKokkos_ENABLE_OPENMP=ON\
  -DKokkos_ARCH_AMPERE86=ON -DKokkos_ENABLE_CUDA=ON\
  -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON\
  -DKokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON\
  -DKokkos_ENABLE_IMPL_DESUL_ATOMICS=OFF\
  -DCMAKE_CXX_STANDARD=17\
  -DCMAKE_INSTALL_PREFIX=/home/bovera/libkokkos
```

If you don't want to use CUDA, just need to specify `Kokkos_ENABLE_OPENMP=ON`, `Kokkos_ENABLE_SERIAL=ON`, `CMAKE_INSTALL_PREFIX=/some/path`(install location, according to yourself) and `CMAKE_CXX_STANDARD=17`. If CUDA is used, `Kokkos_ENABLE_CUDA=ON`, `Kokkos_ENABLE_CUDA_LAMBDA=ON` and `Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON` must be specified. `Kokkos_ENABLE_IMPL_DESUL_ATOMICS=OFF` and `Kokkos_ARCH_XXX=ON` is needed due to the following reasons.

## Bugs Needing Special Attention

### NVCC Default Arch
*This bug has been reported here: [CMake configuration issue](https://github.com/kokkos/kokkos/issues/5868).*

You may need to set the `default_arch` option of `$KOKKOS_SOURCE_DIR$/bin/nvcc_wrapper` manually even you have already choosed your platform with `Kokkos_ARCH_XXX=ON`. It seems that `nvcc_wrapper` won't follow the settings in `KOKKOS_CUDA_OPTIONS` and you will get this error when setting with cmake:
```sh
nvcc fatal   : Value 'sm_35' is not defined for option 'gpu-architecture'
CMake Error at cmake/kokkos_compiler_id.cmake:12 (STRING):
  STRING sub-command REPLACE requires at least four arguments.
```

Though I don't think I have made such mistakes, it is @dalg24 that said:
> Also this is already cherry-picked into release candidate 3.7.2
> 
> Essentially the arch you specified is wrong which triggers the auto detection attempt that has the sm_35 flag that NVCC did not like. If you fix the typo, you won't get the error.

### Initializion Error
*This bug has been reported here: [Kokkos Atomics Bug](https://github.com/parthenon-hpc-lab/parthenon/issues/720). Also mentioned here: [Cannot Enable CUDA_RELOCATABLE_DEVICE](https://github.com/kokkos/kokkos/issues/5922). Thought @dalg24 said he can't reproduce it. Really?*

When using relocatable device code (`Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON`), `Kokkos::initialize` fails immediately. Just add `Kokkos_ENABLE_IMPL_DESUL_ATOMICS=OFF` to solve.

The error is:
```sh
terminate called after throwing an instance of 'std::runtime_error'
  what():  Desul::Error: init_lock_arrays_cuda: post init kernel error(cudaErrorIllegalAddress): an illegal memory access was encountered
```

### CUDA bug
*This bug has been reported here: [parameter packs not expanded with ‘...’](https://github.com/NVlabs/instant-ngp/issues/119)*

When using CUDA<11.6, you will get like this:
```sh
[  1%] Building CUDA object dependencies/tiny-cuda-nn/src/CMakeFiles/tiny-cuda-nn.dir/common.cu.o
/usr/include/c++/11/bits/std_function.h:435:145: error: parameter packs not expanded with ‘...’:
  435 |         function(_Functor&& __f)
      |                                                                                                                                                 ^
/usr/include/c++/11/bits/std_function.h:435:145: note:         ‘_ArgTypes’
```

You just need to download a higher version of CUDA, or you need to add `noexcept` to corresponding `std_function.h`(search the internet for more details).

## Implicit Capture Of 'this'
You may get an warning when compile:
```sh
warning #20178-D: Implicit capture of 'this' in extended lambda expression
```
This is caused by Kokkos' setting of `KOKKOS_LAMBDA` without explicit capture of 'this'. However C++14 don't surport `KOKKOS_CLASS_LAMBDA`(including explicit capture). Use C++17 (by setting`CMAKE_CXX_STANDARD=17`) instead.