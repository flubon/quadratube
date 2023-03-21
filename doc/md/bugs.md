# Bugs

## Can't use CUDA
I can pass CUDA test in `test-inl.h`. However, when I run the full program, I will get

```sh
cudaDeviceSynchronize() error( cudaErrorIllegalAddress): an illegal memory access was encountered ...kokkos/core/src/Cuda/Kokkos_Cuda_Instance.cpp:161

gdb backtrace:
#17 0x00005555555861d2 in ModelInitializer::Initializer::init (
    this=0x7fffffffd7a0, init_para=...)
    at /home/bovera/quadratube/src/model/initializer3.cpp:263
```

If I delete the parallel for in `initializer3.cpp`, I will get some random `cudaErrorLaunchTimeout` while `update()`. It seems there is no bugs when just using openmp. I just give up using CUDA. Can anyone solve this bug?