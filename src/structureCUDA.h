/* CUDA Bayesian clustering functions */

#pragma once
#include "vcfpop.h"

#define TARGETCUDA extern "C"

TARGETCUDA void WaitCUDA();

TARGETCUDA void* Malloc(uint64 size);

TARGETCUDA void* MallocHostCUDA(uint64 size);

TARGETCUDA void FreeHostCUDA(void* addr);

TARGETCUDA void* MallocDeviceCUDA(uint64 size);

TARGETCUDA void FreeCUDA(void* addr);

TARGETCUDA void MemsetCUDA(void* addr, int val, uint64 size);

TARGETCUDA void MemcpyCUDA(void* dst, void* src, uint64 size, bool todevice);

TARGETCUDA void FreeStructureMemory(int devID);

TARGETCUDA void CopyStructureMemory(int devID);

TARGETCUDA int GetDeviceCountCUDA();

TARGETCUDA void AllocMemoryCUDA();

TARGETCUDA void FreeMemoryCUDA();

TARGETCUDA void CreateStreamCUDA(int devID);

TARGETCUDA void DestroyStreamCUDA();

TARGETCUDA void ShowDevicesCUDA();

TARGETCUDA void ResetDeviceCUDA();

TARGETCUDA void Svd64CUDA(double* A, double* U, double* S, double* VT, int64 m, int64 n);

TARGETCUDA void Svd32CUDA(float* A, float* U, float* S, float* VT, int64 m, int64 n);

TARGETCUDA void Eig64CUDA(double* A, double* U, double* V, int64 n);

TARGETCUDA void Eig32CUDA(float* A, float* U, float* V, int64 n);

TARGETCUDA void MatrixMul64CUDA(double* A, double* B, double* res, int64 m, int64 n, int64 p, bool Atrans, bool Btrans);

TARGETCUDA void MatrixMul32CUDA(float * A, float * B, float * res, int64 m, int64 n, int64 p, bool Atrans, bool Btrans);