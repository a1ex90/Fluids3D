
__kernel void DotProd(__global double* mata, __global double* matb, __global double* outMat) {
	outMat[get_global_id(0)] = mata[get_global_id(0)] * matb[get_global_id(0)];
}

__kernel void NumericalReduction(__global double* data, __local double* localData, __global double* outData) {
	size_t globalId = get_global_id(0);
	size_t localSize = get_local_size(0);
	size_t localId = get_local_id(0);

	localData[localId] = data[globalId];
	//make sure everything sofar has been calculated
	barrier(CLK_LOCAL_MEM_FENCE);

	//try to avoid for loops
	for (int i = localSize >> 1; i > 0; i >>= 1) {
		if (localId < i) {
			localData[localId] += localData[localId + i];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	if (localId == 0) {
		outData[get_group_id(0)] = localData[0];
	}

}