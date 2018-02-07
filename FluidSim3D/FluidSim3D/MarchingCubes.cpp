#include "MarchingCubes.h"



namespace MarchingCubes {
	SimUtil::Mesh3D meshData(SimUtil::Mat3Df &grid, std::vector<std::vector<glm::vec3>> &cubeCases, std::vector<std::vector<int>> &cubeIndices, int width, int height, int depth, float tol) {

		int maxGridSize = maxSize(width, height, depth);

		std::vector<glm::vec3> vertices;
		std::vector<glm::vec3> normals;
		std::vector<int> globalIndices;
		int curInd = 0;

		//FOR NORMALS FROM 1 to -2
		for (int i = 0; i < width - 1; i++) {
			for (int j = 0; j < height - 1; j++) {
				for (int k = 0; k < depth - 1; k++) {
					//determine which of the 256 different cases for one square exists
					int selectCase = 0;
					if (grid.get(i, j, k + 1) > tol)
						selectCase += 128;
					if (grid.get(i + 1, j, k + 1) > tol)
						selectCase += 64;
					if (grid.get(i + 1, j + 1, k + 1) > tol)
						selectCase += 32;
					if (grid.get(i, j + 1, k + 1) > tol)
						selectCase += 16;
					if (grid.get(i, j, k) > tol)
						selectCase += 8;
					if (grid.get(i + 1, j, k) > tol)
						selectCase += 4;
					if (grid.get(i + 1, j + 1, k) > tol)
						selectCase += 2;
					if (grid.get(i, j + 1, k) > tol)
						selectCase += 1;
					for (int l = 0; l < cubeIndices[selectCase].size(); l++) {
						globalIndices.push_back(cubeIndices[selectCase][l] + curInd);
					}
					for (int l = 0; l < cubeCases[selectCase].size(); l++) {
						float offsetX = cubeCases[selectCase][l].x;
						float offsetY = cubeCases[selectCase][l].y;
						float offsetZ = cubeCases[selectCase][l].z;

						glm::vec3 normal;
						//offset Interpolation wrong!

						//the offset of 0.5f indicates that this axis needs interpolation
						/*if (offsetX == 0.5f) {
							int offsetJ = (int)offsetY;
							int offsetK = (int)offsetZ;
							offsetX = 1 / (grid.get(i + 1, j + offsetJ, k + offsetK) - grid.get(i, j + offsetJ, k + offsetK)) * (tol - grid.get(i, j + offsetJ, k + offsetK));

							normal.x = offsetX * grid.get(i, j + offsetJ, k + offsetK) - (1.0f - offsetX) * grid.get(i + 1, j + offsetJ, k + offsetK);
							normal.y = offsetX * (grid.get(i, j + offsetJ - 1, k + offsetK) - grid.get(i, j + offsetJ + 1, k + offsetK)) / 2.0 +
								(1.0f - offsetX) * (grid.get(i + 1, j + offsetJ - 1, k + offsetK) - grid.get(i + 1, j + offsetJ + 1, k + offsetK)) / 2.0;
							normal.z = offsetX * (grid.get(i, j + offsetJ, k + offsetK - 1) - grid.get(i, j + offsetJ, k + offsetK + 1)) / 2.0 +
								(1.0f - offsetX) * (grid.get(i + 1, j + offsetJ, k + offsetK - 1) - grid.get(i + 1, j + offsetJ, k + offsetK + 1)) / 2.0;

						}
						else if (offsetY == 0.5f) {
							int offsetI = (int)offsetX;
							int offsetK = (int)offsetZ;
							offsetY = 1 / (grid.get(i + offsetI, j + 1, k + offsetK) - grid.get(i + offsetI, j, k + offsetK)) * (tol - grid.get(i + offsetI, j, k + offsetK));

							normal.y = offsetY * grid.get(i + offsetI, j, k + offsetK) - (1.0f - offsetY) * grid.get(i + offsetI, j + 1, k + offsetK);
							normal.x = offsetY * (grid.get(i + offsetI - 1, j, k + offsetK) - grid.get(i + offsetI + 1, j, k + offsetK)) / 2.0 +
								(1.0f - offsetY) * (grid.get(i + offsetI - 1, j + 1, k + offsetK) - grid.get(i + offsetI + 1, j + 1, k + offsetK)) / 2.0;
							normal.z = offsetY * (grid.get(i + offsetI, j, k + offsetK - 1) - grid.get(i + offsetI, j, k + offsetK + 1)) / 2.0 +
								(1.0f - offsetY) * (grid.get(i + offsetI, j + 1, k + offsetK - 1) - grid.get(i + offsetI, j + 1, k + offsetK + 1)) / 2.0;
						}
						else if (offsetZ == 0.5f) {
							int offsetI = (int)offsetX;
							int offsetJ = (int)offsetY;
							offsetY = 1 / (grid.get(i + offsetI, j + offsetJ, k + 1) - grid.get(i + offsetI, j + offsetJ, k)) * (tol - grid.get(i + offsetI, j + offsetJ, k));

							normal.z = offsetZ * grid.get(i + offsetI, j + offsetJ, k) - (1.0f - offsetZ) * grid.get(i + offsetI, j + offsetJ, k + 1);
							normal.x = offsetZ * (grid.get(i + offsetI - 1, j + offsetJ, k) - grid.get(i + offsetI + 1, j + offsetJ, k)) / 2.0 +
								(1.0f - offsetZ) * (grid.get(i + offsetI - 1, j + offsetJ, k + 1) - grid.get(i + offsetI + 1, j + offsetJ, k + 1)) / 2.0;
							normal.y = offsetZ * (grid.get(i + offsetI, j + offsetJ - 1, k) - grid.get(i + offsetI, j + offsetJ + 1, k)) / 2.0 +
								(1.0f - offsetZ) * (grid.get(i + offsetI, j + offsetJ - 1, k + 1) - grid.get(i + offsetI, j + offsetJ + 1, k + 1)) / 2.0;
						}*/

						float x = 2.0f * (i + offsetX) / (maxGridSize - 1) - 1;
						float y = 2.0f * (j + offsetY) / (maxGridSize - 1) - 1;
						float z = 2.0f * (k + offsetZ) / (maxGridSize - 1) - 1;

						//normalize normals vector
						glm::normalize(normal);

						vertices.push_back(glm::vec3{ x,y,z });
						normals.push_back(normal);

						//increment the current Index
						curInd++;
					}
				}
			}
		}
		return SimUtil::Mesh3D(vertices, normals, globalIndices);
	}

	void initCase(SimUtil::Mat3Df &grid, int caseNo) {
		if (caseNo > 255 || caseNo < 0) {
			return;
		}
		if (caseNo >= 128) {
			grid.set(0, 0, 1, 1);
			caseNo -= 128;
		}
		if (caseNo >= 64) {
			grid.set(1, 0, 1, 1);
			caseNo -= 64;
		}
		if (caseNo >= 32) {
			grid.set(1, 1, 1, 1);
			caseNo -= 32;
		}
		if (caseNo >= 16) {
			grid.set(0, 1, 1, 1);
			caseNo -= 16;
		}
		if (caseNo >= 8) {
			grid.set(0, 0, 0, 1);
			caseNo -= 8;
		}
		if (caseNo >= 4) {
			grid.set(1, 0, 0, 1);
			caseNo -= 4;
		}
		if (caseNo >= 2) {
			grid.set(1, 1, 0, 1);
			caseNo -= 2;
		}
		if (caseNo >= 1) {
			grid.set(0, 1, 0, 1);
			caseNo -= 1;
		}
	}

	void corners(std::vector<glm::vec3> &darkDots, std::vector<glm::vec3> &brightDots, int caseNo) {
		if (caseNo >= 128) {
			darkDots.push_back(glm::vec3{ -1,-1,1 });
			caseNo -= 128;
		}
		else {
			brightDots.push_back(glm::vec3{ -1,-1,1 });
		}
		if (caseNo >= 64) {
			darkDots.push_back(glm::vec3{ 1,-1,1 });
			caseNo -= 64;
		}
		else {
			brightDots.push_back(glm::vec3{ 1,-1,1 });
		}
		if (caseNo >= 32) {
			darkDots.push_back(glm::vec3{ 1,1,1 });
			caseNo -= 32;
		}
		else {
			brightDots.push_back(glm::vec3{ 1,1,1 });
		}
		if (caseNo >= 16) {
			darkDots.push_back(glm::vec3{ -1,1,1 });
			caseNo -= 16;
		}
		else {
			brightDots.push_back(glm::vec3{ -1,1,1 });
		}
		if (caseNo >= 8) {
			darkDots.push_back(glm::vec3{ -1,-1,-1 });
			caseNo -= 8;
		}
		else {
			brightDots.push_back(glm::vec3{ -1,-1,-1 });
		}
		if (caseNo >= 4) {
			darkDots.push_back(glm::vec3{ 1,-1,-1 });
			caseNo -= 4;
		}
		else {
			brightDots.push_back(glm::vec3{ 1,-1,-1 });
		}
		if (caseNo >= 2) {
			darkDots.push_back(glm::vec3{ 1,1,-1 });
			caseNo -= 2;
		}
		else {
			brightDots.push_back(glm::vec3{ 1,1,-1 });
		}
		if (caseNo >= 1) {
			darkDots.push_back(glm::vec3{ -1,1,-1 });
			caseNo -= 1;
		}
		else {
			brightDots.push_back(glm::vec3{ -1,1,-1 });
		}
	}


	void initMarchingCubesCases(std::vector<std::vector<glm::vec3>> &points, std::vector<std::vector<int>> &indices) {
		points.reserve(256);
		indices.reserve(256);
		for (int i = 0; i < 256; i++) {
			points.push_back(std::vector<glm::vec3> {});
			indices.push_back(std::vector<int> {});
		}
		//There are 15 ambigious cases for the 256 total cases
		// ---CASE 0------------
		points[0] = std::vector<glm::vec3>{};
		indices[0] = std::vector<int>{};
		points[255] = std::vector<glm::vec3>{};
		indices[255] = std::vector<int>{};
		// ---CASE 1------------
		points[1] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
		indices[1] = std::vector<int>{ 0,1,2 };
		points[2] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[2] = std::vector<int>{ 0,1,2 };
		points[4] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[4] = std::vector<int>{ 0,1,2 };
		points[8] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[8] = std::vector<int>{ 0,1,2 };
		points[16] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[16] = std::vector<int>{ 0,1,2 };
		points[32] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[32] = std::vector<int>{ 0,1,2 };
		points[64] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[64] = std::vector<int>{ 0,1,2 };
		points[128] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[128] = std::vector<int>{ 0,1,2 };
		// ---CASE 2------------
		points[3] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[3] = std::vector<int>{ 0,1,2,1,2,3 };
		points[6] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[6] = std::vector<int>{ 0,1,2,1,2,3 };
		points[9] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[9] = std::vector<int>{ 0,1,2,1,2,3 };
		points[12] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[12] = std::vector<int>{ 0,1,2,1,2,3 };
		points[17] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[17] = std::vector<int>{ 0,1,2,1,2,3 };
		points[34] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[34] = std::vector<int>{ 0,1,2,1,2,3 };
		points[48] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[48] = std::vector<int>{ 0,1,2,1,2,3 };
		points[68] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[68] = std::vector<int>{ 0,1,2,1,2,3 };
		points[96] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[96] = std::vector<int>{ 0,1,2,1,2,3 };
		points[136] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[136] = std::vector<int>{ 0,1,2,1,2,3 };
		points[144] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[144] = std::vector<int>{ 0,1,2,1,2,3 };
		points[192] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[192] = std::vector<int>{ 0,1,2,1,2,3 };
		// ---CASE 3------------
		points[5] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[5] = std::vector<int>{ 0,1,2,3,4,5 };
		points[10] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[10] = std::vector<int>{ 0,1,2,3,4,5 };
		points[24] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[24] = std::vector<int>{ 0,1,2,3,4,5 };
		points[38] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[38] = std::vector<int>{ 0,1,2,3,4,5 };
		points[66] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[66] = std::vector<int>{ 0,1,2,3,4,5 };
		points[80] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[80] = std::vector<int>{ 0,1,2,3,4,5 };
		points[129] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[129] = std::vector<int>{ 0,1,2,3,4,5 };
		points[160] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[160] = std::vector<int>{ 0,1,2,3,4,5 };
		// ---CASE 4------------
		points[19] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[19] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[50] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[50] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[140] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[140] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[196] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[196] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[98] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[98] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[38] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[38] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };




		// ---CASE 5------------
		points[15] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[15] = std::vector<int>{ 0,1,2,1,3,2 };
		points[51] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[51] = std::vector<int>{ 0,1,2,0,2,3 };
		points[102] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[102] = std::vector<int>{ 0,1,2,1,2,3 };
		points[153] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[153] = std::vector<int>{ 0,1,2,1,2,3 };
		points[204] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[204] = std::vector<int>{ 0,1,2,0,2,3 };
		points[240] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[240] = std::vector<int>{ 0,1,2,1,3,2 };

		// ---CASE 7------------
		points[90] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[90] = std::vector<int>{ 0,1,2,3,4,5,6,7,8,9,10,11 };
		points[165] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[165] = std::vector<int>{ 0,1,2,3,4,5,6,7,8,9,10,11 };


		// ---CASE 10------------
		points[20] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[20] = std::vector<int>{ 0,1,2,3,4,5 };
		points[40] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[40] = std::vector<int>{ 0,1,2,3,4,5 };
		points[65] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[65] = std::vector<int>{ 0,1,2,3,4,5 };
		points[130] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[130] = std::vector<int>{ 0,1,2,3,4,5 };
		// ---CASE 11------------
		points[21] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[21] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[22] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[22] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[28] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[28] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[41] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[41] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[42] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[42] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[44] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[44] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[52] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[52] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[56] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[56] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[67] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[67] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[69] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
		indices[69] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[73] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[73] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[81] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[81] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[84] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[84] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[97] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
		indices[97] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[104] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[104] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[131] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f },  glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[131] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[134] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[134] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[138] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[138] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[146] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[146] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[148] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[148] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[162] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[162] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[168] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[168] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[193] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
		indices[193] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[194] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[194] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		// ---CASE 12------------
		points[37] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[37] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
		points[74] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[74] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
		points[82] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[82] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
		points[88] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[88] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
		points[133] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[133] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
		points[138] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[138] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
		points[161] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[161] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
		points[164] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[164] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };

	}

	int maxSize(int width, int height, int depth) {
		if (width > height) {
			if (width > depth)
				return width;
			else
				return depth;
		}
		else {
			if (height > depth)
				return height;
			else
				return depth;
		}
	}
}
