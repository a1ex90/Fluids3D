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
		// ---CASE 1------------ (TOTAL 16 Cases NOT CHECKED)
		points[1] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
		indices[1] = std::vector<int>{ 0,1,2 };
		points[254] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
		indices[254] = std::vector<int>{ 0,1,2 };

		points[2] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[2] = std::vector<int>{ 0,1,2 };
		points[253] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[253] = std::vector<int>{ 0,1,2 };

		points[4] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[4] = std::vector<int>{ 0,1,2 };
		points[251] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[251] = std::vector<int>{ 0,1,2 };

		points[8] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[8] = std::vector<int>{ 0,1,2 };
		points[247] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[247] = std::vector<int>{ 0,1,2 };

		points[16] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[16] = std::vector<int>{ 0,1,2 };
		points[239] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[239] = std::vector<int>{ 0,1,2 };

		points[32] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[32] = std::vector<int>{ 0,1,2 };
		points[223] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[223] = std::vector<int>{ 0,1,2 };

		points[64] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[64] = std::vector<int>{ 0,1,2 };
		points[191] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[191] = std::vector<int>{ 0,1,2 };

		points[128] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[128] = std::vector<int>{ 0,1,2 };
		points[127] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[127] = std::vector<int>{ 0,1,2 };
		// ---CASE 2------------ (TOTAL 24 Cases NOT CHECKED)
		points[3] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[3] = std::vector<int>{ 0,1,2,1,2,3 };
		points[252] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[252] = std::vector<int>{ 0,1,2,1,2,3 };

		points[6] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f } };
		indices[6] = std::vector<int>{ 0,1,2,1,2,3 };
		points[249] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f } };
		indices[249] = std::vector<int>{ 0,1,2,1,2,3 };

		points[9] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f } };
		indices[9] = std::vector<int>{ 0,1,2,1,2,3 };
		points[246] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f } };
		indices[246] = std::vector<int>{ 0,1,2,1,2,3 };

		points[12] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[12] = std::vector<int>{ 0,1,2,1,2,3 };
		points[243] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[243] = std::vector<int>{ 0,1,2,1,2,3 };

		points[17] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[17] = std::vector<int>{ 0,1,2,1,2,3 };
		points[238] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[238] = std::vector<int>{ 0,1,2,1,2,3 };

		points[34] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[34] = std::vector<int>{ 0,1,2,1,2,3 };
		points[221] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[221] = std::vector<int>{ 0,1,2,1,2,3 };

		points[48] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[48] = std::vector<int>{ 0,1,2,1,2,3 };
		points[207] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[207] = std::vector<int>{ 0,1,2,1,2,3 };

		points[68] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[68] = std::vector<int>{ 0,1,2,1,2,3 };
		points[187] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[187] = std::vector<int>{ 0,1,2,1,2,3 };

		points[96] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[96] = std::vector<int>{ 0,1,2,1,2,3 };
		points[159] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[159] = std::vector<int>{ 0,1,2,1,2,3 };

		points[136] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[136] = std::vector<int>{ 0,1,2,1,2,3 };
		points[119] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[119] = std::vector<int>{ 0,1,2,1,2,3 };

		points[144] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[144] = std::vector<int>{ 0,1,2,1,2,3 };
		points[111] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[111] = std::vector<int>{ 0,1,2,1,2,3 };

		points[192] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[192] = std::vector<int>{ 0,1,2,1,2,3 };
		points[63] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[63] = std::vector<int>{ 0,1,2,1,2,3 };
		// ---CASE 3------------ (TOTAL 24 Cases NOT CHECKED)
		points[5] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[5] = std::vector<int>{ 0,1,2,3,4,5 };
		points[250] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[250] = std::vector<int>{ 0,1,2,3,4,5 };

		points[10] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[10] = std::vector<int>{ 0,1,2,3,4,5 };
		points[245] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[245] = std::vector<int>{ 0,1,2,3,4,5 };

		points[18] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[18] = std::vector<int>{ 0,1,2,3,4,5 };
		points[237] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[237] = std::vector<int>{ 0,1,2,3,4,5 };

		points[24] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[24] = std::vector<int>{ 0,1,2,3,4,5 };
		points[231] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[231] = std::vector<int>{ 0,1,2,3,4,5 };

		points[33] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[33] = std::vector<int>{ 0,1,2,3,4,5 };
		points[222] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[222] = std::vector<int>{ 0,1,2,3,4,5 };

		points[36] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[36] = std::vector<int>{ 0,1,2,3,4,5 };
		points[219] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[219] = std::vector<int>{ 0,1,2,3,4,5 };

		points[66] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[66] = std::vector<int>{ 0,1,2,3,4,5 };
		points[189] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[189] = std::vector<int>{ 0,1,2,3,4,5 };

		points[72] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[72] = std::vector<int>{ 0,1,2,3,4,5 };
		points[183] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[183] = std::vector<int>{ 0,1,2,3,4,5 };

		points[80] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[80] = std::vector<int>{ 0,1,2,3,4,5 };
		points[175] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[175] = std::vector<int>{ 0,1,2,3,4,5 };

		points[129] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[129] = std::vector<int>{ 0,1,2,3,4,5 };
		points[126] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[126] = std::vector<int>{ 0,1,2,3,4,5 };

		points[132] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[132] = std::vector<int>{ 0,1,2,3,4,5 };
		points[123] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[123] = std::vector<int>{ 0,1,2,3,4,5 };

		points[160] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[160] = std::vector<int>{ 0,1,2,3,4,5 };
		points[95] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[95] = std::vector<int>{ 0,1,2,3,4,5 };
		// ---CASE 4------------ (TOTAL 48 Cases NOT CHECKED)
		//8x--Z
		points[19] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[19] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[236] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[236] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[50] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[50] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[205] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[205] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[140] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[140] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[115] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[115] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[196] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[196] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[59] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[59] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[200] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[200] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[55] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[55] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[76] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[76] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[179] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[179] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[49] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[49] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[206] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[206] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[35] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[35] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[220] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[220] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		//8x--Y
		points[38] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[38] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[217] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[217] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[70] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[70] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[185] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[185] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[98] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[98] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[157] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[157] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[100] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[100] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[155] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[155] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[25] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[25] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[230] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[230] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[137] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
		indices[137] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[118] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
		indices[118] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[145] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[145] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[110] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[110] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[152] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[152] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[103] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[103] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		//8x--Z
		points[7] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[7] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[248] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[248] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[11] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[11] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[244] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[244] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[13] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f } };
		indices[13] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[242] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f } };
		indices[242] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[14] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f } };
		indices[14] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[241] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f } };
		indices[241] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[112] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[112] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[143] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[143] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[176] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[176] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[79] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[79] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[208] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f } };
		indices[208] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[47] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f } };
		indices[47] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		points[224] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f } };
		indices[224] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };
		points[31] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f } };
		indices[31] = std::vector<int>{ 0,1,2,1,2,3,1,3,4 };

		// ---CASE 5------------ (TOTAL 6 Cases CHECKED)
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
		
		// ---CASE 6------------ (TOTAL 24 Cases CHECKED)
		//8x--Z
		points[83] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[83] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };
		points[58] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[58] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };

		points[172] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[172] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };
		points[197] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
		indices[197] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };

		points[202] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[202] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };
		points[92] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[92] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };

		points[53] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[53] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };
		points[163] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[163] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };
		//8x--Y
		points[166] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[166] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };
		points[86] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f },  glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[86] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };
		points[106] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[106] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };
		points[101] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
		indices[101] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };

		points[89] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[89] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };
		points[169] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[169] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };
		points[149] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[149] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };
		points[154] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[154] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };
		//8x--Z
		points[135] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[135] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };
		points[75] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[75] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };
		points[45] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[45] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };
		points[30] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[30] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };

		points[120] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[120] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };
		points[180] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[180] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };
		points[210] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[210] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };
		points[225] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
		indices[225] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,5,6,7 };

		// ---CASE 7------------ (TOTAL 2 Cases CHECKED)
		points[90] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[90] = std::vector<int>{ 0,1,2,3,4,5,6,7,8,9,10,11 };
		points[165] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[165] = std::vector<int>{ 0,1,2,3,4,5,6,7,8,9,10,11 };

		// ---CASE 8------------ (TOTAL 8 Cases CHECKED)
		points[216] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{0.0f, 0.5f, 0.0f}, glm::vec3{0.0f, 1.0f, 0.5f} };
		indices[216] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,3,4,5 };
		points[228] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[228] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,3,4,5 };
		points[78] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[78] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,3,4,5 };
		points[141] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f } };
		indices[141] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,3,4,5 };

		points[39] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
		indices[39] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,3,4,5 };
		points[27] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[27] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,3,4,5 };
		points[177] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[177] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,3,4,5 };
		points[114] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f } };
		indices[114] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,3,4,5 };
		// ---CASE 9------------ (TOTAL 12 Cases CHECKED)
		points[212] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[212] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,3,4,5 };
		points[108] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[108] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,3,4,5 };
		points[142] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[142] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,3,4,5 };
		points[201] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[201] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,3,4,5 };

		points[43] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
		indices[43] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,3,4,5 };
		points[147] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[147] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,3,4,5 };
		points[113] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[113] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,3,4,5 };
		points[54] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f } };
		indices[54] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,3,4,5 };

		points[226] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[226] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,3,4,5 };
		points[71] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[71] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,3,4,5 };
		points[29] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f } };
		indices[29] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,3,4,5 };
		points[184] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
		indices[184] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,3,4,5 };
		// ---CASE 10------------ (TOTAL 8 Cases CHECKED)
		points[20] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[20] = std::vector<int>{ 0,1,2,3,4,5 };
		points[235] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[235] = std::vector<int>{ 0,1,2,3,4,5 };

		points[40] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[40] = std::vector<int>{ 0,1,2,3,4,5 };
		points[215] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[215] = std::vector<int>{ 0,1,2,3,4,5 };

		points[65] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[65] = std::vector<int>{ 0,1,2,3,4,5 };
		points[190] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[190] = std::vector<int>{ 0,1,2,3,4,5 };

		points[130] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[130] = std::vector<int>{ 0,1,2,3,4,5 };
		points[125] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[125] = std::vector<int>{ 0,1,2,3,4,5 };
		// ---CASE 11------------ (TOTAL 48 Cases CHECKED)
		points[21] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[21] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[234] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[234] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[22] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[22] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[233] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[233] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[28] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[28] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[227] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[227] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[41] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[41] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[214] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[214] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[42] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[42] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[213] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[213] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[44] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[44] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[211] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[211] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[52] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[52] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[203] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[203] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[56] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[56] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[199] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[199] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[67] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[67] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[188] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[188] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[69] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
		indices[69] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[186] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
		indices[186] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[73] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[73] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[182] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[182] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[81] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[81] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[174] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[174] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[84] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[84] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[171] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[171] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[97] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
		indices[97] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[158] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
		indices[158] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[104] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[104] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[151] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[151] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[131] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f },  glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[131] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[124] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f },  glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[124] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[134] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[134] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[121] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[121] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[138] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[138] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[117] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[117] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[146] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[146] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[109] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[109] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[148] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[148] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[107] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[107] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[162] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[162] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[93] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[93] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[168] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[168] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[87] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[87] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[193] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
		indices[193] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[62] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
		indices[62] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };

		points[194] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[194] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		points[61] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[61] = std::vector<int>{ 0,1,2,1,2,3,4,5,6 };
		// ---CASE 12------------ (TOTAL 16 Cases CHECKED)
		points[26] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[26] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
		points[229] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[229] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };

		points[37] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[37] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
		points[218] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[218] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };

		points[74] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[74] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
		points[181] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[181] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };

		points[82] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[82] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
		points[173] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[173] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };

		points[88] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[88] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
		points[167] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[167] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };

		points[133] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[133] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
		points[122] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[122] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };

		points[161] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[161] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
		points[94] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[94] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };

		points[164] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[164] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
		points[91] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[91] = std::vector<int>{ 0,1,2,3,4,5,6,7,8 };
		// ---CASE 13------------ (TOTAL 6 Cases CHECKED)
		points[195] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[195] = std::vector<int>{ 0,1,2,1,2,3,4,5,6,5,6,7 };

		points[60] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f },  glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[60] = std::vector<int>{ 0,1,2,1,2,3,4,5,6,5,6,7 };

		points[105] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[105] = std::vector<int>{ 0,1,2,1,2,3,4,5,6,5,6,7 };

		points[150] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[150] = std::vector<int>{ 0,1,2,1,2,3,4,5,6,5,6,7 };

		points[170] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[170] = std::vector<int>{ 0,1,2,1,2,3,4,5,6,5,6,7 };

		points[85] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[85] = std::vector<int>{ 0,1,2,1,2,3,4,5,6,5,6,7 };

		// ---CASE 14------------ (TOTAL 12 Cases CHECKED)
		points[232] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f } };
		indices[232] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,2,3,5 };
		points[198] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f } };
		indices[198] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,2,3,5 };
		points[77] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f } };
		indices[77] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,2,3,5 };
		points[156] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f } };
		indices[156] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,2,3,5 };

		points[23] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 0.0f } };
		indices[23] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,2,3,5 };
		points[57] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f } };
		indices[57] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,2,3,5 };
		points[178] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 0.5f, 0.0f, 1.0f } };
		indices[178] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,2,3,5 };
		points[99] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f } };
		indices[99] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,2,3,5 };

		points[209] = std::vector<glm::vec3>{ glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 0.0f }, glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 1.0f } };
		indices[209] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,2,3,5 };
		points[116] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 1.0f }, glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 0.0f } };
		indices[116] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,2,3,5 };
		points[46] = std::vector<glm::vec3>{ glm::vec3{ 1.0f, 0.0f, 0.5f }, glm::vec3{ 1.0f, 0.5f, 1.0f }, glm::vec3{ 0.0f, 0.0f, 0.5f }, glm::vec3{ 0.5f, 1.0f, 0.0f }, glm::vec3{ 0.5f, 1.0f, 1.0f }, glm::vec3{ 0.0f, 0.5f, 0.0f } };
		indices[46] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,2,3,5 };
		points[139] = std::vector<glm::vec3>{ glm::vec3{ 0.5f, 0.0f, 0.0f }, glm::vec3{ 1.0f, 0.5f, 0.0f }, glm::vec3{ 0.5f, 0.0f, 1.0f }, glm::vec3{ 0.0f, 1.0f, 0.5f }, glm::vec3{ 1.0f, 1.0f, 0.5f }, glm::vec3{ 0.0f, 0.5f, 1.0f } };
		indices[139] = std::vector<int>{ 0,1,2,1,2,3,1,3,4,2,3,5 };
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
