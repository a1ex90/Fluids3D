#include "MarchingCubes.h"

namespace MarchingCubes {
	SimUtil::Mesh3D meshData(SimUtil::Mat3Df &grid, int width, int height, int depth, float tol) {

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
					for (int l = 0; l < indiciesLength[selectCase]; l++) {
						globalIndices.push_back(indices[selectCase][l] + curInd);
					}
					for (int l = 0; l < pointsLength[selectCase]; l++) {
						float offsetX = points[selectCase][3 * l] / 2.0f;
						float offsetY = points[selectCase][3 * l + 1] / 2.0f;
						float offsetZ = points[selectCase][3 * l + 2] / 2.0f;

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
			grid.set(0.0f, 0, 1, 1);
			caseNo -= 128;
		}
		if (caseNo >= 64) {
			grid.set(1.0f, 0, 1, 1);
			caseNo -= 64;
		}
		if (caseNo >= 32) {
			grid.set(1.0f, 1, 1, 1);
			caseNo -= 32;
		}
		if (caseNo >= 16) {
			grid.set(0.0f, 1, 1, 1);
			caseNo -= 16;
		}
		if (caseNo >= 8) {
			grid.set(0.0f, 0, 0, 1);
			caseNo -= 8;
		}
		if (caseNo >= 4) {
			grid.set(1.0f, 0, 0, 1);
			caseNo -= 4;
		}
		if (caseNo >= 2) {
			grid.set(1.0f, 1, 0, 1);
			caseNo -= 2;
		}
		if (caseNo >= 1) {
			grid.set(0.0f, 1, 0, 1);
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
