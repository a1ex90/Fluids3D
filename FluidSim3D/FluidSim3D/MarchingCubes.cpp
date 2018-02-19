#include "MarchingCubes.h"

namespace MarchingCubes {
	SimUtil::Mesh3D meshData(SimUtil::Mat3Df &grid, int width, int height, int depth, float tol) {

		int maxGridSize = maxSize(width, height, depth);

		std::vector<glm::vec3> vertices;
		std::vector<glm::vec3> normals;
		std::vector<int> globalIndices;
		int curInd = 0;

		//FOR NORMALS FROM 1 to -2
		for (int i = 1; i < width - 2; i++) {
			for (int j = 1; j < height - 2; j++) {
				for (int k = 1; k < depth - 2; k++) {
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

						int offsetI = (int)offsetX;
						int offsetJ = (int)offsetY;
						int offsetK = (int)offsetZ;

						glm::vec3 normal;

						//the offset of 0.5f indicates that this axis needs interpolation
						if (offsetX == 0.5f) {
							int offsetJ = (int)offsetY;
							int offsetK = (int)offsetZ;
							offsetX = 1 / (grid.get(i + 1, j + offsetJ, k + offsetK) - grid.get(i, j + offsetJ, k + offsetK)) * (tol - grid.get(i, j + offsetJ, k + offsetK));

							normal.x = aboveTol(grid.get(i, j + offsetJ, k + offsetK), tol) - aboveTol(grid.get(i + 1, j + offsetJ, k + offsetK), tol);
							normal.y = 0.5f * (aboveTol(grid.get(i, j + offsetJ - 1, k + offsetK), tol) + aboveTol(grid.get(i + 1, j + offsetJ - 1, k + offsetK), tol)) -
								0.5f * (aboveTol(grid.get(i, j + offsetJ + 1, k + offsetK), tol) + aboveTol(grid.get(i + 1, j + offsetJ + 1, k + offsetK), tol));
							normal.z = 0.5f * (aboveTol(grid.get(i, j + offsetJ, k + offsetK - 1), tol) + aboveTol(grid.get(i + 1, j + offsetJ, k + offsetK - 1), tol)) -
								0.5f * (aboveTol(grid.get(i, j + offsetJ, k + offsetK + 1), tol) + aboveTol(grid.get(i + 1, j + offsetJ, k + offsetK + 1), tol));

						}
						else if (offsetY == 0.5f) {
							int offsetI = (int)offsetX;
							int offsetK = (int)offsetZ;
							offsetY = 1 / (grid.get(i + offsetI, j + 1, k + offsetK) - grid.get(i + offsetI, j, k + offsetK)) * (tol - grid.get(i + offsetI, j, k + offsetK));

							normal.y = aboveTol(grid.get(i + offsetI, j, k + offsetK), tol) - aboveTol(grid.get(i + offsetI, j + 1, k + offsetK), tol);
							normal.x = 0.5f * (aboveTol(grid.get(i + offsetI - 1, j, k + offsetK), tol) + aboveTol(grid.get(i + offsetI - 1, j + 1, k + offsetK), tol)) -
								0.5f * (aboveTol(grid.get(i + offsetI + 1, j, k + offsetK), tol) + aboveTol(grid.get(i + offsetI + 1, j + 1, j + offsetK), tol));
							normal.z = 0.5f * (aboveTol(grid.get(i + offsetI, j, k + offsetK - 1), tol) + aboveTol(grid.get(i + offsetI, j + 1, k + offsetK - 1), tol)) -
								0.5f * (aboveTol(grid.get(i + offsetI, j, k + offsetK + 1), tol) + aboveTol(grid.get(i + offsetI, j + 1, k + offsetK + 1), tol));
						}
						else if (offsetZ == 0.5f) {
							int offsetI = (int)offsetX;
							int offsetJ = (int)offsetY;
							offsetZ = 1 / (grid.get(i + offsetI, j + offsetJ, k + 1) - grid.get(i + offsetI, j + offsetJ, k)) * (tol - grid.get(i + offsetI, j + offsetJ, k));

							normal.z = aboveTol(grid.get(i + offsetI, j + offsetJ, k), tol) - aboveTol(grid.get(i + offsetI, j + offsetJ, k + 1), tol);
							normal.y = 0.5f * (aboveTol(grid.get(i + offsetI, j + offsetJ - 1, k), tol) + aboveTol(grid.get(i + offsetI, j + offsetJ - 1, k + 1), tol)) -
								0.5f * (aboveTol(grid.get(i + offsetI, j + offsetJ + 1, k), tol) + aboveTol(grid.get(i + offsetI, j + offsetJ + 1, k + 1), tol));
							normal.x = 0.5f * (aboveTol(grid.get(i + offsetI - 1, j + offsetJ, k), tol) + aboveTol(grid.get(i + offsetI - 1, j + offsetJ, k + 1), tol)) -
								0.5f * (aboveTol(grid.get(i + offsetI + 1, j + offsetJ, k), tol) + aboveTol(grid.get(i + offsetI + 1, j + offsetJ, k + 1), tol));
						}

						/*float x = 2.0f * (i + offsetX) / (maxGridSize - 1) - 1;
						float y = 2.0f * (j + offsetY) / (maxGridSize - 1) - 1;
						float z = 2.0f * (k + offsetZ) / (maxGridSize - 1) - 1;*/
						//Not correct but better looking with interpolation
						float x = 2.0f * (i + offsetX + 1) / (maxGridSize + 1) - 1;
						float y = 2.0f * (j + offsetY + 1) / (maxGridSize + 1) - 1;
						float z = 2.0f * (k + offsetZ + 1) / (maxGridSize + 1) - 1;

						//normalize normals vector
						normal = glm::normalize(normal);

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

	int aboveTol(float val, float tol) {
		if (val > tol)
			return 1;
		else
			return 0;
	}
}
