#pragma once

#include "SimUtil.h"

namespace MarchingCubes {
	SimUtil::Mesh3D meshData(SimUtil::Mat3Df &grid, std::vector<std::vector<glm::vec3>> &cubeCases, std::vector<std::vector<int>> &cubeIndices, int width, int height, int depth, float tol);
	void initCase(SimUtil::Mat3Df &grid, int caseNo);
	void corners(std::vector<glm::vec3> &darkDots, std::vector<glm::vec3> &brightDots, int caseNo);
	void initMarchingCubesCases(std::vector<std::vector<glm::vec3>> &points, std::vector<std::vector<int>> &indices);
	int maxSize(int width, int height, int depth);
}
