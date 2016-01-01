#pragma once
#include <glm\glm.hpp>
#include <glm\gtx\norm.hpp>
#include <array>

const int NUM_PARTICLES = 1;

class SystemState
{
public:
	SystemState();
	~SystemState();
	
	std::array<glm::dvec3, NUM_PARTICLES> positions;
};

