#pragma once

#include <glm\glm.hpp>
#include <glm\gtx\norm.hpp>
class Nucleus
{
public:
	Nucleus(glm::dvec3 pos, double charge = 1.0);
	~Nucleus();

	glm::dvec3 pos;
	double charge;
};

