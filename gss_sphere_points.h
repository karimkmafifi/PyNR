#ifndef GSS_SPHERE_POINTS_H
#define GSS_SPHERE_POINTS_H
#include "common.h"

struct Gss_sphere_points
{
private:
    std::vector<glm::vec3> n;
public:
    Gss_sphere_points(){}

	void calculate_sphere_points(unsigned int num_sphere_points)
	{
        float golden_angle = pi * (3.0 - std::sqrt(5));
        float	offset = 2.0 / float(num_sphere_points);
		for (int i = 0; i < num_sphere_points; i++)
		{
            float y = i * offset - 1.0 + (offset / 2);
            float r = std::sqrt(1.0 - y*y);
            float phi = i * golden_angle;
            n.push_back(glm::vec3(cos(phi)*r, y, sin(phi)*r));
		}
	}

    inline std::vector<glm::vec3> get_sphere_points()
	{
		return n;
	}
};

#endif
