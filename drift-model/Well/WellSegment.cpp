#include "WellSegment.h"

WellSegment::WellSegment(double l, double theta, double eps, double d):
	length(l),
	tilt_angle(theta),
	relative_roughness(eps),
	diameter(d)
{
}
