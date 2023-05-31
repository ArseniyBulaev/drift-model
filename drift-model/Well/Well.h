#pragma once

#include <vector>

#include"WellSegment.h"

class Well
{
public:
	Well(const std::vector<WellSegment>& segments);
	const std::vector<WellSegment>& getSegments() const;
	const double GetBottomCrossSectionArea() const;
	const double GetLength() const;

private:
	std::vector<WellSegment> _segments;
	double _length = 0;

	double CaculateLength();
};

