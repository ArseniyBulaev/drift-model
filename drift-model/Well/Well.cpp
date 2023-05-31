#include "Well.h"

Well::Well(const std::vector<WellSegment>& segments): _segments(segments)
{
	_length = CaculateLength();
}

const std::vector<WellSegment>& Well::getSegments() const
{
	return _segments;
}

const double Well::GetBottomCrossSectionArea() const
{
	const double pi = 3.1415;
	const double d = _segments.back().diameter; // Диаметр в последнем сегменте трубы
	double s = pi * d * d / 4.0;

	return s;
}

const double Well::GetLength() const
{
	return _length;
}

double Well::CaculateLength()
{
	double total_length = 0;

	for (const WellSegment & segment: _segments)
	{
		total_length += segment.length;
	}

	return total_length;
}
