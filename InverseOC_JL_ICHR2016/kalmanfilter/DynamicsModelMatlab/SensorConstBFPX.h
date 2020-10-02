#pragma once
#include "sensordecorator.h"
class SensorConstBFPX :
	public SensorDecorator
{
public:
	SensorConstBFPX(SensorAbstract * inner);
	~SensorConstBFPX(void);
	//Calls same method of wrapped object
	void getMeasurement(rl::math::Vector &mes);
	//Appends Jacobian Part
	void getObservationJacobian(rl::math::Matrix &J);
	static const unsigned int type_binary = (1 << 10);
};
