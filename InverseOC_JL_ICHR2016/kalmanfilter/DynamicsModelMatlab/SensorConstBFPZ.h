#pragma once
#include "sensordecorator.h"
class SensorConstBFPZ :
	public SensorDecorator
{
public:
	SensorConstBFPZ(SensorAbstract * inner);
	~SensorConstBFPZ(void);
	//Calls same method of wrapped object
	void getMeasurement(rl::math::Vector &mes);
	//Appends Jacobian Part
	void getObservationJacobian(rl::math::Matrix &J);
	static const unsigned int type_binary = (1 << 12);
};
