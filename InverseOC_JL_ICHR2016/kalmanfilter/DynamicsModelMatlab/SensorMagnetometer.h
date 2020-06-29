#pragma once
#include "sensordecorator.h"
class SensorMagnetometer :
	public SensorDecorator
{
public:
	SensorMagnetometer(SensorAbstract * inner);
	~SensorMagnetometer(void);

	//Appends Linear Acceleration to the measurement vector
	void getMeasurement(rl::math::Vector &mes);
	//Appends Jacobian Part
	void getObservationJacobian(rl::math::Matrix &J);
	static const unsigned int type_binary = (1 << 4);
};