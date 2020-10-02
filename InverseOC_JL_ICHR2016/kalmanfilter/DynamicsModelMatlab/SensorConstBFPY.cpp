#include "SensorConstBFPY.h"


SensorConstBFPY::SensorConstBFPY(SensorAbstract * inner) : SensorDecorator(inner)
{
	//x,y,z world position coordinates
	this->m_size = 1;
	this->type = "ConstBFPY";
	this->typeBin = SensorConstBFPY::type_binary;
}


SensorConstBFPY::~SensorConstBFPY(void)
{
}


void SensorConstBFPY::getMeasurement(rl::math::Vector &mes)
{
	//Call Inner
	SensorDecorator::getMeasurement(mes);
	//Append X position
	mes(this->getMeasurementSize() - this->m_size) = getFrame()->t.translation().y();
}

void SensorConstBFPY::getObservationJacobian(rl::math::Matrix &J)
{
	//Call Inner
	SensorDecorator::getObservationJacobian(J);
	//Append Own notice when we get Measurement Size we will getmeasurement size of all wrapped objects plus this wrapper 
	//a.linear is in sensor frame

	//First get the jacobian in base frame
	rl::math::Matrix Jb;
	this->getBaseJacobian(Jb);
	int row = this->getMeasurementSize() - this->m_size;
	//World Frame Y position
	//dw/dq
	J.block(row, 0, 1, Jb.cols()) = Jb.block(Jb.rows() / 2 + 1, 0, 1, Jb.cols());
}