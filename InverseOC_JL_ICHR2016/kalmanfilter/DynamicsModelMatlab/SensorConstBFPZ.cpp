#include "SensorConstBFPZ.h"


SensorConstBFPZ::SensorConstBFPZ(SensorAbstract * inner) : SensorDecorator(inner)
{
	//x,y,z world position coordinates
	this->m_size = 1;
	this->type = "ConstBFPZ";
	this->typeBin = SensorConstBFPZ::type_binary;
}


SensorConstBFPZ::~SensorConstBFPZ(void)
{
}


void SensorConstBFPZ::getMeasurement(rl::math::Vector &mes)
{
	//Call Inner
	SensorDecorator::getMeasurement(mes);
	//Append X position
	mes(this->getMeasurementSize() - this->m_size) = getFrame()->t.translation().z();
}

void SensorConstBFPZ::getObservationJacobian(rl::math::Matrix &J)
{
	//Call Inner
	SensorDecorator::getObservationJacobian(J);
	//Append Own notice when we get Measurement Size we will getmeasurement size of all wrapped objects plus this wrapper 
	//a.linear is in sensor frame

	//First get the jacobian in base frame
	rl::math::Matrix Jb;
	this->getBaseJacobian(Jb);
	int row = this->getMeasurementSize() - this->m_size;
	//World Frame Z position
	//dw/dq
	J.block(row, 0, 1, Jb.cols()) = Jb.block(Jb.rows() / 2 + 2, 0, 1, Jb.cols());
}