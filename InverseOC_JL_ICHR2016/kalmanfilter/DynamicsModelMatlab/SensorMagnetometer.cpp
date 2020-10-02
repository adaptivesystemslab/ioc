#include "SensorMagnetometer.h"


SensorMagnetometer::SensorMagnetometer(SensorAbstract * inner) : SensorDecorator(inner)
{
	//Set measurement size
	this->m_size = 3;
	//Set Type
	this->type = "magnetometer";
	this->typeBin = SensorMagnetometer::type_binary;
}


SensorMagnetometer::~SensorMagnetometer(void)
{
}

void SensorMagnetometer::getMeasurement(rl::math::Vector &mes)
{
	//Call Inner
	SensorDecorator::getMeasurement(mes);
	//Append Own notice when we get Measurement Size we will getmeasurement size of all wrapped objects plus this wrapper 
	rl::mdl::Frame * frame = getFrame().get();

	rl::math::Vector3 m;
	SensorDecorator::getMagneticField(m);
	//Magnetic field to sensor frame
	mes.segment(this->getMeasurementSize() - this->m_size, this->m_size) = frame->t.rotation().transpose()*m;

}

void SensorMagnetometer::getObservationJacobian(rl::math::Matrix &J)
{
	//Call Inner
	SensorDecorator::getObservationJacobian(J);
	//Append Own, notice when we get Measurement Size we will getmeasurement size of all wrapped objects plus this wrapper 
	//a.linear is in sensor frame

	//First get the jacobian in end effector frame
	rl::math::Matrix Jm;
	this->getMagnetometerJacobian(Jm);

	int row = this->getMeasurementSize() - this->m_size;
	//d(M)/dq
	J.block(row, 0, Jm.rows(), Jm.cols()) = Jm;
}
