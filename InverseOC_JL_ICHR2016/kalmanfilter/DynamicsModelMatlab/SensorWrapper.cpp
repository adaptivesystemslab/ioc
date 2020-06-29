#include <rl\math\Vector.h>
#include "SensorWrapper.h"
//Decorated Sensor Stuff
#include "SensorAbstract.h"
#include "SensorCore.h"
#include "SensorGyroscope.h"
#include "SensorAccelerometer.h"
#include "SensorMagnetometer.h"
#include "SensorOrientation.h"
#include "SensorPosition.h"
#include "SensorQuaternion.h"
#include "SensorVelocity.h"
#include "SensorAngularVelocity.h"
#include "SensorYaw.h"
#include "SensorConstBFPX.h"
#include "SensorConstBFPY.h"
#include "SensorConstBFPZ.h"

void SensorWrapper(mxArray *plhs[],const mxArray *prhs[])
{
	int command = (int)mxGetScalar(prhs[0]);

	rl::math::Vector return_vector;
	rl::math::Matrix return_matrix;

	SensorAbstract * sens;

	if(command != 0 && command != 1)
		sens = (*(SensorAbstract**)mxGetData(prhs[1]));

	switch(command)
	{
		case 0:
			SensorWrapper_Usage();
			break;

		case 1:
			//Create new sensor core, this core can then be wrapped by decorators
			{
				//Get the name
				char name[256];
				mxGetString(prhs[1],name,256);
				//Create the sensor
				SensorAbstract * core = new SensorCore(name);
				//Return the pointer
				plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
				*(intptr_t *)mxGetData(plhs[0]) = (intptr_t)core;
			}
			break;
		case 2:
			//Destructor
			delete sens;
			break;
		case 3:
			//Wrap sensor in a decorator
			{
				//Lets get the type of decorator 
				char type[256];
				mxGetString(prhs[2],type,256);

				if(!std::strcmp(type,"gyroscope"))
				{
					SensorAbstract * decorated_sens = new SensorGyroscope(sens);
					plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
					*(intptr_t *)mxGetData(plhs[0]) = (intptr_t)decorated_sens;
				} 
				else if(!std::strcmp(type,"accelerometer"))
				{
					SensorAbstract * decorated_sens = new SensorAccelerometer(sens);
					plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
					*(intptr_t *)mxGetData(plhs[0]) = (intptr_t)decorated_sens;
				} 
				else if (!std::strcmp(type, "magnetometer"))
				{
					SensorAbstract * decorated_sens = new SensorMagnetometer(sens);
					plhs[0] = mxCreateNumericMatrix(1, 1, mxINDEX_CLASS, mxREAL);
					*(intptr_t *)mxGetData(plhs[0]) = (intptr_t)decorated_sens;
				}
				else if(!std::strcmp(type,"orientation"))
				{
					SensorAbstract * decorated_sens = new SensorOrientation(sens);
					plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
					*(intptr_t *)mxGetData(plhs[0]) = (intptr_t)decorated_sens;
				}
				else if(!std::strcmp(type,"position"))
				{
					SensorAbstract * decorated_sens = new SensorPosition(sens);
					plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
					*(intptr_t *)mxGetData(plhs[0]) = (intptr_t)decorated_sens;
				}
				else if(!std::strcmp(type,"quaternion"))
				{
					SensorAbstract * decorated_sens = new SensorQuaternion(sens);
					plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
					*(intptr_t *)mxGetData(plhs[0]) = (intptr_t)decorated_sens;
				}
				else if(!std::strcmp(type,"velocity"))
				{
					SensorAbstract * decorated_sens = new SensorVelocity(sens);
					plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
					*(intptr_t *)mxGetData(plhs[0]) = (intptr_t)decorated_sens;
				}
				else if(!std::strcmp(type,"angularvelocity"))
				{
					SensorAbstract * decorated_sens = new SensorAngularVelocity(sens);
					plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
					*(intptr_t *)mxGetData(plhs[0]) = (intptr_t)decorated_sens;
				}
				else if(!std::strcmp(type,"yaw"))
				{
					SensorAbstract * decorated_sens = new SensorYaw(sens);
					plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
					*(intptr_t *)mxGetData(plhs[0]) = (intptr_t)decorated_sens;
				}
				else if (!std::strcmp(type, "constbfpx"))
				{
					SensorAbstract * decorated_sens = new SensorConstBFPX(sens);
					plhs[0] = mxCreateNumericMatrix(1, 1, mxINDEX_CLASS, mxREAL);
					*(intptr_t *)mxGetData(plhs[0]) = (intptr_t)decorated_sens;
				}
				else if (!std::strcmp(type, "constbfpy"))
				{
					SensorAbstract * decorated_sens = new SensorConstBFPY(sens);
					plhs[0] = mxCreateNumericMatrix(1, 1, mxINDEX_CLASS, mxREAL);
					*(intptr_t *)mxGetData(plhs[0]) = (intptr_t)decorated_sens;
				}
				else if (!std::strcmp(type, "constbfpz"))
				{
					SensorAbstract * decorated_sens = new SensorConstBFPZ(sens);
					plhs[0] = mxCreateNumericMatrix(1, 1, mxINDEX_CLASS, mxREAL);
					*(intptr_t *)mxGetData(plhs[0]) = (intptr_t)decorated_sens;
				}
			}
			break;
		case 4:
			//Get Sensor Velocity [angular linear]
			plhs[0] = mxCreateDoubleMatrix(6,1,mxREAL);
			Eigen::Map<rl::math::Vector>(mxGetPr(plhs[0]),6,1) = sens->getFrame()->v.matrix();
			break;
		case 5:
			//Get Sensor Acceleration [angular linear]
			plhs[0] = mxCreateDoubleMatrix(6,1,mxREAL);
			Eigen::Map<rl::math::Vector>(mxGetPr(plhs[0]),6,1) = sens->getFrame()->a.matrix();
			break;
		case 6:
			//Get Sensor Measurement
			return_vector = rl::math::Vector(sens->getMeasurementSize());
			plhs[0] = mxCreateDoubleMatrix(return_vector.size(),1,mxREAL);
			sens->getMeasurement(return_vector);
			memcpy(mxGetData(plhs[0]),return_vector.data(),return_vector.size()*sizeof(double));
			break;
		case 7:
			//Get Sensor Observation Jacobian
			sens->getJacobian(return_matrix);
			return_matrix = Eigen::Matrix<double,-1,-1,Eigen::ColMajor>(sens->getMeasurementSize(),return_matrix.cols()*3);
			return_matrix.setZero();
			//return_matrix.resize((*(CSensor**)mxGetData(prhs[1]))->getMeasurementSize(),return_matrix.cols()*3);
			sens->getObservationJacobian(return_matrix);
			plhs[0] = mxCreateDoubleMatrix(return_matrix.rows(),return_matrix.cols(),mxREAL);
			memcpy(mxGetData(plhs[0]),return_matrix.data(),return_matrix.size()*sizeof(double));
			break;
		case 8: 
			//Get Sensor Transform
			plhs[0] = mxCreateDoubleMatrix(4,4,mxREAL);
			Eigen::Map<rl::math::Matrix>(mxGetPr(plhs[0]),4,4) = sens->getFrame()->t.matrix();
			break;
		case 9:
			//Get Sensor Type
			plhs[0] = mxCreateString((*(SensorAbstract**)mxGetData(prhs[1]))->getType().c_str());
			break;
		case 10:
			//Get Difference Between predicted observation and actual observation
			plhs[0] = mxCreateDoubleMatrix(sens->getMeasurementSize(),1,mxREAL);
			return_vector = rl::math::Vector(sens->getMeasurementSize());
			sens->makeDZ(Eigen::Map<rl::math::Vector>(mxGetPr(prhs[2]),mxGetM(prhs[2])*mxGetN(prhs[2])),
				Eigen::Map<rl::math::Vector>(mxGetPr(prhs[3]),mxGetM(prhs[3])*mxGetN(prhs[3])),return_vector);
			memcpy(mxGetData(plhs[0]),return_vector.data(),return_vector.size()*sizeof(double));
			break; 
		case 11:
			//Get Jacobian (end effector frame)
			sens->getJacobian(return_matrix);
			plhs[0] = mxCreateDoubleMatrix(return_matrix.rows(),return_matrix.cols(),mxREAL);
			memcpy(mxGetData(plhs[0]),return_matrix.data(),return_matrix.size()*sizeof(double));
			break;
		case 12:
			//Get Jacobian (base frame)
			sens->getBaseJacobian(return_matrix);
			plhs[0] = mxCreateDoubleMatrix(return_matrix.rows(),return_matrix.cols(),mxREAL);
			memcpy(mxGetData(plhs[0]),return_matrix.data(),return_matrix.size()*sizeof(double));
			break;
		case 13:
			//Get Tranformation Derivative base to ee
			//Returns [[dR/dq1 dt/dq1],[dR/dq2 dt/dq2],...] in 4x4xdof matrix
			{
				sens->getJacobian(return_matrix);
				const mwSize dims[3] = {4,4,return_matrix.cols()};
				plhs[0] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
				double * data = (double *)mxGetData(plhs[0]);
				const std::vector<rl::math::Matrix> derivative = sens->getBaseToEETransformDerivative();
				
				for (int i=0;i<dims[2];i++)
				{
					return_matrix = Eigen::Matrix<double,-1,-1,Eigen::ColMajor>(4,4);
					return_matrix = derivative[i];
					memcpy(&data[i*4*4],return_matrix.data(),return_matrix.size()*sizeof(double));		
				}
			}
			break;
		case 14:
			//Get Tranformation Derivative EE to base
			//Returns [[dR/dq1 dt/dq1],[dR/dq2 dt/dq2],...] in 4x4xdof matrix
			{
				sens->getJacobian(return_matrix);
				const mwSize dims[3] = {4,4,return_matrix.cols()};
				plhs[0] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
				double * data = (double *)mxGetData(plhs[0]);
				const std::vector<rl::math::Matrix> derivative = sens->getEEToBaseTransformDerivative();
				
				for (int i=0;i<dims[2];i++)
				{
					return_matrix = Eigen::Matrix<double,-1,-1,Eigen::ColMajor>(4,4);
					return_matrix = derivative[i];
					memcpy(&data[i*4*4],return_matrix.data(),return_matrix.size()*sizeof(double));		
				}
			}
			break;
		case 15:
			//Get the underlying frame
			plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
			*(intptr_t *)mxGetData(plhs[0]) = (intptr_t)sens->getFrame().get();
			break;
		case 16:
			//Get Joints Vector
			//return_vector = sens->getJoints();
			mexPrintf("Get Joint is no longer supported\n");
			plhs[0] = mxCreateDoubleMatrix(return_vector.rows(),return_vector.cols(),mxREAL);
			memcpy(mxGetData(plhs[0]),return_vector.data(),return_vector.size()*sizeof(double));
			break;
		case 17:
			//Get Derivative of W in body frame wrt Q
			sens->getdWdQ(return_matrix);
			plhs[0] = mxCreateDoubleMatrix(return_matrix.rows(), return_matrix.cols(), mxREAL);
			memcpy(mxGetData(plhs[0]), return_matrix.data(), return_matrix.size() * sizeof(double));
			break;
		case 18:
			//Get Sensor Type in Binary Format
			plhs[0] = mxCreateDoubleScalar((double)(*(SensorAbstract**)mxGetData(prhs[1]))->getTypeBin());
			break;
		case 19:
			//Get Base Jacobian Derivative
			sens->getBaseJacobianDerivative(return_matrix); 
			plhs[0] = mxCreateDoubleMatrix(return_matrix.rows(), return_matrix.cols(), mxREAL);
			memcpy(mxGetData(plhs[0]), return_matrix.data(), return_matrix.size() * sizeof(double));
			break;
		case 20:
			//Set the base of the sensor
			sens->setBase((*(rl::mdl::Frame**)mxGetData(prhs[2])));
			break;
		case 21:
			//Get the base of the sensor
			plhs[0] = mxCreateNumericMatrix(1, 1, mxINDEX_CLASS, mxREAL);
			*(intptr_t *)mxGetData(plhs[0]) = (intptr_t)sens->getBase();
			break;
	}
}

void SensorWrapper_Usage()
{

}

