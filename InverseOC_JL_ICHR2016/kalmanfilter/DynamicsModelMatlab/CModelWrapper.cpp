#include "stdafx.h"
#include "CModelWrapper.h"

#include <rl\mdl\XmlFactory.h>
#include <rl\mdl\UrdfFactory.h>
#ifdef _MATLAB
#include <mex.h>
#endif
#include "Model.h"


void CModelWrapper(mxArray *plhs[],const mxArray *prhs[])
{
	int command = (int)mxGetScalar(prhs[0]);

	//The first parameter is always the model pointer unless we are loading model from file
	CModel * model;
	if(command != 1 && command != 0)
		 model = (*(CModel**)mxGetData(prhs[1]));

	switch(command)
	{
		case 0:
			CModelWrapper_Usage();
			break;
		case 1:
			//Create C++ Model object and return the pointer
			{
				char str[256];
				mxGetString(prhs[1],str,256);
				std::string file_path = std::string(str);
				int period = file_path.find_last_of(".");
				std::string ext = file_path.substr(period + 1);
				CModel * model = new CModel();
				if (ext.compare("xml") == 0) {
					rl::mdl::XmlFactory factory;
					factory.load(str, model);
				}
				else {
					rl::mdl::UrdfFactory factory;
					factory.load(str, model);
				}
				plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
				intptr_t * data = (intptr_t *)mxGetData(plhs[0]);
				data[0] = (intptr_t)model;
			}
			break;
		case 2:
			//destructor
			delete model;
			break;
		case 3:
			//addSensor by frame index
			{
				//Get the transformation matrix
				rl::math::Transform t;
				t.matrix() = Eigen::Map<Eigen::Matrix<double,4,4,Eigen::ColMajor>>(mxGetPr(prhs[4]));
				//Get the sensor
				SensorAbstract *sens = (*(SensorAbstract**)mxGetData(prhs[2]));
				//Add sensor to model
				model->addSensor((int)mxGetScalar(prhs[2]),t,sens);
			}
		case 4:
			//addSensor by frame name
			{
				char frame_name[256];
				//Get the required data
				mxGetString(prhs[3],frame_name,256);
				//Get the transformation matrix
				rl::math::Transform t;
				t.matrix() = Eigen::Map<Eigen::Matrix<double,4,4,Eigen::ColMajor>>(mxGetPr(prhs[4]));
				//Get the sensor
				SensorAbstract *sens = (*(SensorAbstract**)mxGetData(prhs[2]));
				//Add sensor to model
				model->addSensor(frame_name,t,sens);
			}
			break;
		case 5:
			//remove sensor
			plhs[0] = mxCreateDoubleScalar(model->removeSensor(*(SensorAbstract**)mxGetData(prhs[2])));
			break;
		case 6:
			//Get sensor by index, returns C++ AbstractSensor *
			plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
			((intptr_t *)mxGetData(plhs[0]))[0] = (intptr_t)model->getSensor((int)mxGetScalar(prhs[2]));
			break;
		case 7:
			//Get Number of Attached Sensors
			plhs[0] = mxCreateDoubleScalar(model->getNumSensors());
			break;
		case 8:
			//Change the base of the model to a different frame based on frame index
			//TODO
			break;
		case 9 :
			//Change the base of the model based on frame name
			char frame_name[256];
			//Get the required data
			mxGetString(prhs[2],frame_name,256);
			model->setBase(frame_name);
			break;
		case 10:
			//Get the base of the model, returns rl::mdl::Frame *
			//TODO
			break;
		case 11:
			//Get Number of transforms
			plhs[0] = mxCreateDoubleScalar((double)model->getTransforms());
			break;
		case 12:
			//Get transfrom between two frames by index
			plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
			((intptr_t *)mxGetData(plhs[0]))[0] = (intptr_t)((rl::mdl::Transform *)model->getTransform((size_t)mxGetScalar(prhs[2])));
			break;
		case 13:
			//Get frame by name, returns rl::mdl::Frame *
			{
				char frame_name[256];
				//Get the required data
				mxGetString(prhs[2],frame_name,256);
				plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
				((intptr_t *)mxGetData(plhs[0]))[0] = (intptr_t)model->getFrame(frame_name);
			}
			break;
		case 14:
			//Calculate Sensor Jacobians
			model->setSensorMeasurements();
			//model->getSensorJacobians(*(SensorAbstract**)mxGetData(prhs[2]));
			break;
		case 15:
			//Get transform by name, returns rl::mdl::Transform *
			{
				char name[256];
				mxGetString(prhs[2], name, 256);
				rl::math::Transform t = (*(CModel**)mxGetData(prhs[1]))->getTransform(name);
				plhs[0] = mxCreateDoubleMatrix(4, 4, mxREAL);
				memcpy(mxGetData(plhs[0]), t.data(), 4 * 4 * sizeof(double));
				break;
			}
		case 16:
			//Get Magnetic Field
			{
				plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
				rl::math::Vector m(3);
				model->getMagField(m);
				memcpy(mxGetPr(plhs[0]), m.data(), 3 * sizeof(double));
			}
			break;
		case 19:
			//Set Magnetic Field
			model->setMagField(Eigen::Map<rl::math::Vector>(mxGetPr(prhs[2]), 3, 1));
			break;
		default:
			mexPrintf("Unknown Command %d\r\n",command);
			CModelWrapper_Usage();
			break;
	}
}

void CModelWrapper_Usage(void)
{


}
