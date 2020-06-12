#include "stdafx.h"
#include <rl\mdl\Model.h>
#include <rl\mdl\XmlFactory.h>


#include "ModelWrapper.h"


void ModelWrapper(mxArray *plhs[],const mxArray *prhs[])
{
	int command = (int)mxGetScalar(prhs[0]);

	//The first parameter is always the model pointer unless we are loading model from file
	rl::mdl::Model * model;
	if(command != 1)
		 model = (*(rl::mdl::Model**)mxGetData(prhs[1]));

	//Return Structures
	rl::math::Matrix return_mat;
	rl::math::Vector return_vec;
	rl::math::MotionVector return_mv;
	rl::math::ForceVector return_fv;
	rl::math::Transform return_t;
	
	switch(command)
	{
		case 0:
			ModelWrapper_Usage();
			break;
		case 5:
			//get acceleration
			plhs[0] = mxCreateDoubleMatrix(model->getDof(),1,mxREAL);
			return_vec = model->getAcceleration();
			Eigen::Map<rl::math::Vector>(mxGetPr(plhs[0]),return_vec.rows(),return_vec.cols()) = return_vec;
			break;
		case 6:
			//Number of bodies
			plhs[0] = mxCreateDoubleScalar((double)model->getBodies());
			break;
		case 7:
			//Get specific Body, returns the C++ pointer
			plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
			((intptr_t *)mxGetData(plhs[0]))[0] = (intptr_t)model->getBody((int)mxGetScalar(prhs[2]));
			break;

		case 8:
			//getDof
			plhs[0] = mxCreateDoubleScalar((double)model->getDof());
			break;
		case 9:
			//getDofPosition
			plhs[0] = mxCreateDoubleScalar((double)model->getDofPosition());
			break;
		case 10:
			//getFrame
			break;
		case 11:
			//getJoint
			plhs[0] = mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
			((intptr_t *)mxGetData(plhs[0]))[0] = (intptr_t)model->getJoint((int)mxGetScalar(prhs[2]));
			break;
		case 12:
			//getJoints
			plhs[0] = mxCreateDoubleScalar((double)model->getJoints());
			break;
		case 13:
			//getOperationalAcceleration
			plhs[0] = mxCreateDoubleMatrix(model->getDof(),1,mxREAL);
			return_mv = model->getOperationalAcceleration((int)mxGetScalar(prhs[2]));
			memcpy(mxGetData(plhs[0]),return_mv.matrix().data(),6*sizeof(double));
			break;
		case 14:
			//getOperationalDof 
			plhs[0] = mxCreateDoubleScalar((double)model->getOperationalDof());
			break;
		case 15:
			//getOperationalForce
			plhs[0] = mxCreateDoubleMatrix(6,1,mxREAL);
			return_fv = model->getOperationalForce((int)mxGetScalar(prhs[2]));
			memcpy(mxGetData(plhs[0]),return_mv.matrix().data(),6*sizeof(double));
			break;
		case 16:
			//getOperationalPosition
			plhs[0] = mxCreateDoubleMatrix(4,4,mxREAL);
			memcpy(mxGetData(plhs[0]),model->getOperationalPosition((int)mxGetScalar(prhs[2])).data(),4*4*sizeof(double));
			break;
		case 17:
			//getOperationalVelocity
			plhs[0] = mxCreateDoubleMatrix(6,1,mxREAL);
			memcpy(mxGetData(plhs[0]),model->getOperationalVelocity((int)mxGetScalar(prhs[2])).matrix().data(),6*sizeof(double));
			break;
		case 18:
			//getManufacturer
			plhs[0] = mxCreateString(model->getManufacturer().c_str());
			break;
		case 19:
			//getMaximum
			plhs[0] = mxCreateDoubleMatrix(model->getDof(),1,mxREAL);
			return_vec = rl::math::Vector(Eigen::Map<rl::math::Vector>(mxGetPr(plhs[0]),model->getDof()));
			return_vec = model->getMaximum();
			break;
		case 20:
			//getMinimum
			plhs[0] = mxCreateDoubleMatrix(model->getDof(),1,mxREAL);
			return_vec = rl::math::Vector(Eigen::Map<rl::math::Vector>(mxGetPr(plhs[0]),model->getDof()));
			return_vec = model->getMinimum();
			break;
		case 21:
			//getName
			plhs[0] = mxCreateString(model->getName().c_str());
			break;
		case 22:
			//getPosition
			plhs[0] = mxCreateDoubleMatrix(model->getDofPosition(),1,mxREAL);
			return_vec = model->getPosition();
			Eigen::Map<rl::math::Vector>(mxGetPr(plhs[0]),return_vec.rows(),return_vec.cols()) = return_vec;
			break;
		case 23:
			//getTorque
			plhs[0] = mxCreateDoubleMatrix(model->getDof(),1,mxREAL);
			return_vec = model->getTorque();
			Eigen::Map<rl::math::Vector>(mxGetPr(plhs[0]),return_vec.rows(),return_vec.cols()) = return_vec;
			break;
		case 24:
			//getSpeed
			plhs[0] = mxCreateDoubleMatrix(model->getDof(),1,mxREAL);
			return_vec = model->getSpeed();
			Eigen::Map<rl::math::Vector>(mxGetPr(plhs[0]),return_vec.rows(),return_vec.cols()) = return_vec;
			break;
		case 25:
			//getVelocity
			plhs[0] = mxCreateDoubleMatrix(model->getDof(),1,mxREAL);
			return_vec = model->getVelocity();
			Eigen::Map<rl::math::Vector>(mxGetPr(plhs[0]),return_vec.rows(),return_vec.cols()) = return_vec;
			break;
		case 26:
			//isColliding
			plhs[0] = mxCreateDoubleScalar(model->isColliding((int)mxGetScalar(prhs[2])));
			break;
		case 27:
			//setAcceleration
			model->setAcceleration(Eigen::Map<rl::math::Vector>(mxGetPr(prhs[2]),model->getDof()));
			break;
		case 28:
			//setManufacturer
			//TODO
			break;
		case 29:
			//setName
			//TODO
			break;
		case 30:
			//setOperationalVelocity
			//TODO
			break;
		case 31:
			//setPosition
			model->setPosition(Eigen::Map<rl::math::Vector>(mxGetPr(prhs[2]),model->getDofPosition()));
			break;
		case 32:
			//setTorque
			model->setTorque(Eigen::Map<rl::math::Vector>(mxGetPr(prhs[2]),model->getDof()));
			break;
		case 33:
			//setVelocity
			model->setVelocity(Eigen::Map<rl::math::Vector>(mxGetPr(prhs[2]),model->getDof()));
			break;
		default:
			mexPrintf("Model Wrapper: Unknown Command %d\r\n",command);
			ModelWrapper_Usage();
			break;
	}

}


void ModelWrapper_Usage()
{
	mexPrintf("------------------Robotics Library rl::mdl::Model Matlab WRAPPER--------------------------------\n\r");
	mexPrintf("Usage: [return] ModelWrapper([command id],params...)\r");
	mexPrintf("ModelWrapper([0] *) --> Print This Usage Message\r");
	
	mexPrintf("----------------------------Function----------------------------------------------\r");
	mexPrintf("[5]  Get Acceleration\r");
	mexPrintf("		Params: rl::mdl::Model *\r");
	mexPrintf("		Returns: Qdd\n\r");
	mexPrintf("[6]  Get Number of Bodies\r");
	mexPrintf("		Params: rl::mdl::Model *\r");
	mexPrintf("		Returns: size_t\n\r");
	mexPrintf("[7]  Get Body by number\r");
	mexPrintf("		Params: rl::mdl::Model *, size_t\r");
	mexPrintf("		Returns: rl::mdl::Body*\n\r");
	mexPrintf("[8]  Get Degrees of Freedom\r");
	mexPrintf("		Params: rl::mdl::Model *\r");
	mexPrintf("		Returns: size_t\n\r");
	mexPrintf("[9]  Get Degrees of Freedom Position\r");
	mexPrintf("		Params: rl::mdl::Model *\r");
	mexPrintf("		Returns: size_t\n\r");
	mexPrintf("[10]  Get Frame\r");
	mexPrintf("		Params: rl::mdl::Model *\r");
	mexPrintf("		Returns: Qdd\n\r");
	mexPrintf("[11]  Get Joint by number\r");
	mexPrintf("		Params: rl::mdl::Model *, size_t\r");
	mexPrintf("		Returns: rl::mdl::Joint *\n\r");
	mexPrintf("[12]  Get Acceleration\r");
	mexPrintf("		Params: rl::mdl::Model *\r");
	mexPrintf("		Returns: Qdd\n\r");
	mexPrintf("[13]  Get Acceleration\r");
	mexPrintf("		Params: rl::mdl::Model *\r");
	mexPrintf("		Returns: Qdd\n\r");
	mexPrintf("[14]  Get Acceleration\r");
	mexPrintf("		Params: rl::mdl::Model *\r");
	mexPrintf("		Returns: Qdd\n\r");
	mexPrintf("[15]  Get Acceleration\r");
	mexPrintf("		Params: rl::mdl::Model *\r");
	mexPrintf("		Returns: Qdd\n\r");
}
