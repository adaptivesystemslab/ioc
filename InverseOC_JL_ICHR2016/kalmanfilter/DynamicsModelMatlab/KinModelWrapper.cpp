#include "stdafx.h"
#include <rl\mdl\Kinematic.h>
#include <rl\mdl\XmlFactory.h>
#include <rl\mdl\UrdfFactory.h>

#include "KinModelWrapper.h"

void KinematicModelWrapper(mxArray *plhs[],const mxArray *prhs[])
{
	int command = (int)mxGetScalar(prhs[0]);

	//The first parameter is always the model pointer unless we are loading model from file
	rl::mdl::Kinematic * model;
	if(command != 1)
		 model = (*(rl::mdl::Kinematic**)mxGetData(prhs[1]));

	//Return Structures
	rl::math::Matrix return_mat;
	rl::math::Vector return_vec;



	switch(command)
	{
		case 0:
			KinematicModelWrapper_Usage();
			break;

		case 1:
			//Load Kinematic Model From File have to break scope but its ok since this will not be called often
			{
				rl::mdl::XmlFactory factory;
				char str[256];
				mxGetString(prhs[1],str,256);

				rl::mdl::Kinematic * model = new rl::mdl::Kinematic();
				factory.load(str,model);
				plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
				intptr_t * data = (intptr_t *)mxGetData(plhs[0]);
				data[0] = (intptr_t)model;
			}
			break;
		case 2:
			//Cleanup 
			delete model;
			break;
		case 3:
			//calculateInversePosition 
			{
				rl::math::Transform t;
				t.matrix() = Eigen::Map<Eigen::Matrix<double,4,4,Eigen::ColMajor>>(mxGetPr(prhs[2]));
				plhs[0] = mxCreateDoubleScalar(model->calculateInversePosition(t,(int)mxGetScalar(prhs[3])));
			}
			break;
		case 4:
			//calculateJacobianDerivative
			model->calculateJacobian();
			break;
		case 5:
			//calculateJacobianDerivative
			model->calculateJacobianDerivative();
			break;
		case 6:
			//calculateJacobianInverse
			model->calculateJacobianInverse();
			break;
		case 7: 
			//calculateManipulabilityMeasure
			plhs[0] = mxCreateDoubleScalar(model->calculateManipulabilityMeasure());
			break;
		case 8:
			//forwardAcceleration
			model->forwardAcceleration();
			break;
		case 9:
			//forwardPosition
			model->forwardPosition();
			break;
		case 10:
			//forwardVelocity
			model->forwardVelocity();
			break;
		case 11:
			//getJacobian
			return_mat = model->getJacobian();
			plhs[0] = mxCreateDoubleMatrix(return_mat.rows(),return_mat.cols(),mxREAL);
			memcpy(mxGetData(plhs[0]),return_mat.data(),return_mat.size()*sizeof(double));
			break;
		case 12:
			//getJacobianDerivative
			return_vec = model->getJacobianDerivative();
			plhs[0] = mxCreateDoubleMatrix(return_vec.size(),1,mxREAL);
			memcpy(mxGetData(plhs[0]),return_vec.data(),return_vec.size()*sizeof(double));
			break;
		case 13:
			//getJacobianInverse
			return_mat = model->getJacobianInverse();
			plhs[0] = mxCreateDoubleMatrix(return_mat.rows(),return_mat.cols(),mxREAL);
			memcpy(mxGetData(plhs[0]),return_mat.data(),return_mat.size()*sizeof(double));
			break;
		default:
			mexPrintf("KinModelWrapper: UNKNOWN COMMAND \r\n");
			KinematicModelWrapper_Usage();
			break;
	}
}



void KinematicModelWrapper_Usage()
{
	mexPrintf("------------------Robotics Library Kinematics Model Matlab WRAPPER--------------------------------\n\r");
	mexPrintf("Usage: [return] KinematicModelWrapper([command id],params...)\r");
	mexPrintf("KinematicModelWrapper([0] *) --> Print This Usage Message\r");
	
	mexPrintf("----------------------------Function----------------------------------------------\r");
	mexPrintf("[1]  Loads Model From File\r");
	mexPrintf("		Params: filename\r");
	mexPrintf("		Returns: rl::mdl::Kinematic *\n\r");
	mexPrintf("[2]  Delete the model c++ ptr, should be called by Matlab Destructor\r");
	mexPrintf("		Params: rl::mdl::Kinematic *\r");
	mexPrintf("		Returns: none\n\r");
	mexPrintf("[3]  Calculate inverse position based on end effector position using a recursive method\r");
	mexPrintf("		Params: rl::mdl::Kinematic *, 4x4 ee T mat, ee number \r");
	mexPrintf("		Returns: none *\n\r");
	mexPrintf("[4]  Calculate Jacobian \r");
	mexPrintf("		Params: rl::mdl::Kinematic *\r");
	mexPrintf("		Returns: none *\n\r");
	mexPrintf("[5]  Calculate Jacobian Derivative JdQd\r");
	mexPrintf("		Params: rl::mdl::Kinematic *\r");
	mexPrintf("		Returns: none *\n\r");
	mexPrintf("[6]  Calculate Jacobian Inverse\r");
	mexPrintf("		Params: rl::mdl::Kinematic *\r");
	mexPrintf("		Returns: none *\n\r");
	mexPrintf("[7]  Calculate Manipulability Measure\r");
	mexPrintf("		Params: rl::mdl::Kinematic *\r");
	mexPrintf("		Returns: none *\n\r");
	mexPrintf("[8]  Forward Acceleration\r");
	mexPrintf("		Params: rl::mdl::Kinematic *\r");
	mexPrintf("		Returns: none *\n\r");
	mexPrintf("[9]  Forward Position\r");
	mexPrintf("		Params: rl::mdl::Kinematic *\r");
	mexPrintf("		Returns: none *\n\r");
	mexPrintf("[10] Forward Velocity\r");
	mexPrintf("		Params: rl::mdl::Kinematic *\r");
	mexPrintf("		Returns: none *\n\r");
	mexPrintf("[11] Get Jacobian\r");
	mexPrintf("		Params: rl::mdl::Kinematic *\r");
	mexPrintf("		Returns: Jacobian Matrix\n\r");
	mexPrintf("[12] Get Jacobian Derivative\r");
	mexPrintf("		Params: rl::mdl::Kinematic *\r");
	mexPrintf("		Returns: Jd*Qd Vector\n\r");
	mexPrintf("[13] Get Jacobian Inverse\r");
	mexPrintf("		Params: rl::mdl::Kinematic *\r");
	mexPrintf("		Returns: inv(J)\n\r");
	mexPrintf("---------------------------------------------------------------------------------\n\r");
}