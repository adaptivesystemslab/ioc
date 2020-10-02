#include "stdafx.h"
#include <rl\mdl\Kinematic.h>
#include <rl\mdl\Dynamic.h>
#include <rl\mdl\XmlFactory.h>

#include "DynModelWrapper.h"

void DynamicModelWrapper(mxArray *plhs[],const mxArray *prhs[])
{
	int command = (int)mxGetScalar(prhs[0]);

	//The first parameter is always the model pointer unless we are loading model from file
	rl::mdl::Dynamic * model;
	if(command != 1 && command != 0)
		 model = (*(rl::mdl::Dynamic**)mxGetData(prhs[1]));

	//Return Structures
	rl::math::Matrix return_mat;
	rl::math::Vector return_vec;
	switch(command)
	{
		case 0:
			DynamicModelWrapper_Usage();
			break;

		case 1:
			//Load Dynamic Model From File have to break scope but its ok since this will not be called often
			{
				rl::mdl::XmlFactory factory;
				char str[256];
				mxGetString(prhs[1],str,256);

				rl::mdl::Kinematic * model = new rl::mdl::Dynamic();
				factory.load(str,model);
				plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
				intptr_t * data = (intptr_t *)mxGetData(plhs[0]);
				data[0] = (intptr_t)model;
			}
			break;
		case 2:
			//Destructor
			delete model;
			break;
		case 3:
			//calculateCentrifugalCoriolis
			model->calculateCentrifugalCoriolis();
			break;
		case 4:
			//calculateGravity
			model->calculateGravity();
			break;
		case 5:
			//calculateMassMatrix
			model->calculateMassMatrix();
			break;
		case 6:
			//calculateMassMatrixInverse
			model->calculateMassMatrixInverse();
			break;
		case 7:
			//calculateOperationalMassMatrixInverse
			model->calculateOperationalMassMatrixInverse();
			break;
		case 8:
			model->eulerCauchy(mxGetScalar(prhs[2]));
			break;
		case 9:
			//forwardDynamics
			model->forwardDynamics();
			break;
		case 10:
			//getCentrifugalCoriolis
			plhs[0] = mxCreateDoubleMatrix(model->getDof(),1,mxREAL);
			Eigen::Map<rl::math::Vector>(mxGetPr(plhs[0]),model->getDof(),1) = model->getCentrifugalCoriolis();
			break;
		case 11:
			//getGravity
			plhs[0] = mxCreateDoubleMatrix(model->getDof(),1,mxREAL);
			Eigen::Map<rl::math::Vector>(mxGetPr(plhs[0]),model->getDof(),1) = model->getGravity();
			break;
		case 12:
			//getMassMatrixInverse
			plhs[0] = mxCreateDoubleMatrix(model->getDof(),model->getDof(),mxREAL);
			Eigen::Map<rl::math::Matrix>(mxGetPr(plhs[0]),model->getDof(),model->getDof()) = model->getMassMatrixInverse();
			break;
		case 13:
			//getMassMatrix
			plhs[0] = mxCreateDoubleMatrix(model->getDof(),model->getDof(),mxREAL);
			Eigen::Map<rl::math::Matrix>(mxGetPr(plhs[0]),model->getDof(),model->getDof()) = model->getMassMatrix();
			break;
		case 14:
			//getOperationalMassMatrixInverse
			plhs[0] = mxCreateDoubleMatrix(model->getOperationalDof()*6,model->getOperationalDof()*6,mxREAL);
			Eigen::Map<rl::math::Matrix>(mxGetPr(plhs[0]),model->getOperationalDof()*6,model->getOperationalDof()*6) = model->getOperationalMassMatrixInverse();
			break;
		case 15:
			//getWorldGravity
			{
				plhs[0] = mxCreateDoubleMatrix(3,1,mxREAL);
				rl::math::Vector g(3);
				model->getWorldGravity(g);
				memcpy(mxGetPr(plhs[0]),g.data(),3*sizeof(double));
			}
			break;
		case 16:
			//inverseDynamics
			model->inverseDynamics();
			break;
		case 17:
			//inverseForce
			model->inverseForce();
			break;
		case 18:
			//rungeKuttaNystrom
			model->rungeKuttaNystrom(mxGetScalar(prhs[2]));
			break;
		case 19:
			//setWorldGravity
			model->setWorldGravity(Eigen::Map<rl::math::Vector>(mxGetPr(prhs[2]),3,1));
			break;
	}
}


void DynamicModelWrapper_Usage(void)
{
	mexPrintf("------------------Robotics Library Dynamic Model Matlab WRAPPER--------------------------------\n\r");
	mexPrintf("Usage: [return] DynamicModelWrapper([command id],params...)\r");
	mexPrintf("DynamicModelWrapper([0] *) --> Print This Usage Message\r");
	
	mexPrintf("----------------------------Function----------------------------------------------\r");
	mexPrintf("[1]  Loads Model From File\r");
	mexPrintf("		Params: filename\r");
	mexPrintf("		Returns: rl::mdl::Dynamic *\n\r");
	mexPrintf("[2]  Delete the model c++ ptr, should be called by Matlab Destructor\r");
	mexPrintf("		Params: rl::mdl::Kinematic *\r");
	mexPrintf("		Returns: none\n\r");
}


