#include "stdafx.h"
#include <rl\mdl\Transform.h>
#include <rl\mdl\Joint.h>
#include <rl\mdl\Revolute.h>
#include <rl\mdl\Prismatic.h>
#include <rl\mdl\Spherical.h>
#include "JointWrapper.h"

void JointWrapper(mxArray *plhs[],const mxArray *prhs[])
{
	int command = (int)mxGetScalar(prhs[0]);

	//The first parameter is always the model pointer unless we are loading model from file
	rl::mdl::Joint * t;
	if(command != 0)
		 t = (*(rl::mdl::Joint**)mxGetData(prhs[1]));

	switch(command)
	{
		case 0:
			JointWrapper_Usage();
			break;
		case 1:
			//Get Type of Joint
			if(dynamic_cast<rl::mdl::Revolute *>(t))
			{
				plhs[0] = mxCreateString("revolute");
			}
			else if(dynamic_cast<rl::mdl::Prismatic *>(t))
			{
				plhs[0] = mxCreateString("prismatic");
			}
			else if(dynamic_cast<rl::mdl::Spherical *>(t))
			{
				plhs[0] = mxCreateString("spherical");
			}
			else
			{
				plhs[0] = mxCreateString("unknown");
			}
			break;
		case 2:
			//Get Joint Position
			plhs[0] = mxCreateDoubleMatrix(t->getDofPosition(),1,mxREAL);
			Eigen::Map<rl::math::Vector>(mxGetPr(plhs[0]),t->getDofPosition(),1) = t->getPosition();
			break;
		case 3:
			//Set Position
			t->setPosition(Eigen::Map<rl::math::Vector>(mxGetPr(prhs[2]),t->getDofPosition()));
		case 4:
			//Get Joint Velocity
			plhs[0] = mxCreateDoubleMatrix(t->getDof(),1,mxREAL);
			Eigen::Map<rl::math::Vector>(mxGetPr(plhs[0]),t->getDof(),1) = t->getVelocity();
			break;
		case 5:
			//Set Joint Velocity
			t->setVelocity(Eigen::Map<rl::math::Vector>(mxGetPr(prhs[2]),t->getDof()));
			break;
		case 6:
			//Get Joint Velocity
			plhs[0] = mxCreateDoubleMatrix(t->getDof(),1,mxREAL);
			Eigen::Map<rl::math::Vector>(mxGetPr(plhs[0]),t->getDof(),1) = t->getAcceleration();
			break;
		case 7:
			//Set Joint Velocity
			t->setAcceleration(Eigen::Map<rl::math::Vector>(mxGetPr(prhs[2]),t->getDof()));
			break;
		default:
			JointWrapper_Usage();
			break;
	}
}

void JointWrapper_Usage(void)
{
	mexPrintf("------------------JOINT WRAPPER--------------------------------\n\r");
	mexPrintf("Usage: [return] JointWrapper([command id],params...)\r");
	mexPrintf("VisualizerWrapper([0] *) --> Print This Usage Message\r");
	
	mexPrintf("----------------------------Function----------------------------------------------\r");
	mexPrintf("[1]  Get Joint Type \r");
	mexPrintf("		Params: none \r");
	mexPrintf("		Returns: String \n\r");
}