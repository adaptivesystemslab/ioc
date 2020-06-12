#include "stdafx.h"
#include <rl\mdl\Body.h>
#include "BodyWrapper.h"

void BodyWrapper(mxArray *plhs[],const mxArray *prhs[])
{
	int command = (int)mxGetScalar(prhs[0]);

	//The first parameter is always the model pointer unless we are loading model from file
	rl::mdl::Body * body;
	if(command != 0)
		 body = (*(rl::mdl::Body**)mxGetData(prhs[1]));

	switch(command)
	{
		case 0:
			BodyWrapper_Usage();
			break;
		case 1:
			//set center of mass
			{
				double * com = mxGetPr(prhs[2]);
				body->setCenterOfMass(com[0],com[1],com[2]);
			}
			break;
		case 2:
			//set Mass
			body->setMass(mxGetScalar(prhs[2]));
			break;
		case 3:
			//Set Inertia
			{
				double * data = mxGetPr(prhs[2]);
				body->setInertia(data[0],data[1],data[2],data[3],data[4],data[5]);
			}
			break;
		case 4:
			//get com
			plhs[0] = mxCreateDoubleMatrix(3,1,mxREAL);
			Eigen::Map<rl::math::Vector>(mxGetPr(plhs[0]),3,1) = body->cm.matrix();
			break;
		case 5:
			//get force
			plhs[0] = mxCreateDoubleMatrix(6,1,mxREAL);
			Eigen::Map<rl::math::Vector>(mxGetPr(plhs[0]),6,1) = body->fX.matrix();
			break;
		case 6:
			//get inertia matrix
			plhs[0] = mxCreateDoubleMatrix(3,3,mxREAL);
			Eigen::Map<rl::math::Matrix>(mxGetPr(plhs[0]),3,3) = body->ic.matrix();
			break;
		case 7:
			//get mass
			plhs[0] = mxCreateDoubleScalar(body->m);
			break;
		case 8:
			//Set fX
			{
				rl::math::ForceVector fX;
				memcpy((void *)fX.matrix().data(),mxGetPr(prhs[2]),6*sizeof(double));
				body->fX = fX;
			}
			break;
	}

}

void BodyWrapper_Usage(void)
{

}