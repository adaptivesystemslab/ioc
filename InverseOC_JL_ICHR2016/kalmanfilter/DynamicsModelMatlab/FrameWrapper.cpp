#include "stdafx.h"
#include <rl\mdl\Frame.h>
#include <rl\mdl\Compound.h>
#include "FrameWrapper.h"

void FrameWrapper(mxArray *plhs[],const mxArray *prhs[])
{
	int command = (int)mxGetScalar(prhs[0]);

	//The first parameter is always the model pointer unless we are loading model from file
	rl::mdl::Frame * frame;
	if(command != 0)
		 frame = (*(rl::mdl::Frame**)mxGetData(prhs[1]));

	//Return Structures
	rl::math::Matrix return_mat;
	rl::math::Vector return_vec;

	switch(command)
	{
		case 0:
			FrameWrapper_Usage();
			break;
		case 1:
			//get name
			plhs[0] = mxCreateString(frame->getName().c_str());
			break;
		case 2:
			//get acceleration
			plhs[0] = mxCreateDoubleMatrix(6,1,mxREAL);
			Eigen::Map<rl::math::Vector>(mxGetPr(plhs[0]),6,1) = frame->a.matrix();
			break;
		case 3:
			//get c
			plhs[0] = mxCreateDoubleMatrix(6,1,mxREAL);
			Eigen::Map<rl::math::Vector>(mxGetPr(plhs[0]),6,1) = frame->c.matrix();
			break;
		case 4:
			//get force
			plhs[0] = mxCreateDoubleMatrix(6,1,mxREAL);
			Eigen::Map<rl::math::Vector>(mxGetPr(plhs[0]),6,1) = frame->f.matrix();
			break;
		case 5:
			//get rigid body inertia
			plhs[0] = mxCreateDoubleMatrix(6,6,mxREAL);
			Eigen::Map<rl::math::Matrix>(mxGetPr(plhs[0]),6,6) = frame->i.matrix();
			break;
		case 6:
			//get articulated body inertia
			plhs[0] = mxCreateDoubleMatrix(6,6,mxREAL);
			Eigen::Map<rl::math::Matrix>(mxGetPr(plhs[0]),6,6) = frame->iA.matrix();
			break;
		case 7:
			//get pA
			plhs[0] = mxCreateDoubleMatrix(6,1,mxREAL);
			Eigen::Map<rl::math::Vector>(mxGetPr(plhs[0]),6,1) = frame->pA.matrix();
			break;
		case 8:
			//get frame transform
			plhs[0] = mxCreateDoubleMatrix(4,4,mxREAL);
			Eigen::Map<rl::math::Matrix>(mxGetPr(plhs[0]),4,4) = frame->t.matrix();
			break;
		case 9:
			//get frame velocity
			plhs[0] = mxCreateDoubleMatrix(6,1,mxREAL);
			Eigen::Map<rl::math::Vector>(mxGetPr(plhs[0]),6,1) = frame->v.matrix();
			break;
		case 10:
			//get frame pleuker transform
			plhs[0] = mxCreateDoubleMatrix(6,6,mxREAL);
			//return_mat = frame->x.matrixMotion().matrix();
			Eigen::Map<rl::math::Matrix>(mxGetPr(plhs[0]),6,6) = frame->x.matrixMotion().matrix();
			break;

		//TODO: MAYBE SOMEDAY MOVE THIS INTO COMPOUND SO Frame < Element < Compound & handle
		case 11:
			//Get input transform
			plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
			if (frame->compound)
			{
				((intptr_t *)mxGetData(plhs[0]))[0] = (intptr_t)frame->compound->inTransform;
			}
			else
			{
				((intptr_t *)mxGetData(plhs[0]))[0] = 0;
			}

			break;
		case 12:
			//Get output transform
			plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
			if (frame->compound)
			{
				((intptr_t *)mxGetData(plhs[0]))[0] = (intptr_t)frame->compound->outTransform;
			}
			else
			{
				((intptr_t *)mxGetData(plhs[0]))[0] = 0;
			}

			break;
	}

}

void FrameWrapper_Usage(void)
{

}