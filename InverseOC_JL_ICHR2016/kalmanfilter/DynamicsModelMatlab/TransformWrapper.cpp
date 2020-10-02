#include "stdafx.h"
#include <rl\mdl\Transform.h>
#include "TransformWrapper.h"

void TransformWrapper(mxArray *plhs[],const mxArray *prhs[])
{
	int command = (int)mxGetScalar(prhs[0]);

	//The first parameter is always the model pointer unless we are loading model from file
	rl::mdl::Transform * t;
	if(command != 0)
		 t = (*(rl::mdl::Transform**)mxGetData(prhs[1]));

	switch(command)
	{
		case 0:
			TransformWrapper_Usage();
			break;
		case 1:
			//Get input frame
			plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
			((intptr_t *)mxGetData(plhs[0]))[0] = (intptr_t)t->in;
			break;
		case 2:
			//Get output frame
			plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
			((intptr_t *)mxGetData(plhs[0]))[0] = (intptr_t)t->out;
			break;
		case 3:
			//Get transformation Matrix
			plhs[0] = mxCreateDoubleMatrix(4,4,mxREAL);
			Eigen::Map<rl::math::Matrix>(mxGetPr(plhs[0]),4,4) = t->t.matrix();
			break;
		case 4:
			//Get Pleuker Transform
			plhs[0] = mxCreateDoubleMatrix(4,4,mxREAL);
			Eigen::Map<rl::math::Matrix>(mxGetPr(plhs[0]),4,4) = t->x.transform().matrix();
			break;
		case 5:
			//Get Name
			plhs[0] = mxCreateString(t->getName().c_str());
			break;
		case 6:
			//Set Transform
			t->t.matrix() = Eigen::Map<Eigen::Matrix<double,4,4,Eigen::ColMajor>>(mxGetPr(prhs[2]));
			t->x.rotation() = t->t.linear().transpose();
			t->x.translation() = t->t.translation();
			break;
	}
}

void TransformWrapper_Usage(void)
{

}