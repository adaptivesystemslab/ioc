// DynamicsModelMatlab.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include <matrix.h>
#include <mex.h>

#include "DynamicsModelMatlab.h"

//Include Robotics Library Wrappers
#include "ModelWrapper.h"
#include "KinModelWrapper.h"
#include "DynModelWrapper.h"
#include "CModelWrapper.h"
#include "VisualizerWrapper.h"
#include "SensorWrapper.h"
#include "FrameWrapper.h"
#include "BodyWrapper.h"
#include "TransformWrapper.h"
#include "JointWrapper.h"
#include "CCameraWrapper.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, 
  const mxArray *prhs[])
{
	int class_id = (int)mxGetScalar(prhs[0]); 



	switch(class_id)
	{
		case 0:
			usage();
			break;
		case 1:
			//Jump into the rl::mdl::Model Wrapper
			ModelWrapper(plhs,&prhs[1]);
			break;
		case 2:
			//Jump into the Kinematic Model Wrapper
			KinematicModelWrapper(plhs,&prhs[1]);
			break;
		case 3:
			//Jump into the Dynamic Model Wrapper
			DynamicModelWrapper(plhs,&prhs[1]);
			break;
		case 4:
			//Jump into the CModel Model Wrapper
			CModelWrapper(plhs,&prhs[1]);
			break;
		case 5:
			//Jump into the Visualizer Wrapper
			VisualizerWrapper(plhs,&prhs[1]);
			break;
		case 6:
			//Sensor and their decorator Wrapper 
			SensorWrapper(plhs,&prhs[1]);
			break;
		case 7:
			//Frame Wrapper
			FrameWrapper(plhs,&prhs[1]);
			break;
		case 8:
			//Body Wrapper
			BodyWrapper(plhs,&prhs[1]);
			break;
		case 9:
			//Body Wrapper
			TransformWrapper(plhs,&prhs[1]);
			break;
		case 10:
			//Joint Wrapper
			JointWrapper(plhs,&prhs[1]);
			break;

		case 11:
			CCameraWrapper(plhs, &prhs[1]);
			break;

		default:
			mexPrintf("Unknown Class ID: %d\r\n",class_id);
			usage();
			break;
	}
}

