#include <stdio.h>
#include <mex.h>
#include <matrix.h>
#ifndef DYNAMICS_MODEL_MATLAB_H
#define DYNAMICS_MODEL_MATLAB_H

void usage(void)
{
	mexPrintf("------------------Robotics Library Matlab WRAPPER--------------------------------\n\r");
	mexPrintf("Usage: [return] DynamicsModelMatlab([class id][command id],params...)\r");
	mexPrintf("DynamicsModelMatlab([0] *) --> Print This Usage Message\r");
	mexPrintf("DynamicsModelMatlab([class id] [0] *) --> Print Usage For Specific Class Wrapper\n\r");
	mexPrintf("----------------------------CLASSES----------------------------------------------\r");
	mexPrintf("[1] Kinematic Model\r");
	mexPrintf("[2] Dynamic Model\r");
	mexPrintf("[3] Joint\r");
	mexPrintf("[4] Frame\r");
	mexPrintf("[5] Body\r");
	mexPrintf("[6] Transform\r");
	mexPrintf("[11] Camera Object\r");
	mexPrintf("---------------------------------------------------------------------------------\n\r");
}
#endif