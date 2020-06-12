#ifndef KINEMATIC_MODEL_WRAPPER_H
#define KINEMATIC_MODEL_WRAPPER_H
#include <stdio.h>
#include <mex.h>

void KinematicModelWrapper(mxArray *plhs[],const mxArray *prhs[]);

void KinematicModelWrapper_Usage();

#endif