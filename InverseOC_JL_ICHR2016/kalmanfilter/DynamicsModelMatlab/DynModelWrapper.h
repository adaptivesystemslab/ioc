#ifndef DYNAMIC_MODEL_WRAPPER_H
#define DYNAMIC_MODEL_WRAPPER_H
#include <stdio.h>
#include <mex.h>

void DynamicModelWrapper(mxArray *plhs[],const mxArray *prhs[]);

void DynamicModelWrapper_Usage();

#endif