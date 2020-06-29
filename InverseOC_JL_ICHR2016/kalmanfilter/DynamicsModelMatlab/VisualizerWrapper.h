#ifndef VISUALIZER_WRAPPER_H
#define VISUALIZER_WRAPPER_H
#include <stdio.h>
#include <mex.h>

void VisualizerWrapper(mxArray *plhs[],const mxArray *prhs[]);

void VisualizerWrapper_Usage();

#endif