#include "CCameraWrapper.h"
#include <mrpt/math.h>
#include <mrpt/utils.h>
#include "CameraSensor.h"

void CCameraWrapper(mxArray *plhs[],const mxArray *prhs[])
{
	int command = (int)mxGetScalar(prhs[0]);

	//The first parameter is always the model pointer unless we are loading model from file
	CameraSensor * cam;
	if(command != 1 && command != 0)
		cam = (*(CameraSensor**)mxGetData(prhs[1]));
	switch(command)
	{
		case 0:
			CCameraWrapper_Usage();
			break;
		case 1:
			{
				plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
				intptr_t * ptr = (intptr_t *)mxGetData(plhs[0]);
	
				char name[256];
				mxGetString(prhs[1],name,256);

				cam = new CameraSensor(std::string(name));
				*ptr = (intptr_t)cam;
			}
			break;

		case 2:
			delete cam;
			break;

		case 3:
			//Set Camera Pose
			{
				using MapMatrix44d = Eigen::Map<Eigen::Matrix<double, 4, 4>>;
				double * posePtr = mxGetPr(prhs[2]);
				MapMatrix44d pose(posePtr);
				cam->pose(pose);
			}
			break;

		case 4:
			//Get Camera Pose
			{				
				using MapMatrix44d = Eigen::Map<Eigen::Matrix<double, 4, 4>>;
				plhs[0] = mxCreateDoubleMatrix(4, 4, mxREAL);
				double* data = mxGetPr(plhs[0]);
				MapMatrix44d outPose(data);
				outPose = cam->pose();
			}
			break;

		case 5:
			{
				using MapArray33d = Eigen::Map<Eigen::Array33d>;
				double* matCalibration = mxGetPr(prhs[2]);
				cam->calibration(MapArray33d(matCalibration));
			}
			break;

		case 6:
			{
				using MapArray33d = Eigen::Map<Eigen::Array33d>;
				plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
				MapArray33d(mxGetPr(plhs[0])) = cam->calibration();
			}
			break;

		case 7:
			{
				using MapMatrixXd = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>;
				bool waitForNewImage = (int)mxGetScalar(prhs[2]) ? true : false;

				// Create the Output Image Matrix for Matlab and map it!
				size_t imgDims[3] = { cam->height, cam->width, 3 };
				plhs[0] = mxCreateNumericArray(3, imgDims, mxDOUBLE_CLASS, mxREAL);
				auto outData = mxGetPr(plhs[0]);
				
				// Map the output matrix...
				MapMatrixXd outR(outData, imgDims[0], imgDims[1]),
							outG(outData + (imgDims[0]*imgDims[1]), imgDims[0], imgDims[1]),
							outB(outData + 2 * (imgDims[0]*imgDims[1]), imgDims[0], imgDims[1]);

				cam->observation(waitForNewImage, outR, outG, outB);
			}
			break;

		default:
			mexPrintf("Unknown CCamera Command %d\r\n",command);
			CCameraWrapper_Usage();
			break;
	}
}


void CCameraWrapper_Usage(void)
{
	mexPrintf("------------------CCamera WRAPPER--------------------------------\n\r");
	mexPrintf("Usage: [return] CCameraWrapper([command id],params...)\r");
	mexPrintf("CCameraWrapper([0] *) --> Print This Usage Message\r");
	
	mexPrintf("----------------------------Function----------------------------------------------\r");
	mexPrintf("[1]  Creates CCamera \r");
	mexPrintf("		Params: string name, int width, int height\r");
	mexPrintf("		Returns: CCameraPtr * *\n\r");
	mexPrintf("[2]  Delete the CCamera c++ ptr, should be called by Matlab Destructor\r");
	mexPrintf("		Params: CCameraPtr *\r");
	mexPrintf("		Returns: none\n\r");
	mexPrintf("[3]  Sets the Camera Pose \r");
	mexPrintf("		Params: CCameraPtr *, mxArray (6x1)\r");
	mexPrintf("		Returns: none\n\r");
	mexPrintf("[4]  Returns the Camera Pose \r");
	mexPrintf("		Params: CCameraPtr *\r");
	mexPrintf("		Returns: mxArray(6x1)\n\r");
	mexPrintf(".\r");
	mexPrintf(".\r");
	mexPrintf(". please check switch statement in cpp file for other functionality \r\n");
}