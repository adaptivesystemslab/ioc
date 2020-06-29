#include "VisualizerWrapper.h"
#include "Visualizer.h"
#include <mrpt/opengl.h>
#include <Windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include "Model.h"
#include <mrpt/synch/CCriticalSection.h>
#include <Eigen\src\Core\Matrix.h>
#include <mrpt/opengl.h>
#include "CCameraWrapper.h"

void VisualizerWrapper(mxArray *plhs[],const mxArray *prhs[])
{
	int command = (int)mxGetScalar(prhs[0]);

	//mexPrintf("Size of int: %d\n Size of intptr_t: %d\n Size of void*: %d", sizeof(int), sizeof(intptr_t), sizeof(void*));

	//The first parameter is always the model pointer unless we are loading model from file
	CVisualizer * vis;
	if(command != 1 && command != 0)
		 vis = (*(CVisualizer**)mxGetData(prhs[1]));

	switch(command)
	{
		case 0:
			VisualizerWrapper_Usage();
			break;
		case 1:
			//Create C++ Model object and return the pointer
			{
				plhs[0] = mxCreateNumericMatrix(1,1,mxINDEX_CLASS,mxREAL);
				intptr_t * ptr = (intptr_t *)mxGetData(plhs[0]);
	
				char name[256];
				mxGetString(prhs[1],name,256);

				CVisualizer *vis = new CVisualizer(name,mxGetScalar(prhs[2]),mxGetScalar(prhs[3]));
				vis->init();
				*ptr = (intptr_t)vis;
			}
			break;
		case 2:
			//destructor
			{
				mrpt::opengl::COpenGLScenePtr &theScene = vis->get3DSceneAndLock();
				theScene->clear();
				vis->unlockAccess3DScene();
				// Update window, if required
				vis->forceRepaint();
				delete vis;
			}
			break;
		case 3:
			//addModel
			vis->addModel(*(CModel**)mxGetData(prhs[2]));
			break;
		case 4:
			//update
			vis->update();
			break;
		case 5:
			//addMarker
			{				
				using MapVector3d = Eigen::Map<Eigen::Vector3d>;
				using MapVector4d = Eigen::Map<Eigen::Vector4d>;

				char name[256];
				mxGetString(prhs[2],name,256);

				double* position = mxGetPr(prhs[3]);
				double* colour = mxGetPr(prhs[4]);

				(*(CVisualizer **)mxGetData(prhs[1]))->addMarker(std::string(name), 
																 MapVector3d(position),
																 MapVector4d(colour),
																 (bool)mxGetScalar(prhs[5]));
				break;
			}

		case 6:
		//Visualize Sensor Covariance
			(*(CVisualizer **)mxGetData(prhs[1]))->setSensorCovar((*(CModel**)mxGetData(prhs[2])),(*(SensorAbstract**)mxGetData(prhs[3])),
				mrpt::math::CMatrixDouble::Map(mxGetPr(prhs[4]),3,3),(bool)mxGetScalar(prhs[5]));
		break;
		case 7:
		//Change Sensor Color
		(*(CVisualizer **)mxGetData(prhs[1]))->setSensorColor((*(CModel**)mxGetData(prhs[2])),(*(SensorAbstract**)mxGetData(prhs[3])),
				mrpt::math::CMatrixDouble41::Map(mxGetPr(prhs[4])));
		break;
		case 8:
			//load scene into visualizer
			{
				char name[256];
				mxGetString(prhs[2],name,256);
				(*(CVisualizer **)mxGetData(prhs[1]))->addScene(std::string(name));
			}
			break;
		case 9:
			//load AssimpModel into visualizer
			{
				char name[256];
				mxGetString(prhs[2],name,256);
				double * posePtr = mxGetPr(prhs[3]);
				mrpt::math::TPose3D pose(posePtr[0],posePtr[1],posePtr[2],posePtr[3],posePtr[4],posePtr[5]);
				(*(CVisualizer **)mxGetData(prhs[1]))->add3DModel(std::string(name),pose);
			}
			break;
		case 10:
			//Add camera to visualizer
			{
				auto camera = *(CameraSensor**)mxGetData(prhs[2]);
				vis->addCamera(camera);
			}
			break;

		case 11:		
			// Set Marker Covariance
			{
				using MapVector4d = Eigen::Map<Eigen::Vector4d>;
				using MapMatrix3d = Eigen::Map<Eigen::Matrix3d>;

				char name[256];
				mxGetString(prhs[2], name, 256);

				double* covar = mxGetPr(prhs[3]);
				double* colour = mxGetPr(prhs[4]);

				//Visualize Sensor Covariance
				vis->setMarkerCovar(std::string(name), MapMatrix3d(covar), MapVector4d(colour), (bool)mxGetScalar(prhs[5]));
				break;
			}

		case 12:
			// Set Background Colour for Visualizer
			{
				using MapVector3d = Eigen::Map<Eigen::Vector3d>;

				double* colour = mxGetPr(prhs[2]);

				auto viewport = (*(CVisualizer **)mxGetData(prhs[1]))->getDefaultViewport();
				viewport->setCustomBackgroundColor(mrpt::utils::TColorf(colour[0], colour[1], colour[2]));
				break;
			}

		case 13:
			// Capture Image in the Visualizer.
			{				
				using MapMatrixXd = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>;
				 
				auto vis = (*(CVisualizer **)mxGetData(prhs[1]));
				mrpt::utils::CImage img;
				if (vis->getLastWindowImage(img))
				{
					// Create the Output Image Matrix for Matlab and map it!
					size_t imgDims[3] = { img.getHeight(), img.getWidth(), 3 };
					plhs[1] = mxCreateNumericArray(3, imgDims, mxDOUBLE_CLASS, mxREAL);
					auto outData = mxGetPr(plhs[1]);

					// Map the output matrix...
					MapMatrixXd outR(outData, imgDims[0], imgDims[1]),
								outG(outData + (imgDims[0] * imgDims[1]), imgDims[0], imgDims[1]),
								outB(outData + 2 * (imgDims[0] * imgDims[1]), imgDims[0], imgDims[1]);

					mrpt::math::CMatrixFloat r, g, b;
					img.getAsRGBMatrices(b, g, r);
					
					outR = r.cast<double>();
					outG = g.cast<double>();
					outB = b.cast<double>();

					plhs[0] = mxCreateLogicalScalar(true);
				}
				else 
				{
					plhs[0] = mxCreateLogicalScalar(false);
				}
				break;
			}
			
		case 14:
			//Get Visualizer Main Camera Pose
		{
			//Create mx array 
			plhs[0] = mxCreateDoubleMatrix(6, 1, mxREAL);
			double * data = mxGetPr(plhs[0]);
			float x, y, z, Az,El,Zo;
			//Populate
			(*(CVisualizer **)mxGetData(prhs[1]))->getCameraPointingToPoint(x, y, z);
			Az = (*(CVisualizer **)mxGetData(prhs[1]))->getCameraAzimuthDeg();
			El = (*(CVisualizer **)mxGetData(prhs[1]))->getCameraElevationDeg();
			Zo = (*(CVisualizer **)mxGetData(prhs[1]))->getCameraZoom();
			data[0] = x;data[1] = y;data[2] = z;
			data[3] = Az;data[4] = El;data[5] = Zo;  
			break;
		}
		case 15:
			//Set Visualizer Pose
		{
			//The pose should be a 6x1 vector in the format [pointing at (x, y, z) Azimuth Elevation Zoom]'
			double * data = mxGetPr(prhs[2]);
			(*(CVisualizer **)mxGetData(prhs[1]))->setCameraPointingToPoint(data[0], data[1], data[2]);
			(*(CVisualizer **)mxGetData(prhs[1]))->setCameraAzimuthDeg(data[3]);
			(*(CVisualizer **)mxGetData(prhs[1]))->setCameraElevationDeg(data[4]);
			(*(CVisualizer **)mxGetData(prhs[1]))->setCameraZoom(data[5]);
			(*(CVisualizer **)mxGetData(prhs[1]))->forceRepaint();
			break;
		}
		case 16:
			//Set transform color
			(*(CVisualizer **)mxGetData(prhs[1]))->setTransformColor((*(rl::mdl::Model**)mxGetData(prhs[2])), (*(rl::mdl::Transform**)mxGetData(prhs[3])),
				mrpt::math::CMatrixDouble41::Map(mxGetPr(prhs[4])));
			break;

		default:
			mexPrintf("Unknown Visualizer Command %d\r\n",command);
			VisualizerWrapper_Usage();
			break;
	}
}

void VisualizerWrapper_Usage(void)
{
	mexPrintf("------------------CModel Visualizer WRAPPER--------------------------------\n\r");
	mexPrintf("Usage: [return] VisualizerWrapper([command id],params...)\r");
	mexPrintf("VisualizerWrapper([0] *) --> Print This Usage Message\r");
	
	mexPrintf("----------------------------Function----------------------------------------------\r");
	mexPrintf("[1]  Creates Visualizer Window\r");
	mexPrintf("		Params: string name, int width, int height\r");
	mexPrintf("		Returns: CVisualizer *\n\r");
	mexPrintf("[2]  Delete the Visualizer c++ ptr, should be called by Matlab Destructor\r");
	mexPrintf("		Params: CVisualizer *\r");
	mexPrintf("		Returns: none\n\r");
	mexPrintf(".\r");
	mexPrintf(".\r");
	mexPrintf(". please check switch statement in cpp file for other functionality \r\n");
}