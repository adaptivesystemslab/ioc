# Inverse optimal control

## To add an additional IOC model and trials

- Make a new file in `Data_json/[new_folder]/[new_trial.json]`. These json files sets up the run configurations and specific files to load. So it applies "runParamGlobal" to all the settings, then gets overwritten with "runTemplates", then gets overwritten with "Files". You can see the final output in `inverseIOC/ioc_main.m` line 58, in `trialInfo`. This sets the parameters for what to load. 
- In `InverseIOC/IOCRun.m`, you see that the line 5 loads model. This getmodel file is in `logic/getmodel.m`. You need to add an entry that matches `trialInfo.baseModel` (the experiment type) and `trialInfo.model` (any sub models you may have) to the switch function. 
- Then make a copy of https://github.com/adaptivesystemslab/ioc/blob/master/Logic/ModelRL_Template.m and update this appropriately. It's mainly making sure the joint angle array are loaded properly, as well as any kinematic/dynamic models. 

- To add new cost functions, you need to name it in `common\featuresEnum`. This name is what `trialInfo` is looking for. Then just add a corresponding entry in `featuresCalc`. If you look in `featuresCalc.initCf`, you will see that each enum has a set of basis fucntions (bf), and the actual cost function (cf) itself. The basis function is the fundamental feature required. So how to go from the model kinematics/dynamics to the cost function. It's listed here so it gets calculated at each iteration, to save computational effort. You can either just use a inline function (ie featuresEnums.cartVeloSumSqu), or write your own function if the cf calculation has a few steps (ie featuresEnums.cartCurvatureSumSqu)
- The gradient should be calculated

Once this is done, you should be able to run the function, and graph the outputs using `analysis\iocAnalysis00_plot_demo'
