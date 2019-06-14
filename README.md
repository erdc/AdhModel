# AdhModel
| All | Windows | Mac |
| --------- | --------- | ------------- |
|[![Build Status](https://dev.azure.com/kimberlycpevey/kimberlycpevey/_apis/build/status/kcpevey.AdhModel?branchName=master)](https://dev.azure.com/kimberlycpevey/kimberlycpevey/_build/latest?definitionId=4&branchName=master)|[![Build Status](https://dev.azure.com/kimberlycpevey/kimberlycpevey/_apis/build/status/kcpevey.AdhModel?branchName=master&jobName=vs2017-win2016)](https://dev.azure.com/kimberlycpevey/kimberlycpevey/_build/latest?definitionId=4&branchName=master)|[![Build Status](https://dev.azure.com/kimberlycpevey/kimberlycpevey/_apis/build/status/kcpevey.AdhModel?branchName=master&jobName=macOS-10.13)](https://dev.azure.com/kimberlycpevey/kimberlycpevey/_build/latest?definitionId=4&branchName=master)

  The AdhModel object contains all the data necessary to describe an AdH model. Each AdhModel contains an AdhMesh object and each mesh object contains one or more AdhSimulation objects. 
    
  The AdhModel object is not directly tied to any visualization package. While it does contain Holoviews objects, these may be rendered with either Bokeh or Matplotlib.
  
### Structure
* AdhModel
  * AdhMesh
    * AdhSimulation
     
  
### Inheritance
AdhModel inherits from the GeneSIs base Model class. 
AdhMesh inherits from the GeneSIs base Mesh classes. 
AdhSimulation inherits from the GeneSIs base Simulation class. 
  
#### AdhMesh
Each model object contains one mesh object. A mesh object may have one or more simulations. While this package is not the proper location for visualizations, the mesh class does contain several convenience functions for viewing the mesh. 

#### AdhSimulation
The simulation object is the container for the boundary conditions and the initial conditions for a simulation. Additionally, the results, if available, are also located in the AdHSimulation object. Each mesh may contain one or more simulations. For example, for a model which is analyzing base versus plan, there will be two simulation objects - one each for base and plan.
  
### Developer Notes
To opt out of a Azure Pipelines CI build, add [skip ci] or [ci skip] to the commit message. 