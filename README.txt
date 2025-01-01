This folder contains: 

Folders: 
* cases -> Test cases reported in the manuscript
* fiberOrientation -> OpenFOAM code with the functionObject solving the evolution of the 2nd-order orientation tensor 
* myPimpleFoam -> Custom version of pimpleFoam

Scripts:
* Allwmake -> Compile functionObject and myPimpleFoam
* Allwclean -> Clean functionObject and myPimpleFoam

* to assure all the requeired files to run are executable make the file makeExecutable executable by running the following commands:
>> chmod +x makeExecutable
>> ./makeExecutable