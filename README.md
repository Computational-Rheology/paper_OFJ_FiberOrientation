Repository associated to the paper entitled **An Open-Source Framework for Modeling the Evolution of Fiber Orientation** by Ramoa *et al.*
* doi: https://doi.org/10.51560/ofj.v5.131


## The contents of the repository folders are the following: 
 * cases &rarr; Test cases reported in the paper
 * fiberOrientation &rarr; OpenFOAM code with the *functionObject* developed for solving the evolution of the 2nd-order orientation tensor 
 * fiberOrientationPython &rarr; Python scripts developed to generate the reference solutions
 * fibersClosureSymbolic &rarr;  Python scripts developed to generate the C-code code required to calculate the contractions between fourth- and second-order tensors required for the different closure models employed in the paper
 * myPimpleFoam &rarr; Custom version of pimpleFoam used to run the paper case s
 
## The scripts available in the repository's root folder have the following functionalities:
* Allwmake &rarr; Compile functionObject and myPimpleFoam
* Allwclean &rarr; Clean functionObject and myPimpleFoam


## How to run the cases?
### Make all the root folder scripts executable
```
chmod +x Allw*
```
### To compile the solvers and utilities run the Allwmake script
```
./Allwmake
```
### After compilation change to the *cases* folder check the respetive *README.md* files