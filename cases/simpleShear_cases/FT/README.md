## The contents of this folder are the following:

### Scripts:
* Allrun &rarr; run the cases 
* Allclean &rarr; clean the directory

### Folders: 
* PostProcessingData &rarr; python data and postProcessing scripts
* base &rarr; base case 

## How to run the cases, generate and check the results:
* To run the case execute the *Allrun* script
```
./Allrun
```
* After running change to the *PostProcessData* folder and run the *runBatch* script, which is created as symbolik link to *...simpleShear_cases/pythonScriptToPlotData_simpleShear/runBatch_simpleShear* script:
```
./runBatch FT
```
* The postprocessing Python script requires the following packages: *pandas*, *numpy* and *matplotlib*. Refer to the Python documentation to install these packages before running the script
* The graphs and sampling points are stored, respectively, in "svg" and "csv" formats, in folders .... */simpleShear_cases/FT/`<closure>`/mesh`<i>`/postProcessing/sampling/`<finalTime>`*
