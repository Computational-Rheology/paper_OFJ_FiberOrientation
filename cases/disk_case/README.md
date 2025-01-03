## The contents of this folder are the following

### Scripts:
* Allrun &rarr; run the cases 
* Allclean &rarr; clean the directory

### Folders: 
* PostProcessingData &rarr; python data and postProcessing scripts
* base &rarr; base case 

## How to run the cases and generate the results:
* To run the case execute the *Allrun* script
```
./Allrun
```
* After running move to the PostProcessData folder and run the runBatch_disk script (an internet connection is required to copy the theorectical results from the paper repository):
```
./runBatch_disk
```
* The graphs and sampling points are stored in "svg" and "csv" formats, respectively, in folders .... *disk_case/IBOF/mesh`<i>`/postProcessing/sampling/`<finalTime>`*

## *Notes* :
* The developed Python script requires the packages pandas, numpy and matplotlib (refer to the Python scripts to )



