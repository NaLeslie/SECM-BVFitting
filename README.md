# SECM-BVFitting
Program for fitting SECM images pairs with Butler-Volmer Kinetics

## Usage
Unzipping `if_252mV-3.zip` (~47MB uncompressed) gives an example instruction file which the given SecmBV.java is set-up to fit.  
Linux commands shown here. COMSOL is called from the console slightly differently on Windows.
- Edit any of the physical parameters defined at the top of the SecmBV class.
- Edit `control_file = ""` in the main method to the desired instruction file (or leave alone if fitting the file in `if_252mV-3.zip`).
- Compile: `comsol compile SecmBV.java`
- Run: `comsol batch -inputfile SecmBV.class -batchlog batlog.log > fitlog.log`

Note: you will need to enable the following in COMSOL's `Preferences > Security > Methods and Java Libraries` menu:
- Allow access to system properties
- File system access: All files

### Instruction file
Instruction files have the following format:  
```
##ENCODING: csv
##X-Sampling: StartIndex,StepSize,NumberOfSteps
#22,4,17
##Y-Sampling: StartIndex,StepSize,NumberOfSteps
#9,4,17
##xindex,yindex,switch,xcoord/m,ycoord/m,OXcurrent/A,REDcurrent/A
0,0,0,0.0,0.0,9.228357295097126E-10,-9.228357295097126E-10
0,1,0,0.0,2.0E-6,9.22923353766823E-10,-9.228357295097126E-10
 ...
 ```
These files can be generated in part by my [SECM-SEM-Align](https://github.com/NaLeslie/SECM-SEM-Align) project.
Switch defines whether the pixel at the given set of coordinates is initially 'reactive' or not, and the OX and RED currents are the microelectrode currents at that position for an oxidising and reducing microelectrode respectively.

## Publication
This work is associated with the publication:  
Leslie, N.; Mauzeroll, J. Quantification of Butler-Volmer Electron Transfer Kinetics from Scanning Electrochemical Microscopy Images. ACS Electrochemistry 
which is currently under peer review.
