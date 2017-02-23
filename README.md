For a complete overview and instructions for script execution for both R and Matlab, please visit the project site [https://maciverlab.github.io/projects/bigeye/bigeye-explanation.html](https://maciverlab.github.io/projects/bigeye/bigeye-explanation.html) 

##Preparing your system and running our R code

1) Install R.

Go to https://www.r-project.org/ and download the version appropriate
for your operating system. The latest R version we checked our code
with is R 3.3.2.

Note for Mac users: we are not using any X11 graphics devices in
our code, but if you plan to use R with X11 for other pojects, make
sure to install XQuartz (https://www.xquartz.org/).


2) Ensure your system can execute Rscript from the command line or
terminal.

For your convenience, we have written the code so you can run all
the scripts simply from command line (Windows) or terminal (Mac).
We'll call it terminal from here on to be less wordy.

When working with Windows operating systems it is necessary to add
the location of the Rscript executable file (Rscript.exe) to the
PATH, an environment variable that specifies the directory of
executable programs. The first step is to find the directory of the
Rscript.exe file using the Explorer. A typical location may look
like this:

C:\Program Files\R\R-3.3.2\bin\

In the Control Panel, navigate to System => Advanced System Settings
and click on Environment Variables. That will open a new window
displaying user and system variables boxes. In the system variables
box, highlight the row with the Path variable and click Edit. In
Windows 7 and 8 the PATH is a simple text string. Add the R executable
directory at the very end of that string, preceded by a semicolon:

;C:\Program Files\R\R-3.3.2\bin\

and click OK after which you are finished. In Windows 10, a new
dialogue box will appear that lists all PATH directories on separate
rows in a box. Click New and paste the Rscript.exe file path (without
semicolon!) into the new row. Click OK and you are done.

When working with Macs, you do not have to worry about adding the
directory to the PATH as it should already be included. You can
check this by typing

which Rscript

in the terminal and hitting the return key. The output may look
something like this:

/usr/local/bin/Rscript

If nothing is returned, you normally can fix the problem by
reinstalling R. This may become necessary after a major OS X upgrade.


3) Installing R packages

Our code relies on a number of R packages (collections of existing
functions). Instead of having you manually install each package one
by one, we have written an installer script that will take care of
this for you. Open the terminal in the main folder of our code
repository, type

Rscript installer.R -cwd

and hit return.


4) Executing code

Once you are done with steps 1)-3) you are ready to run our R
scripts. We have created a unique folder for each figure, or in
some cases specific panels of a composite figure. For example, if
you want to recreate Figure 2, simply open the terminal in the
folder for this figure and execute the R file contained within that
folder:

Rscript fig02.R -cwd

The output, in this case Figure 2 in *.pdf format, will be stored
in the same folder.

We have also provided code to carry out the entire selective regime
analysis implemented in bayou. You will find two R scripts in the
bayou folder and you can execute the code as outlined above. Both
scripts perform the analysis over the entire tree set, one calculating
residuals with a BM correlation struture (presented in the paper),
the other with an OU correlation structure. Even though the code
is written to use all cores of your processor, the run time will
likely be several days. Results will be stored in form of *.csv
files in the bayou folder at the very end of the computation. We
supply our bayou-results in /data/paleo/bayou_output.


5) In the unlikely event our code will not work with future R
releases you can find and install previous releases of R here:
https://cran.r-project.org/bin/windows/base/old/ The latest R version
we checked our code with is R 3.3.2. If it turns out you need to
roll back to this R version you must also roll back the packages
required for our code, as packages are continuously developed
alongside the new R updates. Old packages may not work with current
versions of R and vice versa. You can roll back to old packages by
executing an alternate installer-file, but you really should only
do this if you have run into problems related to R and package
updates/new versions.

Rscript rollbackinstaller.R -cwd

Windows users may have to install Rtools first (there will be a
prompt), which requires a subsequent restart.

## Running Matlab

Matlab functions require curve fitting and image processing toolboxes

Before running the Matlab code run the startup code [bigeye_startup](https://github.com/maciverlab/bigeye/blob/master/bigeye_startup.m), to link all files and folders. 

The data files from HydroLight (Sequoia Scientific Inc., Bellevue WA) simulations can be found under [figs/data/vision/hydrolight](https://github.com/maciverlab/bigeye/tree/master/figs/data/vision/hydrolight). Baseline water conditions that are used to generate Figure 4 in the main paper are located in folders base_moon, base_stars, and base_sun. The folder base_sun is also used as the baseline water condition for simulating the effect of water conditions on visual range, as depicted in Supplementary Appendix Figure 6. 


The function *getBugContrast.m* in [figs/figExt07_contrast/image_contrast](https://github.com/maciverlab/bigeye/blob/master/figs/figExt07_contrast/image_contrast/getBugContrast.m) will ask the user to select the object of interest, do not close the figure window until all the images have been exhausted.

## Literature
All cited literature can be found in the cited_literature folder; for details and instructions please refer to the README in that directory.

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">This repository for code and data to reproduce results within ”Massive increase in visual range preceded the origin of terrestrial vertebrates”</span> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
