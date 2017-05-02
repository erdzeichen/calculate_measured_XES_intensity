This is a little python module to predict the measured flux in a emission spectroscopy setup. At the moment it is written for a crystal based setup, however with a few tricks it can be used for microcalorimeters too(drop me a mail to get explanations). The Knut och Alice Wallenbergs Stiftelse financed me during the stay at the ESRF during which this code was written and I am very grateful for their support.

Simply double click the "gui_generated.py".  each field has its explanations if you hover over it. The basis of this program is the xraylib library from the ESRF. I however compiled a offline version of the tool, so it can be run without the library installed.



#------------In the gui-------------------------

The gui is only the frontend and dumps all information first in three pandas.Dataframe and then as Ascii in three files on the disc. After this it calls the functions from "python_calculation.py"
By using the functions from "python_calculation.py" in the command line the information can be used further. "python_calculation.py" runs if not otherwise set the "standard"experiment which is the last run.
So i mostly use the gui to set the values and then go to the command line

(floating the mouse over a field gives a text which might help to use the tool)
pressing enter in any text field executes the program and gives the expected countrate for the settings the same does clicking on "run single". The "run loop" function loops over different parameters.
I programmed two loops. the inner loop is designed to be run over many values and typically is used for energy scans and the like.
The outer loop is run over any of the parameters but takes much longer for large sets of data
if the loop is run the results are shown in a plot.

#--------command line usage :-----------------

#----general usage---------- (examples below) 
all functions nees something close to:
get_AtomicLevelWidth(Atoms='Ar',Shells='K')
get_AtomicLevelWidth(Atoms=['Ar','Ag'],Shells=['K','L3'])
get_LineWidth(Atoms='Ar',Lines='K-M3')
compounds()  returns a long list with compopunds the program knows already (can be easily edited)

usefull functions include:
get_LineEnergy    
get_EdgeEnergy
get_absorb        #returns a DataFrame with the material absorption
get_LineWidth
get_AtomicLevelWidth
get_partial_XRF_cross     #returns  a DataFrame with fluoresence in a specific line
get_total_cross           #returns  a DataFrame with the total fluoresence
relation_wrapper          #neatly combining several materials and calculates solutions! very usefull for work with liquids

emission_detection
vary_something
vary_a_second

#---emission----
emission_detection() runs the tool once and gives the expected emission intensity using the data from the "Default directory"
This is also the data that was last run , so I recomment do run the tool in a separate console to generate this comfortably there are also a nubmer of other usefull tools:
run_loop_from_save()  runs the loops

#---Absorption---
calculates the absorption as function of Energy. This can be a single value or a iterable. 

get_absorb(compound="Kapton",Energy=range(1000,2000,10)) returns a nice Dataframe with the compund data, if the compounds are not known then you have to enter something like this instead:
get_absorb(compound="FeC10H10",density=6,Energy=range(1000,2000,10))



#---------comments-------------------------------

I tried to comment each function extensively, if you have questions however don't hesitate to drop me an email or leave a comment.
I'm actively working on the tool so there might still be some bugs and i'm working on more functionality. Just let me know about bugs or wishes and i will try to fix/include them.

Hint beside the standard python modules you will need the python "wx" and "pandas" packages. They come with most installations based on Pythonxy or EPD (see https://sites.google.com/site/quickwrapupoffreesoftware/home/programming). Under Linux you will have to install these packages. The first start might take a few seconds, since some modules will be precompiled.

I will spend some time later and include a module to calculate the absorption spectra and a special module for microcalorimeter and broadband excitations.
for now:

Broadband calculation: make a plot with emission intensity as function of incoming photon energy and a fixed number of incoming photons. Multiply this afterwards with the source spectrum
For microcalorimeter: set the fudge factor 1, the reflectivity 1, one crystal and the diameter corresponding to the effective area of the detector. After the emission is calculated multiply with the spectral sensitivity of the detector, et voila

The compund library was snatched from Bruce Ravels Hephaestus
