This project provides two program parts
Gui_generated launches a gui that allows you to easily enter all your
parameters. these are then saved as plain text files in the experiments
under "standard" and whatever name you enter. you can always reload an
saved experiment. the "run" button or the "run loop" button simply
starts the second code "python_calculations"
This code does load whatever is written in the experiment\standard and
executes the program from there.
that means you can also use this in the command line.

The whole program is pretty much a more fancy wrapper for the xraylib
project

'''this tool is to calculate the total fluescence yield measured for
solid and liquids
the data used is either generated from the installed xraylib
http://ftp.esrf.eu/pub/scisoft/xraylib/readme.html
or taken from a database that was previously generated from this source
all set values are handled in three panda Dataframes and explained in
them
The frames are stored in a folder that can be either specificly given or
the standard folder is used
I can highly recomment to use the xraylib life since additional
functions like the compound parser will be available
see the import values function for more details
ALL VALUES USED FOR CALCULATION ARE IN SI!!!!
developed by
Jens Uhlig 2013
The compund library was snatched from Bruce ravels Hephaestus