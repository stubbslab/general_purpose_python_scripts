# useful-python-tools
A repository of useful python code to accomplish general tasks. 

Motivation: often, users of Python want a set of generally useful functions/operations that do not, by default, exist. 
We hope that members of the Stubbs research group will share with each other, via this repository, generally useful python functionality.  

Usage: The user can download some or all of the Python scripts in this directory for their personnal use.  
       Most scripts require third-part Python libraries (such as numpy and astropy).  
       Python should inform you which libraries you must install when you attempt to import this Python code by throwing an error like:
       ModuleNotFoundError: No module named 'numpy' 
       To get around this error, you need to install numpy.  This can be done simply through, for example, pip: 
       $ pip install numpy 
       
Interdependence: Some of these python code is interdependent.  For example, the IntuitiveFitsFileManagement.py uses some of the functions in cantrips.py. 
                 To use the former, you need the latter.  
                 If you attempt to use IntuitiveFitsFileManagement.py without also downloading cantrips.py, you will get an error. like: 
                 ModuleNotFoundError: No module named 'cantrips'
                 To circumvent this error, you must also download the cantrips.py file and place it either in the same directory as IntuitiveFitsFileManagement.py or somewhere in the Python path. 
                 
                 
                 
       
