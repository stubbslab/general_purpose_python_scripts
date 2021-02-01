#This script executes the three commands necessary to make a function updated after making changes to the original script.
#This assumes that the name of the function and the script are the same

def qRL(funct_to_reload):
    import funct_to_reload
    reload(funct_to_reload)
    from funct_to_reload import funct_to_relaod
    

