import numpy as np
import sextractorObject_py3 as SO 

def generateWCSFromSexCat(target_cat, target_dir):

    CatObject = SO.PySex(target_dir + target_cat)
    CatObject.generateFileForAstrometry()
    
    return 1
