import numpy as np
import time
import multiprocessing as mp
import math as m
from netCDF4 import Dataset

####################################################################################
# GENERATING THE INTERACTION FUNCTIONS                                             #
####################################################################################

def interaction_maker( names, interaction ):
    f = open("interaction_functions.py","w")
    contain_all = import_statements()
    intp = interpreter( names, interaction )

    ind_holder = []

    for i in np.arange(len(intp)):
        f_str = "def f%d( conc_00, n1 ):\n" % i
        f_str += space() + "output = "
        ind_holder.append(intp[i][0])

        for j in np.arange(1,len(intp[i])):

            if j != 1:
                f_str += " + "

            for k in np.arange(1, len(intp[i][j])):
                
                if   k % 2 == 1:
                    f_str += "np.power( conc_00[n1,%d,:,:,:] , " % intp[i][j][k]

                elif k % 2 == 0:
                    f_str += "%d)" % intp[i][j][k]

                if (k != len(intp[i][j])-1) and (k %2 == 0):
                    f_str += " * "
                elif (k == len(intp[i][j])-1):
                    f_str += " * (%f)" % intp[i][j][0]

        f_str += "\n"+returnStr() + "\n\n"
        contain_all += f_str

    f.write(contain_all)
    return ind_holder
        
        


        
       

def import_statements():
    return "import numpy as np\nimport time\nimport multiprocessing as mp\nimport math as m\nfrom netCDF4 import Dataset\n\n"

def space():
    return "    "

def returnStr():
    return "    return output"

####################################################################################
# INTERPRETING INTERACTION LIST                                                    #
####################################################################################
def interpreter( names, interaction ):
  interpreted_interaction = []
  for i in np.arange(len(interaction)):
    obj = [ names.index(interaction[i][0]) ]
    for j in np.arange(1,len(interaction[i])):

      newLine = []

      for k in np.arange(0,len(interaction[i][j])):

        if k % 2 == 1:
          newLine.append( names.index( interaction[i][j][k] ) )
        elif k % 2 == 0:
          newLine.append( interaction[i][j][k] )

      obj.append(newLine)


    interpreted_interaction.append(obj)

  return interpreted_interaction





      
