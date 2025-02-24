from cosmosis.runtime.pipeline import LikelihoodPipeline
from cosmosis.runtime.config import Inifile


#################################################
# Read in the cosmosis configuration file pipeline.ini 
inifile = 'config/pipeline.ini'
ini = Inifile(inifile) 

output_section_name = ini.get('scale_cuts','output_section_name')
#################################################

##################################################
# (optional) You can change things or read from ini: 

# set the value of MY_PATH in the [DEFAULT] section to /Users/username/Documents/repos/
# ini.set("DEFAULT","MY_PATH","/Users/username/Documents/repos/")

# get the value of modules in the [pipeline] section
# ini.get("pipeline","modules")

# ini.set("pipeline","modules","sample_S8 sigma8toAs camb")

ini.set("DEFAULT","PROJECT_PATH","/Users/masgari/Documents/CosmicShear/emulator_project")
####################################################

####################################################
# setup the pipeline and give it the ini values
pipeline = LikelihoodPipeline(ini) 
# Get the fiducial values


####################################################
# Now execute the program. When you vary the parameters, you don't need to do the setup again.

params_fiducial = pipeline.start_vector()
# run the pipeline for the given parameters. All results are saved in block
new_params =LHCS[0]
import time
t1=time.time()
block    = pipeline.run_parameters(params_fiducial)
t2=time.time()
print(t2-t1)
theory_vector = block[output_section_name,'theory']

# theory_vector should have 75 values: 5 modes x 15 redshift bin combinations (we have 5 redshift bins which create 5*6/2 Unique combinations)
# You want to create a theory_vector for different values of parameters. Currently it will only do this for the fiducial values. 

# Other useful commands:
# list of varied parameter names
param_names_varied = pipeline.varied_params

# prints the varied paramaters as set by the values.ini file
for param in pipeline.varied_params:
 	print(param)

# # You can modify which parameters you vary or fix like this
# pipeline.set_varied("cosmological_parameters", "omch2", 0.2, 0.4)
# pipeline.set_fixed("cosmological_parameters", "h0", 0.72)

# import numpy as np
# Or you can vary them by creating a new params_values array
# params_values = np.array([0.11,  0.0225,  0.6898,  0.969 ,  0.8128,  7.8   ,  1.    ,
        # 0.    , -0.181 , -1.11  , -1.395 ,  1.265])
# Then run the pipeline again
# block    = pipeline.run_parameters(params_values)