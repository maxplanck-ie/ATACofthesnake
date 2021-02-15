import os
import yaml
from ATACofthesnake import misc

# Read / set variables.
with open('Parameters.yaml') as f:
	paramDic = yaml.load(f, Loader=yaml.FullLoader)
Conditions = list(paramDic['Cond'])