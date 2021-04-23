#!/usr/bin/env python3
from rich import print
import sys
import subprocess
from ATACofthesnake import misc
import yaml
import os

# Read / set variables.
with open('Parameters.yaml') as f:
	paramDic = yaml.load(f, Loader=yaml.FullLoader)

