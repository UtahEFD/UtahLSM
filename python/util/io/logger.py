# 
# UtahLSM
# 
# Copyright (c) 2017–2025 Jeremy A. Gibbs
# Copyright (c) 2017–2025 Rob Stoll
# Copyright (c) 2017–2025 Eric Pardyjak
# Copyright (c) 2017–2025 Pete Willemsen
# 
# This file is part of UtahLSM.
# 
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
# 

import json
import time
import netCDF4 as nc
import numpy as np
from util import constants

class Logger(object):
	
	def print_double(x,section=False,label=False):
		
		if label:
			print('%s%s: %0.17g'%(section,label,x))
		else:
			print('LOG -> %0.17g'%x)
	
	def print_hex(x,section=False,label=False):
		
		if label:
			print('%s%s: '%(section,label)+float.hex(x))
		else:
			print('LOG -> '+float.hex(x))