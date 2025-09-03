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

import logging
import numpy as np

# local logger
logger = logging.getLogger(__name__)

# Solve tridiagonal matrix using the Thomas algorithm
def tridiagonal(a,b,c,r,u):

	# Local variables
	n   = len(a)
	gam = np.zeros(n)
	bet = b[0]
	
	# Make sure diagonal band is not zero
	if (b[0] == 0.0):
		logger.error("Error 1 in tridiag")
		raise SystemExit(1)

	# Initialize first element of solution vector
	u[0] = r[0]/(bet)
	
	# Forward sweep 
	for j in range(1,n):
		gam[j] = c[j-1]/bet
		bet    = b[j]-a[j]*gam[j]
		if bet == 0.0:
			logger.error("Error 2 in tridiag")
			raise SystemExit(1)
		u[j]=(r[j]-a[j]*u[j-1])/bet
	
	# Backward sweep
	for j in range(n-2,-1,-1):
		u[j] -= gam[j+1]*u[j+1]
