# 
# UtahLSM
# 
# Copyright (c) 2017–2023 Jeremy A. Gibbs
# Copyright (c) 2017–2023 Rob Stoll
# Copyright (c) 2017–2023 Eric Pardyjak
# Copyright (c) 2017–2023 Pete Willemsen
# 
# This file is part of UtahLSM.
# 
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
# 

import sys
import numpy as np

# Solve tridiagonal matrix using the Thomas algorithm
def tridiagonal(a,b,c,r,u):
	print()        
	print("tridiag---------")
	# Local variables
	n   = len(a)
	gam = np.zeros(n)
	bet = b[0]
		
	print('n: %d'%n)
	print("----------------")
	
	# Make sure diagonal band is not zero
	if (b[0] == 0.0):
		print("Error 1 in tridiag")
		sys.exit()

	# Initialize first element of solution vector
	u[0] = r[0]/(bet)
	print("uj (0): %.17f"%(u[0]))
	# Forward sweep 
	for j in range(1,n):
		print("----------------")
		gam[j] = c[j-1]/bet
		bet    = b[j]-a[j]*gam[j]
		if bet == 0.0:
			print("Error 2 in tridiag")
			sys.exit()
		u[j]=(r[j]-a[j]*u[j-1])/bet
		print("gm (%d): %.17f"%(j,gam[j]))
		print("bt (%d): %.17f"%(j,bet))
		print("uj (%d): %.17f"%(j,u[j]))
	print("----------------")
	# Backward sweep
	for j in range(n-2,-1,-1):
		u[j] -= gam[j+1]*u[j+1]