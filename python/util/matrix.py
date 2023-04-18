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
	
	# Local variables
	n   = len(a)
	gam = np.zeros(n)
	bet = b[0]
	print("----MATRIX----")
	for ii in range(n):
		print('{:.5f}'.format(u[ii]))
	print("--------------")
	
	# Make sure diagonal band is not zero
	if (bet == 0.0):
		print("Error 1 in tridiag")
		sys.exit()
	
	# Initialize first element of solution vector
	u[0] = r[0]/(bet)
	print('%.5f %.5f %.5f %.5f'%(a[0],b[0],c[0],r[0]))
	# Forward sweep 
	for j in range(1,n):
		gam[j] = c[j-1]/bet
		bet    = b[j]-a[j]*gam[j]
		if bet == 0.0:
			print("Error 2 in tridiag")
			sys.exit()
		u[j]=(r[j]-a[j]*u[j-1])/bet

	# Backward sweep
	for j in range(n-2,-1,-1):
		u[j] -= gam[j+1]*u[j+1]
	
	print("----MATRIX 2----")
	for ii in range(n):
		print('{:.5f}'.format(u[ii]))
	print("--------------")