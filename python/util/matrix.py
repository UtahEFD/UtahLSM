#!/usr/bin/env python
import sys
import numpy as np

# Solve tridiagonal matrix using the Thomas algorithm
def tridiagonal(a,b,c,r,u):
	
	# Local variables
	n   = len(a)
	gam = np.zeros(n)
	bet = b[0]
	
	# Make sure diagonal band is not zero
	if (bet == 0.0):
		print("Error 1 in tridiag")
		sys.exit()
	
	# Initialize first element of solution vector
	u[0] = r[0]/(bet)
	
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