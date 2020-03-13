import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['figure.figsize'] = [10, 10]
matplotlib.rcParams['font.size'] = 20

Twist_angles = np.array([-15,-10,-5,-3,-1,0,1,3,5,10,15])
layer_number = 10

for alpha in Twist_angles:
	ax = plt.subplot(111, projection='polar')
	ax.grid(False)
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	ax.spines['polar'].set_visible(False)
	ax.set_rlim(0,2*layer_number)

	k = 1.-alpha/360.
	sqrt_3 = np.sqrt(3.)


	beta0_0 = (4*np.pi/6.)/k
	beta0_1 = (7*np.pi/6.)/k
	beta0_2 = (np.arctan(sqrt_3*(-1)/5)+2*np.pi)/k


	r0_1 = np.arange(0, 1., 0.01)
	theta0_1 = np.full(r0_1.shape, beta0_0)
	
	theta0_2 = np.arange(beta0_0, beta0_1, 0.01)
	r0_2 = sqrt_3*(0.5)*k/np.cos(k*theta0_2-5*np.pi/6.)
	
	theta0_3 = np.arange(beta0_1, beta0_2, 0.01)
	r0_3 = sqrt_3*(0.5)*k/np.cos(k*theta0_3-9*np.pi/6.)
	
	ax.plot(theta0_1, r0_1,color = "black")
	ax.plot(theta0_1+np.pi, r0_1,color = "red")
	
	ax.plot(theta0_2, r0_2,color = "black")
	ax.plot(theta0_2+np.pi, r0_2,color = "red")
	
	ax.plot(theta0_3, r0_3,color = "black")
	ax.plot(theta0_3+np.pi, r0_3,color = "red")

	for n in range(1,layer_number):
		two_n_pi = 2*n*np.pi

		beta0 = (np.arctan(sqrt_3*(-1)/(4*n+1))+two_n_pi)/k
		beta1 = (np.arctan(sqrt_3*(2*n+1)/(2*n-1))+two_n_pi)/k
		beta2 = (2*np.pi/3.+two_n_pi)/k
		beta3 = (np.arctan(sqrt_3*(1)/(4*n+3))+np.pi+two_n_pi)/k
		beta4 = (np.arctan(sqrt_3*(2*n+1)/(2*n+3))+np.pi+two_n_pi)/k
		beta5 = (np.arctan(sqrt_3*(-1)*(2*n+1)/(2*n+5))+2*np.pi+two_n_pi)/k
		beta6 = (np.arctan(sqrt_3*(-1)/(4*n+5))+2*np.pi+two_n_pi)/k

		theta1 = np.arange(beta0, beta1, 0.01)
		r1 = sqrt_3*n*k/np.cos(k*theta1-1*np.pi/6.+two_n_pi)
		theta2 = np.arange(beta1, beta2, 0.01)
		r2 = sqrt_3*(n+0.5)*k/np.cos(k*theta2-3*np.pi/6.+two_n_pi)
		theta3 = np.arange(beta2, beta3, 0.01)
		r3 = sqrt_3*(n+0.5)*k/np.cos(k*theta3-5*np.pi/6.+two_n_pi)
		theta4 = np.arange(beta3, beta4, 0.01)
		r4 = sqrt_3*(n+1)*k/np.cos(k*theta4-7*np.pi/6.+two_n_pi)
		theta5 = np.arange(beta4, beta5, 0.01)
		r5 = sqrt_3*(n+0.5)*k/np.cos(k*theta5-9*np.pi/6.+two_n_pi)
		theta6 = np.arange(beta5, beta6, 0.01)
		r6 = sqrt_3*(n+1.5)*k/np.cos(k*theta6-11*np.pi/6.+two_n_pi)

		ax.plot(theta1, r1,color = "black")
		ax.plot(theta1+np.pi, r1,color = "red")
		ax.plot(theta2, r2,color = "black")
		ax.plot(theta2+np.pi, r2,color = "red")
		ax.plot(theta3, r3,color = "black")
		ax.plot(theta3+np.pi, r3,color = "red")		
		ax.plot(theta4, r4,color = "black")
		ax.plot(theta4+np.pi, r4,color = "red")		
		ax.plot(theta5, r5,color = "black")
		ax.plot(theta5+np.pi, r5,color = "red")		
		ax.plot(theta6, r6,color = "black")
		ax.plot(theta6+np.pi, r6,color = "red")
	

	ax.set_title("Twist Angles = "+str(alpha), va='bottom')
	plt.savefig("Hex_Twist_angles_10_"+str(alpha) + ".pdf", format='pdf', dpi=150, transparent=True)
	plt.clf()
	plt.close()
