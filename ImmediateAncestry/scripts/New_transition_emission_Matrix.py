import numpy

def transformation_TransitionMatrix(TransMat_old,num_recombs):

	n,m=numpy.shape(TransMat_old)
	mat_ppal=numpy.eye(n)
	mat_sec=numpy.fliplr(mat_ppal)
	mat_double=numpy.multiply((mat_ppal*-1)+1,(mat_sec*-1)+1)
	
	T1=numpy.diagflat(numpy.diag(TransMat_old))
	T2=numpy.multiply(TransMat_old,mat_double)
	T3=numpy.multiply(TransMat_old,mat_sec)


	new_matrix=[]

	if num_recombs == 1:	# toy example. 1 recombs allowed
		return numpy.bmat([[T1,T2],[numpy.zeros((n,m)),T1]])
		
	elif num_recombs == 2:	# toy example. 2 recombs allowed
		return numpy.bmat([[T1,T2,T3],[numpy.zeros((n,m)),T1,T2],[numpy.zeros((n,m)),numpy.zeros((n,m)),T1]])
		
	elif num_recombs >2 :	# More than 2 recombs

		for i in range (num_recombs+1):	# num_recombs+1 = num matrix rows, i= currently row number
			fila=[]

			if i == 0:	#First row i=0
				#fila.append(T1)
				#fila.append(T2)
				#fila.append(T3)
				#print ("T1,T2,T3")
				fila=numpy.bmat([[T1,T2,T3]])

				for j in range ( num_recombs-i-2):	#cerosdesp
					fila=numpy.bmat([[fila,numpy.zeros((n,m))]])
					#fila.append(numpy.zeros((n,m)))
					#print ("0")

				new_matrix=fila

				
			elif (i >0 and i < num_recombs-1): # Rows with T1,T2,T3 
				
				fila=numpy.zeros((n,m))
				for j in range (i-1): #ceros antes

					fila=numpy.bmat([[numpy.zeros((n,m)),fila]])
					#fila.append(numpy.zeros((n,m)))
					
				#fila.append(T1)
				#fila.append(T2)
				#fila.append(T3)
				fila=numpy.bmat([[fila,T1,T2,T3]])

				if ((num_recombs-i-2) >=0):
				
					for j in range ( num_recombs-i-2):	#cerosdesp
						#fila.append(numpy.zeros((n,m)))
						fila=numpy.bmat([[fila,numpy.zeros((n,m))]])
				
			
				new_matrix=numpy.bmat([[new_matrix],[fila]] )
				
			else:	#Rows with T2 or T3 missing
				
				fila=numpy.zeros((n,m))
				if i == num_recombs-1:

					for j in range (i-1): #ceros antes
						fila=numpy.bmat([[numpy.zeros((n,m)),fila]])
						#fila.append(numpy.zeros((n,m)))
						
					#fila.append(T1)
					#fila.append(T2)
					fila=numpy.bmat([[fila,T1,T2]])
					
				else:
					for j in range (i-1): #ceros antes
						fila=numpy.bmat([[numpy.zeros((n,m)),fila]])
						#fila.append(numpy.zeros((n,m)))
					fila=numpy.bmat([[fila,T1]])	
					#fila.append(T1)
	
				new_matrix=numpy.bmat([[new_matrix],[fila]])
				

	else:	
		#Less than 2 recombs allowed, change the matrix_generation function to make the OldtransitinMat
		print(" num recombs doesn't match with the choosen Transsition Matrix type")
	
	return new_matrix
				


def transformation_EmissionMatrix(EmMatr_old,num_recombs):
	em_mat=EmMatr_old
	for i in range (num_recombs):
		
		em_mat=numpy.bmat([[em_mat],[EmMatr_old]])

	return em_mat

def initial_matrix(gen,max_numberof_recombs,num_recombs):
	M=(2**(gen-1))**2
	M=[1.0/M]*M
	leng=len(M)
	if  max_numberof_recombs != -1:
		for i in range (num_recombs):
	
			for j in range (leng):
				M.append(0)			
	return M
		


##################
    
#TransMat_old=transition_matrix(0.03, 2)
#TransMat_new=transformation_TransitionMatrix(TransMat_old,4)

#print ("\nFinal\n",TransMat_new[0])

#EmsMat_old=generate_emission_matrix(,2)
#EmsMat_new=transformation_EmissionMatrix(EmMat_old,4)
