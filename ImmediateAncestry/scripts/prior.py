from ggrandparents_model import find_smallet_equivalence_class
from math import factorial
from itertools import combinations, combinations_with_replacement
import sys, math
from time import time, localtime, strftime
from itertools import product
import re 


lista= ['e', 'v', 't', 's']
generation=[1,2,3,4]
priors={}		#Configuraciones representativas "Smallest class"
archivo = open("archivo2.txt","w")


for gen in generation:
	nindiv=(2**gen)

	for config in product(lista, repeat=nindiv):
		small=find_smallet_equivalence_class(config)

		if small in priors:
			priors[small]=int(priors[small])+1

		else:
			priors[small]=1



for key in priors:
#for key in configs:
	archivo.write(str((key)))
	archivo.write(str(("\t")))
	archivo.write(str((priors[key])))
	archivo.write(str(("\n")))

archivo.close()



def prior (x):
	small = find_smallet_equivalence_class(x)




#_____________________________________________________________

#prior("evsttsve")

#permutacion(["e","v","s","t","t","s","v","e"],)

x=prior (['e', 's', 's', 's', 's', 't', 't', 't'])
print('x=', x)







