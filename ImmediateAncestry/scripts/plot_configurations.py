import sys, math
import re 
import os
import subprocess
import numpy

from os import scandir, getcwd


# ## # ## # ## # ## # ## # ## # ## # 

def ls(ruta = getcwd()):
    return [arch.name for arch in scandir(ruta) if arch.is_file()]



def take_config (filepath,gen):
	cols=2*2**gen
	config=[]
	file_object  = open(filepath, "r")
	for conf in file_object.read(cols):
		if conf != " ":
			config.append(conf)
	return config



def read_plot_files (path, generations):
	lista=ls(path)

	configs={}
	configs_rho={}

	for fil in lista:

		compl_path=path+'/'+fil
		config="".join(take_config(compl_path,generations))

		if 'rhoinf_out' not in fil:


			if config in configs:

				num=configs[config]
				configs[config]=num+1
			else:

				configs[config]=int(1)

		else:

			if config in configs_rho:

				num=configs_rho[config]
				configs_rho[config]=num+1
			else:

				configs_rho[config]=int(1)


	filpath=path+"/rho_ConfigsCount.txt"
	f = open(filpath, 'w')

	for key in configs:
		f.write(key)
		f.write("\t")
		f.write(str(configs[key]))
		f.write("\n")
	#	print ("config",key, " : ", configs[key])

	f.close()

	ffilpath=path+"/rhoInfinity_ConfigsCount.txt"
	ff = open(ffilpath, 'w')

	for key in configs_rho:
		ff.write(key)
		ff.write("\t")
		ff.write(str(configs_rho[key]))
		ff.write("\n")
	#	print ("config",key, " : ", configs_rho[key])

	ff.close()


def read_pure_files (path,gen):
	pures=['i231_Cindy', 'i232_Mirinda', 'i233_Alfred', 'i234_Ula', 'i235_Lara', 'i236_Luky', 'i237_Gamin', 'i238_Brigitte', 'i239_Vaillant', 'i240_Doris', 'i241_Julie', 'i242_Clara', 'i243_Noemie', 'i244_Yogui', 'i245_Tibe', 'i246_Blanquita', 'i247_Negrita', 'i248_Marlin', 'i249_Bosco', 'i250_Donald', 'i251_Jimmie', 'i252_Berta', 'i253_Annie', 'i254_Mike', 'i255_SeppToni', 'i256_Linda', 'i257_Clint', 'i258_McVean', 'i259_Cindy', 'i260_Alice', 'i261_Koby', 'i204_Akwaya_Jean', 'i205_Banyo', 'i206_Basho', 'i207_Damian', 'i208_Julie', 'i209_Kopongo', 'i210_Koto', 'i211_Paquita', 'i212_Taweh', 'i213_Tobi', 'i201_Diana', 'i202_Ikuru', 'i203_Trixie', 'i214_Vincent', 'i215_Andromeda', 'i216_Harriet', 'i217_Bwambale', 'i218_Kidongo', 'i219_Nakuu', 'i220_Padda', 'i221_Cindy', 'i222_Frederike', 'i223_Washu', 'i224_Athanga', 'i225_Coco', 'i226_Mgbadolite', 'i227_Tongo', 'i228_Cleo', 'i229_Bihati', 'i230_Maya']
	i=0
	filpath=path+"/Pures.txt"
	f = open(filpath, 'w')
	ffilpath=path+"/Pures_rhoInfinity.txt"
	ff = open(ffilpath, 'w')

	for fil in pures:

		compl_path=path+'/'+fil+'_2.out'
		compl_path_rhoinf=path+'/'+fil+'_2.rhoinf_out'

		config="".join(take_config(compl_path,gen))
		config_rhoinf="".join(take_config(compl_path_rhoinf,gen))

		i=i+1
	
		if i == 1 :
			f.write('Ptt\n')
			ff.write('Ptt\n')
		elif i == 19 :
			f.write('\nPtv\n')
			ff.write('\nPtv\n')
		elif i ==  32:
			f.write('\nPte\n')
			ff.write('\nPte\n')
		elif i == 42 :
			f.write('\nPts\n')
			ff.write('\nPts\n')
		
		f.write(fil)
		f.write("\t")
		f.write(str(config))
		f.write("\n")
		
	
		ff.write(fil)
		ff.write("\t")
		ff.write(str(config))
		ff.write("\n")
	f.close()
	ff.close()
	
	
	
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## D a t a ##
path='/home/nat/NATALIA/ImmediateAncestry/ImmediateAncestry/scripts/SIMULATIONS/LikelihoodSimulation/4Gen_1RecombMax'
generations=4


read_plot_files (path,generations)
read_pure_files (path,generations)
