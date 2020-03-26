# -*- coding: utf-8 -*-
# Name: parsejason_v2.py
# Autor: Ivan
# Fecha: 06.08.2019
# Uso:  parsejason_v2.py -f inFile.json
# Descripción: parsea un archivo en formato json
from __future__ import division
import argparse
import re
import phyphy


####################################################################
##### lee argumento: archivo de entrada en linea de comandos    ####
####################################################################
parseoArgs = argparse.ArgumentParser()
parseoArgs.add_argument(
		"-u", "--usage",
			help = "Usage  parsejason_v2.py [-f] InFile.json",
			action="store_true")

parseoArgs.add_argument(
		"-f", "--file",
			help="Archivo de entrada en formato json", required=True)


my_args = parseoArgs.parse_args()


####################################################################
##### Main                                                      ####
####################################################################

data = re.search(r"([\w_.-]+)\.json", my_args.file)
if(data != None):
	### Define a FEL Extractor, for example
	e = phyphy.Extractor(my_args.file)
	e.extract_csv(data.group(1)+".csv")  ## save to fel.csv

	### tab-delimited output, as fel.tsv
	e.extract_csv(data.group(1)+".tsv", delim = "\t")
