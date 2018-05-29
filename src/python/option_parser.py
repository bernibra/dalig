#!/usr/bin/python

###########################################

# system modules
from optparse import OptionParser
import sys

###########################################

def parse_cl_options():
	#usage = "usage: %prog [-r ROLE_FILENAME] [-s ROLE_FILENAME] NETWORK_FILE NETWORK_FILE"
	usage = "usage: %prog [OPTION] FIRST_NETWORK_FILE SECOND_NETWORK_FILE ALIGNMENT"
	parser = OptionParser(usage)

	parser.add_option("-k", "--degree",
					  action="store", dest="degree", type="int",
					  help="degree of alignment to conduct [default: %default]",
					  default=0,
					 )
	parser.add_option("-l", "--cost_function",
					  action="store", dest="cost_function", type="int",
					  help="Euclidean distance (0), Pearson's correlation coeficient (1) or Chi-squared test (2) [default: %default]",
					  default=2,
					 )
	parser.add_option("-n", "--normalization",
					  action="store", dest="normalization", type="int",
					  help="Remove all null nodes and normalize (0); the structure of aligned nodes and the neigbours (1); the structure of aligned nodes (2) [default: %default]",
					  default=0,
					 )
	parser.add_option("-a", "--alignment",
					  action="store_true", dest="palignment",
					  help="Print alignment (it might be wrong... I need to doublecheck) [default: %default]",
					  default=False,
					 )
	parser.add_option("-s", "--neighbours",
					  action="store_true", dest="neighbours",
					  help="Print number of neighbours of each species [default: %default]",
					  default=False,
					 )
	parser.add_option("-b", "--bipartite",
					  action="store_true", dest="bipartite",
					  help="consider the networks as bipartite [default: %default]",
					  default=False,
					 )

	parser.set_defaults(roles1=None, roles2=None)

	(options, args) = parser.parse_args()

	if len(args) != 3:
		parser.print_help()
		sys.exit()
	else:
		return (options, args)
