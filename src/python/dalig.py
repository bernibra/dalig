
import os
import random
import subprocess
import sys
import tempfile
import copy

from pymfinder import motif_roles,print_roles
from option_parser import parse_cl_options

def read_network(filename, nodes, normalization=0, red=False, a_map=None):
	network = []
	predators=set([])
	prey=set([])
	unipartite=False
	new_nodes=set([])
	inFile = open(filename,'r')
	for line in inFile:
		interaction=line.strip().split()
		if interaction[0] not in nodes or interaction[1] not in nodes:
			sys.stderr.write("The nodes in the alignment file don't agree with the network files.\n")
			sys.exit()
		if unipartite==False and interaction[0] not in prey and interaction[1] not in predators:
			predators.add(interaction[0])
			prey.add(interaction[1])
		else:
			unipartite=True
		if red:
			if (a_map[interaction[0]]=="NULL" or a_map[interaction[1]]=="NULL") and normalization==2:
				continue
			if (a_map[interaction[0]]=="NULL" and a_map[interaction[1]]=="NULL") and normalization==1:
				continue
		network.append(interaction)
		new_nodes.add(interaction[0])
		new_nodes.add(interaction[1])

	inFile.close()

	##########THIS IS SPECIFICALLY THE SUPERSHIT PART.... BUT IT WORKS!
	if red:
		size=len(interaction)
		nodes=nodes.difference(new_nodes)
		for k in nodes:
			if a_map[k]!="NULL":
				if size==3:
					interact=[str(k),str(k),str(1)]
				elif size==2:
					interact=[k,k]
				else:
					sys.stderr.write("There is something weird in one of your networks.\n")
				network.append(interact)
				new_nodes.add(k)
	####################################################################


	return network, unipartite, new_nodes

def read_alignment(alignment_file):
	inFile = open(alignment_file, 'r')
	n1=set([])
	n2=set([])
	alignment=[]
	for line in inFile:
		pairs=line.strip().split()
		if len(pairs)!=2:
			sys.stderr.write("There is something peculiar in the alignment file.\n")
			sys.exit()
		if (pairs[0]!="NULL" and pairs[0] in n1) or (pairs[1]!="NULL" and pairs[1] in n2):
			sys.stderr.write("There are nodes that appear more than once in the alignment file.\n")
			sys.exit()
		n1.add(pairs[0])
		n2.add(pairs[1])
		alignment.append(pairs)
	inFile.close()
	return alignment, n1, n2

def check_alignment(alignment):
	n1=set([])
	n2=set([])
	newalignment=[]
	for line in alignment:
		if len(line)!=2:
			sys.stderr.write("There is something peculiar in the alignment file.\n")
			sys.exit()
		if (str(line[0])!="NULL" and str(line[0]) in n1) or (str(line[1])!="NULL" and str(line[1]) in n2):
			sys.stderr.write("There are nodes that appear more than once in the alignment file.\n")
			sys.exit()
		n1.add(str(line[0]))
		n2.add(str(line[1]))
		newalignment.append([str(line[0]), str(line[1])])

	return newalignment, n1, n2

def rewrite_alignment(alignment, n1, n2, net1, net2):
	newalignment=[]
	for line in alignment:
		if (line[0] in n1 or line[0]=="NULL") and (line[1] in n2 or line[1]=="NULL"):
			newalignment.append([str(line[0]), str(line[1])])
	return newalignment


def dalig_input(net1, net2, net1_roles, net2_roles, alignment):
	all_roles = net1_roles[net1_roles.keys()[0]].keys()

	input = []

	# write out the links for network 1
	input.append('\n'.join([' '.join(link) for link in net1]))
	input.append('###')

	# write out the roles for network 1
	input.append('\n'.join([' '.join(map(str,[n]+[net1_roles[n][r] for r in all_roles])) for n in net1_roles]))
	input.append('///')

	# write out the links for network 2
	input.append('\n'.join([' '.join(link) for link in net2]))
	input.append('###')

	# write out the roles for network 2
	input.append('\n'.join([' '.join(map(str,[n]+[net2_roles[n][r] for r in all_roles])) for n in net2_roles]))

	# write out the alignment
	input.append('&&&')
	input.append('\n'.join([' '.join(link) for link in alignment]))

	return '\n'.join(input)

def dalig(options, args, r_alignment=True):

	#if options.degree==0 and options.normalization>0:
	#	sys.stderr.write("Not implemented yet the combination degree==0 and normalization>0.\n")

	if r_alignment:
		alignment, n1, n2 = read_alignment(args[2])
	else:
		alignment, n1, n2 = check_alignment(args[2])

	try:
		n1.remove("NULL")
	except KeyError:
		pass
	try:
		n2.remove("NULL")
	except KeyError:
		pass

	##########THIS IS KIND OF SUPERSHIT.... BUT IT WORKS!
	if options.normalization>0:
		if len(n1)<=len(n2):
			reduce_n1=False
			#a_map={key:value for (value, key) in alignment}
			a_map={}
			for (value, key) in alignment:
				a_map[key]=value
		else:
			reduce_n1=True
			#a_map={key:value for (key, value) in alignment}
			a_map={}
			for (key, value) in alignment:
				a_map[key]=value

		net1, net1type, n1 = read_network(args[0], n1, options.normalization, red=reduce_n1, a_map=a_map)
		net2, net2type, n2 = read_network(args[1], n2, options.normalization, red=reduce_n1==False, a_map=a_map)

		alignment = rewrite_alignment(alignment, n1, n2, net1, net2)
	else:
		net1, net1type, n1 = read_network(args[0], n1)
		net2, net2type, n2 = read_network(args[1], n2)
	##########################################################

	if alignment==[]:
		return 'optimal = [ (NULL,NULL:0,0:1)] \noverlap = (0,0)\nEnergy = nan\nNormalized energy = 1'

	#Are the networks unipartite or bipartite?
	if options.bipartite:
		if net1type==False and net2type==False:
			unipartite=net1type
		else:
			print "You are comparing a unipartite network with a bipartite one. Both will be considered as unipartite"
			unipartite=True
	else:
		unipartite=True

	if unipartite:
		net1_roles = motif_roles(copy.deepcopy(net1),motifsize=2,)
		net1_roles2 = motif_roles(copy.deepcopy(net1),motifsize=3,)
		for i in net1_roles:
			net1_roles[i].update(net1_roles2[i])
	else:
		net1_roles = motif_roles(copy.deepcopy(net1),motifsize=2, networktype = "bipartite",)
		for k in range(3,5):
			net1_roles2 = motif_roles(copy.deepcopy(net1),motifsize=k, networktype = "bipartite",)
			for i in net1_roles:
				net1_roles[i].update(net1_roles2[i])


	if unipartite:
		net2_roles = motif_roles(copy.deepcopy(net2),motifsize=2,)
		net2_roles2 = motif_roles(copy.deepcopy(net2),motifsize=3,)
		for i in net2_roles:
			net2_roles[i].update(net2_roles2[i])
	else:
		net2_roles = motif_roles(copy.deepcopy(net2),motifsize=2, networktype = "bipartite",)
		for k in range(3,5):
			net2_roles2 = motif_roles(copy.deepcopy(net2),motifsize=k, networktype = "bipartite",)
			for i in net2_roles:
				net2_roles[i].update(net2_roles2[i])

	dalig_in = dalig_input(net1, net2, net1_roles, net2_roles, alignment)
	# get a random seed
	rnd_seed = random.randint(0,sys.maxint)

	# do we want to start from a random initial alignment?
	if options.palignment:
		aflag = "-a"
	else:
		aflag = ""

	# do we want to start from a random initial alignment?
	if options.neighbours:
		sflag = "-s"
	else:
		sflag = ""

	# call the dalig alignment code
	command = "GSL_RNG_SEED=%s PATH=$PATH:/home/bbr36/.local/bin dalig.x -k %s -l %s -n %s %s %s" % (rnd_seed,
								options.degree,
								options.cost_function,
								options.normalization,
								aflag,
								sflag)

	#dalig_out = tempfile.TemporaryFile()
	process = subprocess.Popen(command,
					bufsize=0,
					stdin=subprocess.PIPE,
					stdout=subprocess.PIPE,
					stderr=subprocess.PIPE,
					shell=True,
					)

	# write the network and role data to dalig
	process.stdin.write(dalig_in)
	#process.stdin.close()

	output = process.communicate()

	# print out and store the dalig stdout line by line as it comes
	#output=""
	#for line in iter(poutput[0].readline, ''):
	#	output+=line

	#process.wait()
		
	# get rid of the GSL seed info and any empty lines from stderr
	#dalig_stderr = [i for i in process.stderr.readlines() if "GSL_RNG_SEED" not in i and i != '']
	#if dalig_stderr:
	#	sys.stderr.write(''.join(dalig_stderr))

	return output[0]

def pydalig(network1=None, network2=None,
		alignment=None, normalization=0,
		cost_function=1, bipartite=False,
		degree=1, palignment=False, 
		neighbours=False):

	if type(network1)!=str or type(network2)!=str:
		sys.stderr.write("There is something peculiar in your input files.\n")
		sys.exit()

	class optin:
		def __init__(self, roles1 = None, cost_function =cost_function,
						bipartite = bipartite, normalization=normalization,
						degree = degree, roles2 = None,
						palignment=palignment, neighbours=neighbours):
			self.roles1 = roles1
			self.cost_function = cost_function
			self.normalization = normalization
			self.bipartite = bipartite
			self.neighbours = neighbours
			self.palignment = palignment
			self.degree = degree
			self.roles2 = roles2

	options=optin()

	if type(alignment)==str:
		r_alignment=True
	elif type(alignment)==type([(1,2),(3,4)]):
		r_alignment=False
	else:
		sys.stderr.write("There is something peculiar in your input files.\n")
		sys.exit()

	args=(network1, network2, alignment)

	output=dalig(options, args, r_alignment)
	return output


########################################
########################################
########################################
########################################

def main():
	options, args = parse_cl_options()
	output=dalig(options, args)
	print(output)

if __name__ == "__main__":
	main()
