import os,sys
import argparse
import subprocess
import numpy as np
parser = argparse.ArgumentParser(description='funbiased')
subparsers = parser.add_subparsers(help='sub-command help')


parser_f2 = subparsers.add_parser('F2', help='compute F2')
parser_f2.add_argument('p1vcf', help='vcf file population 1')
parser_f2.add_argument('p2vcf', help='vcf file population 2')
parser_f2.add_argument('p1kin', help='kinship matrix population 1')
parser_f2.add_argument('p2kin', help='kinship matrix population 2')
parser_f2.set_defaults(mode='F2')


parser_f3 = subparsers.add_parser('F3', help='compute F3')
parser_f3.add_argument('p1vcf', help='vcf file population 1')
parser_f3.add_argument('p2vcf', help='vcf file population 2')
parser_f3.add_argument('p3vcf', help='vcf file population 3')
parser_f3.add_argument('p1kin', help='kinship matrix population 1')
parser_f3.add_argument('p2kin', help='kinship matrix population 2')
parser_f3.add_argument('p3kin', help='kinship matrix population 3')
parser_f3.set_defaults(mode='F3')


parser_f3n = subparsers.add_parser('F3norm', help='compute normalized F3')
parser_f3n.add_argument('p1vcf', help='vcf file population 1')
parser_f3n.add_argument('p2vcf', help='vcf file population 2')
parser_f3n.add_argument('p3vcf', help='vcf file population 3')
parser_f3n.add_argument('p1kin', help='kinship matrix population 1')
parser_f3n.add_argument('p2kin', help='kinship matrix population 2')
parser_f3n.add_argument('p3kin', help='kinship matrix population 3')
parser_f3n.set_defaults(mode='F3norm')


parser_f4n = subparsers.add_parser('F4norm', help='compute normalized F4')
parser_f4n.add_argument('p1vcf', help='vcf file population 1')
parser_f4n.add_argument('p2vcf', help='vcf file population 2')
parser_f4n.add_argument('p3vcf', help='vcf file population 3')
parser_f4n.add_argument('p4vcf', help='vcf file population 4')
parser_f4n.add_argument('p1kin', help='kinship matrix population 1')
parser_f4n.add_argument('p2kin', help='kinship matrix population 2')
parser_f4n.add_argument('p3kin', help='kinship matrix population 3')
parser_f4n.add_argument('p4kin', help='kinship matrix population 4')
parser_f4n.add_argument('popP', help='population P for normalization')
parser_f4n.set_defaults(mode='F4norm')


ins = parser.parse_args()
allins = vars(ins)


def readkin(kin):
	kincoeff=0.00
	with open(kin,'r') as kinmat:
		for eachind in kinmat:
			for eacheachkin in eachind.rstrip('\n').split():
				weight=1/float(len(eachind.split()))
				kincoeff+=weight*weight*float(eacheachkin)
	return(kincoeff)


def readvcf(vcf):
	plist=[]
	misslist=[]
	with open(vcf,'r') as popdata:
		for indv in popdata:
			if indv[0]!='#':
				indnum = len([x for x in indv.split() if '|' in x or '/' in x ])
				missing=[]
				if ".|." in indv or "./." in indv:
					for ind in range(len(indv.split())):
						if indv.split()[ind]=='.|.' or indv.split()[ind]=='./.':
							missing.append(len(indv.split())-ind)
					misslist.append(missing)
				else:
					misslist.append('No')
				if len(missing)==indnum:
					plist.append(np.nan)
                                if len(missing)!=indnum:
					pc0=0.00
					pc1=0.00
					for eachind in indv.split():
						if "|" in eachind or "/" in eachind:
							pc0+=eachind.count('0')
							pc1+=eachind.count('1')
					plist.append(pc0/(pc0+pc1))
					
	return(plist,misslist)


def miskin(kin,sites):
	kincoeff=0.00
	with open(kin,'r') as kinmat:
		kindat=kinmat.readlines()
		for eachind in range(len(kindat)):
			for eachindind in range(len(kindat[eachind].split())):
				if eachindind != sites[0]-len(kindat[eachind]) and eachind!= sites[0]-len(kindat[eachind]):
					weight=1/float(len(kindat[eachind].split()))
					#print kindat[eachind].split()
					kincoeff+=weight*weight*float(kindat[eachind].split()[eachindind])
	return(kincoeff)


if allins['mode'] == 'F2':

	p1 = allins['p1vcf']
	p2 = allins['p2vcf']
	k1 = allins['p1kin']
	k2 = allins['p2kin']
	kincoeffa=readkin(k1)
	kincoeffb=readkin(k2)
	p1list,p1mis=readvcf(p1)
	p2list,p2mis=readvcf(p2)
	countsite=0	
	ftot=0.0
	ftotunb=0.0
	for r in range(len(p1list)):
		if p1list[r] is np.nan or p2list[r] is np.nan:
			continue
		if p1mis[r]!='No':
			kincoeffa=miskin(k1,p1mis[r])
		countsite+=1
		a,b= p1list[r],p2list[r]
		biased,unb= (a-b)**2, ((a-b)**2)-(kincoeffa*((a*(1-a))/kincoeffa))-(kincoeffb*((b*(1-b))/kincoeffb))
		ftot+=biased
		ftotunb+=unb
	print ftot/countsite,ftotunb/countsite


if allins['mode'] == 'F3':

	p1 = allins['p1vcf']
	p2 = allins['p2vcf']
	p3 = allins['p3vcf']
	k1 = allins['p1kin']
	k2 = allins['p2kin']
	k3 = allins['p3kin']
	kincoeffa=readkin(k1)
	kincoeffb=readkin(k2)
	kincoeffc=readkin(k3)
	p1list,p1mis=readvcf(p1)
	p2list,p2mis=readvcf(p2)
	p3list,p3mis=readvcf(p3)
	countsite=0
	ftot=0.0
	ftotunb=0.0
	for r in range(len(p1list)):
		if p1list[r] is np.nan or p2list[r] is np.nan or p3list[r] is np.nan:
			continue
		if p1mis[r]!='No':
			kincoeffa=miskin(k1,p1mis[r])
		countsite+=1
		a,b,c= p1list[r],p2list[r],p3list[r]
		biased,unb= (a-b)*(a-c), ((a-b)*(a-c))-(kincoeffa*((a*(1-a))/kincoeffa))
		ftot+=biased
		ftotunb+=unb
	print ftot/countsite,ftotunb/countsite


if allins['mode'] == 'F3norm':

	p1 = allins['p1vcf']
	p2 = allins['p2vcf']
	p3 = allins['p3vcf']
	k1 = allins['p1kin']
	k2 = allins['p2kin']
	k3 = allins['p3kin']
	kincoeffa=readkin(k1)
	kincoeffb=readkin(k2)
	kincoeffc=readkin(k3)
	p1list,p1mis=readvcf(p1)
	p2list,p2mis=readvcf(p2)
	p3list,p3mis=readvcf(p3)
	countsite=0
	ftot=0.0
	ftotunb=0.0
	fden=0.0
	fdenunb=0.0
	for r in range(len(p1list)):
		if p1list[r] is np.nan or p2list[r] is np.nan or p3list[r] is np.nan:
			continue
		if p1mis[r]!='No':
			kincoeffa=miskin(k1,p1mis[r])
		countsite+=1
		a,b,c= p1list[r],p2list[r],p3list[r]
		biased,unb,denom,unbdenom= (a-b)*(a-c), ((a-b)*(a-c))-(kincoeffa*((a*(1-a))/kincoeffa)),a*(1-a),((a*(1-a))/kincoeffa)
		ftot+=biased
		ftotunb+=unb
		fden+=denom
		fdenunb+=unbdenom
	print (ftot/countsite)/(2*(fden/countsite)),(ftotunb/countsite)/(2*(fdenunb/countsite))


if allins['mode'] == 'F4norm':

	p1 = allins['p1vcf']
	p2 = allins['p2vcf']
	p3 = allins['p3vcf']
	p4 = allins['p4vcf']
	k1 = allins['p1kin']
	k2 = allins['p2kin']
	k3 = allins['p3kin']
	k4 = allins['p4kin']
	normpop = allins['popP']
	p1list,p1mis=readvcf(p1)
	p2list,p2mis=readvcf(p2)
	p3list,p3mis=readvcf(p3)
	p4list,p4mis=readvcf(p4)
	if normpop=='A':
		plist=p1list
		kincoeffp=readkin(k1)
	elif normpop=='B':
		plist=p2list
		kincoeffp=readkin(k2)
	elif normpop=='C':
		plist=p3list
		kincoeffp=readkin(k3)
	elif normpop=='D':
		plist=p4list
		kincoeffp=readkin(k4)
	countsite=0
	ftot=0.0
	ftotunb=0.0
	fden=0.0
	fdenunb=0.0
	for r in range(len(p1list)):
		if p1list[r] is np.nan or p2list[r] is np.nan or p3list[r] is np.nan or p4list[r] is np.nan:
			continue
		if p1mis[r]!='No':
			kincoeffp=miskin(k1,p1mis[r])
		countsite+=1
		a,b,c,d,p= p1list[r],p2list[r],p3list[r],p4list[r],plist[r]
		biased,denom,unbdenom= (a-b)*(c-d),p*(1-p),((p*(1-p))/kincoeffp)
		ftot+=biased
		fden+=denom
		fdenunb+=unbdenom
	print (ftot/countsite)/(fden/countsite),(ftot/countsite)/(fdenunb/countsite)
