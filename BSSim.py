#!/usr/bin/env python2.6
#!/home/local/bin/python2.6
#-*- coding:utf-8 -*-

'''
Description:

BSSim: bisulfite sequencing simulator for next-generation sequencing

Copyright (c) 2012, Luke You <xigyou@gmail.com>.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License.

@Compile Date: Friday, June 29, 2012
@version: 1.2
@author:  Luke You
@Author Affiliations: Institute of Genomic Medicine, Wenzhou Medical College, Wenzhou 325035, China
@contact: xigyou@gmail.com
'''

import os
import sys
import getopt
import string
import time
from collections import defaultdict
import pyfasta
import re
import math
import random
import numpy
import bisect
import cPickle
from multiprocessing import *


def read_fa(faFile):
	#Read user's input fata file.
	from pyfasta import Fasta
	fa=Fasta(faFile)
	return fa
	del fa

def random_snp(fa,p1,p2,DS,Polyploid):
	if Polyploid == 1:
		p2=0
	p0=1-p1-p2					#don't have SNP
	tmp="".join(fa)
	tmp=tmp.upper()
	ref=['A']*len(tmp)
	for i in range(0,len(tmp)):
		ref[i]=tmp[i]
		ref[i]=re.sub(r'N' , random.choice('ATCG'),tmp[i])									#use another random deoxyribonucleotide replace 'N'
	del tmp
	if not DS:
		for i in range(0,len(ref)):
			if random.random() > p0:									#have SNP
				if random.random() <= p1/(p1+p2):					#homozygote
					if ref[i] == 'A':
						if random.random() <= 4/6.0:
							ref[i]='G'
						else:
							if random.random() <= 0.5:
								ref[i]='C'
							else:
								ref[i]='T'
					if ref[i] == 'G':
						if random.random() <= 4/6.0:
							ref[i]='A'
						else:
							if random.random() <= 0.5:
								ref[i]='C'
							else:
								ref[i]='T'
					if ref[i] == 'C':
						if random.random() <= 4/6.0:
							ref[i]='T'
						else:
							if random.random() <= 0.5:
								ref[i]='A'
							else:
								ref[i]='G'
					if ref[i] == 'T':
						if random.random() <= 4/6.0:
							ref[i]='C'
						else:
							if random.random() <= 0.5:
								ref[i]='A'
							else:
								ref[i]='G'
				else:												#heterozygote
					# print "333333333\theterozygote\t%d\t%s" % (i,ref[i])
					ref[i]='N'
	ref="".join(ref)
	return ref
	del ref

def random_snp_diploid(fa,p1,p2,DS):
	p0=1-p1-p2					#don't have SNP
	tmp="".join(fa)
	tmp=tmp.upper()
	ref=['A']*len(tmp)
	for i in range(0,len(tmp)):
		ref[i]=tmp[i]
		ref[i]=re.sub(r'N' , random.choice('ATCG'),tmp[i])									#use another random deoxyribonucleotide replace 'N'
	del tmp
	ref="".join(ref)
	if not DS:
		ref_A=[0]*len(ref)
		ref_B=[0]*len(ref)
		for i in range(0,len(ref)):
			if random.random() <= p0:							#don't have SNP
				ref_A[i]=ref[i]
				ref_B[i]=ref[i]
				# print "111111111\tdon't have SNP\tref[i]"
				continue
			elif random.random() <= p1/(p1+p2):					#homozygote
				# print "222222222\thomozygote\t%d\t%s" % (i,ref[i])
				if ref[i] == 'A':
					if random.random() <= 4/6.0:
						ref_A[i]=ref_B[i]='G'
					else:
						if random.random() <= 0.5:
							ref_A[i]=ref_B[i]='C'
						else:
							ref_A[i]=ref_B[i]='T'
				if ref[i] == 'G':
					if random.random() <= 4/6.0:
						ref_A[i]=ref_B[i]='A'
					else:
						if random.random() <= 0.5:
							ref_A[i]=ref_B[i]='C'
						else:
							ref_A[i]=ref_B[i]='T'
				if ref[i] == 'C':
					if random.random() <= 4/6.0:
						ref_A[i]=ref_B[i]='T'
					else:
						if random.random() <= 0.5:
							ref_A[i]=ref_B[i]='A'
						else:
							ref_A[i]=ref_B[i]='G'
				if ref[i] == 'T':
					if random.random() <= 4/6.0:
						ref_A[i]=ref_B[i]='C'
					else:
						if random.random() <= 0.5:
							ref_A[i]=ref_B[i]='A'
						else:
							ref_A[i]=ref_B[i]='G'
			else:												#heterozygote
				# print "333333333\theterozygote\t%d\t%s" % (i,ref[i])
				if ref[i] == 'A':
					if random.random() <= 4/6.0:
						ref_A[i]='G'
						ref_B[i]='A'
					else:
						if random.random() <= 0.5:
							ref_A[i]='C'
							ref_B[i]='A'
						else:
							ref_A[i]='T'
							ref_B[i]='A'
				if ref[i] == 'G':
					if random.random() <= 4/6.0:
						ref_A[i]='A'
						ref_B[i]='G'
					else:
						if random.random() <= 0.5:
							ref_A[i]='C'
							ref_B[i]='G'
						else:
							ref_A[i]='T'
							ref_B[i]='G'
				if ref[i] == 'C':
					if random.random() <= 4/6.0:
						ref_A[i]='T'
						ref_B[i]='C'
					else:
						if random.random() <= 0.5:
							ref_A[i]='A'
							ref_B[i]='C'
						else:
							ref_A[i]='G'
							ref_B[i]='C'
				if ref[i] == 'T':
					if random.random() <= 4/6.0:
						ref_A[i]='C'
						ref_B[i]='T'
					else:
						if random.random() <= 0.5:
							ref_A[i]='A'
							ref_B[i]='T'
						else:
							ref_A[i]='G'
							ref_B[i]='T'
		ref_A="".join(ref_A)
		ref_B="".join(ref_B)
	else:
		ref_A=ref_B=ref
	return ref_A,ref_B
	del ref_A,ref_B

def read_dbsnp(dbsnpfile):
	#read dbsnpfile and create database.
	f=open(dbsnpfile,'r')
	max_len=0
	dbsnp=defaultdict(int)
	while (f):
		line=f.readline().rstrip("\n")
		if len(line) == 0: # Zero length indicates EOF
			break
		line=line.split('\t')
		chrom = line[0]
		k=str(line[1])
		if line[2] == "-":
			dbsnp[chrom,k]=line[4]+"\t"+line[3]+"\t"+line[6]+"\t"+line[5]
		else:
			dbsnp[chrom,k]=line[3]+"\t"+line[4]+"\t"+line[5]+"\t"+line[6]
	f.close() # close the file
	return dbsnp
	del dbsnp


def methyl(ref,conversion_rate,CG_conversion_rate,CHG_conversion_rate,CHH_conversion_rate,mC_rate,mCG_rate,mCHG_rate,mCHH_rate,CG_beta_distribution,mCG_mu,mCG_sigma,CHG_beta_distribution,mCHG_mu,mCHG_sigma,CHH_beta_distribution,mCHH_mu,mCHH_sigma):
	ref_rate=[0]*len(ref)
	ref_pattern=[""]*len(ref)
	for j in range(0,len(ref)):
		if ref[j] == 'C' :
			if (len(ref)-j)<3:
				ref_rate[j]=conversion_rate
				continue
			if ref[j+1] == 'G':															#CG
				ref_pattern[j]="CG"
				if random.random() <= mCG_rate:
					ref_rate[j]=(1-ref_methyl_rate('CG',CG_beta_distribution,mCG_mu,mCG_sigma))*CG_conversion_rate
				else:
					ref_rate[j]=CG_conversion_rate
			else:
				if ref[j+2] == 'G':														#CHG
					ref_pattern[j]="CHG"
					if random.random() <= mCHG_rate:
						ref_rate[j]=(1-ref_methyl_rate('CHG',CHG_beta_distribution,mCHG_mu,mCHG_sigma))*CHG_conversion_rate
					else:
						ref_rate[j]=CHG_conversion_rate
				else:
					ref_pattern[j]="CHH"
					if random.random() <= mCHH_rate:				#CHH
						ref_rate[j]=(1-ref_methyl_rate('CHH',CHH_beta_distribution,mCHH_mu,mCHH_sigma))*CHH_conversion_rate
					else:
						ref_rate[j]=CHH_conversion_rate
		elif ref[j] == 'G':
			if j<3:
				ref_rate[j]=conversion_rate
				continue
			if ref[j-1] == 'C':															#CG
				ref_pattern[j]="CG"
				if random.random() <= mCG_rate:
					ref_rate[j]=(1-ref_methyl_rate('CG',CG_beta_distribution,mCG_mu,mCG_sigma))*CG_conversion_rate
				else:
					ref_rate[j]=CG_conversion_rate
			else:
				if ref[j-2] == 'C' :													#CHG
					ref_pattern[j]="CHG"
					if random.random() <= mCHG_rate:
						ref_rate[j]=(1-ref_methyl_rate('CHG',CHG_beta_distribution,mCHG_mu,mCHG_sigma))*CHG_conversion_rate
					else:
						ref_rate[j]=CHG_conversion_rate
				else:
					ref_pattern[j]="CHH"
					if random.random() <= mCHH_rate:				#CHH
						ref_rate[j]=(1-ref_methyl_rate('CHH',CHH_beta_distribution,mCHH_mu,mCHH_sigma))*CHH_conversion_rate
					else:
						ref_rate[j]=CHH_conversion_rate
		else:
			ref_rate[j]=0
			ref_pattern[j]=""
		if (ref[j] == 'A' or ref[j] == 'T') and ref_rate[j] > 0:
			print "%d\t%s\t%s" %(j,ref[j],ref_rate[j])
	return ref_rate,ref_pattern
	del ref_rate,ref_pattern


def methyl_for_random_SNP_of_diploid(ref,ref2,conversion_rate,CG_conversion_rate,CHG_conversion_rate,CHH_conversion_rate,mC_rate,mCG_rate,mCHG_rate,mCHH_rate,CG_beta_distribution,mCG_mu,mCG_sigma,CHG_beta_distribution,mCHG_mu,mCHG_sigma,CHH_beta_distribution,mCHH_mu,mCHH_sigma):
	ref_rate=[0]*len(ref)
	ref_pattern=[""]*len(ref)
	for j in range(0,len(ref)):
		if ref[j] == 'C' :
			if (len(ref)-j)<3:
				ref_rate[j]=conversion_rate
				continue
			if ref[j+1] == 'G':															#CG
				ref_pattern[j]="CG"
				if random.random() <= mCG_rate*mCG_mu**(1/2.0):
					ref_rate[j]=(1-ref_methyl_rate('CG',CG_beta_distribution,mCG_mu**(1/2.0),mCG_sigma))*CG_conversion_rate
				else:
					ref_rate[j]=CG_conversion_rate
			else:
				if ref[j+2] == 'G':														#CHG
					ref_pattern[j]="CHG"
					if random.random() <= mCHG_rate*mCHG_mu**(1/2.0):
						ref_rate[j]=(1-ref_methyl_rate('CHG',CHG_beta_distribution,mCHG_mu**(1/2.0),mCHG_sigma))*CHG_conversion_rate
					else:
						ref_rate[j]=CHG_conversion_rate
				else:
					ref_pattern[j]="CHH"
					if random.random() <= mCHH_rate*mCHH_mu**(1/2.0):				#CHH
						ref_rate[j]=(1-ref_methyl_rate('CHH',CHH_beta_distribution,mCHH_mu**(1/2.0),mCHH_sigma))*CHH_conversion_rate
					else:
						ref_rate[j]=CHH_conversion_rate
		elif ref[j] == 'G':
			if j<3:
				ref_rate[j]=conversion_rate
				continue
			if ref[j-1] == 'C':															#CG
				ref_pattern[j]="CG"
				if random.random() <= mCG_rate*mCG_mu**(1/2.0):
					ref_rate[j]=(1-ref_methyl_rate('CG',CG_beta_distribution,mCG_mu**(1/2.0),mCG_sigma))*CG_conversion_rate
				else:
					ref_rate[j]=CG_conversion_rate
			else:
				if ref[j-2] == 'C' :													#CHG
					ref_pattern[j]="CHG"
					if random.random() <= mCHG_rate*mCHG_mu**(1/2.0):
						ref_rate[j]=(1-ref_methyl_rate('CHG',CHG_beta_distribution,mCHG_mu**(1/2.0),mCHG_sigma))*CHG_conversion_rate
					else:
						ref_rate[j]=CHG_conversion_rate
				else:
					ref_pattern[j]="CHH"
					if random.random() <= mCHH_rate*mCHH_mu**(1/2.0):				#CHH
						ref_rate[j]=(1-ref_methyl_rate('CHH',CHH_beta_distribution,mCHH_mu**(1/2.0),mCHH_sigma))*CHH_conversion_rate
					else:
						ref_rate[j]=CHH_conversion_rate
		else:
			ref_rate[j]=0
			ref_pattern[j]=""
		if (ref[j] == 'A' or ref[j] == 'T') and ref_rate[j] > 0:
			print "%d\t%s\t%s" %(j,ref[j],ref_rate[j])
	ref_rate2=[0]*len(ref2)
	ref_pattern2=[""]*len(ref2)
	for j in range(0,len(ref2)):
		if ref2[j] == 'C' :
			if (len(ref2)-j)<3:
				ref_rate2[j]=conversion_rate
				continue
			if ref2[j+1] == 'G':															#CG
				ref_pattern2[j]="CG"
				if (ref_rate[j] < CG_conversion_rate or ( ref_rate[j] == CG_conversion_rate and random.random() <= ((mCG_rate*(1-mCG_mu**(1/2.0)))/(1-mCG_rate*mCG_mu**(1/2.0))))) and random.random() <= mCG_mu**(1/2.0):
					ref_rate2[j]=(1-ref_methyl_rate('CG',CG_beta_distribution,mCG_mu**(1/2.0),mCG_sigma))*CG_conversion_rate
				else:
					ref_rate2[j]=CG_conversion_rate
			else:
				if ref2[j+2] == 'G':														#CHG
					ref_pattern2[j]="CHG"
					if (ref_rate[j] < CHG_conversion_rate or ( ref_rate[j] == CHG_conversion_rate and random.random() <= ((mCHG_rate*(1-mCHG_mu**(1/2.0)))/(1-mCHG_rate*mCHG_mu**(1/2.0))))) and random.random() <= mCHG_mu**(1/2.0):
						ref_rate2[j]=(1-ref_methyl_rate('CHG',CHG_beta_distribution,mCHG_mu**(1/2.0),mCHG_sigma))*CHG_conversion_rate
					else:
						ref_rate2[j]=CHG_conversion_rate
				else:
					ref_pattern2[j]="CHH"
					if (ref_rate[j] < CHH_conversion_rate or ( ref_rate[j] == CHH_conversion_rate and random.random() <= ((mCHH_rate*(1-mCHH_mu**(1/2.0)))/(1-mCHH_rate*mCHH_mu**(1/2.0))))) and random.random() <= mCHH_mu**(1/2.0):
						ref_rate2[j]=(1-ref_methyl_rate('CHH',CHH_beta_distribution,mCHH_mu**(1/2.0),mCHH_sigma))*CHH_conversion_rate
					else:
						ref_rate2[j]=CHH_conversion_rate
		elif ref2[j] == 'G':
			if j<3:
				ref_rate2[j]=conversion_rate
				continue
			if ref2[j-1] == 'C':															#CG
				ref_pattern2[j]="CG"
				if (ref_rate[j] < CG_conversion_rate or ( ref_rate[j] == CG_conversion_rate and random.random() <= ((mCG_rate*(1-mCG_mu**(1/2.0)))/(1-mCG_rate*mCG_mu**(1/2.0))))) and random.random() <= mCG_mu**(1/2.0):
					ref_rate2[j]=(1-ref_methyl_rate('CG',CG_beta_distribution,mCG_mu**(1/2.0),mCG_sigma))*CG_conversion_rate
				else:
					ref_rate2[j]=CG_conversion_rate
			else:
				if ref2[j-2] == 'C' :													#CHG
					ref_pattern2[j]="CHG"
					if (ref_rate[j] < CHG_conversion_rate or ( ref_rate[j] == CHG_conversion_rate and random.random() <= ((mCHG_rate*(1-mCHG_mu**(1/2.0)))/(1-mCHG_rate*mCHG_mu**(1/2.0))))) and random.random() <= mCHG_mu**(1/2.0):
						ref_rate2[j]=(1-ref_methyl_rate('CHG',CHG_beta_distribution,mCHG_mu**(1/2.0),mCHG_sigma))*CHG_conversion_rate
					else:
						ref_rate2[j]=CHG_conversion_rate
				else:
					ref_pattern2[j]="CHH"
					if (ref_rate[j] < CHH_conversion_rate or ( ref_rate[j] == CHH_conversion_rate and random.random() <= ((mCHH_rate*(1-mCHH_mu**(1/2.0)))/(1-mCHH_rate*mCHH_mu**(1/2.0))))) and random.random() <= mCHH_mu**(1/2.0):
						ref_rate2[j]=(1-ref_methyl_rate('CHH',CHH_beta_distribution,mCHH_mu**(1/2.0),mCHH_sigma))*CHH_conversion_rate
					else:
						ref_rate2[j]=CHH_conversion_rate
		else:
			ref_rate2[j]=0
			ref_pattern2[j]=""
		if (ref2[j] == 'A' or ref2[j] == 'T') and ref_rate2[j] > 0:
			print "%d\t%s\t%s" %(j,ref2[j],ref_rate2[j])
	return ref_rate,ref_pattern,ref_rate2,ref_pattern2
	del ref_rate,ref_pattern,ref_rate2,ref_pattern2


def ref_methyl_rate(mC,beta_distribution,mu,sigma):
	i=1
	if mC == 'CG':
		while(i):
			if beta_distribution:
				k=bata_fun(mu,sigma)
			else:
				k=random.gauss(1.696,0.3761)
			if k<=1 and k>0 :
				return k
				i=i-1
	elif mC == 'CHG':
		while(i):
			if beta_distribution:
				k=bata_fun(mu,sigma)
			else:
				k=random.gauss(-0.3255,0.1995)
			if k<=1 and k>0 :
				return k
				i=i-1
	elif mC == 'CHH':
		while(i):
			if beta_distribution:
				k=bata_fun(mu,sigma)
				if k<=1 and k>0 :
					return k
					i=i-1
			else:
				if random.random() <= 0.9091:
					k=random.gauss(-0.3994,0.2396)
					if k<=0.75 and k>0 :
						return k
						i=i-1
					else:
						k=random.gauss(2.669,0.6518)
						if k<=1 and k>0.75 :
							return k
							i=i-1


def create_Reads(PE,directional,ref_start,Ref,rate,start,reads_length,max_err):
	reads1T=[0]*reads_length
	reads2T=[0]*reads_length
	reads1B=[0]*reads_length
	reads2B=[0]*reads_length
	strand1T='+.W'								#OT
	start1T=start+1
	end1T=start+reads_length
	strand2T='-.W'								#CTOT
	start2T=start+len(Ref)-reads_length+1
	end2T=start+len(Ref)
	strand2B='+.C'								#CTOB
	start2B=start+1
	end2B=start+reads_length
	strand1B='-.C'								#OB
	start1B=start+len(Ref)-reads_length+1
	end1B=start+len(Ref)
	if max_err == 0:
		error=0
	else:
		error=random.randrange(0,max_err)
	for i in range(0,reads_length):
		if (reads_length-i)<=error:
			reads1T[i]=random.choice('ATCG')
		if rate[i] == 0:
			reads1T[i]=Ref[i]
		else:
			if Ref[i] == 'C':
				if random.random() <= rate[i]:
					reads1T[i]='T'
				else:
					reads1T[i]='C'
			elif Ref[i] == 'G':
				reads1T[i]='G'
			else:
				print "11111\t%d\t%s\t%s" % (i,Ref[i],rate[i])
	if PE or (not PE and not directional):
		for i in range((len(Ref)-1),(len(Ref)-reads_length-1),-1):
			if (i-(len(Ref)-reads_length-1))<=error:
				reads1B[len(Ref)-1-i]=random.choice('ATCG')
			if rate[i] == 0:
				reads1B[len(Ref)-1-i]=reverse_base(Ref[i])
			else:
				if Ref[i] == 'G':
					if random.random() <= rate[i]:
						reads1B[len(Ref)-1-i]='T'
					else:
						reads1B[len(Ref)-1-i]='C'
				elif Ref[i] == 'C':
					reads1B[len(Ref)-1-i]='G'
				else:
					print "44444\t%d\t%s\t%s" % (i,Ref[i],rate[i])
	if directional :
		reads1="".join(reads1T)
		strand1=strand1T
		start1=start1T
		end1=end1T
		if PE :
			reads2="".join(reads1B)
			strand2=strand1B
			start2=start1B
			end2=end1B
	else:
		if random.random() <= 0.5:
			reads1="".join(reads1T)
			strand1=strand1T
			start1=start1T
			end1=end1T
			if PE :
				for i in range((len(Ref)-1),(len(Ref)-reads_length-1),-1):
					if (i-(len(Ref)-reads_length-1))<=error:
						reads2T[len(Ref)-1-i]=random.choice('ATCG')
					if rate[i] == 0:
						reads2T[len(Ref)-1-i]=reverse_base(Ref[i])
					else:
						if Ref[i] == 'C':
							if random.random() <= rate[i]:
								reads2T[len(Ref)-1-i]='A'
							else:
								reads2T[len(Ref)-1-i]='G'
						elif Ref[i] == 'G':
							reads2T[len(Ref)-1-i]='C'
						else:
							print "22222\t%d\t%s\t%s" % (i,Ref[i],rate[i])
				reads2="".join(reads2T)
				strand2=strand2T
				start2=start2T
				end2=end2T
		else:
			reads1="".join(reads1B)
			strand1=strand1B
			start1=start1B
			end1=end1B
			if PE :
				for i in range(0,reads_length):
					if (reads_length-i)<=error:
						reads2B[i]=random.choice('ATCG')
					if rate[i] == 0:
						reads2B[i]=Ref[i]
					else:
						if Ref[i] == 'G':
							if random.random() <= rate[i]:
								reads2B[i]='A'
							else:
								reads2B[i]='G'
						elif Ref[i] == 'C':
							reads2B[i]='C'
						else:
							print "333333\t%d\t%s\t%s" % (i,Ref[i],rate[i])
				reads2="".join(reads2B)
				strand2=strand2B
				start2=start2B
				end2=end2B
	# print "%d\t%d" % (start1 , start2)
	if PE :
		return reads1,strand1,start1+ref_start,end1+ref_start,reads2,strand2,start2+ref_start,end2+ref_start
	else:
		return reads1,strand1,start1+ref_start,end1+ref_start


def reverse_base(orin):
	if orin == 'A':
		return 'T'
	if orin == 'T':
		return 'A'
	if orin == 'G':
		return 'C'
	if orin == 'C':
		return 'G'


def bata_fun(mu,sigma):
	k=mu*(1-mu)/sigma-1
	alpha=mu*k
	beta=(1-mu)*k
	# print "mu%f\tsigma%f\talpha%f\tbeta%f\n" % (mu,sigma,alpha,beta)
	return random.betavariate(alpha,beta)


def dynamic_quality(reads,quality,random_sequencing_errors,front_point,front_point_sigma,front_qual,front_qual_sigma,end_point,end_point_sigma,end_qual,end_qual_sigma):
	qual=[quality]*len(reads)
	reads_out=[""]*len(reads)
	front=len(reads)*random.normalvariate(front_point,front_point_sigma)
	end=len(reads)*random.normalvariate(end_point,end_point_sigma)
	for i in range(0,len(reads)):
		if i<front:
			qual[i]=bata_fun(quality-front_qual,front_qual_sigma)
			if random_sequencing_errors and random.random<sequencing_error_posibility(qual[i]):
				reads_out[i]=random.choice('ATCG')
			else:
				reads_out[i]=reads[i]
		if front<=i<=end:
			if random_sequencing_errors and random.random<sequencing_error_posibility(qual[i]):
				reads_out[i]=random.choice('ATCG')
			else:
				reads_out[i]=reads[i]
		if i>end:
			qual[i]=bata_fun(quality-end_qual,end_qual_sigma)
			if random_sequencing_errors and random.random<sequencing_error_posibility(qual[i]):
				reads_out[i]=random.choice('ATCG')
			else:
				reads_out[i]=reads[i]
	reads="".join(reads_out)
	return reads,qual


def sequencing_error_posibility(qual):
	return 10**(-(qual*36+2)/10.0)


def output(SAM,SAMW,SAMC,PE,FQ1,FQ2,READS_R1,QUALS_R1,READS_R2,QUALS_R2,technology,chrom,reads1,strand1,start1,end1,reads2,strand2,start2,end2,out,position,index,qual_mu,qual_sigma,random_sequencing_errors,front_point,front_point_sigma,front_qual,front_qual_sigma,end_point,end_point_sigma,end_qual,end_qual_sigma,Dynamic_qual):
	if Dynamic_qual:
		quality=bata_fun(qual_mu,qual_sigma)
		if quality<0.5:
			quality=0.5
		(reads1,qual1)=dynamic_quality(reads1,quality,random_sequencing_errors,front_point,front_point_sigma,front_qual,front_qual_sigma,end_point,end_point_sigma,end_qual,end_qual_sigma)
		(reads2,qual2)=dynamic_quality(reads2,quality,random_sequencing_errors,front_point,front_point_sigma,front_qual,front_qual_sigma,end_point,end_point_sigma,end_qual,end_qual_sigma)
	else:
		qual1=[qual_mu]*len(reads1)
		qual2=[qual_mu]*len(reads2)
	if technology == 'Solexa':
		output_Solexa(SAM,SAMW,SAMC,PE,FQ1,qual1,FQ2,qual2,chrom,reads1,strand1,start1,end1,reads2,strand2,start2,end2,out,position,index) 
	elif technology == 'SOLiD':
		output_SOLiD(SAM,SAMW,SAMC,PE,READS_R1,QUALS_R1,qual1,READS_R2,QUALS_R2,qual2,chrom,reads1,strand1,start1,end1,reads2,strand2,start2,end2,out,position,index) 
	elif technology == '454':
		output_454(SAM,SAMW,SAMC,PE,READS_R1,QUALS_R1,qual1,READS_R2,QUALS_R2,qual2,chrom,reads1,strand1,start1,end1,reads2,strand2,start2,end2,out,position,index)


def output_Solexa(SAM,SAMW,SAMC,PE,FQ1,qual1,FQ2,qual2,chrom,reads1,strand1,start1,end1,reads2,strand2,start2,end2,out,position,index):
	basic=':1'+':17'+':'+str(random.randint(1000,99999))+':'+str(random.randint(1000,99999))+'#0'+index
	quality1=[""]*len(qual1)
	for i in range(0,len(qual1)):
		quality1[i]=chr(int(qual1[i]*40.499)+64)
	quality1="".join(quality1)
	if position:
		first=chrom+'.'+str(start1)+'-'+str(end1)+'.'+strand1+'.'+str(start2)+'-'+str(end2)+'.'+strand2+basic+'/1'
	else:
		first='FC61FL8AAXX'+basic+'/1'
	FQ1.write('@'+first+'\n'+reads1+'\n+\n'+quality1+'\n')
	first2=quality2=''
	if reads2:
		quality2=[""]*len(qual2)
		for i in range(0,len(qual2)):
			quality2[i]=chr(int(qual2[i]*40.499)+64)
		quality2="".join(quality2)
		if position:
			first2=chrom+'.'+str(start1)+'-'+str(end1)+'.'+strand1+'.'+str(start2)+'-'+str(end2)+'.'+strand2+basic+'/2'
		else:
			first2='FC61FL8AAXX'+basic+'/2'
		FQ2.write('@'+first2+'\n'+reads2+'\n+\n'+quality2+'\n')
	if SAM:
		output_SAM(SAM,SAMW,SAMC,chrom,first,reads1,strand1,start1,end1,qual1,first2,reads2,strand2,start2,end2,qual2)




def output_SOLiD(SAM,SAMW,SAMC,PE,READS_R1,QUALS_R1,qual1,READS_R2,QUALS_R2,qual2,chrom,reads1,strand1,start1,end1,reads2,strand2,start2,end2,out,position,index):
	basic=' 17'+'_'+str(random.randint(1000,99999))+'_'+str(random.randint(1000,99999))+'#0'+index
	cspace_read1=reads_to_colorreads(reads1)
	if reads2:
		cspace_read2=reads_to_colorreads(reads2)
	quality1=[""]*(len(qual1)-1)
	for i in range(0,(len(qual1)-1)):
		quality1[i]=str(int(qual1[i]*43.499))
	quality1=" ".join(quality1)
	if position:
		first='SRR062657'+chrom+'.'+str(start1)+'-'+str(end1)+'.'+strand1+'.'+str(start2)+'-'+str(end2)+'.'+strand2+basic+'/1'
	else:
		first='SRR062657'+basic+'/1'
	READS_R1.write('>'+first+'\nT'+cspace_read1+'\n')
	QUALS_R1.write('>'+first+'\n'+quality1+'\n')
	first2=quality2=''
	if reads2:
		# fq2=open(out+'_2.fq','w')
		quality2=[""]*(len(qual2)-1)
		for i in range(0,(len(qual2)-1)):
			quality2[i]=str(int(qual2[i]*43.499))
		quality2=" ".join(quality2)
		if position:
			first2='SRR062657'+chrom+'.'+str(start1)+'-'+str(end1)+'.'+strand1+'.'+str(start2)+'-'+str(end2)+'.'+strand2+basic+'/2'
		else:
			first2='SRR062657'+basic+'/2'
		READS_R2.write('>'+first2+'\nT'+cspace_read2+'\n')
		QUALS_R2.write('>'+first2+'\n'+quality2+'\n')
	if SAM:
		output_SAM(SAM,SAMW,SAMC,chrom,first,reads1,strand1,start1,end1,qual1,first2,reads2,strand2,start2,end2,qual2)



def reads_to_colorreads(reads):
	color_reads=[""]*len(reads)
	for i in range(0,len(reads)-1):
		if reads[i] == 'A':
			if reads[i+1] == 'A':
				color_reads[i]= str(0)
			if reads[i+1] == 'C':
				color_reads[i]= str(1)
			if reads[i+1] == 'G':
				color_reads[i]= str(2)
			if reads[i+1] == 'T':
				color_reads[i]= str(3)
		if reads[i] == 'C':
			if reads[i+1] == 'A':
				color_reads[i]= str(1)
			if reads[i+1] == 'C':
				color_reads[i]= str(0)
			if reads[i+1] == 'G':
				color_reads[i]= str(3)
			if reads[i+1] == 'T':
				color_reads[i]= str(2)
		if reads[i] == 'G':
			if reads[i+1] == 'A':
				color_reads[i]= str(2)
			if reads[i+1] == 'C':
				color_reads[i]= str(3)
			if reads[i+1] == 'G':
				color_reads[i]= str(0)
			if reads[i+1] == 'T':
				color_reads[i]= str(1)
		if reads[i] == 'T':
			if reads[i+1] == 'A':
				color_reads[i]= str(3)
			if reads[i+1] == 'C':
				color_reads[i]= str(2)
			if reads[i+1] == 'G':
				color_reads[i]= str(1)
			if reads[i+1] == 'T':
				color_reads[i]= str(0)
	color_reads="".join(color_reads)
	return color_reads


def output_454(SAM,SAMW,SAMC,PE,READS_R1,QUALS_R1,qual1,READS_R2,QUALS_R2,qual2,chrom,reads1,strand1,start1,end1,reads2,strand2,start2,end2,out,position,index):
	basic=':1'+':17'+':'+str(random.randint(1000,99999))+':'+str(random.randint(1000,99999))+'#0'+index
	quality1=[""]*len(qual1)
	for i in range(0,len(qual1)):
		quality1[i]=str(int(qual1[i]*39.499))
	quality1=" ".join(quality1)
	if position:
		first='EM7LVYS01DSWC4.'+chrom+'.'+str(start1)+'-'+str(end1)+'.'+strand1+'.'+str(start2)+'-'+str(end2)+'.'+strand2+basic+'/1'
	else:
		first='EM7LVYS01DSWC4'+basic+'/1'
	READS_R1.write('>'+first+'\n'+reads1+'\n')
	QUALS_R1.write('>'+first+'\n'+quality1+'\n')
	first2=quality2=''
	if reads2:
		# fq2=open(out+'_2.fq','w')
		quality2=[""]*len(qual2)
		for i in range(0,len(qual2)):
			quality2[i]=str(int(qual2[i]*39.499))
		quality2=" ".join(quality2)
		if position:
			first2='EM7LVYS02FUAUF.'+chrom+'.'+str(start1)+'-'+str(end1)+'.'+strand1+'.'+str(start2)+'-'+str(end2)+'.'+strand2+basic+'/2'
		else:
			first2='EM7LVYS02FUAUF'+basic+'/2'
		READS_R2.write('>'+first2+'\n'+reads2+'\n')
		QUALS_R2.write('>'+first2+'\n'+quality2+'\n')
	if SAM:
		output_SAM(SAM,SAMW,SAMC,chrom,first,reads1,strand1,start1,end1,qual1,first2,reads2,strand2,start2,end2,qual2)



def output_SAM(SAM,samw,samc,chrom,first,reads1,strand1,start1,end1,qual1,first2,reads2,strand2,start2,end2,qual2):
	if reads2:
		if strand1 =='+.W' and strand2 =='-.C':						#PE and directional
			quality1=quality_sanger(qual1,'')
			samw.write(first+'\t0\t'+chrom+'\t'+str(start1)+'\t0\t'+str(len(reads1))+'M\t=\t0\t0\t'+reads1+'\t'+quality1+'\n')
			reads2=reverse_sequence(reads2)
			quality2=quality_sanger(qual2,True)
			samc.write(first2+'\t0\t'+chrom+'\t'+str(start2)+'\t0\t'+str(len(reads2))+'M\t=\t0\t0\t'+reads2+'\t'+quality2+'\n')
		if strand1 =='+.W' and strand2 =='-.W':						#PE and Watson
			quality1=quality_sanger(qual1,'')
			samw.write(first+'\t0\t'+chrom+'\t'+str(start1)+'\t0\t'+str(len(reads1))+'M\t=\t'+str(start2)+'\t'+str(start2-start1+len(reads2))+'\t'+reads1+'\t'+quality1+'\n')
			reads2=reverse_sequence(reads2)
			quality2=quality_sanger(qual2,True)
			samw.write(first2+'\t0\t'+chrom+'\t'+str(start2)+'\t0\t'+str(len(reads2))+'M\t=\t'+str(start1)+'\t-'+str(start2-start1+len(reads2))+'\t'+reads2+'\t'+quality2+'\n')
		if strand1 =='-.C' and strand2 =='+.C':						#PE and Crick
			quality2=quality_sanger(qual2,'')
			samc.write(first2+'\t0\t'+chrom+'\t'+str(start2)+'\t0\t'+str(len(reads2))+'M\t=\t'+str(start1)+'\t'+str(start1-start2+len(reads2))+'\t'+reads2+'\t'+quality2+'\n')
			reads1=reverse_sequence(reads1)
			quality1=quality_sanger(qual1,True)
			samc.write(first+'\t0\t'+chrom+'\t'+str(start1)+'\t0\t'+str(len(reads1))+'M\t=\t'+str(start2)+'\t-'+str(start1-start2+len(reads2))+'\t'+reads1+'\t'+quality1+'\n')
	else:
		if strand1 =='+.W':											#SE and +.Watson
			quality1=quality_sanger(qual1,'')
			samw.write(first+'\t0\t'+chrom+'\t'+str(start1)+'\t0\t'+str(len(reads1))+'M\t=\t0\t0\t'+reads1+'\t'+quality1+'\n')
		if strand1 =='-.W': 										#SE and -.Watson
			reads1=reverse_sequence(reads1)
			quality1=quality_sanger(qual1,True)
			samw.write(first2+'\t0\t'+chrom+'\t'+str(start1)+'\t0\t'+str(len(reads1))+'M\t=\t0\t0\t'+reads1+'\t'+quality1+'\n')
		if strand1 =='-.C':											#SE and -.Crick
			reads1=reverse_sequence(reads1)
			quality1=quality_sanger(qual1,True)
			samc.write(first+'\t0\t'+chrom+'\t'+str(start1)+'\t0\t'+str(len(reads1))+'M\t=\t0\t0\t'+reads1+'\t'+quality1+'\n')
		if strand1 =='+.C':											#SE and +.Crick
			quality1=quality_sanger(qual1,'')
			samc.write(first2+'\t0\t'+chrom+'\t'+str(start1)+'\t0\t'+str(len(reads1))+'M\t=\t0\t0\t'+reads1+'\t'+quality1+'\n')



def reverse_sequence(orin):
	l=len(orin)
	new=[""]*l
	for i in range(l,0,-1):
		new[l-i]=reverse_base(orin[i-1])
	new="".join(new)
	return new


def quality_sanger(qual1,reverse):
	quality1=[""]*len(qual1)
	if reverse:
		for i in range(len(qual1),0,-1):
			quality1[len(qual1)-i]=chr(int(qual1[i-1]*60.499)+33)
		quality1="".join(quality1)
	else:
		for i in range(0,len(qual1)):
			quality1[i]=chr(int(qual1[i]*60.499)+33)
		quality1="".join(quality1)
	return quality1


def cat_tmp(samw,samc,fq1,fq2,reads_r1,quals_r1,reads_r2,quals_r2):
	if fq1:
		os.system('cat '+fq1+'.tmp.* >>'+fq1+' && rm '+fq1+'.tmp.*')
	if fq2:
		os.system('cat '+fq2+'.tmp.* >>'+fq2+' && rm '+fq2+'.tmp.*')
	if samw:
		os.system('cat '+samw+'.tmp.* >>'+samw+' && rm '+samw+'.tmp.*')
	if samc:
		os.system('cat '+samc+'.tmp.* >>'+samc+' && rm '+samc+'.tmp.*')
	if reads_r1:
		os.system('cat '+reads_r1+'.tmp.* >>'+reads_r1+' && rm '+reads_r1+'.tmp.*')
	if reads_r2:
		os.system('cat '+reads_r2+'.tmp.* >>'+reads_r2+' && rm '+reads_r2+'.tmp.*')
	if quals_r1:
		os.system('cat '+quals_r1+'.tmp.* >>'+quals_r1+' && rm '+quals_r1+'.tmp.*')
	if quals_r2:
		os.system('cat '+quals_r2+'.tmp.* >>'+quals_r2+' && rm '+quals_r2+'.tmp.*')


#cerate reads for user input SNP
def cerate_reads_for_input_SNP(dbsnp,cpu,number_tmp,ref_start,ref,directional,Rate,max_err,output_ref,Fragment_length,Fragment_length_sigma,reads_length_mu,reads_length_sigma,SAM,samw,samc,PE,fq1,fq2,reads_r1,quals_r1,reads_r2,quals_r2,technology,chrom,out,position,index,qual_mu,qual_sigma,random_sequencing_errors,front_point,front_point_sigma,front_qual,front_qual_sigma,end_point,end_point_sigma,end_qual,end_qual_sigma,Dynamic_qual):
	i=0
	# open_output_file_as_a+
	FQ1=FQ2=READS_R1=READS_R2=QUALS_R1=QUALS_R2=SAMW=SAMC=''
	if technology == 'Solexa':
		FQ1=open(fq1+'.tmp.'+str(cpu),'w')
		if PE:
			FQ2=open(fq2+'.tmp.'+str(cpu),'w')
	if technology == 'SOLiD':
		READS_R1=open(reads_r1+'.tmp.'+str(cpu),'w')
		QUALS_R1=open(quals_r1+'.tmp.'+str(cpu),'w')
		if PE:
			READS_R2=open(reads_r2+'.tmp.'+str(cpu),'w')
			QUALS_R2=open(quals_r2+'.tmp.'+str(cpu),'w')
	if technology == '454':
		READS_R1=open(reads_r1+'.tmp.'+str(cpu),'w')
		QUALS_R1=open(quals_r1+'.tmp.'+str(cpu),'w')
		if PE:
			READS_R2=open(reads_r2+'.tmp.'+str(cpu),'w')
			QUALS_R2=open(quals_r2+'.tmp.'+str(cpu),'w')
	if SAM:
		SAMW=open(samw+'.tmp.'+str(cpu),'w')
		SAMC=open(samc+'.tmp.'+str(cpu),'w')
	
	while (i<number_tmp):
		fragment_length=int(random.normalvariate(Fragment_length,Fragment_length_sigma))
		reads_length=int(random.normalvariate(reads_length_mu,reads_length_sigma))
		if reads_length>fragment_length:
			reads_length=fragment_length-20
		start=random.randint(0,len(ref)-fragment_length)
		seq=[0]*fragment_length
		rate=[0]*fragment_length
		for j in range(start,start+fragment_length):
			# print "%d\t%d" % (number[chrom] , start)
			if not dbsnp[chrom,j]:
			# if srt(chrom,str(j)) in dbsnp
				seq[j-start]=ref[j]
				rate[j-start]=Rate[j]
			else:
				snp=dbsnp[chrom,j].split('\t')
				r=random.random()
				if r<=snp[0]:
					seq[j-start]="A"
					rate[j-start]=0
				if snp[0]<r<=(snp[0]+snp[1]):
					seq[j-start]="T"
					rate[j-start]=0
				if (snp[0]+snp[1])<r<=(snp[0]+snp[1]+snp[2]):
					seq[j-start]="C"
					if ref[j] == "C":
						rate[j-start]=rate[j]
					else:
						rate[j-start]=0
				if (snp[0]+snp[1]+snp[2])<r<=(snp[0]+snp[1]+snp[2]+snp[3]):
					seq[j-start]="G"
					if ref[j] == "G":
						rate[j-start]=rate[j]
					else:
						rate[j-start]=0
		seq="".join(seq)
		if PE :
			(reads1,strand1,start1,end1,reads2,strand2,start2,end2)=create_Reads(PE,directional,ref_start,seq,rate,start,reads_length,max_err)
		else:
			(reads1,strand1,start1,end1)=create_Reads(PE,directional,ref_start,seq,rate,start,reads_length,max_err)
		# output reads to file
		if PE :
			output(SAM,SAMW,SAMC,PE,FQ1,FQ2,READS_R1,QUALS_R1,READS_R2,QUALS_R2,technology,chrom,reads1,strand1,start1,end1,reads2,strand2,start2,end2,out,position,index,qual_mu,qual_sigma,random_sequencing_errors,front_point,front_point_sigma,front_qual,front_qual_sigma,end_point,end_point_sigma,end_qual,end_qual_sigma,Dynamic_qual)
		else:
			output(SAM,SAMW,SAMC,PE,FQ1,FQ2,READS_R1,QUALS_R1,READS_R2,QUALS_R2,technology,chrom,reads1,strand1,start1,end1,'','','','',out,position,index,qual_mu,qual_sigma,random_sequencing_errors,front_point,front_point_sigma,front_qual,front_qual_sigma,end_point,end_point_sigma,end_qual,end_qual_sigma,Dynamic_qual)
		i=i+1


# cerate reads for random SNP of diploid
def cerate_reads_for_random_SNP_of_diploid(cpu,number_tmp,ref_start,ref_A,ref_B,directional,rate_A,rate_B,max_err,output_ref,Fragment_length,Fragment_length_sigma,reads_length_mu,reads_length_sigma,SAM,samw,samc,PE,fq1,fq2,reads_r1,quals_r1,reads_r2,quals_r2,technology,chrom,out,position,index,qual_mu,qual_sigma,random_sequencing_errors,front_point,front_point_sigma,front_qual,front_qual_sigma,end_point,end_point_sigma,end_qual,end_qual_sigma,Dynamic_qual):
	i=0
	# open_output_file_as_a+
	FQ1=FQ2=READS_R1=READS_R2=QUALS_R1=QUALS_R2=SAMW=SAMC=''
	if technology == 'Solexa':
		FQ1=open(fq1+'.tmp.'+str(cpu),'w')
		if PE:
			FQ2=open(fq2+'.tmp.'+str(cpu),'w')
	if technology == 'SOLiD':
		READS_R1=open(reads_r1+'.tmp.'+str(cpu),'w')
		QUALS_R1=open(quals_r1+'.tmp.'+str(cpu),'w')
		if PE:
			READS_R2=open(reads_r2+'.tmp.'+str(cpu),'w')
			QUALS_R2=open(quals_r2+'.tmp.'+str(cpu),'w')
	if technology == '454':
		READS_R1=open(reads_r1+'.tmp.'+str(cpu),'w')
		QUALS_R1=open(quals_r1+'.tmp.'+str(cpu),'w')
		if PE:
			READS_R2=open(reads_r2+'.tmp.'+str(cpu),'w')
			QUALS_R2=open(quals_r2+'.tmp.'+str(cpu),'w')
	if SAM:
		SAMW=open(samw+'.tmp.'+str(cpu),'w')
		SAMC=open(samc+'.tmp.'+str(cpu),'w')

	while (i<number_tmp):
		fragment_length=int(random.normalvariate(Fragment_length,Fragment_length_sigma))
		reads_length=int(random.normalvariate(reads_length_mu,reads_length_sigma))
		if reads_length>fragment_length:
			reads_length=fragment_length-20
		if random.random()<0.5:
			start=random.randint(0,len(ref_A)-fragment_length)
			seq=ref_A[(start):(start+fragment_length)]
			rate=rate_A[(start):(start+fragment_length)]
		else:
			start=random.randint(0,len(ref_B)-fragment_length)
			seq=ref_B[(start):(start+fragment_length)]
			rate=rate_B[(start):(start+fragment_length)]
		if PE :
			(reads1,strand1,start1,end1,reads2,strand2,start2,end2)=create_Reads(PE,directional,ref_start,seq,rate,start,reads_length,max_err)
		else:
			(reads1,strand1,start1,end1)=create_Reads(PE,directional,ref_start,seq,rate,start,reads_length,max_err)
		# output reads to file
		if PE :
			output(SAM,SAMW,SAMC,PE,FQ1,FQ2,READS_R1,QUALS_R1,READS_R2,QUALS_R2,technology,chrom,reads1,strand1,start1,end1,reads2,strand2,start2,end2,out,position,index,qual_mu,qual_sigma,random_sequencing_errors,front_point,front_point_sigma,front_qual,front_qual_sigma,end_point,end_point_sigma,end_qual,end_qual_sigma,Dynamic_qual)
		else:
			output(SAM,SAMW,SAMC,PE,FQ1,FQ2,READS_R1,QUALS_R1,READS_R2,QUALS_R2,technology,chrom,reads1,strand1,start1,end1,'','','','',out,position,index,qual_mu,qual_sigma,random_sequencing_errors,front_point,front_point_sigma,front_qual,front_qual_sigma,end_point,end_point_sigma,end_qual,end_qual_sigma,Dynamic_qual)
		i=i+1


def cerate_reads_for_random_SNP_of_haploid_and_polyploid(cpu,number_tmp,ref_start,ref,directional,Rate,max_err,output_ref,Fragment_length,Fragment_length_sigma,reads_length_mu,reads_length_sigma,SAM,samw,samc,PE,fq1,fq2,reads_r1,quals_r1,reads_r2,quals_r2,technology,chrom,out,position,index,qual_mu,qual_sigma,random_sequencing_errors,front_point,front_point_sigma,front_qual,front_qual_sigma,end_point,end_point_sigma,end_qual,end_qual_sigma,Dynamic_qual):
	i=0
	# open_output_file_as_a+
	FQ1=FQ2=READS_R1=READS_R2=QUALS_R1=QUALS_R2=SAMW=SAMC=''
	if technology == 'Solexa':
		FQ1=open(fq1+'.tmp.'+str(cpu),'w')
		if PE:
			FQ2=open(fq2+'.tmp.'+str(cpu),'w')
	if technology == 'SOLiD':
		READS_R1=open(reads_r1+'.tmp.'+str(cpu),'w')
		QUALS_R1=open(quals_r1+'.tmp.'+str(cpu),'w')
		if PE:
			READS_R2=open(reads_r2+'.tmp.'+str(cpu),'w')
			QUALS_R2=open(quals_r2+'.tmp.'+str(cpu),'w')
	if technology == '454':
		READS_R1=open(reads_r1+'.tmp.'+str(cpu),'w')
		QUALS_R1=open(quals_r1+'.tmp.'+str(cpu),'w')
		if PE:
			READS_R2=open(reads_r2+'.tmp.'+str(cpu),'w')
			QUALS_R2=open(quals_r2+'.tmp.'+str(cpu),'w')
	if SAM:
		SAMW=open(samw+'.tmp.'+str(cpu),'w')
		SAMC=open(samc+'.tmp.'+str(cpu),'w')

	while (i<number_tmp):
		fragment_length=int(random.normalvariate(Fragment_length,Fragment_length_sigma))
		reads_length=int(random.normalvariate(reads_length_mu,reads_length_sigma))
		if reads_length>fragment_length:
			reads_length=fragment_length-20
		start=random.randint(0,len(ref)-fragment_length)
		seq=[0]*fragment_length
		rate=[0]*fragment_length
		for j in range(start,start+fragment_length):
			if ref[j] == "N":
				seq[j-start]=random.choice('ATCG')
			else:
				seq[j-start]=ref[j]
			rate[j-start]=Rate[j]
		seq="".join(seq)
		if PE :
			(reads1,strand1,start1,end1,reads2,strand2,start2,end2)=create_Reads(PE,directional,ref_start,seq,rate,start,reads_length,max_err)
		else:
			(reads1,strand1,start1,end1)=create_Reads(PE,directional,ref_start,seq,rate,start,reads_length,max_err)
		# output reads to file
		if PE :
			output(SAM,SAMW,SAMC,PE,FQ1,FQ2,READS_R1,QUALS_R1,READS_R2,QUALS_R2,technology,chrom,reads1,strand1,start1,end1,reads2,strand2,start2,end2,out,position,index,qual_mu,qual_sigma,random_sequencing_errors,front_point,front_point_sigma,front_qual,front_qual_sigma,end_point,end_point_sigma,end_qual,end_qual_sigma,Dynamic_qual)
		else:
			output(SAM,SAMW,SAMC,PE,FQ1,FQ2,READS_R1,QUALS_R1,READS_R2,QUALS_R2,technology,chrom,reads1,strand1,start1,end1,'','','','',out,position,index,qual_mu,qual_sigma,random_sequencing_errors,front_point,front_point_sigma,front_qual,front_qual_sigma,end_point,end_point_sigma,end_qual,end_qual_sigma,Dynamic_qual)
		i=i+1



def usage():
	print '''\n\n##############################################################################
#    BSSim: bisulfite sequencing simulator for next-generation sequencing    #
##############################################################################\n
	\n	BSSim is implemented in the Python language and run in an operating system-independent manner. It can allow users to mimic various methylation level (total methylation level of cytosines, percentage of cytosines that methylated and methylation level of total methylcytosines) and bisulfite conversion rate in CpG, CHG and CHH context, respectively. It can also simulate genetic variations that are divergent from the reference sequence along with the sequencing error and quality distributions. In the output, both directional/non-directional, various read length, single/paired-end reads and alignment data in the SAM format can be generated. BSSim is a cross-platform BS-seq simulator offers output read datasets not only suitable for Illumina's Solexa, but also for Roche's 454 and Applied Biosystems' SOLiD.
	\nUsage: ./BSSim.py [options]
	\nOptions:
	\nGeneral
	-h	Help.
	-i	Input reference sequence in the fasta format.
	-d	Sequencing depth (>0). Default: 30.
	-U	The max number of processes core (>0). Default: 2.
	-l	Read length (>0). Default: 90 bp.
	-s	Single-end pattern.
	-p	Paired-end pattern (default).
	-t	Sequencing platform: Solexa/SOLiD/454. Default: Solexa.
	-f	Fragment length (library size) (>0). Default: 300 bp.
	--FR	Standard deviation of -f (0~(f/2.58)). Default: 20 bp.
	-n	Number of reads to be generated (>0).
	-q	Quality score (mean value of quality score) [0~1]. Default: 0.95 (95% of highest score).
	-e	Number of max error base at the end of a reads [0~l]. Default: 0.
	-N	The number of bases to be read into RAM one time (>0). Default is 1000000.
	-o	Prefix of output file. Default is set by name of input file.
	-D	Directional: reads 1 is same direction with reference sequencing (Watson strand) and read 2 is from Crick strand. Default: non-directional.
	-P	Output position information into the output file. Default is not.
	-A	Output alignment result in SAM format. Default is not.
	-V	Version information.
	-R	Output the reference methylation information. Default is not.
			format is:\nchromosome	position	ref_genome	ref_A	methylation pattern	default methylation rate(ignore it if ref is A,T)	ref_B(homologous chromosome)	methylation pattern	default methylation rate




	\nDNA methylation
	--ML	Total methylation level of cytosines (overall DNA methylation level) (0~1). Default: 0.0612.
		--CL	CG methylation level (0~1). Default: 0.8653.
		--GL	CHG methylation level (0~1). Default: 0.0203.
		--HL	CHH methylationlevel (0~1). Default: 0.0227.

	--MR	All mC/C rate (the ratio of total methylcytosines/total cytosines) (0~1). Default: 0.073.
		--CR	mCG/CG rate (0~1). Default: 0.852.
		--GR	mCHG/CHG rate (0~1). Default: 0.019.
		--HR	mCHH/CHH rate (0~1). Default: 0.025.

	--MM	Methylation level of total methylcytosines. Default: 0.8529.
		--CM	mCG methylation level (0~1). Default: 0.8529.
		--GM	mCHG methylation level (0~1). Default: 0.0837.
		--HM	mCHH methylation level (0~1). Default: 0.9091*0.0994+(1-0.9091)*0.8965.

	--MCS	Standard deviation of --MM (0~(1-MM)*MM). Default: 0.01.
		--CS	Standard deviation of --CM (0~(1-CM)*CM). Default: (1-CM)*CM/2.0.
		--GS	Standard deviation of --GM (0~(1-GM)*GM). Default: (1-GM)*GM/2.0.
		--HS	Standard deviation of --HM (0~(1-HM)*HM). Default: (1-HM)*HM/2.0.

	--BC	All cytosines' bisulfite conversion rate [0~1]. Defualt is 0.998.
		--CC	CG conversion rate [0~1]. Default: 0.998.
		--GC	CHG conversion rate [0~1]. Default: 0.998.
		--HC	CHH conversion rate [0~1]. Default: 0.998.




	\nSNP
	-S	SNP file with SNP information, specifying location and frequency of SNPs.
		format is:\nChromosome	position	strand	frequency of A	frequency of T	frequency of C	frequency of G
		chr10	1	+	0	0.4	0	0.6
		chr10	2	+	0.3	0.2	0.1	0.4
	--DS	Do not add SNP. Default is add (based on prior probability).
	-G	Polyploid type of reference sequencing (>0). Default: 2.
	-Y	The frequency of homozygous SNPs [0~(1-Z)]. Default: 0.0005.
	-Z	The frequency of heterozygous SNPs [0~(1-Y)]. Default: 0.001.




	\nRead quality
	-q	Quality score (mean value of quality score). Default: 0.95 (95% of highest score).
		--DQ	Randomly introduce quality value. Default: uniform quality score.
		--RE	Randomly introduce sequencing errors by sequencing quality value (Q =-10*log10(e),Q is the sequencing quality value (Phred score), e is the error rate, Massingham, et al., 2012). Default is not.
		--QS	Standard deviation of -q (0~(1-q)*q). Default: (1-q)*q/2.

	BSSim can also allow users to split every read into three part (The head part, The end part and interval) to add different quality value along the read.


		--FP		(Lengh of the head part)/(total read length) (0~1). Default: 0.01 (1% of total read length).
		--FPS		Standard deviation of --FP (0~(1-FP)*FP). Default: (1-FP)*FP/2.
		--FQ		The mean quality value of the head part less than -q (0~q). Default: 0.1
		--FS		Standard deviation of --FQ (0~(1-(q-FQ))*(q-FQ)). Default: (1-(q-FQ))*(q-FQ)/8.
		--EP		(Lengh of the end part)/(total read length) (0~1). Default: 0.8 (80% of total read length).
		--EPS		Standard deviation of --EP (0~(1-EP)*EP). Default: (1-EP)*EP/2.
		--EQ		The mean quality value of the end part less than -q (0~q). Default: 0.2
		--ES		Standard deviation of --EQ (0~(1-(q-EQ))*(q-EQ)). Default: (1-(q-EQ))*(q-EQ)/4.

 \nExample:
	./BSSim.py -i test.fa
	./BSSim.py -i test.fa -d 10 -t 454 -U 4 -G 1 -s -A -R -o out
	./BSSim.py -i test.fa -d 10 -t 454 -U 4 -G 1 -s -A -R --CR 0.9 --CM 0.4 --HM 0.3 --HC 0.7 -o out\n\n'''


def main(argv):
	print(time.strftime('BSsim start work at: %Y-%m-%d %H:%M:%S'))
	starttime = time.time()
	refFile=''
	depth=30
	CPUS=2
	reads_length_mu=90
	reads_length_sigma=0
	abund=''
	number=''
	PE=True
	technology='Solexa'
	Fragment_length=300
	Fragment_length_sigma=20
	number_of_seqs=''
	qual_mu=0.95
	qual_sigma=0.005
	Dynamic_qual=False
	random_sequencing_errors=False
	front_point=0.01
	front_point_sigma=(1-front_point)*front_point/2.0
	front_qual=0.1
	front_qual_sigma=(1-(qual_mu-front_qual))*(qual_mu-front_qual)/8.0
	end_point=0.8
	end_point_sigma=(1-end_point)*end_point/2.0
	end_qual=0.2
	end_qual_sigma=(1-(qual_mu-end_qual))*(qual_mu-end_qual)/4.0
	base_num=1000000
	max_err=0
	mode=''
	index=''
	CG_conversion_rate = CHG_conversion_rate = CHH_conversion_rate = conversion_rate = 0.998
	CG_methylation_level = CHG_methylation_level = CHH_methylation_level = C_methylation_level = 0
	mC_rate=0.06122667								#Date from Yanhuang
	mCG_rate=0.8653									#Date from Yanhuang
	mCHG_rate=0.0203								#Date from Yanhuang
	mCHH_rate=0.0227								#Date from Yanhuang
	mCG_mu=0.8529
	mCHG_mu=0.0837
	mCHH_mu=0.9091*0.0994+(1-0.9091)*0.8965
	mCG_sigma=mCHG_sigma=mCHH_sigma=mC_sigma=0.01
	models=''
	directional=False
	circular=False
	out=''
	position=False
	DS=False
	p1=0.0005					#homozygote
	p2=0.001					#heterozygote
	CG_beta_distribution=CHG_beta_distribution=CHH_beta_distribution=beta_distribution=False
	output_ref=False
	dbsnpfile=''
	Polyploid=2
	meta=False
	complementary = 1
	SAM = False
	adapter = False
	adapter1 = "GATCGGAAGAGCACACGTCT"
	adapter2 = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	try:
		opts, args=getopt.getopt(sys.argv[1:], "Rhi:d:l:spt:f:n:q:e:o:AVPS:G:N:U:Y:Z:DB",["CB","GB","HB","DS","FR=","BC=","CC=","HC=","GC=","ML=","CL=","GL=","HL=","MR=","CR=","GR=","HR=","MM=","CM=","GM=","HM=","MCS=","CS=","GS=","HS=","DQ","RE","QS=","FP=","FPS=","FQ=","FS=","MQ=","MS=","EP=","EPS=","EQ=","ES="])
	except getopt.GetoptError as e:
		usage()
		print e.msg
		print "Please check your options. maybe you should use \"--CL\", not \"-CL\""
		sys.exit()
	for opt, arg in opts:
		if opt == '-h':
			usage()
			sys.exit()
		elif opt == '-i':
			refFile=arg
		elif opt == '-d':
			depth=int(arg)
		elif opt == '-U':
			CPUS=int(arg)
		elif opt == '-s':
			SE=True
			PE=False
		elif opt == '-p':
			PE=True
		elif opt == '-t':
			technology=arg
			if technology == 'SOLiD':
				reads_length_mu=50
				qual_mu=0.98
			elif technology == '454':
				reads_length_mu=250
				reads_length_sigma=100
				qual_mu=0.99
				end_qual=0.03
				end_qual_sigma=0.01
			elif technology == 'Solexa':
				tmp=1
			else:
				print "Sequencing platform must be Solexa/SOLiD/454.\n"
				sys.exit()
		elif opt == '-l':
			reads_length_mu=int(arg)
		elif opt == '-f':
			Fragment_length=int(arg)
		elif opt == '--FR':
			Fragment_length_sigma=int(arg)
			if Fragment_length_sigma<(Fragment_length/2.58):
				tmp=1
			else:
				print "Standard deviation of -f must at range of (0~(f/2.58)).\n"
				sys.exit()
		elif opt == '-n':
			number_of_seqs=int(arg)
		elif opt == '-q':
			qual_mu=float(arg)
			if qual_mu>=0 and qual_mu<=1:
				tmp=1
			else:
				print "Quality score (mean value of quality score) must at range of [0~1].\n"
				sys.exit()
		elif opt == '-e':
			max_err=int(arg)
		elif opt == '-o':
			out=arg
		elif opt == '-a':
			adapter=arg
		elif opt == '-P':
			position=True
		elif opt == '-D':
			directional=True
			print "Simulate the directional reads!"
		elif opt == '-B':
			CG_beta_distribution=CHG_beta_distribution=CHH_beta_distribution=beta_distribution=True
		elif opt == '--CB':
			CG_beta_distribution=True
		elif opt == '--GB':
			CHG_beta_distribution=True
		elif opt == '--HB':
			CHH_beta_distribution=True
		elif opt == '--DS':
			DS=True
		elif opt == '-S':
			dbsnpfile=arg
		elif opt == '-G':
			Polyploid=int(arg)
		elif opt == '-N':
			base_num=int(arg)
		elif opt == '-Y':
			p1=float(arg)
		elif opt == '-Z':
			p2=float(arg)
		elif opt == '-R':
			output_ref=True
		elif opt == '-A':
			SAM=True
		elif opt == '-V':
			print '''Program: BSSim\nCompile Date: Friday, June 29, 2012\nAuthor Affiliations: Institute of Genomic Medicine, Wenzhou Medical College, Wenzhou 325035, China\nVersion: 1.2\nContact: xigyou@gmail.com'''
			sys.exit()
		elif opt == '--BC':
			conversion_rate=float(arg)
			if conversion_rate>=0 and conversion_rate<=1:
				CHH_conversion_rate=CHG_conversion_rate=CG_conversion_rate=conversion_rate
			else:
				print "All cytosines' bisulfite conversion rate must at range of [0~1].\n"
				sys.exit()
		elif opt == '--CC':
			CG_conversion_rate=float(arg)
			if CG_conversion_rate>=0 and CG_conversion_rate<=1:
				tmp=1
			else:
				print "CG conversion rate must at range of [0~1].\n"
				sys.exit()
		elif opt == '--GC':
			CHG_conversion_rate=float(arg)
			if CHG_conversion_rate>=0 and CHG_conversion_rate<=1:
				tmp=1
			else:
				print "CHG conversion rate must at range of [0~1].\n"
				sys.exit()
		elif opt == '--HC':
			CHH_conversion_rate=float(arg)
			if CHH_conversion_rate>=0 and CHH_conversion_rate<=1:
				tmp=1
			else:
				print "CHH conversion rate must at range of [0~1].\n"
				sys.exit()
		elif opt == '--ML':
			C_methylation_level=float(arg)
			if C_methylation_level>0 and C_methylation_level<1:
				CHH_methylation_level=CHG_methylation_level=CG_methylation_level=C_methylation_level
			else:
				print "Total methylation level of cytosines must at range of (0~1).\n"
				sys.exit()
		elif opt == '--CL':
			CG_methylation_level=float(arg)
			if CG_methylation_level>0 and CG_methylation_level<1:
				tmp=1
			else:
				print "CG methylation level of cytosines must at range of (0~1).\n"
				sys.exit()
		elif opt == '--GL':
			CHG_methylation_level=float(arg)
			if CHG_methylation_level>0 and CHG_methylation_level<1:
				tmp=1
			else:
				print "CHG methylation level of cytosines must at range of (0~1).\n"
				sys.exit()
		elif opt == '--HL':
			CHH_methylation_level=float(arg)
			if CHH_methylation_level>0 and CHH_methylation_level<1:
				tmp=1
			else:
				print "CHH methylation level of cytosines must at range of (0~1).\n"
				sys.exit()
		elif opt == '--MR':
			mC_rate=float(arg)
			if mC_rate>0 and mC_rate<1:
				mCHH_rate=mCHG_rate=mCG_rate=mC_rate
			else:
				print "The ratio of total methylcytosines/total cytosines must at range of (0~1).\n"
				sys.exit()
		elif opt == '--CR':
			mCG_rate=float(arg)
			if mCG_rate>0 and mCG_rate<1:
				tmp=1
			else:
				print "mCG/CG rate must at range of (0~1).\n"
				sys.exit()
		elif opt == '--GR':
			mCHG_rate=float(arg)
			if mCHG_rate>0 and mCHG_rate<1:
				tmp=1
			else:
				print "mCHG/CHG rate must at range of (0~1).\n"
				sys.exit()
		elif opt == '--HR':
			mCHH_rate=float(arg)
			if mCHH_rate>0 and mCHH_rate<1:
				tmp=1
			else:
				print "mCHH/CHH rate must at range of (0~1).\n"
				sys.exit()
		elif opt == '--MM':
			mC_mu=float(arg)
			if mC_mu>0 and mC_mu<1:
				mCG_mu=mCHG_mu=mCHH_mu=mC_mu
				mCG_sigma=mCHG_sigma=mCHH_sigma=(1-mC_mu)*(mC_mu)/2.0
				beta_distribution=True
			else:
				print "Methylation level of total methylcytosines must at range of (0~1).\n"
				sys.exit()
		elif opt == '--CM':
			mCG_mu=float(arg)
			if mCG_mu>0 and mCG_mu<1:
				mCG_sigma=(1-mCG_mu)*(mCG_mu)/2.0
				CG_beta_distribution=True
			else:
				print "mCG methylation level must at range of (0~1).\n"
				sys.exit()
		elif opt == '--GM':
			mCHG_mu=float(arg)
			if mCHG_mu>0 and mCHG_mu<1:
				mCHG_sigma=(1-mCHG_mu)*(mCHG_mu)/2.0
				CHG_beta_distribution=True
			else:
				print "mCHG methylation level must at range of (0~1).\n"
				sys.exit()
		elif opt == '--HM':
			mCHH_mu=float(arg)
			if mCHH_mu>0 and mCHH_mu<1:
				mCHH_sigma=(1-mCHH_mu)*(mCHH_mu)/2.0
				CHH_beta_distribution=True
			else:
				print "mCHH methylation level must at range of (0~1).\n"
				sys.exit()
		elif opt == '--MCS':
			mC_sigma=float(arg)
			if mC_sigma>0 and mC_sigma<((1-mC_mu)*(mC_mu)):
				mCG_sigma=mCHG_sigma=mCHH_sigma=mC_sigma
				beta_distribution=True
			else:
				print "Standard deviation of --MM must at range of (0~(1-MM)*MM).\n"
				sys.exit()
		elif opt == '--CS':
			mCG_sigma=float(arg)
			if mCG_sigma>0 and mCG_sigma<((1-mCG_mu)*(mCG_mu)):
				CG_beta_distribution=True
			else:
				print "Standard deviation of --CM must at range of (0~(1-CM)*CM).\n"
				sys.exit()
		elif opt == '--GS':
			mCHG_sigma=float(arg)
			if mCHG_sigma>0 and mCHG_sigma<((1-mCHG_mu)*(mCHG_mu)):
				CHG_beta_distribution=True
			else:
				print "Standard deviation of --GM must at range of (0~(1-GM)*GM).\n"
				sys.exit()
		elif opt == '--HS':
			mCHH_sigma=float(arg)
			if mCHH_sigma>0 and mCHH_sigma<((1-mCHH_mu)*(mCHH_mu)):
				CHH_beta_distribution=True
			else:
				print "Standard deviation of --HM must at range of (0~(1-HM)*HM).\n"
				sys.exit()
		elif opt == '--DQ':
			Dynamic_qual=True
			print "Boot the randomly introduce quality value!\n"
		elif opt == '--RE':
			random_sequencing_errors=True
			print "Randomly introduce sequencing errors by sequencing quality value!\n"
		elif opt == '--QS':
			qual_sigma=float(arg)
			Dynamic_qual=True
			print "Boot the randomly introduce quality value!\n"
		elif opt == '--FP':
			front_point=float(arg)
			if front_point>0 and front_point<1:
				front_point_sigma=(1-front_point)*(front_point)/2.0
			else:
				print "(Lengh of the head part)/(total read length) must at range of (0~1).\n"
				sys.exit()
		elif opt == '--FPS':
			front_point_sigma=float(arg)
			if front_point_sigma>0 and front_point_sigma<((1-front_point)*(front_point)):
				tmp=1
			else:
				print "Standard deviation of --FP must at range of (0~(1-FP)*FP).\n"
				sys.exit()
		elif opt == '--FQ':
			front_qual=float(arg)
			if front_qual>0 and front_qual<qual_mu:
				front_qual_sigma=(1-qual_mu+front_qual)*(qual_mu-front_qual)/8.0
			else:
				print "The mean quality value of the head part less than -q must at range of (0~q).\n"
				sys.exit()
		elif opt == '--FS':
			front_qual_sigma=float(arg)
			if front_qual_sigma>0 and front_qual_sigma<((1-qual_mu+front_qual)*(qual_mu-front_qual)):
				tmp=1
			else:
				print "Standard deviation of --FQ must at range of (0~(1-(q-FQ))*(q-FQ)).\n"
				sys.exit()
		elif opt == '--EP':
			end_point=float(arg)
			if end_point>0 and end_point<1:
				end_point_sigma=(1-end_point)*(end_point)/2.0
			else:
				print "(Lengh of the end part)/(total read length) must at range of (0~1).\n"
				sys.exit()
		elif opt == '--EPS':
			end_point_sigma=float(arg)
			if end_point_sigma>0 and end_point_sigma<((1-end_point)*(end_point)):
				tmp=1
			else:
				print "Standard deviation of --EP must at range of (0~(1-EP)*EP).\n"
				sys.exit()
		elif opt == '--EQ':
			end_qual=float(arg)
			if end_qual>0 and end_qual<qual_mu:
				end_qual_sigma=(1-qual_mu+end_qual)*(qual_mu-end_qual)/4.0
			else:
				print "The mean quality value of the end part less than -q must at range of (0~q).\n"
				sys.exit()
		elif opt == '--ES':
			end_qual_sigma=float(arg)
			if end_qual_sigma>0 and end_qual_sigma<((1-qual_mu+end_qual)*(qual_mu-end_qual)):
				tmp=1
			else:
				print "Standard deviation of --EQ must at range of (0~(1-(q-EQ))*(q-EQ)).\n"
				sys.exit()
		

	if refFile == '':
		usage()
		print "Input a fa file as reference, please!"
		sys.exit()

	if CG_methylation_level > 0 and not CG_beta_distribution:
		if CG_methylation_level > mCG_mu:
			mCG_mu=CG_methylation_level/mCG_rate
			if mCG_mu<1:
				# print "The CG methylation level is higher than default mCG_mu, auto changed, new mCG_mu is %s." % (mCG_mu)
				CG_beta_distribution=True
			if mCG_mu>1:
				mCG_mu=mCG_rate=CG_methylation_level**(1/2.0)
				CG_beta_distribution=True
		else:
			mCG_rate=CG_methylation_level/mCG_mu
	else:
		CG_methylation_level=mCG_rate*mCG_mu

	if CHG_methylation_level > 0 and not CHG_beta_distribution:
		if CHG_methylation_level > mCHG_mu:
			mCHG_mu=CHG_methylation_level/mCHG_rate
			if mCHG_mu<1:
				# print "The CHG_methylation_level is higher than default mCHG_mu, auto changed, new mCHG_mu is %s." % (mCHG_mu)
				CHG_beta_distribution=True
			if mCHG_mu>1:
				mCHG_mu=mCHG_rate=CHG_methylation_level**(1/2.0)
				CHG_beta_distribution=True
				# print "mCHG_mu is larger than 1, you need to set two of the follow.(mCHG_mu=CHG_methylation_level/mCHG_rate)\n"
		else:
			mCHG_rate=CHG_methylation_level/mCHG_mu
	else:
		CHG_methylation_level=mCHG_rate*mCHG_mu

	if CHH_methylation_level > 0 and not CHH_beta_distribution:
		if CHH_methylation_level > mCHH_mu:
			mCHH_mu=CHH_methylation_level/mCHH_rate
			if mCHH_mu<1:
				# print "The CHH_methylation_level is higher than default mCHH_mu, auto changed, new mCHH_mu is %s." % (mCHH_mu)
				CHH_beta_distribution=True
			if mCHH_mu>1:
				mCHH_mu=mCHH_rate=CHH_methylation_level**(1/2.0)
				CHH_beta_distribution=True
				# print "mCHH_mu is larger than 1, you need to set two of the follow.(mCHH_mu=CHH_methylation_level/mCHH_rate)\n"
		else:
			mCHH_rate=CHH_methylation_level/mCHH_mu
	else:
		CHH_methylation_level=mCHH_rate*mCHH_mu

# ### OPENING OUTPUT FILEHANDLES
	fq1=fq2=reads_r1=reads_r2=quals_r1=quals_r2=samw=samc=''
	if out == '':
		out=re.sub( r'\.fa' , '' , refFile )
		out=out+'-Fragment_length-'+str(Fragment_length)+'-depth-'+str(depth)+'-methylation_level_CG-'+str(CG_methylation_level)+'-CHG-'+str(CHG_methylation_level)+'-CHH-'+str(CHH_methylation_level)
	if technology == 'Solexa':
		fq1=out+'.1.fq'
		FQ1=open(fq1,'w')
		FQ1.close()
		if PE:
			fq2=out+'.2.fq'
			FQ2=open(fq2,'w')
			FQ1.close()
	if technology == 'SOLiD':
		reads_r1=out+'.1.csfasta'
		READS_R1=open(out+'.1.csfasta','w')
		READS_R1.close()
		quals_r1=out+'.1_QV.qual'
		QUALS_R1=open(quals_r1,'w')
		QUALS_R1.close()
		if PE:
			reads_r2=out+'.2.csfasta'
			READS_R2=open(reads_r2,'w')
			READS_R2.close()
			quals_r2=out+'.2_QV.qual'
			QUALS_R2=open(quals_r2,'w')
			QUALS_R2.close()
	if technology == '454':
		reads_r1=out+'.1.fna'
		READS_R1=open(reads_r1,'w')
		READS_R1.close()
		quals_r1=out+'.1.qual'
		QUALS_R1=open(quals_r1,'w')
		QUALS_R1.close()
		if PE:
			reads_r2=out+'.2.fna'
			READS_R2=open(reads_r2,'w')
			READS_R2.close()
			quals_r2=out+'.2.qual'
			QUALS_R2=open(quals_r2,'w')
			QUALS_R2.close()
	if output_ref:
		Ref=out+'.ref'
		REF=open(Ref,'w')

	print(time.strftime('start read fa at: %Y-%m-%d %H:%M:%S'))
	fa=read_fa(refFile)

	if SAM :
		samw=out+'.Watson.sam'
		SAMW=open(samw,'w')
		samc=out+'.Crick.sam'
		SAMC=open(samc,'w')
		SAMW.write('@HD\tVN:1.0\tSO:unsorted\tYou can get methylcytosine points on the Watson strand, the cytosines\' methylation information on reference.\n')
		SAMC.write('@HD\tVN:1.0\tSO:unsorted\tYou can get methylcytosine points on the Crick strand, the guanines\' methylation information on reference.\n')
		for chrom in sorted(fa.keys()):
			SAMW.write('@SQ\tSN:'+chrom+'\tLN:'+str(len(fa[chrom]))+'\n')
			SAMC.write('@SQ\tSN:'+chrom+'\tLN:'+str(len(fa[chrom]))+'\n')
		SAMW.close()
		SAMC.close()


	#calculate reads number
	if number_of_seqs:
		total=0
		for chrom in sorted(fa.keys()):
			total=total+len(fa[chrom])
		depth=2*reads_length_mu*number_of_seqs/total
	number={}
	for chrom in sorted(fa.keys()):
		number[chrom]=int(depth*len(fa[chrom])/reads_length_mu)
		if PE:
			number[chrom]=int(number[chrom]/2)


	# print(time.strftime('start add SNP at: %Y-%m-%d %H:%M:%S'))
# create SNP database
#user input SNP
	if dbsnpfile:
		dbsnp=read_dbsnp(dbsnpfile)
		for chrom in sorted(fa.keys()):
			print "\nsimulate reads from: %s" % (chrom)
			print(time.strftime('start at: %Y-%m-%d %H:%M:%S'))
			DS=True
			to_len=len(fa[chrom])
			for n in range(0,to_len/base_num+1):
				ref_start=n*base_num
				print "%d of %d" % (n+1,to_len/base_num+1)
				if to_len<base_num :
					ref_tmp=fa[chrom][0:to_len]
					number_tmp=number[chrom]*(to_len-ref_start)/to_len
				elif n==(to_len/base_num) :
					ref_tmp=fa[chrom][ref_start:to_len]
					number_tmp=number[chrom]*(to_len-ref_start)/to_len
				else:
					ref_tmp=fa[chrom][ref_start:(n+1)*base_num]
					number_tmp=number[chrom]*base_num/to_len
				ref=random_snp(ref_tmp,p1,p2,DS,Polyploid)
			#identify methylation point
				(Rate,ref_pattern)=methyl(ref,conversion_rate,CG_conversion_rate,CHG_conversion_rate,CHH_conversion_rate,mC_rate,mCG_rate,mCHG_rate,mCHH_rate,CG_beta_distribution,mCG_mu,mCG_sigma,CHG_beta_distribution,mCHG_mu,mCHG_sigma,CHH_beta_distribution,mCHH_mu,mCHH_sigma)
			#out put ref
				if output_ref:
					for j in range(0,len(ref)):
						if dbsnp[chrom,str(j+1)]!=0:
							REF.write(chrom+'\t'+str(ref_start+j+1)+'\t'+fa[chrom][ref_start+j]+'\t'+ref[j]+'\t'+str(ref_pattern[j])+'\t'+str(1-Rate[j])+'\t'+dbsnp[chrom,str(j+1)]+'\n')
						else:
							REF.write(chrom+'\t'+str(ref_start+j+1)+'\t'+fa[chrom][ref_start+j]+'\t'+ref[j]+'\t'+str(ref_pattern[j])+'\t'+str(1-Rate[j])+'\n')
			#cerate reads
				number_tmp=int(number_tmp/CPUS+1)
				process = []
				pool = Pool(processes=CPUS)	# set the processes max number
				for cpu in range(CPUS):
					# print "Proc(%d) Start..."%cpu
					result = pool.apply_async(cerate_reads_for_input_SNP,(dbsnp,cpu,number_tmp,ref_start,ref,directional,Rate,max_err,output_ref,Fragment_length,Fragment_length_sigma,reads_length_mu,reads_length_sigma,SAM,samw,samc,PE,fq1,fq2,reads_r1,quals_r1,reads_r2,quals_r2,technology,chrom,out,position,index,qual_mu,qual_sigma,random_sequencing_errors,front_point,front_point_sigma,front_qual,front_qual_sigma,end_point,end_point_sigma,end_qual,end_qual_sigma,Dynamic_qual))
				pool.close()
				pool.join()
				if result.successful():
					print '%d CPUS run successful!' %CPUS
					cat_tmp(samw,samc,fq1,fq2,reads_r1,quals_r1,reads_r2,quals_r2)
				else:
					print 'Some CPUS not run successful!!'
				# cpu=1
				# cerate_reads_for_input_SNP(dbsnp,cpu,number_tmp,ref_start,ref,directional,Rate,max_err,output_ref,Fragment_length,Fragment_length_sigma,reads_length_mu,reads_length_sigma,SAM,samw,samc,PE,fq1,fq2,reads_r1,quals_r1,reads_r2,quals_r2,technology,chrom,out,position,index,qual_mu,qual_sigma,random_sequencing_errors,front_point,front_point_sigma,front_qual,front_qual_sigma,end_point,end_point_sigma,end_qual,end_qual_sigma,Dynamic_qual)
				print(time.strftime('\tend at: %Y-%m-%d %H:%M:%S'))

#random SNP
	else:
	#random SNP of diploid
		if Polyploid == 2:
			CG_beta_distribution=CHG_beta_distribution=CHH_beta_distribution=beta_distribution=True
			for chrom in sorted(fa.keys()):
				print "\nsimulate reads from: %s" % (chrom)
				print(time.strftime('start at: %Y-%m-%d %H:%M:%S'))
				to_len=len(fa[chrom])
				for n in range(0,to_len/base_num+1):
					ref_start=n*base_num
					print "%d of %d" % (n+1,to_len/base_num+1)
					if to_len<base_num :
						ref_tmp=fa[chrom][0:to_len]
						number_tmp=number[chrom]*(to_len-ref_start)/to_len
					elif n==(to_len/base_num) :
						ref_tmp=fa[chrom][ref_start:to_len]
						number_tmp=number[chrom]*(to_len-ref_start)/to_len
					else:
						ref_tmp=fa[chrom][ref_start:(n+1)*base_num]
						number_tmp=number[chrom]*base_num/to_len
					[ref_A,ref_B]=random_snp_diploid(ref_tmp,p1,p2,DS)
				#identify methylation point
					(rate_A,ref_pattern_A,rate_B,ref_pattern_B)=methyl_for_random_SNP_of_diploid(ref_A,ref_B,conversion_rate,CG_conversion_rate,CHG_conversion_rate,CHH_conversion_rate,mC_rate,mCG_rate,mCHG_rate,mCHH_rate,CG_beta_distribution,mCG_mu,mCG_sigma,CHG_beta_distribution,mCHG_mu,mCHG_sigma,CHH_beta_distribution,mCHH_mu,mCHH_sigma)
				#output ref
					if output_ref:
						for j in range(0,len(ref_A)):
							REF.write(chrom+'\t'+str(ref_start+j+1)+'\t'+fa[chrom][ref_start+j]+'\t'+ref_A[j]+'\t'+str(ref_pattern_A[j])+'\t'+str(1-rate_A[j])+'\t'+ref_B[j]+'\t'+str(ref_pattern_B[j])+'\t'+str(1-rate_B[j])+'\n')
				# cerate reads for diploid
					number_tmp=int(number_tmp/CPUS+1)
					process = []
					pool = Pool(processes=CPUS)	# set the processes max number
					for cpu in range(CPUS):
						# print "Proc(%d) Start..."%cpu
						result = pool.apply_async(cerate_reads_for_random_SNP_of_diploid,(cpu,number_tmp,ref_start,ref_A,ref_B,directional,rate_A,rate_B,max_err,output_ref,Fragment_length,Fragment_length_sigma,reads_length_mu,reads_length_sigma,SAM,samw,samc,PE,fq1,fq2,reads_r1,quals_r1,reads_r2,quals_r2,technology,chrom,out,position,index,qual_mu,qual_sigma,random_sequencing_errors,front_point,front_point_sigma,front_qual,front_qual_sigma,end_point,end_point_sigma,end_qual,end_qual_sigma,Dynamic_qual))
					pool.close()
					pool.join()
					if result.successful():
						print '%d CPUS run successful!' %CPUS
						cat_tmp(samw,samc,fq1,fq2,reads_r1,quals_r1,reads_r2,quals_r2)
					else:
						print 'Some CPUS not run successful!!'
					# cpu=1
					# cerate_reads_for_random_SNP_of_diploid(cpu,number_tmp,ref_A,ref_B,directional,rate_A,rate_B,max_err,output_ref,Fragment_length,Fragment_length_sigma,reads_length_mu,reads_length_sigma,SAM,samw,samc,PE,fq1,fq2,reads_r1,quals_r1,reads_r2,quals_r2,technology,chrom,out,position,index,qual_mu,qual_sigma,random_sequencing_errors,front_point,front_point_sigma,front_qual,front_qual_sigma,end_point,end_point_sigma,end_qual,end_qual_sigma,Dynamic_qual)
					print(time.strftime('This part finish at: %Y-%m-%d %H:%M:%S'))
	#random SNP of haploid and polyploid
		else:
			for chrom in sorted(fa.keys()):
				print "\nsimulate reads from: %s" % (chrom)
				print(time.strftime('start at: %Y-%m-%d %H:%M:%S'))
				to_len=len(fa[chrom])
				for n in range(0,to_len/base_num+1):
					ref_start=n*base_num
					print "%d of %d" % (n+1,to_len/base_num+1)
					if to_len<base_num :
						ref_tmp=fa[chrom][0:to_len]
						number_tmp=number[chrom]*(to_len-ref_start)/to_len
					elif n==(to_len/base_num) :
						ref_tmp=fa[chrom][ref_start:to_len]
						number_tmp=number[chrom]*(to_len-ref_start)/to_len
					else:
						ref_tmp=fa[chrom][ref_start:(n+1)*base_num]
						number_tmp=number[chrom]*base_num/to_len
					ref=random_snp(ref_tmp,p1,p2,DS,Polyploid)
				#identify methylation point
					(Rate,ref_pattern)=methyl(ref,conversion_rate,CG_conversion_rate,CHG_conversion_rate,CHH_conversion_rate,mC_rate,mCG_rate,mCHG_rate,mCHH_rate,CG_beta_distribution,mCG_mu,mCG_sigma,CHG_beta_distribution,mCHG_mu,mCHG_sigma,CHH_beta_distribution,mCHH_mu,mCHH_sigma)
				#out put ref
					if output_ref:
						for j in range(0,len(ref)):
							if ref[j] == "N":
								REF.write(chrom+'\t'+str(ref_start+j+1)+'\t'+fa[chrom][ref_start+j]+'\t'+ref[j]+'\t'+str(ref_pattern[j])+'\t'+str(1-Rate[j])+'\t'+"0.25\t0.25\t0.25\t0.25"+'\n')
							else:
								REF.write(chrom+'\t'+str(ref_start+j+1)+'\t'+fa[chrom][ref_start+j]+'\t'+ref[j]+'\t'+str(ref_pattern[j])+'\t'+str(1-Rate[j])+'\n')
				#cerate reads
					number_tmp=int(number_tmp/CPUS+1)
					process = []
					pool = Pool(processes=CPUS)	# set the processes max number
					for cpu in range(CPUS):
						# print "Proc(%d) Start..."%cpu
						result = pool.apply_async(cerate_reads_for_random_SNP_of_haploid_and_polyploid,(cpu,number_tmp,ref_start,ref,directional,Rate,max_err,output_ref,Fragment_length,Fragment_length_sigma,reads_length_mu,reads_length_sigma,SAM,samw,samc,PE,fq1,fq2,reads_r1,quals_r1,reads_r2,quals_r2,technology,chrom,out,position,index,qual_mu,qual_sigma,random_sequencing_errors,front_point,front_point_sigma,front_qual,front_qual_sigma,end_point,end_point_sigma,end_qual,end_qual_sigma,Dynamic_qual))
					pool.close()
					pool.join()
					if result.successful():
						print '%d CPUS run successful!' %CPUS
						cat_tmp(samw,samc,fq1,fq2,reads_r1,quals_r1,reads_r2,quals_r2)
					else:
						print 'Some CPUS not run successful!!'
					# cpu=1
					# cerate_reads_for_random_SNP_of_haploid_and_polyploid(cpu,number_tmp,ref,directional,Rate,max_err,output_ref,Fragment_length,Fragment_length_sigma,reads_length_mu,reads_length_sigma,SAM,samw,samc,PE,fq1,fq2,reads_r1,quals_r1,reads_r2,quals_r2,technology,chrom,out,position,index,qual_mu,qual_sigma,random_sequencing_errors,front_point,front_point_sigma,front_qual,front_qual_sigma,end_point,end_point_sigma,end_qual,end_qual_sigma,Dynamic_qual)
					print(time.strftime('\tend at: %Y-%m-%d %H:%M:%S'))
	elapsed = (time.time() - starttime)/3600.0
	print "Total time used: %f h" %elapsed


if __name__ == "__main__":
	main(sys.argv[1:])
