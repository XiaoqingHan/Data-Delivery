# -*- coding: utf-8 -*-

"""
Author: hanxiaoqing
Created in: 2022-05
"""

### packages
import os
import numpy as np
import pandas as pd
import json
import time
from datetime import datetime,date

### replace the unit to num
def remove_unit(dic,s):

	strings = "total_reads,mapped_reads,MID_filter_reads"
	if s in strings:
		alpha = ''.join([x for x in dic[s] if x.isalpha()])
		index = dic[s].find(alpha)
		num = dic[s][:index]
		if alpha == "K":
			dic[s] = int(float(num)*1000)
		elif alpha == "M":
			dic[s] = int(float(num)*1000000)
		elif alpha == "G":
			dic[s] = int(float(num)*1000000000)
		else:
			dic[s] = float(dic[s].split('(')[0])

		return dic[s]

### get info2
def para_info(para2):

	df2 = pd.DataFrame()
	res2 = {"TotalReads":[],"valid_CID_ratio(%)":[],"bin200_Mean_MID_per_spot":[],"Fraction_MID_in_spots_under_tissue(%)":[]}
	for i in para2:
		if isinstance(i,str) and os.path.exists(i):
			res_path = ''
			web = i.split('/')[-1][:-7]
			res_path += i+'/result/'+web+"/new_final_result.json"
			if os.path.exists(res_path):
				f3 = open(res_path,'r')
				f3 = json.load(f3)
				part1 = f3["1.Filter_and_Map"]["1.1.Adapter_Filter"]
				total,mapped,filter_ = 0,0,0
				for j in part1:
					j["total_reads"] = remove_unit(j,"total_reads")
					total += j["total_reads"]
					j["mapped_reads"] = remove_unit(j,"mapped_reads")
					mapped += j["mapped_reads"]
					j["MID_filter_reads"] = remove_unit(j,"MID_filter_reads")
					filter_ += j["MID_filter_reads"]
				cid_reads = mapped-filter_
				cid_ratio = '{:.2%}'.format(cid_reads/total)
				if "4.TissueCut" in f3.keys():
					part2 = f3["4.TissueCut"]["4.2.TissueCut_Bin_stat"][4]
					mid = part2["Mean_MID_per_spot"]
					part3 = f3["4.TissueCut"]["4.1.TissueCut_Total_Stat"][0]
					if "Fraction_MID_in_spots_under_tissue" in part3.keys():
						fraction = part3["Fraction_MID_in_spots_under_tissue"]				
					elif "Fraction_MID_in_Spots_Under_Tissue" in part3.keys():
						fraction = part3["Fraction_MID_in_Spots_Under_Tissue"]
				else:
					mid = np.nan
					fraction = np.nan
				"""
				# suitable for the situation that the fraction value is null: fraction = MID_counts/UNIQUE_READS
				elif "MID_counts" in part3.keys():
					umi = f3["2.Alignment"]["2.1.Input_Read"][5]
					unique = f3["2.Alignment"]["2.6.Filter_And_Deduplication"]["UNIQUE_READS"]
					print(umi,unique)					
				"""
				res2["TotalReads"] = total
				res2["valid_CID_ratio(%)"] = cid_ratio
				res2["bin200_Mean_MID_per_spot"] = mid
				res2["Fraction_MID_in_spots_under_tissue(%)"] = fraction
				df2 = df2.append(res2,ignore_index=True)
			else:
				res2["TotalReads"] = np.nan
				res2["valid_CID_ratio(%)"] = np.nan
				res2["bin200_Mean_MID_per_spot"] = np.nan
				res2["Fraction_MID_in_spots_under_tissue(%)"] = np.nan
				df2 = df2.append(res2,ignore_index=True)
	
		else:
			res2["TotalReads"] = np.nan
			res2["valid_CID_ratio(%)"] = np.nan
			res2["bin200_Mean_MID_per_spot"] = np.nan
			res2["Fraction_MID_in_spots_under_tissue(%)"] = np.nan
			df2 = df2.append(res2,ignore_index=True)
	
	return df2

### get info4
def error_log(para4):

	res4 = []
	for i in para4:
		if isinstance(i,str):
			res_path = ''
			res_path += i+'/logs/stderr'
			if os.path.exists(res_path):
				f4 = open(res_path)
				iter_f = iter(f4)
				text = ""
				for line in iter_f:
					text += line
				res4.append(text)
			else:
				res4.append(np.nan)
		else:
			res4.append(np.nan)

	return res4

### CID
def cids(cid):

	if isinstance(cid,str):
		cid = float(cid)
		if cid < 10:
			return "疑似SN号错误"
		elif 10 <= cid < 30:
			return "疑似混样"
	else:
		return np.nan

### fraction and meanMID
def frac_bin200(frac,bin200):

	if isinstance(frac,str) and isinstance(bin200,str):
		frac = float(frac)
		bin200 = bin200.replace('K','*1000')	
		bin200 = float(eval(str(bin200)))
		if bin200 < 5000 and frac >= 50:
			return "bin200 Mean MID per spot < 5k"
		elif bin200 >= 5000 and frac < 50:
			return "Fraction MID in spots under tissue < 50%"
		elif bin200 < 5000 and frac < 50:
			return "bin200 Mean MID < 5k & Fraction MID < 50%"
	else:
		return np.nan

### update info

def update(f_o,f_n):

	strings = ["errorType","Treatment","Remarks"]
	for i in f_n["ANALYSZ_ID"].values:
		if i in f_o["ANALYSZ_ID"].values:
			for j in strings:
				f_n.loc[f_n["ANALYSZ_ID"]==i,j] = f_o.loc[f_o["ANALYSZ_ID"]==i,j]
		else:
			f_n.loc[f_n["ANALYSZ_ID"]==i,"errorType"] = f_n.loc[f_n["ANALYSZ_ID"]==i,"errorType"]
		
	return f_n
