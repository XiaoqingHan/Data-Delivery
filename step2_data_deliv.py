# -*- coding: utf-8 -*-

"""
Author: hanxiaoqing
Created in: 2022-05
"""

import function
from function import *

### time monitor
start = time.time()

### get original info
#f1 = pd.read_table("test",delimiter='\t')
f1 = pd.read_table("out.txt",delimiter='\t')
#print(f1.head())
f1 = f1[['ANALYSZ_ID','sampleid','correctSampleid','autoName','rawdataParameter','specificParameter','owner','analysisStat','analysisDir'
,'createTime','finishTime']]
f2 = f1[f1['autoName'].isin(["spatialRNAvisualization_v2","spatialRNAvisualization_v3","spatialRNAvisualization_v3tmp","spatialRNAvisualization_cell","spatialRNAvisualization_develop"])]
#print(f2.head())
#print(f1.shape,f2.shape)

### info1
df1 = pd.DataFrame()
para1 = f2["specificParameter"]
for x in para1:
	res1 = {"Tissue":[],"Reference":[],"SlideArea":[],"ImagePath":''}
	y = list(x.split("|"))
	for i in range(len(y)):
		if y[i].split("=")[0] == "Tissue":
			res1["Tissue"] = y[i].split("=")[1]
		elif y[i].split("=")[0] == "referenceIndex":
			res1["Reference"] = y[i].split("=")[1]
		elif y[i].split("=")[0] == "SlideArea":
			res1["SlideArea"] = y[i].split("=")[1]
		elif y[i].split("=")[0] == "Image":
			if len(y[i].split("=")[1]) != 0:
				res1["ImagePath"] = y[i].split("=")[1]
			else:
				res1["ImagePath"] = np.nan
	df1 = df1.append(res1,ignore_index=True)
#print(df1.shape)

### info2
#print(para_info(f2["analysisDir"]).shape)

### info3
df3 = pd.DataFrame()
res3 = {"DeltaTime(h)":[]}
para3 = f2["createTime"]
para4 = f2["finishTime"]
for (t1,t2) in zip(para3,para4):
	if isinstance(t1,str) and isinstance(t2,str):
		t1,t2 = datetime.strptime(str(t1),"%Y-%m-%d %H:%M:%S"),datetime.strptime(str(t2),"%Y-%m-%d %H:%M:%S")
		hour1 = (t2-t1).seconds/3600
		hour2 = (t2-t1).days*24
		delta_time = round(hour1+hour2,5)
		res3["DeltaTime(h)"] = delta_time
		df3 = df3.append(res3,ignore_index=True)
	else:
		res3["DeltaTime(h)"] = np.nan
		df3 = df3.append(res3,ignore_index=True)
#print(df3.shape)

### merge all info into one dataframe
df4 = pd.DataFrame(error_log(f2["analysisDir"]))
df4.columns = ["errorLog"]
f4 = df1.join(para_info(f2["analysisDir"])).join(df3).join(df4)
f2.index = range(len(f2))
df = pd.concat([f2,f4],axis=1)
#print(df.shape)

### output of dataframe 
df["valid_CID_ratio(%)"] = df["valid_CID_ratio(%)"].astype(str).str.strip("%")
df["errorType"] = df["valid_CID_ratio(%)"].map(cids)
df["Fraction_MID_in_spots_under_tissue(%)"] = df["Fraction_MID_in_spots_under_tissue(%)"].astype(str).str.strip("%")
df["Warning"] = df.apply(lambda x: frac_bin200(x["Fraction_MID_in_spots_under_tissue(%)"],x["bin200_Mean_MID_per_spot"]),axis=1)
df["Treatment"] = np.nan
df["Remarks"] = np.nan
df["Monitor"] = np.nan
df["createTime"] = pd.to_datetime(df["createTime"])
df["finishTime"] = pd.to_datetime(df["finishTime"])
for i in range(len(df)):
	# total reads can be considerred as a factor: total reads <= 10g
	if (df.loc[i,"analysisStat"]=="running") & (df.loc[i,"DeltaTime(h)"]>48) & (int(df.loc[i,"SlideArea"])<=72):
		df.loc[i,"Monitor"] = "已运行超48h"
	else:
		df.loc[i,"Monitor"] = ""
for i in range(len(df)):
	if isinstance(df.loc[i,"errorLog"],str):
		if ".fq.gz' doesn't exist, quit now" in df.loc[i,"errorLog"]:
			df.loc[i,"errorType"] = "参数错误：命中无效路径"
		elif "Error code: SAW-A10125" in df.loc[i,"errorLog"]:
			df.loc[i,"errorType"] = "参数错误：UMI_Start_Pos错误"
		elif "Error code: SAW-A10183" in df.loc[i,"errorLog"]:
			df.loc[i,"errorType"] = "计算资源不足"
		elif "Killed" in df.loc[i,"errorLog"]:
			df.loc[i,"errorType"] = "计算资源不足"
		elif "Disk quota exceeded" in df.loc[i,"errorLog"]:
			df.loc[i,"errorType"] = "磁盘存储不足"
		elif "Error code: SAW-A10180" in df.loc[i,"errorLog"]:
			df.loc[i,"errorType"] = "磁盘存储不足或IO问题"
		elif "Error code: SAW-A10400" in df.loc[i,"errorLog"]:
			df.loc[i,"errorType"] = "外部信号引发退出"
		elif "Error code: SAW-A10150" in df.loc[i,"errorLog"]:
			df.loc[i,"errorType"] = "参数错误：fq文件错误"
		elif "Error code: SAW-A10169" in df.loc[i,"errorLog"]:
			df.loc[i,"errorType"] = "参数错误：fq文件重复"
		elif "Error code: SAW-A10171" in df.loc[i,"errorLog"]:
			df.loc[i,"errorType"] = "参数错误：mask文件错误"
	else:
		df.loc[i,"errorType"] = df.loc[i,"errorType"]
order = ["ANALYSZ_ID","sampleid","correctSampleid","autoName","owner","analysisStat","createTime","finishTime","DeltaTime(h)",\
        "Warning","errorType","Treatment","Remarks","Monitor","valid_CID_ratio(%)","Fraction_MID_in_spots_under_tissue(%)","bin200_Mean_MID_per_spot",\
        "errorLog","Reference","SlideArea","Tissue","TotalReads","rawdataParameter","specificParameter","analysisDir","ImagePath"]
df = df[order]
df = df.replace('nan','')
#print(df.shape)
df.to_excel("daily_new.xlsx",index=False)

### generate final file
if os.path.exists("daily_old.xlsx"):
	file_old = pd.read_excel("daily_old.xlsx",index=False)
	file_new = pd.read_excel("daily_new.xlsx",index=False)
	update(file_old,file_new).to_excel("daily_merge.xlsx",index=False)	
	today = date.today()
	print("The daily-update task of %s has done!" %today)
else:
	print("The output of spatial tasks until yesterday has done!")

### the end
end = time.time()
print("The consuming time is:",end-start)

