#!/usr/bin/python3
# coding: utf-8

#20210614プログレスバー追加
import os
import glob
import csv
import pandas as pd
import numpy as np
import math
import tqdm
import time

d_name="\\\\172.16.17.153\\rdc-share\\☆研開共有\\波形分析\\AI原因推定\\DAS_test\\"
o_name="\\\\172.16.17.153\\rdc-share\\☆研開共有\\波形分析\\AI原因推定\\CSV_test\\"


skip=len(d_name)
file_name=glob.glob(d_name+"*.das",recursive=True)
#print(file_name)
for j in tqdm.tqdm(range(len(file_name))):
  #  time.sleep(1.01)
  #  print(file_name[j])
    chk=file_name[j]
    chkd=chk[skip:skip+1]
    if chkd=="N" or chkd=="I":
        f = open(file_name[j],"rb")
    else:
        print("utf8")
        f = open(file_name[j],"r",encoding="utf-8")
    das_file=chk[skip:int(len(chk)-4)]
    #print(das_file)
    out_file_name=o_name+das_file+".csv"
    io=[]
    vo=[]
    data_num=16384
    header_len=1024
        
    
    hxnum=f.read()


    if hxnum[0:2]==b"IR":
        header_len=2048
        v_cof=0.125
        i_cof=0.03125
        
        for i in range(data_num):
        
            i_data=hxnum[header_len+i*12:header_len+i*12+2]
            v_data=hxnum[header_len+i*12+2:header_len+i*12+4]
            v_data=bin(int.from_bytes(v_data,byteorder='big'))
            i_data=bin(int.from_bytes(i_data,byteorder='big'))
            polarity_v = int(bin(int(v_data,0)>> 15),0)
                         
            if polarity_v == 1:
                v_data = (int(bin(int(v_data,0) ^ 0b1111111111111111), 0)+1) * (-1)
            else:  #1
                v_data= int(v_data,base=2)
            polarity_i = int(bin(int(i_data,0)>> 15),0)
            if polarity_i == 1:
                i_data =(int(bin(int(i_data,0) ^ 0b1111111111111111),0)+1)*(-1)
            else: #2
                i_data= int(i_data,base=2)


   
            io.append(i_data*i_cof)
            vo.append(v_data*v_cof)
    elif hxnum[0:2]=="AA" or hxnum[0:2]=="GS":

        data_len=4  
        v_cof=0.25
        i_cof=0.0625
      
        for i in range(data_num):
            v_data=hxnum[header_len+i*4:header_len+i*4+4]
  
            polarity_v = int(bin(int("0x"+v_data, 0) >> 15), 0)
            if polarity_v == 1:
                #print(v_data)
                v_data = (int(bin(int("0x"+v_data, 0) ^ 0b1111111111111111), 0)+1) * (-1)
            else:
                v_data= int("0x"+v_data,0)
       
            i_start=header_len+data_num*data_len
            i_data=hxnum[i_start+i*4:i_start+i*4+4]
            polarity_i = int(bin(int("0x"+i_data, 0) >> 15), 0)
            if polarity_i == 1:
                i_data = (int(bin(int("0x"+i_data, 0) ^ 0b1111111111111111), 0)+1) * (-1)
                   # i_data = (int(bin(int("0x"+i_data, 0) ^ 0b1111111111111111), 0)) * (-1)
            else:
                i_data= int("0x"+i_data,0)
            io.append(i_data*i_cof)
            vo.append(v_data*v_cof)
    else:  #DAS

        data_len=2  
        v_cof=2.0
        i_cof=0.5
     
        for i in range(data_num):
            v_data=hxnum[header_len+i*4:header_len+i*4+2]
            i_data=hxnum[header_len+i*4+2:header_len+i*4+4]
            #print("hx",v_data,i_data)
            v_data=bin(int.from_bytes(v_data,byteorder='little'))
            i_data=bin(int.from_bytes(i_data,byteorder='little'))
            polarity_v = int(bin(int(v_data,0)>> 15),0)
                         
            if polarity_v == 1:
                v_data = (int(bin(int(v_data,0) ^ 0b1111111111111111), 0)+1) * (-1)
            else:  #1
                v_data= int(v_data,base=2)
            polarity_i = int(bin(int(i_data,0)>> 15),0)
            if polarity_i == 1:
                i_data =(int(bin(int(i_data,0) ^ 0b1111111111111111),0)+1)*(-1)
            else: #2
                i_data= int(i_data,base=2)
      
             
            io.append(i_data*i_cof)
            vo.append(v_data*v_cof)
     
            
                
        
    outlist=[io,vo]
    name=["Ｉｏ(mA)","Ｖｏ(V)"]
    outframe=pd.DataFrame(outlist,index=name)
    outframe=outframe.T
    outframe.to_csv(out_file_name,encoding="utf-8_sig")
    f.close()
pass
  
