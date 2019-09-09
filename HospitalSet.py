#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 15:30:38 2018

@author: sidd
"""
import pandas as pd
import math
from sklearn import preprocessing
from sklearn.linear_model import LinearRegression
#import psycopg2 as ps
import pandasql as ps
from sklearn.linear_model import LogisticRegression
from sklearn.cross_validation import train_test_split
from sklearn import metrics
from sklearn.metrics import confusion_matrix

Main = pd.read_csv('Patient_data_kims_2017.csv')

Main = Main.dropna()
AgeR=[-math.inf,25,30,34,38,41,44,47,50,52,54,56,59,61,64,66,69,72,76,(math.inf)]
mhaemoR = [-math.inf,21,23.2,24.5,25.3,25.8,26.4,26.7,27.1,27.4,27.7,28,28.3,28.6,28.9,29.3,29.7,30.2,math.inf]
monocytesR = [-math.inf,5,6,7,8,9,10,12,math.inf]
wbc_countR=[-math.inf,5600,6500,7100,7600,8100,8500,9000,9500,10000,10600,11400,12100,
13000,14000,15300,16900,18700,21200,26200,math.inf]


def Ranges(mdata,flag):
    mdata= mdata.sort_values([str(mdata.columns[1])])
    if (flag==1):
        Dist_Col = pd.DataFrame((mdata.iloc[:,1]).unique().tolist())
        Dist_Col.columns = ['Dist_' + str(mdata.columns[1])]
        Dist_Col['TotalCount']=0
        Dist_Col['TargetCount']=0
        for x in range(0,len(Dist_Col)): 
            global targetsum        
            searchterm = Dist_Col[str(Dist_Col.columns[0])][x]
            hold= mdata[mdata.iloc[:,1] ==searchterm]
            Dist_Col['TotalCount'][x]=len(hold)
            Dist_Col['TargetCount'][x]=hold.iloc[:, 0].sum()
            print(x)
    if (flag==2):
        Dist_Col=pd.DataFrame((mdata.iloc[:,[0,1,2]]))
        Dist_Col=Dist_Col.sort_index()
    Dist_Col['Probability']=0
    Dist_Col['Percentage_Count']=0
    Dist_Col['Cumm_Perc']=0
    Dist_Col['Odds_Ratio']=0

    #print(Dist_Col.columns[0])
    global targetsum 
    Dist_Col['Percentage_Count']=((Dist_Col['TargetCount'])/(Dist_Col.TargetCount.sum()))*100
    Dist_Col['Probability']=(Dist_Col['TargetCount'])/(Dist_Col['TotalCount'])*100  
    Dist_Col['Odds_Ratio']=((Dist_Col['Probability']))/(100-(Dist_Col['Probability']))
    Dist_Col['Cumm_Perc'] = Dist_Col.Percentage_Count.cumsum()

    return (Dist_Col)
tempAge= Ranges((Main.iloc[:,[0,2]]),1)
tempMono= Ranges((Main.iloc[:,[0,4]]),1)
tempWBC = Ranges((Main.iloc[:,[0,5]]),1)
tempMH = Ranges((Main.iloc[:,[0,3]]),1)

def Binning(df,Ranges,flag):
    BinData = pd.DataFrame(columns = ['Bins','TargetCount','TotalCount','Probability','Odds_Ratio'])
    if (flag==1):
        for r in range(0,len(Ranges)-1):
            LR=Ranges[r]
            UR=Ranges[r+1]
            temp2= df.query(df.columns[0]+'>='+str(LR) +'&' + df.columns[0]+ '<='+str(UR))
            temp2 = temp2.reset_index(drop=True)
            rlabel=str(LR)+' - '+str(UR)
            TCount = temp2.TargetCount.sum()
            TOCount= temp2.TotalCount.sum()
            BinData=BinData.append({'Bins':rlabel,'TotalCount' : TOCount,'TargetCount':TCount,'Probability':(TCount/TOCount) , 'Odds_Ratio':((TCount/TOCount)/(1-(TCount/TOCount))) } , ignore_index=True)
    return (BinData)
AgesBand= Binning(tempAge,AgeR,1)
WBCBand= Binning(tempWBC,wbc_countR,1 )
MonoBand= Binning(tempMono,monocytesR,1)
MHBand= Binning(tempMH,mhaemoR,1 )

def FlagWithin(Col):
    for x in range(0,len(Col)):
        Col.iloc[:,0] = (Col.iloc[:,0]).replace([Col.iloc[x,0]], x)
    Col.columns =['Bin', 'TargetCount', 'TotalCount','Probs','Odds_Ratio']    
FlagWithin(AgesBand)
FlagWithin(WBCBand)
FlagWithin(MonoBand)
FlagWithin(MHBand)

def LinReg(df,dependent, header):
    lr = LinearRegression()
    X =df.loc[:,['Bin']]
    y =df.loc[:,dependent]
    lr.fit(X,y)
    df[header]=0
    df.loc[:,header] = lr.predict(X)
LinReg(AgesBand,'Odds_Ratio', 'P_Odds')
LinReg(WBCBand,'Odds_Ratio','P_Odds')
LinReg(MonoBand,'Odds_Ratio', 'P_Odds')
LinReg(MHBand,'Odds_Ratio','P_Odds')

#BinData = pd.DataFrame(columns = ['target_100', 'umrno', 'age', 'min_comp_bld_pic_p_mean_cell_haemo',
#       'max_comp_bld_pic_p_monocytes', 'max_comp_bld_pic_p_wbc_count']
def AssBinM(df,ranges,num,name,flag):
    if (flag==1):
        for x in range(0,len(ranges)-1):        
            LR = ranges[x]+.0001
            UR = ranges[x+1]+.0001
            temp2 = None
            temp2= Main.query(Main.columns[num]+'>='+str(LR) +'&' + Main.columns[num]+ '<='+str(UR))
            temp2 = temp2.reset_index(drop=True)
            temp2[name]=0
            temp2[name]=x
            if (x==0):
                Main2=temp2
            else:
                Main2 = Main2.append(temp2)
    if (flag==2):
        for x in range(0,len(ranges)-1):        
            LR = ranges[x]+.0001
            UR = ranges[x+1]+.0001
            temp2 = None
            temp2= df.query(df.columns[num]+'>='+str(LR) +'&' + df.columns[num]+ '<='+str(UR))
            temp2 = temp2.reset_index(drop=True)
            temp2[name]=0
            temp2[name]=x
            if (x==0):
                Main2=temp2
            else:
                Main2 = Main2.append(temp2)
    return (Main2)
Ass = AssBinM(Main,AgeR, 2, 'AgeBin',1)
Bass = AssBinM(Ass,mhaemoR, 3, 'MHBin',2)
Cass = AssBinM(Bass,monocytesR, 4, 'MonoBin',2)
Dass = AssBinM(Cass, wbc_countR, 5, 'WBCBin',2)

Dass.reset_index(drop=True)


#q1 = " SELECT P_Odds FROM AgesBand "
#Ages_Odds = ps.sqldf(q1, locals())
#
#q2 = " SELECT P_Odds FROM MHBand "
#MH_Odds = ps.sqldf(q2, locals())
#
#q3 = " SELECT P_Odds FROM MonoBand "
#Mono_Odds = ps.sqldf(q3, locals())
#
#q4 = " SELECT P_Odds FROM WBCBand "
#WBC_Odds = ps.sqldf(q4, locals())


Man = ps.sqldf("Select A.*, B.P_Odds As Age_Odds From Dass As A Left Join AgesBand As B on A.AgeBin=B.Bin;", locals())
Ban = ps.sqldf("Select A.*, B.P_Odds As MH_Odds From Man As A Left Join MHBand As B on A.MHBin=B.Bin;", locals())
Can = ps.sqldf("Select A.*, B.P_Odds As Mono_Odds From Ban As A Left Join MonoBand As B on A.MonoBin=B.Bin;", locals())
FinalMain = ps.sqldf("Select A.*, B.P_Odds As WBCB_Odds From Can As A Left Join WBCBand As B on A.WBCBin=B.Bin;", locals())
FinalMain = FinalMain.reset_index(drop=True)


#####Logarithmic Regression

independent = ['Age_Odds', 'MH_Odds','Mono_Odds', 'WBCB_Odds']
dependent = 'target_100'


train_x, test_x, train_y, test_y = train_test_split(FinalMain[independent], FinalMain[dependent], train_size=0.5)


print ("train_x size :: ", train_x.shape)
print ("train_y size :: ", train_y.shape)
 
print ("test_x size :: ", test_x.shape)
print ("test_y size :: ", test_y.shape)


logmodel = LogisticRegression()

Mod =  logmodel.fit(train_x,train_y)



def accuracy(trmodel, indep, dep):

    acc = trmodel.score(indep, dep)
    return acc


acc = accuracy(logmodel, train_x, train_y)
 
print ("Train Accuracy :: ", acc)

TR_pred = logmodel.predict(test_x)
Real = test_y

Confusion = metrics.confusion_matrix(Real, TR_pred)
print(Confusion)



test_y['target_100'].sum()
train_y['target_100'].sum()








