# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 16:52:46 2017

@author: Marcin Glowacki
"""
#Redshift estimation for multiple radio sources. 
#Provide a data file in same directory with WISE W1, W2, W3 magnitudes, corresponding errors and source name columns.

import os
import sys
from math import *
from numpy import *
from pylab import *
from matplotlib import rc, rcParams
rc('text', usetex=False)
rc('font',**{'family':'serif','serif':['serif'],'size':10})
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage,AnnotationBbox
from matplotlib.cbook import get_sample_data
from matplotlib._png import read_png
from scipy.stats import kde
from scipy import stats
import pandas as pd
rcParams.update({'font.size': 10})
import warnings
warnings.filterwarnings("ignore")

##################################################
#variables - change was desired. E.g. filename!

#INPUT FILE
#filename of sample to change
samplefile = 'test.csv'
#separation in the sample file. , for CSV, space for ASCISS, etc.
sep = ","
#fluxlim - only change if considering radio-bright objects!
###IMPORTANT - please label columns as W1, W1e, etc; and have a 'name' column for source name.
###If different column names are used, please update lines below (61-67)

#VARIABLES IN CODE
fluxlim = 0
#bin sizes. Minimise to increase computation time, but be aware of loss of accuracy
bin_no = 150;
step = (18.-8.)/bin_no;
bin_switch = True

redshift_min = 0.4 #in case you want to know the chance of object at a lower redshift than this value
redshift_max = 1.0 #in case you want to know the chance of object at a greater redshift
plot_on = True #False if you don't want to plot anything. Files saved to an output folder.
ioff()
##################################################
#important numbers! change limits as you see fit
infile = pd.read_csv(
        filepath_or_buffer=samplefile,    
        sep=sep)
df = pd.read_csv(
    filepath_or_buffer='jhyc_wise_spec.txt',
    sep=' ')

#removes invalide values.
infile = infile.loc[infile['W1'] > 0]
#enable following lines if you have 'null' values in your file
#infile = infile.loc[infile['W1'] != 'null']
#infile = infile.loc[infile['W2'] != 'null']
#infile = infile.loc[infile['W3'] != 'null']

W1_num = infile.columns.get_loc('W1')
try:
    W1e_num = infile.columns.get_loc('W1e')
except:
    try:
        W1e_num = infile.columns.get_loc('W1err')    
    except:
        "Error; please provide errors to the WISE magnitudes."
W2_num = infile.columns.get_loc('W2')
try:
    W2e_num = infile.columns.get_loc('W2e')
except:
    try:
        W2e_num = infile.columns.get_loc('W2err')    
    except:
        "Error; please provide errors to the WISE magnitudes."
W3_num = infile.columns.get_loc('W3')
try:
    W3e_num = infile.columns.get_loc('W3e')
except:
    try:
        W3e_num = infile.columns.get_loc('W3err')    
    except:
        "Error; please provide errors to the WISE magnitudes."
try:
    namenum = infile.columns.get_loc('name')
except:
    try:
        namenum = infile.columns.get_loc('Name')
    except:
        namenum = 0 #assume first column is the name.

W1 = infile.ix[:,W1_num].values
W1e = infile.ix[:,W1e_num].values
W2 = infile.ix[:,W2_num].values
W2e = infile.ix[:,W2e_num].values
W3 = infile.ix[:,W3_num].values
W3e = infile.ix[:,W3e_num].values
sname = infile.ix[:,namenum].values

chance_L = []
chance_H = []
chance_Q = []
highest = []

med_T = []
interval_T = []
med_T2 = []
interval_T2 = []
backup = []
backup2 = []
#for the chance of being at lower/higher redshift 
chance_minz = [] #chance at being less than redshift_min
chance_between = [] #chance at being between redshiftmin/max
chance_maxz = []#chance at being more than redshift_max
chance_minz2 = [] #chance at being less than redshift_min W2
chance_between2 = [] #chance at being between redshiftmin/max W2
chance_maxz2 = []#chance at being more than redshift_max W2

#LARGESS file
w1 = df.ix[:,27].values
w2 = df.ix[:,29].values
w3 = df.ix[:,31].values
w1e = df.ix[:,28].values
w2e = df.ix[:,30].values
w3e = df.ix[:,32].values
sclass = df.ix[:,20].values
flux_f = df.ix[:,11].values
flux_n = df.ix[:,13].values          

##################################################
# stage one - set up to determine the chance of being a LERG, HERG or QSO
nbins = 100;
(x,y) = (w2-w3,w1-w2)
cc = 0;
name = 'LERG'
x_sub = [];
y_sub = [];
#df.loc[(df['BESTCLASS']=='LERG') & ((df['FIRST_TOTINTFLUX'] > fluxlim) | (df['NVSS_TOTINTFLUX'] > fluxlim)),['BESTCLASS']]
for i in sclass:
    if i == name and (float(flux_f[cc]) > fluxlim or float(flux_n[cc]) > fluxlim):
        x_sub.append(x[cc])
        y_sub.append(y[cc])
    cc += 1;
#scatter(x_sub,y_sub,facecolor='r',marker='o',s=18,zorder=101,alpha=1.0,linewidth=0.5)
k_L = kde.gaussian_kde([asarray(x_sub),asarray(y_sub)])
a = asarray(x_sub)
#xi, yi = np.mgrid[min(x_sub):max(a[a<10]):nbins*1j, min(y_sub):max(y_sub):nbins*1j]#not max(x_sub) here, too large
xi, yi = np.mgrid[-1:6:nbins*1j, -1:3:nbins*1j]
zi_L = k_L(np.vstack([xi.flatten(), yi.flatten()]))
tot_l = len(x_sub)

name = 'HERG'
cc = 0
x_sub = [];
y_sub = [];
for i in sclass:
    if i == name and (float(flux_f[cc]) > fluxlim or float(flux_n[cc]) > fluxlim):
        x_sub.append(x[cc])
        y_sub.append(y[cc])
    cc += 1;
#scatter(x_sub,y_sub,facecolor='b',marker='o',s=18,zorder=101,alpha=1.0,linewidth=0.5)
k_H = kde.gaussian_kde([asarray(x_sub),asarray(y_sub)])
a = asarray(x_sub)
#xi, yi = np.mgrid[min(x_sub):max(a[a<10]):nbins*1j, min(y_sub):max(y_sub):nbins*1j]#not max(x_sub) here, too large
xi, yi = np.mgrid[-1:6:nbins*1j, -1:3:nbins*1j]
zi_H = k_H(np.vstack([xi.flatten(), yi.flatten()]))
tot_h = len(x_sub)

name = 'AeB'
x_sub = [];
y_sub = [];
cc = 0
for i in sclass:
    if i == name and (float(flux_f[cc]) > fluxlim or float(flux_n[cc]) > fluxlim):
        x_sub.append(x[cc])
        y_sub.append(y[cc])
    cc += 1;
#scatter(x_sub,y_sub,facecolor='g',marker='o',s=18,zorder=101,alpha=1.0,linewidth=0.5)
k_Q = kde.gaussian_kde([asarray(x_sub),asarray(y_sub)])
a = asarray(x_sub)
#xi, yi = np.mgrid[min(x_sub):max(a[a<10]):nbins*1j, min(y_sub):max(y_sub):nbins*1j]#not max(x_sub) here, too large
xi, yi = np.mgrid[-1:6:nbins*1j, -1:3:nbins*1j]
zi_Q = k_Q(np.vstack([xi.flatten(), yi.flatten()]))
tot_q = len(x_sub)

#loop through each source
count = 0;
for (W1_in,W1_er,W2_in,W2_er,W3_in,W3_er,n) in zip(W1,W1e,W2,W2e,W3,W3e,sname):
    #THE LOOPING BEGINS
    count +=1;
    #firstly, adjust errors. If they are non-sensical or 0.0 (aka upper limit), adjust to some sensible number
    W1_in = float(W1_in)
    W2_in = float(W2_in)
    W3_in = float(W3_in)
#    try:
#        W1_er = float(W1_er)*3
#    except:
#        W1_er = 0.3
#    try:
#        W2_er = float(W2_er)*3
#    except:
#        W2_er = 0.3
#    try:
#        W3_er = float(W3_er)*3
#    except:
#        W3_er = 0.3
    if W1_er > 0.5:
        W1_er = 0.5
    if W2_er > 0.5:
        W2_er = 0.5
    pt = [round(W2_in-W3_in,4),round(W1_in-W2_in,4)]
    ##################################################

    if plot_on == True:
        fig, axes = subplots(ncols=2, nrows=2, sharex=True, sharey=True)
    
    if plot_on == True:
        axes[0,0].pcolormesh(xi, yi, zi_L.reshape(xi.shape))
        axes[1,1].pcolormesh(xi, yi, zi_L.reshape(xi.shape),alpha=0.7)
    percent_L = float(k_L(pt)/max(zi_L))*100
                     
    
    if plot_on == True:
        axes[1,0].pcolormesh(xi, yi, zi_H.reshape(xi.shape))
        axes[1,1].pcolormesh(xi, yi, zi_H.reshape(xi.shape),alpha=0.5)
    percent_H = float(k_H(pt)/max(zi_H))*100

    
    if plot_on == True:
        axes[0,1].pcolormesh(xi, yi, zi_Q.reshape(xi.shape))    
        axes[1,1].pcolormesh(xi, yi, zi_Q.reshape(xi.shape),alpha=0.3)
    
    percent_Q = float(k_Q(pt)/max(zi_Q))*100
    
    if plot_on == True:
        ylim(-1,3)
        xlim(-1,6)
        axes[0,0].set_title('LERG')
        axes[1,0].set_title('HERG')
        axes[0,1].set_title('AeB/QSO')
        axes[1,1].set_title('Combined')
        
        axes[0,0].scatter(pt[0],pt[1], facecolor='r',marker='D',s=30,linewidth=0.5)
        axes[1,0].scatter(pt[0],pt[1], facecolor='r',marker='D',s=30,linewidth=0.5)
        axes[0,1].scatter(pt[0],pt[1], facecolor='r',marker='D',s=30,linewidth=0.5)
        axes[1,1].scatter(pt[0],pt[1], facecolor='r',marker='D',s=30,linewidth=0.5)
        axes[1,0].set_xlabel('W2-W3 (mag)')
        axes[1,0].set_ylabel('W1-W2 (mag)')
        axes[1,1].set_xlabel('W2-W3 (mag)')
        axes[0,0].set_ylabel('W1-W2 (mag)')
        fname = 'output/'+str(n)+'_WISE.png'
        fig.savefig(fname,bbox_inches='tight')
        
    #following is for rescaling of problematic positions and/or cases where you can't really tell by wise colour!
    if percent_L < 0.0001:
        percent_L = 0
    if percent_H < 0.0001:
        percent_H = 0
    if percent_Q < 0.0001:
        percent_Q = 0
    #so unlikely it's practically 0. 
    if percent_L+percent_H+percent_Q > 0.1:
        tot = percent_L+percent_H+percent_Q
        prob_L = round(percent_L/tot*100,2);
        prob_H = round(percent_H/tot*100,2);
        prob_Q = round(percent_Q/tot*100,2);
        if max(prob_L,prob_H,prob_Q) == prob_L:
            mostlikely = 'LERG'
        if max(prob_L,prob_H,prob_Q) == prob_H:
            mostlikely = 'HERG'
        if max(prob_L,prob_H,prob_Q) == prob_Q:
            mostlikely = 'QSO'
            
    else:
        tot = 0;
        prob_L = 33.333;
        prob_H = 33.333;
        prob_Q = 33.333;
        mostlikely = 'N/A'
        #Very unusual WISE colours results in this.
        
    chance_L.append(prob_L)
    chance_H.append(prob_H)
    chance_Q.append(prob_Q)
    highest.append(mostlikely)
    #as QSOs span a wider redshift distribution, we can encforce more bins in plotting.
    #bin_switch = True necessary. 
    if (prob_Q > 10 or W1_in > 16.0) and bin_switch == True:
        bin_no = 150;
        step = (18.-8.)/bin_no;
        #make the switch off if you do not want this!
    ##################################################
    #generating the redshift probability distributions
    
    z_LERG = (df.loc[(df['BESTCLASS']=='LERG') & 
                     ((df['FIRST_TOTINTFLUX'] > fluxlim) | (df['NVSS_TOTINTFLUX'] > fluxlim)),['Z']])
    z_HERG = (df.loc[(df['BESTCLASS']=='HERG') & 
                     ((df['FIRST_TOTINTFLUX'] > fluxlim) | (df['NVSS_TOTINTFLUX'] > fluxlim)),['Z']])
    z_QSO = (df.loc[(df['BESTCLASS']=='AeB') & 
                    ((df['FIRST_TOTINTFLUX'] > fluxlim) | (df['NVSS_TOTINTFLUX'] > fluxlim)),['Z']])
    
    binmax = 0.975;
    if prob_Q > 10.0 or W1_in > 16.0:
        binmax = max(asarray(z_QSO)) + 0.025;
    bins2 = linspace(0.00,binmax,bin_no)
    binwidth = bins2[1]-bins2[0]
    WISErange = arange(8,18,step)
    probW1 = stats.norm.pdf(WISErange,W1_in,W1_er)*step
    probW2 = stats.norm.pdf(WISErange,W2_in,W2_er)*step
    W1diff = abs(WISErange-W1_in)
    W2diff = abs(WISErange-W2_in)
    
    #setup for extrapolation
    P_L = [0]*bin_no
    P_H = P_L
    P_Q = P_L
    PW2_L = P_L
    PW2_H = P_L
    PW2_Q = P_L
    
    means_LW1 = []
    stds_LW1 = []
    means_LW2 = []
    stds_LW2 = []
    means_HW1 = []
    stds_HW1 = []
    means_HW2 = []
    stds_HW2 = []
    means_QW1 = []
    stds_QW1 = []
    means_QW2 = []
    stds_QW2 = []
    for i in WISErange:
        #LERGs, W1 first
        zdistW1 = z_LERG.loc[(abs(df['W1']-i) < (step/2)),['Z']]
        mm,ss = mean(asarray(zdistW1)),std(asarray(zdistW1))
        if (i > 13.0 and i < 15.5) and isnan(mm) == False:
            means_LW1.append(mm)
            stds_LW1.append(ss)
        #W2
        zdistW2 = z_LERG.loc[(abs(df['W2']-i) < (step/2)),['Z']]
        mm,ss = mean(asarray(zdistW2)),std(asarray(zdistW2))
        if (i > 13.0 and i < 15.5) and isnan(mm) == False:
            means_LW2.append(mm)
            stds_LW2.append(ss)
        #HERGs:
        zdistW1 = z_HERG.loc[(abs(df['W1']-i) < (step/2)),['Z']]
        mm,ss = mean(asarray(zdistW1)),std(asarray(zdistW1))
        if (i > 13.0 and i < 15.5) and isnan(mm) == False:
            means_HW1.append(mm)
            stds_HW1.append(ss)
        zdistW2 = z_HERG.loc[(abs(df['W2']-i) < (step/2)),['Z']]
        mm,ss = mean(asarray(zdistW2)),std(asarray(zdistW2))
        if (i > 13.0 and i < 15.5) and isnan(mm) == False:
            means_HW2.append(mm)
            stds_HW2.append(ss)
        #QSOs
        zdistW1 = z_QSO.loc[(abs(df['W1']-i) < (step/2)),['Z']]
        mm,ss = mean(asarray(zdistW1)),std(asarray(zdistW1))
        if (i > 13.0 and i < 15.5) and isnan(mm) == False:
            means_QW1.append(mm)
            stds_QW1.append(ss)
        zdistW2 = z_QSO.loc[(abs(df['W2']-i) < (step/2)),['Z']]
        mm,ss = mean(asarray(zdistW2)),std(asarray(zdistW2))
        if (i > 13.0 and i < 15.5) and isnan(mm) == False:
            means_QW2.append(mm)
            stds_QW2.append(ss)
    #now to fit to these ranges:
    xx = WISErange[where(abs(WISErange - 14.25)<1.25)]
    fitLW1 = polyfit(xx,means_LW1,1)
    fitLW2 = polyfit(xx,means_LW2,1)
    fitHW1 = polyfit(xx,means_HW1,1)
    fitHW2 = polyfit(xx,means_HW2,1)
    fitQW1 = polyfit(xx,means_QW1,1)
    fitQW2 = polyfit(xx,means_QW2,1)
    sfitLW1 = polyfit(xx,stds_LW1,1)
    sfitLW2 = polyfit(xx,stds_LW2,1)
    sfitHW1 = polyfit(xx,stds_HW1,1)
    sfitHW2 = polyfit(xx,stds_HW2,1)
    sfitQW1 = polyfit(xx,stds_QW1,1)
    sfitQW2 = polyfit(xx,stds_QW2,1)
    
    pdfsW1 = []
    pdfsW2 = []
    #LERGS FIRST
    for i in WISErange:
        #i = WISE bin; remember step = 0,5, so 0.25 width
        zdistW1 = z_LERG.loc[(abs(df['W1']-i) < (step/2)),['Z']]
        mm,ss = mean(asarray(zdistW1)),std(asarray(zdistW1))
        if (i > 15.5):
            mm = fitLW1[0]*i+fitLW1[1]
            ss = sfitLW1[0]*i+sfitLW1[1]
        if i < 13.0 and isnan(mm) == True:
            #for the very bright case - as the fits will either curve back up
            #or for linear case go to negative values of redshift (nonsense)
            #we set to sensible values instead
            #numbers based on typical values observed for LARGESS survey
            #larger scatter to account for lack of information
            mm = uniform(0.003,0.02)
            ss = uniform(0.01,0.03)
        if ss == 0.0:
            mm = uniform(0.003,0.02)
            ss = uniform(0.01,0.03)
        if ss > mm:
            ss = mm
            #to reduce chance of a distribution suggesting a negative redshift
        pdf = stats.norm.pdf(bins2,mm,ss)*binwidth
        if isnan(pdf[0]) == True:
            pdf = [0]*bin_no
        if sum(pdf)>0:
            pdf = pdf/sum(pdf)
        pdfsW1.append(pdf)
        #same but for W2
        zdistW2 = z_LERG.loc[(abs(df['W2']-i) < (step/2)),['Z']]
        mm,ss = mean(asarray(zdistW2)),std(asarray(zdistW2))
        if (i > 15.5):
            mm = fitLW2[0]*i+fitLW2[1]
            ss = sfitLW2[0]*i+sfitLW2[1]
        if i < 13.0 and isnan(mm) == True:
            mm = uniform(0.003,0.02)
            ss = uniform(0.01,0.03)
        if ss == 0.0:
            mm = uniform(0.003,0.02)
            ss = uniform(0.01,0.03)
        if ss > mm:
            ss = mm
        pdf = stats.norm.pdf(bins2,mm,ss)*binwidth
        if isnan(pdf[0]) == True:
            pdf = [0]*bin_no
        if sum(pdf)>0:
            pdf = pdf/sum(pdf)
        pdfsW2.append(pdf)
    
    cc = 0
    for i in pdfsW1:
        try:
            P_L = P_L + i*probW1[cc]
        except:
            P_L = P_L
        #print aaa
        cc+=1
    cc = 0
    for i in pdfsW2:
        try:
            PW2_L = PW2_L + i*probW2[cc]
        except:
            PW2_L = PW2_L
        #print aaa
        cc+=1
    #SAME BUT FOR HERGS
    pdfsW1 = []
    pdfsW2 = []
    for i in WISErange:
        #i = WISE bin; remember step = 0,5, so 0.25 width
        zdistW1 = z_HERG.loc[(abs(df['W1']-i) < (step/2)),['Z']]
        mm,ss = mean(asarray(zdistW1)),std(asarray(zdistW1))
        if (i > 15.5):
            mm = fitHW1[0]*i+fitHW1[1]
            ss = sfitHW1[0]*i+sfitHW1[1]
        if i < 13.0 and isnan(mm) == True:
            mm = uniform(0.02,0.05)
            ss = uniform(0.01,0.2)
        if ss == 0.0:
            mm = uniform(0.02,0.05)
            ss = uniform(0.01,0.2)
        if ss > mm:
            ss = mm
        pdf = stats.norm.pdf(bins2,mm,ss)*binwidth
        if isnan(pdf[0]) == True:
            pdf = [0]*bin_no
        if sum(pdf)>0:
            pdf = pdf/sum(pdf)
        pdfsW1.append(pdf)
        #same but for W2
        zdistW2 = z_HERG.loc[(abs(df['W2']-i) < (step/2)),['Z']]
        mm,ss = mean(asarray(zdistW2)),std(asarray(zdistW2))
        if (i > 15.5):
            mm = fitHW2[0]*i+fitHW2[1]
            ss = sfitHW2[0]*i+sfitHW2[1]
        if i < 13.0 and isnan(mm) == True:
            mm = uniform(0.02,0.05)
            ss = uniform(0.01,0.2)
        if ss == 0.0:
            mm = uniform(0.02,0.05)
            ss = uniform(0.01,0.2)
        if ss > mm:
            ss = mm
        pdf = stats.norm.pdf(bins2,mm,ss)*binwidth
        if isnan(pdf[0]) == True:
            pdf = [0]*bin_no
        if sum(pdf)>0:
            pdf = pdf/sum(pdf)
        pdfsW2.append(pdf)
    
    cc = 0
    for i in pdfsW1:
        try:
            P_H = P_H + i*probW1[cc]
        except:
            P_H = P_H
        #print aaa
        cc+=1
    cc = 0
    for i in pdfsW2:
        try:
            PW2_H = PW2_H + i*probW2[cc]
        except:
            PW2_H = PW2_H
        #print aaa
        cc+=1
        
    pdfsW1 = []
    pdfsW2 = []
    for i in WISErange:
        #i = WISE bin; remember step = 0,5, so 0.25 width
        zdistW1 = z_QSO.loc[(abs(df['W1']-i) < (step/2)),['Z']]
        mm,ss = mean(asarray(zdistW1)),std(asarray(zdistW1))
        if (i > 15.5 and (isnan(mm) == True or mm == 0 or ss == 0)):
            mm = fitQW1[0]*i+fitQW1[1]
            ss = sfitQW1[0]*i+sfitQW1[1]
        if i < 13.0 and isnan(mm) == True:
            mm = uniform(0.02,0.1)
            ss = uniform(0.01,0.2)
        if ss == 0.0:
            mm = uniform(0.02,0.1)
            ss = uniform(0.01,0.2)
        if ss > mm:
            ss = mm
        pdf = stats.norm.pdf(bins2,mm,ss)*binwidth
        if isnan(pdf[0]) == True:
            pdf = [0]*bin_no
        if sum(pdf)>0:
            pdf = pdf/sum(pdf)
        pdfsW1.append(pdf)
        #same but for W2
        zdistW2 = z_QSO.loc[(abs(df['W2']-i) < (step/2)),['Z']]
        mm,ss = mean(asarray(zdistW2)),std(asarray(zdistW2))
        if (i > 15.5 and (isnan(mm) == True or mm == 0 or ss == 0)):
            mm = fitQW2[0]*i+fitQW2[1]
            ss = sfitQW2[0]*i+sfitQW2[1]   
        if i < 13.0 and isnan(mm) == True:
            mm = uniform(0.02,0.1)
            ss = uniform(0.01,0.2)
        if ss == 0.0:
            mm = uniform(0.02,0.1)
            ss = uniform(0.01,0.2)
        if ss > mm:
            ss = mm
        pdf = stats.norm.pdf(bins2,mm,ss)*binwidth
        if isnan(pdf[0]) == True:
            pdf = [0]*bin_no
        if sum(pdf)>0:
            pdf = pdf/sum(pdf)
        pdfsW2.append(pdf)
    
    cc = 0
    for i in pdfsW1:
        try:
            P_Q = P_Q + i*probW1[cc]
        except:
            P_Q = P_Q
        #print aaa
        cc+=1
    cc = 0
    for i in pdfsW2:
        try:
            PW2_Q = PW2_Q + i*probW2[cc]
        except:
            PW2_Q = PW2_Q
        #print aaa
        cc+=1
    succ = 1
    succ2 = 1
    try:
        P_Total = P_L*prob_L/100+P_H*prob_H/100+P_Q*prob_Q/100
        if sum(P_Total)<0.25:
            if mostlikely == 'LERG':
                z_LERG_w1 = z_LERG.loc[(abs(df['W1']-W1_in) < W1_er),['Z']]
                m_LW1,sig_LW1 = mean(asarray(z_LERG_w1)), std(asarray(z_LERG_w1))
                P_Total = stats.norm.pdf(bins2, m_LW1, sig_LW1)*binwidth
                succ = 0
                if sum(P_Total)<0.25:
                    P_Total = [0]*bin_no
            if mostlikely == 'HERG':
                z_HERG_w1 = z_HERG.loc[(abs(df['W1']-W1_in) < W1_er),['Z']]
                m_HW1,sig_HW1 = mean(asarray(z_HERG_w1)), std(asarray(z_HERG_w1))
                P_Total = stats.norm.pdf(bins2, m_HW1, sig_HW1)*binwidth
                succ = 0
                if sum(P_Total)<0.25:
                    P_Total = [0]*bin_no
                
            if mostlikely == 'QSO':
                z_QSO_w1 = z_QSO.loc[(abs(df['W1']-W1_in) < W1_er),['Z']]
                m_QW1,sig_QW1 = mean(asarray(z_QSO_w1)), std(asarray(z_QSO_w1))
                P_Total = stats.norm.pdf(bins2, m_QW1, sig_QW1)*binwidth
                succ = 0
                if sum(P_Total)<0.25:
                    P_Total = [0]*bin_no
            if mostlikely == 'N/A':
                P_Total = [0]*bin_no
                succ = 0

    except:
        try:
            succ = 0
            if mostlikely == 'LERG':
                z_LERG_w1 = z_LERG.loc[(abs(df['W1']-W1_in) < W1_er),['Z']]
                m_LW1,sig_LW1 = mean(asarray(z_LERG_w1)), std(asarray(z_LERG_w1))
                P_Total = stats.norm.pdf(bins2, m_LW1, sig_LW1)*binwidth
                if sum(P_Total)<0.25:
                    P_Total = [0]*bin_no
            if mostlikely == 'HERG':
                z_HERG_w1 = z_HERG.loc[(abs(df['W1']-W1_in) < W1_er),['Z']]
                m_HW1,sig_HW1 = mean(asarray(z_HERG_w1)), std(asarray(z_HERG_w1))
                P_Total = stats.norm.pdf(bins2, m_HW1, sig_HW1)*binwidth
                if sum(P_Total)<0.25:
                    P_Total = [0]*bin_no                
            if mostlikely == 'QSO':
                z_QSO_w1 = z_QSO.loc[(abs(df['W1']-W1_in) < W1_er),['Z']]
                m_QW1,sig_QW1 = mean(asarray(z_QSO_w1)), std(asarray(z_QSO_w1))
                P_Total = stats.norm.pdf(bins2, m_QW1, sig_QW1)*binwidth
                if sum(P_Total)<0.25:
                    P_Total = [0]*bin_no
            if mostlikely == 'N/A':
                P_Total = [0]*bin_no
                succ = 0 
        except:
            P_Total = [0]*bin_no
            succ = 0
                      
    #same but for W2 now
    try:
        PW2_Total = PW2_L*prob_L/100+PW2_H*prob_H/100+PW2_Q*prob_Q/100
        if sum(PW2_Total)<0.25:
            if mostlikely == 'LERG':
                z_LERG_w2 = z_LERG.loc[(abs(df['W2']-W2_in) < W2_er),['Z']]
                m_LW2,sig_LW2 = mean(asarray(z_LERG_w2)), std(asarray(z_LERG_w2))
                PW2_Total = stats.norm.pdf(bins2, m_LW2, sig_LW2)*binwidth
                succ = 0
                if sum(PW2_Total)<0.25:
                    PW2_Total = [0]*bin_no
            if mostlikely == 'HERG':
                z_HERG_w2 = z_HERG.loc[(abs(df['W2']-W2_in) < W2_er),['Z']]
                m_HW2,sig_HW2 = mean(asarray(z_HERG_w2)), std(asarray(z_HERG_w2))
                PW2_Total = stats.norm.pdf(bins2, m_HW2, sig_HW2)*binwidth
                succ = 0
                if sum(PW2_Total)<0.25:
                    PW2_Total = [0]*bin_no
                
            if mostlikely == 'QSO':
                z_QSO_w2 = z_QSO.loc[(abs(df['W2']-W2_in) < W2_er),['Z']]
                m_QW2,sig_QW2 = mean(asarray(z_QSO_w2)), std(asarray(z_QSO_w2))
                PW2_Total = stats.norm.pdf(bins2, m_QW2, sig_QW2)*binwidth
                succ = 0
                if sum(PW2_Total)<0.25:
                    PW2_Total = [0]*bin_no
            if mostlikely == 'N/A':
                PW2_Total = [0]*bin_no
                succ = 0

    except:
        try:
            succ = 0
            if mostlikely == 'LERG':
                z_LERG_w2 = z_LERG.loc[(abs(df['W2']-W2_in) < W2_er),['Z']]
                m_LW2,sig_LW2 = mean(asarray(z_LERG_w2)), std(asarray(z_LERG_w2))
                PW2_Total = stats.norm.pdf(bins2, m_LW2, sig_LW2)*binwidth
                if sum(PW2_Total)<0.25:
                    PW2_Total = [0]*bin_no
            if mostlikely == 'HERG':
                z_HERG_w2 = z_HERG.loc[(abs(df['W2']-W2_in) < W2_er),['Z']]
                m_HW2,sig_HW2 = mean(asarray(z_HERG_w2)), std(asarray(z_HERG_w2))
                PW2_Total = stats.norm.pdf(bins2, m_HW2, sig_HW2)*binwidth
                if sum(PW2_Total)<0.25:
                    PW2_Total = [0]*bin_no                
            if mostlikely == 'QSO':
                z_QSO_w2 = z_QSO.loc[(abs(df['W2']-W2_in) < W2_er),['Z']]
                m_QW2,sig_QW2 = mean(asarray(z_QSO_w2)), std(asarray(z_QSO_w2))
                PW2_Total = stats.norm.pdf(bins2, m_QW2, sig_QW2)*binwidth
                if sum(PW2_Total)<0.25:
                    PW2_Total = [0]*bin_no
            if mostlikely == 'N/A':
                PW2_Total = [0]*bin_no
                succ = 0 
        except:
            PW2_Total = [0]*bin_no
            succ = 0
    
    #calculation of median, confidence interval values. 
    try:
        med_T.append(round(mean([min(bins2[where(cumsum(P_Total)>=(0.5))]),max(bins2[where(cumsum(P_Total)<=(0.5))])]),4))
    except:
        med_T.append(-999)
        #something went wrong for the source; this value indicates this in the output file.
    #here the distributions are usually skewed, so assuming a normal distribution would be foolish
    #we just calculate the intervals for the 68/95/99.7% CIs ourselves via cumsum.
    try:
        sig1_low = round(mean([min(bins2[where(cumsum(P_Total)>=(0.5-0.3413))]),max(bins2[where(cumsum(P_Total)<=(0.5-0.3413))])]),4)
    except:
        sig1_low = round(min(bins2),4)
    try:
        sig2_low = round(mean([min(bins2[where(cumsum(P_Total)>=(0.5-0.3413-0.136))]),max(bins2[where(cumsum(P_Total)<=(0.5-0.3413-0.136))])]),4)
    except:
        sig2_low = round(min(bins2),4)            
    try:
        sig3_low = round(mean([min(bins2[where(cumsum(P_Total)>=(0.5-0.3413-0.136-0.0214))]),max(bins2[where(cumsum(P_Total)<=(0.5-0.3413-0.136-0.0214))])]),4)
    except:
        sig3_low = round(min(bins2),4)          
    try:
        sig1_high = round(mean([min(bins2[where(cumsum(P_Total)>=(0.5+0.3413))]),max(bins2[where(cumsum(P_Total)<=(0.5+0.3413))])]),4)
    except:
        sig1_hig = round(max(bins2),4)    
    try:
        sig2_high = round(mean([min(bins2[where(cumsum(P_Total)>=(0.5+0.3413+0.136))]),max(bins2[where(cumsum(P_Total)<=(0.5+0.3413+0.136))])]),4)
    except:
        sig2_high = round(max(bins2),4) 
    try:
        sig3_high = round(mean([min(bins2[where(cumsum(P_Total)>=(0.5+0.3413+0.136+0.0214))]),max(bins2[where(cumsum(P_Total)<=(0.5+0.3413+0.136+0.0214))])]),4)
    except:
        sig3_high = round(max(bins2),4)
    store = ([sig1_low,sig2_low,sig3_low],[sig1_high,sig2_high,sig3_high])
    interval_T.append(store)
    
    #now same but for W2!
    try:
        med_T2.append(round(mean([min(bins2[where(cumsum(PW2_Total)>=(0.5))]),max(bins2[where(cumsum(PW2_Total)<=(0.5))])]),4))
    except:
        med_T2.append(-999)
        #this suggests that we have empty sets! too high/low a W1 (or W2) magnitude to be useable.
        #can extend this code later
    #here the distributions are usually skewed, so assuming a normal distribution would be foolish
    #he just calculate the intervals for 1,2,3 sigma ourselves via cumsum.
    try:
        sig1_low = round(mean([min(bins2[where(cumsum(PW2_Total)>=(0.5-0.3413))]),max(bins2[where(cumsum(PW2_Total)<=(0.5-0.3413))])]),4)
    except:
        sig1_low = round(min(bins2),4)    
    try:
        sig2_low = round(mean([min(bins2[where(cumsum(PW2_Total)>=(0.5-0.3413-0.136))]),max(bins2[where(cumsum(PW2_Total)<=(0.5-0.3413-0.136))])]),4)
    except:
        sig2_low = round(min(bins2),4)           
    try:
        sig3_low = round(mean([min(bins2[where(cumsum(PW2_Total)>=(0.5-0.3413-0.136-0.0214))]),max(bins2[where(cumsum(PW2_Total)<=(0.5-0.3413-0.136-0.0214))])]),4)
    except:
        sig3_low = round(min(bins2),4)             
    try:
        sig1_high = round(mean([min(bins2[where(cumsum(PW2_Total)>=(0.5+0.3413))]),max(bins2[where(cumsum(PW2_Total)<=(0.5+0.3413))])]),4)
    except:
        sig1_high = round(max(bins2),4)
    try:
        sig2_high = round(mean([min(bins2[where(cumsum(PW2_Total)>=(0.5+0.3413+0.136))]),max(bins2[where(cumsum(PW2_Total)<=(0.5+0.3413+0.136))])]),4)
    except:
        sig2_high = round(max(bins2),4)
    try:
        sig3_high = round(mean([min(bins2[where(cumsum(PW2_Total)>=(0.5+0.3413+0.136+0.0214))]),max(bins2[where(cumsum(PW2_Total)<=(0.5+0.3413+0.136+0.0214))])]),4)
    except:
        sig3_high = round(max(bins2),4)
    store2 = ([sig1_low,sig2_low,sig3_low],[sig1_high,sig2_high,sig3_high])
    interval_T2.append(store2)
    if sum(P_Total)>0.25:
        chance_minz.append(round(sum(P_Total[where(bins2 < redshift_min)])/sum(P_Total)*100,2))
        chance_maxz.append(round(sum(P_Total[where(bins2 > redshift_max)])/sum(P_Total)*100,2))
        chance_between.append(round(sum(P_Total[where((bins2 > redshift_min) & (bins2 < redshift_max))])/sum(P_Total)*100,2))
    else:
        chance_minz.append(0)
        chance_maxz.append(0)
        chance_between.append(0)
    if sum(PW2_Total)>0.25:
        chance_minz2.append(round(sum(PW2_Total[where(bins2 < redshift_min)])/sum(PW2_Total)*100,2))
        chance_maxz2.append(round(sum(PW2_Total[where(bins2 > redshift_max)])/sum(PW2_Total)*100,2))
        chance_between2.append(round(sum(PW2_Total[where((bins2 > redshift_min) & (bins2 < redshift_max))])/sum(PW2_Total)*100,2))
    else:
        chance_minz2.append(0)
        chance_maxz2.append(0)
        chance_between2.append(0)
    
    backup.append(succ)
    backup2.append(succ2)
    
    #plotting the probability distributions if desirced.
    if plot_on == True:
        figW1 = plt.figure() #figW1 for obvious reasons
        axx = figW1.add_subplot(111)
        
        axx.plot(bins2,P_L*prob_L,'r--',label='P(LERG)')
        axx.plot(bins2,P_H*prob_H,'b--',label='P(HERG)')
        axx.plot(bins2,P_Q*prob_Q,'g--',label='P(QSO)')
        axx.plot(bins2,P_Total*100,'k',label='P(Total)')
        #axx.plot([0.44,0.44],[0,6.2],'k--',label='z = 0.44 (HI)',alpha=0.5)
        axx.legend()
        #axx.set_xlim(0,1.1)
        axx.set_xlabel('Redshift',size=14)
        axx.set_ylabel('Probability (%)',size=14)
        #axx.set_title('p(z) scaled for LERG/HERG/QSO',size=14)
        fname = 'output/'+str(n)+'_W1_pdf.png'
        figW1.savefig(fname,bbox_inches='tight')

        figW2 = plt.figure() #figW2 for obvious reasons
        axxW2 = figW2.add_subplot(111)
        
        axxW2.plot(bins2,PW2_L*prob_L,'r--',label='P(LERG)')
        axxW2.plot(bins2,PW2_H*prob_H,'b--',label='P(HERG)')
        axxW2.plot(bins2,PW2_Q*prob_Q,'g--',label='P(QSO)')
        
        axxW2.plot(bins2,PW2_Total*100,'k',label='P(Total)')
        #axxW2.plot([0.54,0.54],[0,2.1],'k--',label='z = 0.54 (HI)',alpha=0.5)
        axxW2.legend()
        #axxW2.set_xlim(0,1.1)

        axxW2.set_xlabel('Redshift',size=14)
        axxW2.set_ylabel('Probability (%)',size=14)
        #axxW2.set_title('p(z) scaled for LERG/HERG/QSO',size=14)
        fname = 'output/'+str(n)+'_W2_pdf.png'
        figW2.savefig(fname,bbox_inches='tight')

#Write results to new file

infile['Chance_LERG'] = chance_L
infile['Chance_HERG'] = chance_H
infile['Chance_QSO'] = chance_Q
infile['most_likely'] = highest
infile['median_W1_Total'] = med_T
infile['int_W1_Total_1sig'] = [str(round(max(0.000,item[0][0]),4))+'-'+str(round(max(0.000,item[1][0]),4)) for item in interval_T]
infile['int_W1_Total_2sig'] = [str(round(max(0.000,item[0][1]),4))+'-'+str(round(max(0.000,item[1][1]),4)) for item in interval_T]
infile['int_W1_Total_3sig'] = [str(round(max(0.000,item[0][2]),4))+'-'+str(round(max(0.000,item[1][2]),4)) for item in interval_T]
infile['median_W2_Total'] = med_T2
infile['int_W2_Total_1sig'] = [str(round(max(0.000,item[0][0]),4))+'-'+str(round(max(0.000,item[1][0]),4)) for item in interval_T2]
infile['int_W2_Total_2sig'] = [str(round(max(0.000,item[0][1]),4))+'-'+str(round(max(0.000,item[1][1]),4)) for item in interval_T2]
infile['int_W2_Total_3sig'] = [str(round(max(0.000,item[0][2]),4))+'-'+str(round(max(0.000,item[1][2]),4)) for item in interval_T2]
infile['combine_succ'] = backup
infile['combine_succ_W2'] = backup2
infile['prob_lessthan_minz'] = chance_minz
infile['prob_between_zrange'] = chance_between
infile['prob_morethan_maxz'] = chance_maxz
infile['prob_lessthan_minz_W2'] = chance_minz2
infile['prob_between_zrange_W2'] = chance_between2
infile['prob_morethan_maxz_W2'] = chance_maxz2
infile.to_csv('results.csv', sep = ',')