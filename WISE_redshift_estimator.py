# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 16:52:46 2017

@author: Marcin
"""
#Redshift estimation for a single radio source.
#Input the WISE magnitudes and error for your source.

import os
import sys
from math import *
from numpy import *
from pylab import *
from matplotlib import rc, rcParams
from random import *
rc('text', usetex=False)
rc('font',**{'family':'serif','serif':['serif'],'size':10})
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage,AnnotationBbox
from matplotlib.cbook import get_sample_data
from matplotlib._png import read_png
from astropy.io import ascii
from scipy.stats import kde
from scipy import stats
import pandas as pd
from operator import itemgetter
import warnings
warnings.filterwarnings("ignore")

##################################################
#important numbers! change limits as you see fit
nbins = 100
W1_in = 16.0;
W1_er = 0.1;
W2_in = 15.6;
W2_er = 0.1
W3_in = 13.3;
W3_er = 0.3

pt = [round(W2_in-W3_in,4),round(W1_in-W2_in,4)]
fluxlim = 0; #50 mjy
bin_no = 150;
step = (18.-8.)/bin_no;
#18-minimum for low chance of QSO, 18-30 for decent chance of QSO? 
#Remember redshift ranges are different!
bin_switch = True
#if true, then if chance of QSO is high bin_no is switched to 30. Can turn off.
redshift_min = 0.4; #in case you want to know the chance of object at less than this redshift
plot_on = True #False if you don't want to plot anything to fine
print "WISE colours of radio source are "+str(pt)+"."
##################################################
df = pd.read_csv(
    filepath_or_buffer='jhyc_wise_spec.txt',
    sep=' ')
df = df[df['W1']>0];
df = df[df['W2']>0];
df = df[df['W3']>0];
w1 = df.ix[:,27].values
w2 = df.ix[:,29].values
w3 = df.ix[:,31].values
w1e = df.ix[:,28].values
w2e = df.ix[:,30].values
w3e = df.ix[:,32].values
sclass = df.ix[:,20].values
flux_f = df.ix[:,11].values
flux_n = df.ix[:,13].values          
(x,y) = (w2-w3,w1-w2)             

##################################################
# stage one - determine the chance of being a LERG, etc
if plot_on == True:
    (x,y) = (w2-w3,w1-w2)
    fig, axes = subplots(ncols=2, nrows=2, sharex=True, sharey=True)

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
k = kde.gaussian_kde([asarray(x_sub),asarray(y_sub)])
a = asarray(x_sub)

#xi, yi = np.mgrid[min(x_sub):max(a[a<10]):nbins*1j, min(y_sub):max(y_sub):nbins*1j]#not max(x_sub) here, too large
xi, yi = np.mgrid[-1:6:nbins*1j, -1:3:nbins*1j]
zi = k(np.vstack([xi.flatten(), yi.flatten()]))
ziL = zi
if plot_on == True:
    axes[0,0].pcolormesh(xi, yi, zi.reshape(xi.shape))
    #axes[1,1].pcolormesh(xi, yi, zi.reshape(xi.shape),alpha=0.7)
percent_L = float(k(pt)/max(zi))*100
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
k = kde.gaussian_kde([asarray(x_sub),asarray(y_sub)])
a = asarray(x_sub)
#xi, yi = np.mgrid[min(x_sub):max(a[a<10]):nbins*1j, min(y_sub):max(y_sub):nbins*1j]#not max(x_sub) here, too large
xi, yi = np.mgrid[-1:6:nbins*1j, -1:3:nbins*1j]
zi = k(np.vstack([xi.flatten(), yi.flatten()]))
ziH = zi
if plot_on == True:
    axes[1,0].pcolormesh(xi, yi, zi.reshape(xi.shape))
    #axes[1,1].pcolormesh(xi, yi, zi.reshape(xi.shape),alpha=0.7)
percent_H = float(k(pt)/max(zi))*100
tot_h = len(x_sub)
x_sub = [];
y_sub = [];
name = 'AeB'
cc = 0
for i in sclass:
    if i == name and (float(flux_f[cc]) > fluxlim or float(flux_n[cc]) > fluxlim):
        x_sub.append(x[cc])
        y_sub.append(y[cc])
    cc += 1;
#scatter(x_sub,y_sub,facecolor='g',marker='o',s=18,zorder=101,alpha=1.0,linewidth=0.5)
k = kde.gaussian_kde([asarray(x_sub),asarray(y_sub)])
a = asarray(x_sub)
#xi, yi = np.mgrid[min(x_sub):max(a[a<10]):nbins*1j, min(y_sub):max(y_sub):nbins*1j]#not max(x_sub) here, too large
xi, yi = np.mgrid[-1:6:nbins*1j, -1:3:nbins*1j]
zi = k(np.vstack([xi.flatten(), yi.flatten()]))
ziQ = zi
tot_q = len(x_sub)
tott = tot_l+tot_h+tot_q
ziT = ziL*(float(tott)/tot_l)+ziH*(float(tott)/tot_h)+ziQ*(float(tott)/tot_q) #scaling for appearence
if plot_on == True:
    axes[0,1].pcolormesh(xi, yi, zi.reshape(xi.shape))    
    axes[1,1].pcolormesh(xi, yi, ziT.reshape(xi.shape),alpha=1.0)

percent_Q = float(k(pt)/max(zi))*100


if plot_on == True:
    ylim(-1,3)
    xlim(-1,6)
    axes[0,0].set_title('LERG')
    axes[1,0].set_title('HERG')
    axes[0,1].set_title('QSO')
    axes[1,1].set_title('Combined')
    
    axes[0,0].scatter(pt[0],pt[1], facecolor='r',marker='D',s=30,linewidth=0.5)
    axes[1,0].scatter(pt[0],pt[1], facecolor='r',marker='D',s=30,linewidth=0.5)
    axes[0,1].scatter(pt[0],pt[1], facecolor='r',marker='D',s=30,linewidth=0.5)
    axes[1,1].scatter(pt[0],pt[1], facecolor='r',marker='D',s=30,linewidth=0.5)
    axes[1,0].set_xlabel('W2-W3 (mag)')
    axes[1,0].set_ylabel('W1-W2 (mag)')
    axes[1,1].set_xlabel('W2-W3 (mag)')
    axes[0,0].set_ylabel('W1-W2 (mag)')
    fname = 'output/TestSource_WISE_Colours.png'
    fig.savefig(fname,bbox_inches='tight')

tot = percent_L+percent_H+percent_Q
prob_L = round(percent_L/tot*100,2);
prob_H = round(percent_H/tot*100,2);
prob_Q = round(percent_Q/tot*100,2);
print("LERG probability is "+str(prob_L)+"%")
print("HERG probability is "+str(prob_H)+"%")
print("QSO probability is "+str(prob_Q)+"%")

#as QSOs span a wider redshift distribution, we can encforce more bins in plotting.
#bin_switch = True necessary. 
if (prob_Q > 10 or W1_in > 16.0) and bin_switch == True:
    bin_no = 150;
    step = (18.-8.)/bin_no;
    #make the switch off if you do not want this!

##################################################

z_LERG = (df.loc[(df['BESTCLASS']=='LERG') & 
                 ((df['FIRST_TOTINTFLUX'] > fluxlim) | (df['NVSS_TOTINTFLUX'] > fluxlim)),['Z']])
z_HERG = (df.loc[(df['BESTCLASS']=='HERG') & 
                 ((df['FIRST_TOTINTFLUX'] > fluxlim) | (df['NVSS_TOTINTFLUX'] > fluxlim)),['Z']])
z_QSO = (df.loc[(df['BESTCLASS']=='AeB') & 
                ((df['FIRST_TOTINTFLUX'] > fluxlim) | (df['NVSS_TOTINTFLUX'] > fluxlim)),['Z']])

binmax = 1.2#0.975;
if prob_Q > 10.0 or W1_in >= 16.0:
    binmax = max(asarray(z_QSO)) + 0.025;
bins2 = linspace(0.0,binmax,bin_no)
bins3 = linspace(0.0,max(asarray(z_QSO))-0.025,bin_no)
binwidth = bins2[1]-bins2[0];
binwidth3 = bins3[1]-bins3[0];
WISErange = arange(8,18,step)
W1weight = []
probW1 = stats.norm.pdf(WISErange,W1_in,W1_er)*step
probW2 = stats.norm.pdf(WISErange,W2_in,W2_er)*step
W1diff = abs(WISErange-W1_in)
W2diff = abs(WISErange-W2_in)
#for i in sub_LERG['W1diff']:
#    minlist = [abs(j - i) for j in WISEdiff]
#    W1weight.append(probW1[min(enumerate(minlist), key=itemgetter(1))[0]])
#sub_LERG['W1Weight'] = W1weight

P_L = [0]*bin_no
P_H = P_L
P_Q = P_L
PW2_L = P_L
PW2_H = P_L
PW2_Q = P_L

pdfsW1 = []
pdfsW2 = []
#LERGS FIRST
cc = 0;

#here we calculate the mean and std expected for bright/faint W1, W2 magnitudes
#we loop over the LERGs and HERGs, W1 and W2 
#take the middle, well populated parts
#and extrapolate for higher/lower values
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
#xx_L = WISErange[where(abs(WISErange - 13.25)<2.25)]
mmm = [];
sss = []
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
    #figure(cc)
    #plot(bins2,pdf*len(zdistW2))
    #hist(asarray(zdistW2),bins=bins2)
    #title(i)
    cc += 1
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
#    
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

P_Total = P_L*prob_L/100+P_H*prob_H/100+P_Q*prob_Q/100
PW2_Total = PW2_L*prob_L/100+PW2_H*prob_H/100+PW2_Q*prob_Q/100

try:
    W1_median = (round(mean([min(bins2[where(cumsum(P_Total)>=(0.5))]),max(bins2[where(cumsum(P_Total)<=(0.5))])]),2))
except:
    W1_median = (-999)
    #something went wrong for the source; this value indicates this in the output file.
#here the distributions are usually skewed, so assuming a normal distribution would be foolish
#we just calculate the intervals for the 68/95/99.7% CIs ourselves via cumsum.
try:
    sig1_low = round(mean([min(bins2[where(cumsum(P_Total)>=(0.5-0.3413))]),max(bins2[where(cumsum(P_Total)<=(0.5-0.3413))])]),2)
except:
    sig1_low = round(min(bins2),2)
try:
    sig2_low = round(mean([min(bins2[where(cumsum(P_Total)>=(0.5-0.3413-0.136))]),max(bins2[where(cumsum(P_Total)<=(0.5-0.3413-0.136))])]),2)
except:
    sig2_low = round(min(bins2),2)            
try:
    sig3_low = round(mean([min(bins2[where(cumsum(P_Total)>=(0.5-0.3413-0.136-0.0214))]),max(bins2[where(cumsum(P_Total)<=(0.5-0.3413-0.136-0.0214))])]),2)
except:
    sig3_low = round(min(bins2),2)          
try:
    sig1_high = round(mean([min(bins2[where(cumsum(P_Total)>=(0.5+0.3413))]),max(bins2[where(cumsum(P_Total)<=(0.5+0.3413))])]),2)
except:
    sig1_hig = round(max(bins2),2)    
try:
    sig2_high = round(mean([min(bins2[where(cumsum(P_Total)>=(0.5+0.3413+0.136))]),max(bins2[where(cumsum(P_Total)<=(0.5+0.3413+0.136))])]),2)
except:
    sig2_high = round(max(bins2),2) 
try:
    sig3_high = round(mean([min(bins2[where(cumsum(P_Total)>=(0.5+0.3413+0.136+0.0214))]),max(bins2[where(cumsum(P_Total)<=(0.5+0.3413+0.136+0.0214))])]),2)
except:
    sig3_high = round(max(bins2),2)
print "W1: Median redshift of "+str(W1_median)+", 68\% redshift CI of "+str(sig1_low)+" - "+str(sig1_high)

#now same but for W2!
try:
    W2_median = (round(mean([min(bins2[where(cumsum(PW2_Total)>=(0.5))]),max(bins2[where(cumsum(PW2_Total)<=(0.5))])]),2))
except:
    W2_median = (-999)
    #this suggests that we have empty sets! too high/low a W1 (or W2) magnitude to be useable.
    #can extend this code later
#here the distributions are usually skewed, so assuming a normal distribution would be foolish
#he just calculate the intervals for 1,2,3 sigma ourselves via cumsum.
try:
    sig1_low = round(mean([min(bins2[where(cumsum(PW2_Total)>=(0.5-0.3413))]),max(bins2[where(cumsum(PW2_Total)<=(0.5-0.3413))])]),2)
except:
    sig1_low = round(min(bins2),2)    
try:
    sig2_low = round(mean([min(bins2[where(cumsum(PW2_Total)>=(0.5-0.3413-0.136))]),max(bins2[where(cumsum(PW2_Total)<=(0.5-0.3413-0.136))])]),2)
except:
    sig2_low = round(min(bins2),2)           
try:
    sig3_low = round(mean([min(bins2[where(cumsum(PW2_Total)>=(0.5-0.3413-0.136-0.0214))]),max(bins2[where(cumsum(PW2_Total)<=(0.5-0.3413-0.136-0.0214))])]),2)
except:
    sig3_low = round(min(bins2),2)             
try:
    sig1_high = round(mean([min(bins2[where(cumsum(PW2_Total)>=(0.5+0.3413))]),max(bins2[where(cumsum(PW2_Total)<=(0.5+0.3413))])]),2)
except:
    sig1_high = round(max(bins2),2)
try:
    sig2_high = round(mean([min(bins2[where(cumsum(PW2_Total)>=(0.5+0.3413+0.136))]),max(bins2[where(cumsum(PW2_Total)<=(0.5+0.3413+0.136))])]),2)
except:
    sig2_high = round(max(bins2),2)
try:
    sig3_high = round(mean([min(bins2[where(cumsum(PW2_Total)>=(0.5+0.3413+0.136+0.0214))]),max(bins2[where(cumsum(PW2_Total)<=(0.5+0.3413+0.136+0.0214))])]),2)
except:
    sig3_high = round(max(bins2),2)
print "W2: Median redshift of "+str(W2_median)+", 68\% redshift CI of "+str(sig1_low)+" - "+str(sig1_high)

if plot_on == True:
        figW1 = plt.figure() #figW1 for obvious reasons
        axx = figW1.add_subplot(111)
        
        axx.plot(bins2,P_L*prob_L,'r--',label='p(LERG)')
        axx.plot(bins2,P_H*prob_H,'b--',label='p(HERG)')
        axx.plot(bins2,P_Q*prob_Q,'g--',label='p(QSO)')
        axx.plot(bins2,P_Total*100,'k',label='p(Total)')
        axx.legend()
        axx.set_xlabel('z',size=14)
        axx.set_ylabel('Probability density (%)',size=14)
        fname = 'output/TestSource_W1_pdf.png'
        figW1.savefig(fname,bbox_inches='tight')

        
        figW2 = plt.figure() #figW2 for obvious reasons
        axxW2 = figW2.add_subplot(111)
        
        axxW2.plot(bins2,PW2_L*prob_L,'r--',label='p(LERG)')
        axxW2.plot(bins2,PW2_H*prob_H,'b--',label='p(HERG)')
        axxW2.plot(bins2,PW2_Q*prob_Q,'g--',label='p(QSO)')
        axxW2.plot(bins2,PW2_Total*100,'k',label='p(Total)')
        axxW2.legend()
        axxW2.set_xlabel('z',size=14)
        axxW2.set_ylabel('Probability density (%)',size=14)
        fname = 'output/TestSource_W2_pdf.png'
        figW2.savefig(fname,bbox_inches='tight')
