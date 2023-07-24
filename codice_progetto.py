####Cavazzini Lorenzo
####Matricola:0000916454
####ID:61

##########import libraries ###########################
import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats.stats import pearsonr
import random

dir_origin= os.getcwd()
dir_project=dir_origin+'/Project'
dir_data=dir_project+'/Data'
dir_output=dir_project+'/Output'
dir_plot=dir_project+'/Plot'
dir_step2=dir_plot+'/step_2/'
dir_subsample=dir_step2+'/subsample/'
dir_parent_sample=dir_step2+'/parent_sample/'
dir_step3=dir_plot+'/step_3/'
dir_step4=dir_plot+'/step_4/'
dir_step2_out=dir_output+'/step_2/'
dir_step3_out=dir_output+'/step_3/'
dir_step4_out=dir_output+'/step_4/'

#Directory creation ##################################
list_dir=[dir_project,
          dir_data,
          dir_output,
          dir_plot,
          dir_step2,
          dir_subsample,
          dir_parent_sample,
          dir_step3,
          dir_step4,
          dir_step2_out,
          dir_step3_out,
          dir_step4_out]

for i in list_dir :
    if os.path.exists(i) :
        print('Directory ',i,'already existing')
    else :
        os.mkdir(i)
        print('Directory',i,'created')

#cerca se il file esiste oppure lo sposta ############
os.chdir(dir_data)

filename='data_SDSS_Info.fit'

if os.path.isfile(filename):
    print('file already existing')
else:
    move='mv '+dir_origin+'/'+filename+' .'
    print('moving file')
    os.system(move)

hdul=fits.open(filename)
cols=hdul[1].columns
values=hdul[1].data
ID=values["ID"]
my_ID= ID==61
my_data=values[my_ID]
hdul.close()

#Creazione delle liste con le sole quantita necessarie allo svolgimento del progetto
data_list=[values['petroMag_u'],
           values['petroMag_g'],
           values['petroMag_r'],
           values['petroMag_i'],
           values['petroMag_z'],
           values['h_alpha_flux'],
           values['h_beta_flux'],
           values['oiii_5007_flux'],
           values['nii_6584_flux'],
           values['lgm_tot_p50'],
           values['sfr_tot_p50'],
           values['absMagU'],
           values['absMagG'],
           values['absMagR'],
           values['absMagI'],
           values['absMagZ'],
           values['z']]

my_data_list=[my_data['petroMag_u'],
              my_data['petroMag_g'],
              my_data['petroMag_r'],
              my_data['petroMag_i'],
              my_data['petroMag_z'],
              my_data['h_alpha_flux'],
              my_data['h_beta_flux'],
              my_data['oiii_5007_flux'],
              my_data['nii_6584_flux'],
              my_data['lgm_tot_p50'],
              my_data['sfr_tot_p50'],
              my_data['absMagU'],
              my_data['absMagG'],
              my_data['absMagR'],
              my_data['absMagI'],
              my_data['absMagZ'],
              my_data['z']]

#list with name of properties
name_list=['petroMag_u',
           'petroMag_g',
           'petroMag_r',
           'petroMag_i',
           'petroMag_z',
           'h_alpha_flux',
           'h_beta_flux',
           'oiii_5007_flux',
           'nii_6584_flux',
           'lgm_tot_p50',
           'sfr_tot_p50',
           'absMagU',
           'absMagG',
           'absMagR',
           'absMagI',
           'absMagZ',
           'z']

print('\ndati estratti dal file.')

#############################################################################################################
print('\nINIZIO SECONDO STEP \n')

#pulizia delle liste di dati
#getting rid of non readable data (non-float values)                        
def nan_inf(data):
    check_nan_inf=[]
    for i in range(len(data)):
        check_nan_inf.append(data[i][np.isfinite(data[i])])
    return check_nan_inf

data_list1=nan_inf(data_list)
my_data_list1=nan_inf(my_data_list)

#getting rid of outliers (data outside 4-standard deviation)
def sigma_clipping(data):
    sigma_clipped=[]
    for i in range(len(data)):
        x=data[i]
        x=x[np.logical_and(x>x.mean()-4*x.std(),x<x.mean()+4*x.std())]
        x=x[np.logical_and(x>x.mean()-4*x.std(),x<x.mean()+4*x.std())]
        x=x[np.logical_and(x>x.mean()-4*x.std(),x<x.mean()+4*x.std())]
        x=x[np.logical_and(x>x.mean()-4*x.std(),x<x.mean()+4*x.std())]
        sigma_clipped.append(x)
    return sigma_clipped

data_list2=sigma_clipping(data_list1)                           #final list of cleaned data
my_data_list2=sigma_clipping(my_data_list1)                     #final list of my cleaned data
print('dati puliti da valori non leggibili e outliers')

#gaussian function for plotting
def gaussian(bins,mean,std) :
    x=np.zeros(len(bins)-1)
    for i in range (len(x)) :
        x[i]=(bins[i]+bins[i+1])/2
    y=1/(std*np.sqrt(2*np.pi))*np.exp(-(x-mean)**2/(2*std**2))
    return x, y

#histogram function for step2/step4 
def histogram(data,name) :
    i=0
    plt.figure(8, figsize=(13,8))           
    gs=gridspec.GridSpec(1,3)
    while i<2 :
        ax=plt.subplot(gs[i])
        counts, bins, ignored=ax.hist(data, 25, histtype='bar', density=True, align='mid')  
        xm,mod=gaussian(bins,np.mean(data),np.std(data))                                    
        ax.plot(xm, mod, c='#7CFC00', ls='-.', lw=3) 
        if i==0 :                                                               #mean plot
            ax.axvline(np.mean(data), c='m', ls='--',lw=3)
            ax.axvspan(np.mean(data)-np.std(data), np.mean(data)+np.std(data), color='#FFD700', alpha=0.3) 
            ax.set_xlabel(name,fontsize=15)
            ax.set_title('mean and error',fontsize=14)
        else :                                                                  #median plot
            ax.axvline(np.median(data), c='m', ls='--',lw=3)
            ax.axvspan(np.percentile(data,16), np.percentile(data,84), color='#CD7F32', alpha=0.3)
            ax.set_xlabel(name,fontsize=15)
            ax.set_title('median and 16th to 84th percentiles',fontsize=12)
        i=i+1       
    axr=plt.subplot(gs[i])                                                      
    res=(counts-mod)                            
    axr.scatter(xm,res)
    axr.axhline(np.mean(res), c='b', ls='--')
    axr.set_xlabel(name,fontsize=15)
    axr.set_title('Residuals',fontsize=14)                                                      
    plt.savefig(name+'.png')
    
#mean function for step2/step3/step4
def output_mean(data,filename,name):
    i=0
    mean=[]
    std=[]
    while i<len(data) :
        mean.append(np.mean(data[i]))
        std.append(np.std(data[i]))
        i=i+1
    m=np.array(mean)                    #trasformo le liste in numpy array per poterne fare le operazioni
    s=np.array(std)
    n=len(name)
    r=np.zeros(n, dtype=[('var1', 'U29'), ('var2', float),('var3', float), ('var4', float)])
    r['var1']= name
    r['var2']= m-s
    r['var3']= m
    r['var4']= m+s
    fileout=open(filename+'.dat','w')
    np.savetxt(fileout, r, delimiter="", fmt='%s\t\t\t %f\t\t\t %f\t\t\t %f', newline='\n', header='Name\t\t\t\t mean-std\t\t\t  mean\t\t\t\t mean+std')
    fileout.close()
    
#median function for step2/step3/step4
def output_median(data,filename,name) :
    i=0
    median=[]
    perc_16=[]
    perc_84=[]
    while i<len(data) :
        median.append(np.median(data[i]))
        perc_16.append(np.percentile(data[i],16))
        perc_84.append(np.percentile(data[i],84))
        i=i+1
    n=len(name)
    r=np.zeros(n, dtype=[('var1', 'U29'), ('var2', float),('var3', float), ('var4', float)])
    r['var1']= name
    r['var2']= perc_16
    r['var3']= median
    r['var4']= perc_84
    fileout=open(filename+'.dat','w')
    np.savetxt(fileout ,r ,delimiter="",fmt='%s\t\t\t %f\t\t\t %f\t\t\t %f', newline='\n', header='Name\t\t\t   16th_percentile\t\t\t   median\t\t\t  84th_percentile')
    fileout.close()

#plotting histograms ######################################################################################################
os.chdir(dir_subsample)
k=0
while k<len(my_data_list2) :
    subsample=my_data_list2[k]
    histogram(subsample,name_list[k])
    k=k+1

print('Finito di plottare subsample')
os.chdir(dir_parent_sample)

k=0
while k<len(data_list2) :
    parent_sample=data_list2[k]
    histogram(parent_sample,name_list[k])
    k=k+1

print('Finito di plottare parent_sample')

#calling output functions for mean and median values
os.chdir(dir_step2_out)

output_mean(data_list2,'mean_std',name_list)
output_mean(my_data_list2,'my_mean_std',name_list)
output_median(data_list2,'median_percentiles',name_list)
output_median(my_data_list2,'my_median_percentiles',name_list)

print('Dati salvati in output \n')

#calculating mean and error manually for one case
redshift=my_data_list2[16]
j=0
h=0
summ=0
summ2=0
while j<len(redshift) :
        summ=summ+redshift[j]
        j=j+1
mean_redshift=summ/len(redshift)

while h<len(redshift) :
        summ2=summ2+((redshift[h]-mean_redshift)**2)
        h=h+1
std_redshift=np.sqrt(summ2/(len(redshift)-1))
m=round(mean_redshift,5)
ms1=round(mean_redshift-std_redshift,5)
ms2=round(mean_redshift+std_redshift,5)

print('con il calcolo manuale la media del (mio) redshift è: ',m)
print('mean-std: ',ms1,' mean+std: ',ms2,'\n')

##############################################################################################################################
print ('INIZIO TERZO STEP \n')

def step_3(x,y,name) :                      
    res=np.polyfit(x,y,1)                  
    x_bf=np.arange(x.min(), x.max(), 0.0005)
    y_bf=np.polyval(res,x_bf)              
    r=pearsonr(x,y)   #returns a correation coefficient between x and y values
    
    plt.figure(figsize=(9,7))       
    plt.scatter(x,y,s=18 , c='#77DD77', marker='o', label=name+'(redshift)', edgecolor='#50C878')        
    plt.plot(x_bf,y_bf, c='#E30B5C', ls='--', lw=2, label='best_fit', zorder=1)                          
    plt.legend(loc='best',fontsize='medium')
    plt.axis(xmin=0)
    plt.ylabel(name,fontsize=12)
    plt.xlabel('redshift', fontsize=12)
    plt.title('correlation value='+ str(round(r[0],3)))
    plt.savefig(name+'.png')
    
    if r[0]>0.5 or r[0]<-0.5 :                   #se la correlazione è forte:
        print('rilevata correlazione tra redshift e '+name)
        z=[0.00,0.02,0.04,0.06,0.08,0.10] 
        bin_interv=['0.00-0.02','0.02-0.04','0.04-0.06','0.06-0.08','0.08-0.10']
        z_bins=[]
        i=0
        while i<len(z)-1 :                              
            bin_i=y[np.logical_and(x>=z[i], x<z[i+1])]
            z_bins.append(bin_i)                            
            i=i+1
            
        plt.figure(1, figsize=(12,6))                       #fà istogrammi per correlazione forte (media e std)
        gs=gridspec.GridSpec(1,5)
        gs.update(wspace=0.5, hspace=0.5)
        for j in range(len(z_bins)) :
            a=plt.subplot(gs[j])
            counts, bins, ign=a.hist(z_bins[j], 25, histtype='bar', density=True, align='mid') 
            xm,mod=gaussian(bins,np.mean(z_bins[j]),np.std(z_bins[j]))
            a.plot(xm, mod, c='#7CFC00', ls='-.', lw=3)                          
            a.axvline(np.mean(z_bins[j]), c='m', ls='--',lw=3)
            a.axvspan(np.mean(z_bins[j])-np.std(z_bins[j]), np.mean(z_bins[j])+np.std(z_bins[j]), color='#FFBF00', alpha=0.25) 
            a.set_title('redshift tra '+bin_interv[j], fontsize='small')
        plt.xlabel(name+'(correlation value='+ str(round(r[0],3))+')',fontsize=12)
        plt.savefig('redshift trend for '+name+'(mean and std)',dpi=150)
        plt.close()
        
        plt.figure(1, figsize=(12,6))                       #fa istogrammmi per correlazione forte (mediana e percentili)
        gs=gridspec.GridSpec(1,5)
        gs.update(wspace=0.5, hspace=0.5)
        for k in range(len(z_bins)) :
            b=plt.subplot(gs[k])
            counts, bins, ign=b.hist(z_bins[k], 25, histtype='bar', density=True, align='mid') 
            xm,mod=gaussian(bins,np.mean(z_bins[k]),np.std(z_bins[k]))
            b.plot(xm, mod, c='#7CFC00', ls='-.', lw=3)
            b.axvline(np.median(z_bins[k]), c='m', ls='--',lw=3)
            b.axvspan(np.percentile(z_bins[k],16), np.percentile(z_bins[k],84), color='#CD7F32', alpha=0.25)
            b.set_title('redshift tra '+bin_interv[k] , fontsize='small')
        plt.xlabel(name+'(correlation value='+ str(round(r[0],3))+')',fontsize=10, loc='center')    
        plt.savefig('redshift trend for '+name+'(median and percentiles)',dpi=150)
        plt.close()

        os.chdir(dir_step3_out)                             #gli output della media e mediana per corrleazione
        output_mean(z_bins,'correlation_mean_std_'+name, bin_interv)
        output_median(z_bins,'correlation_median_percentile_'+name, bin_interv)
        os.chdir(dir_step3)
        print('I dati della correlazione sono stati salvati.')

#dobbiamo fare in modo di avere array di dati della stessa lunghezza del redshift per poterli plottare uno in funzione dell'altro
name_step3=['petroMag_z','h_alpha_flux','lgm_tot_p50','Sfr_tot_p50','absMagZ'] 
my_data_list3=[my_data_list[4],
               my_data_list[5],
               my_data_list[9],
               my_data_list[10],
               my_data_list[15]]
os.chdir(dir_step3)
i=0
while i<len(my_data_list3) :    
    x=my_data_list[16]
    y=my_data_list3[i]
    y=y[np.isfinite(x)]         #making redshift and the other properties of the same size
    x=x[np.isfinite(x)]
    x=x[np.isfinite(y)]
    y=y[np.isfinite(y)]
    k=0
    while k<4 :                 #sigma-clipping
        x=x[np.logical_and(y>y.mean()-4*y.std(), y<y.mean()+4*y.std())]
        y=y[np.logical_and(y>y.mean()-4*y.std(), y<y.mean()+4*y.std())]
        k=k+1
    step_3(x,y,name_step3[i])   #calling the function
    print('Finito di plottare '+name_step3[i]+' in funzione di redshift')
    i=i+1

############################################################################################################################
print ('\nINIZIO QUARTO STEP \n')

def plot_4(x,y,z,x_rel,y_rel) :
    x_label=['log(nIII/H_alpha)','log(M/Mo)','log(M/Mo)']
    y_label=['log(oIII/H_beta)','u-r','SFR']
    name=['BPT_diagram','Color-mass_diagram','SFR-mass_diagram']
    color_maps=['viridis','cividis','magma']
    
    plt.figure(figsize=(11,9))
    scat=plt.scatter(x,y,c=z,cmap=color_maps[i],vmin=min(z), vmax=max(z))
    cb=plt.colorbar(scat)
    plt.plot( x_rel,y_rel, c='b', linestyle='--',label='theoretical relation')
    plt.legend(loc='upper left',fontsize=11)
    cb.set_label(r'redshift',fontsize=14)
    plt.xlabel(x_label[i],fontsize=14)
    plt.ylabel(y_label[i],fontsize=14)
    plt.savefig(name[i]+'.png',)

def sigma_clipping4(x,y,z) :
    k=0
    while k<4 :
        y=y[np.logical_and(x>x.mean()-4*x.std(),x<x.mean()+4*x.std())]
        z=z[np.logical_and(x>x.mean()-4*x.std(),x<x.mean()+4*x.std())]
        x=x[np.logical_and(x>x.mean()-4*x.std(),x<x.mean()+4*x.std())]
        k=k+1
    return x,y,z

#estrazione proprietà utili allo svolgimento dello step4 #################################################
u=my_data_list[0]
r=my_data_list[2]
h_a=my_data_list[5]
h_b=my_data_list[6]
oiii=my_data_list[7]
nii=my_data_list[8]
mass=my_data_list[9]
sfr=my_data_list[10]
z=my_data_list[16]

#BPT #####################################################################################################
oiii=oiii[h_a!=0]
nii=nii[h_a!=0]
h_b=h_b[h_a!=0]
z1=z[h_a!=0]
h_a=h_a[h_a!=0]

oiii=oiii[h_b!=0]
nii=nii[h_b!=0]
h_a=h_a[h_b!=0]
z1=z1[h_b!=0]
h_b=h_b[h_b!=0]

x1=nii/h_a
y1=oiii/h_b

y1=y1[x1>0]
z1=z1[x1>0]
x1=x1[x1>0]

x1=x1[y1>0]
z1=z1[y1>0]
y1=y1[y1>0]

log_x=np.log(x1)
log_y=np.log(y1)

log_y,log_x,z1=sigma_clipping4(log_y,log_x,z1)  
log_x,log_y,z1=sigma_clipping4(log_x,log_y,z1)  #dati da usare per i plot: log(nII/H_a), log(oIII/H_b), z

#above and below theorethical relation
log_x_a=[]
log_x_b=[]
log_y_a=[]
log_y_b=[]
for k in range(len(log_y)) :
    if log_y[k] > 0.61/(log_x[k]-0.05)+1.3 :
        log_y_a.append(log_y[k])
    else :
        log_y_b.append(log_y[k])

for k in range(len(log_x)) :
    if log_x[k] > 0.61/(log_y[k]-1.3)+0.05 :
        log_x_a.append(log_x[k])
    else :
        log_x_b.append(log_x[k])
        
#relation for thoretica line in diagram
count,bins1=np.histogram(log_x,100)                         
l= np.zeros(len(bins1)-1)
for i in range(len(l)):
	l[i]=(bins1[i]+bins1[i+1])/2

x1_rel=l
y1_rel=0.61/(l-0.05)+1.3
y1_rel=y1_rel[x1_rel<0]
x1_rel=x1_rel[x1_rel<0]
x1_rel=x1_rel[y1_rel>-5]
y1_rel=y1_rel[y1_rel>-5]

#color-mass ##############################################################################################
y2=u-r
mass2,y2,z2=sigma_clipping4(mass,y2,z)
y2,mass2,z2=sigma_clipping4(y2,mass2,z2)        #dati da usare per i plot: u-r, stellar mass, z

#above and below theorethical relation
y2_a=[]
y2_b=[]
x2_a=[]
x2_b=[]
for k in range(len(y2)) :
    if y2[k] > -0.495+0.25*(mass2[k]) :
        y2_a.append(y2[k])
    else:
        y2_b.append(y2[k])

for k in range(len(mass2)) :
    if mass2[k] > (0.495+y2[k])/0.25 :
        x2_a.append(mass2[k])
    else :
        x2_b.append(mass2[k])

#theoretical relation
y2_rel=[]
for j in range(len(mass2)) :
    rel2=-0.495+0.25*(mass2[j])
    y2_rel.append(rel2)
x2_rel=mass2.copy()

#SFR-mass ################################################################################################

mass3,y3,z3=sigma_clipping4(mass,sfr,z)
y3,mass3,z3=sigma_clipping4(y3,mass3,z3)        #dati da usare per i plot: sfr, stellar mass, z

#above and below theorethical relation
y3_a=[]
y3_b=[]
x3_a=[]
x3_b=[]
for k in range(len(y3)) :
    if y3[k] > -8.64+0.76*(mass3[k]) :
        y3_a.append(y3[k])
    else :
        y3_b.append(y3[k])

for k in range(len(mass3)):
    if mass3[k] > (y3[k]+8.64)/0.76 :
        x3_a.append(mass3[k])
    else :
        x3_b.append(mass3[k])

#theorethical relation
y3_rel=[]
for j in range(len(mass3)) :
    rel3=-8.64+0.76*(mass3[j])
    y3_rel.append(rel3)
x3_rel=mass3.copy()

#chiamata delle funzioni plot ################################################################################
os.chdir(dir_step4)
i=0
while i<3 :
    if i==0 :
        plot_4(log_x,log_y,z1,x1_rel,y1_rel)
        print('Finito di plottare BPT_diagram')
        histogram(log_x_a,'above_BPT_x')
        histogram(log_x_b,'below_BPT_x')
        histogram(log_y_a,'above_BPT_y')
        histogram(log_y_b,'below_BPT_y')
        print('Finito di plottare gli istogrammi di BPT')
    elif i==1 :
        plot_4(mass2,y2,z2,x2_rel,y2_rel)
        print('Finito di plottare color-mass_diagram')
        histogram(x2_a,'above_color-mass_x')
        histogram(x2_b,'below_color-mass_x')
        histogram(y2_a,'above_color-mass_y')
        histogram(y2_b,'below_color-mass_y')
        print('Finito di plottare gli istogrammi di color-mass')
    elif i==2 :
        plot_4(mass3,y3,z3,x3_rel,y3_rel)
        print('Finito di plottare SFR-mass_diagram')
        histogram(x3_a,'above_SFR_x')
        histogram(x3_b,'below_SFR_x')
        histogram(y3_a,'above_SFR_y')
        histogram(y3_b,'below_SFR_y')
        print('Finito di plottare gli istogrammi di SFR-mass')
    i=i+1

#chiamata delle funzioni output 
out_list1=[log_x,log_y,log_x_a,log_x_b,log_y_a,log_y_b]
out_list2=[mass2,y2,x2_a,x2_b,y2_a,y2_b]
out_list3=[mass3,y3,x3_a,x3_b,y3_a,y3_b]
names=['   x   ','   y   ','x_above','x_below','y_above','y_below',]

os.chdir(dir_step4_out)
output_mean(out_list1,'BPT_mean_std_',names)
output_median(out_list1,'BPT_median_percentile',names)
output_mean(out_list2,'color-mass_mean_std_',names)
output_median(out_list2,'color-mass_median_percentile',names)
output_mean(out_list3,'SFR_mean_std_',names)
output_median(out_list3,'SFR_median_percentile',names)
print('Dati salvati in output.')
