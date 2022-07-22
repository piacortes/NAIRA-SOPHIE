# Hi there :-), I'm Neda Heidari. I'm not a real auther but I update this program with new python package and clean it from repeatition and so on! just as a tip, if you are using this program, I guess you should check B-V of your star with literature and be sure about it or maybe at some point update it!

import os, glob
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import subprocess
import numpy as np
from fitpoly import *
from function_make_Halpha import *
import argparse

# set up - directories and files

# Create the parser
parser = argparse.ArgumentParser(description='Target name')

# Add the arguments
parser.add_argument('target', type=str, help='the name of the star you want to create the obslists')

# Execute the parse_args() method
args = parser.parse_args()

target = args.target

#target = 'Gl205'

reduced_dir = "/net/GSP/spirou/sophie/SOPHIE_data/reduced/"
save_dir= "/net/GSP/users/pcortes/SOPHIE"
#night_folders= "/home/nheidari/moon_contamination/day_kepler16.txt"
obslist_dir = '/net/GSP/users/pcortes/SOPHIE/obslist_'+target+'.txt'


d = { "G2": "","M5": "","K5": "","F0": "","K0": "","M4": "", "_ccf_": "_e2ds"}

#because there are a lot of star that their B-V in header is not written and it is known, the best way is you extract the date of oberved star from Sophie archival data the load it hear and extract the activity. 

# read list of folders (generted with SOPHIE archival data)

f = open(obslist_dir, "r")
lines = f.readlines()
night_list = []
for x in lines:
    night_list.append(x.strip('\n')[-10:])

f.close()
night_list = list(dict.fromkeys(night_list))
print night_list 
#night_list=['2021-06-16','2021-06-21','2021-06-22']
#night_list=['2016-07-08','2016-07-11','2016-07-12','2016-07-17','2016-07-21','2016-07-25','2016-08-03','2016-08-06','2016-08-11','2016-08-22','2016-08-25','2016-09-02']
#berv_test=[3.8255179243411179,3.1960131367117648,3.113221205274427,2.5265561773819329,1.922567305438569,1.2658548101592808,-0.059565503919885611,-0.45043277723021824,-1.3426152368399848,-2.8514953527346463,-3.3028799513929519,-4.4706031333105178]
#rv_test=[-22.7627,-20.3745,-20.4507,-23.9428,-29.2530,-35.3483,-46.8899,-47.6483,-40.2453,-20.3868,-21.8696,-32.0693]

HAL=[]
Sig_Hal=[]
Prot=[]
R_hk=[]
Sig_r_hk=[]
S_mw=[]
Sig_s=[]
bjd=[]
Mode=[]
Sn1=[]
Rv=[]
DPRtype=[]
berv=[]
bv=[]
rhk_header=[]
ic_extmeanzone=9
#for night in night_list:
for kk,night in enumerate(night_list):
    print(night)
    p=open(save_dir+'aan_'+target+'.txt','w')
    p.write('ccf\n')
    p.close()
    try:
        os.system('cd ' +str(reduced_dir+night)+ ' ; grep -l -r '+target+' -R | grep ccf| grep _A.fits | grep -v newccf >> '+ save_dir+'aan_'+target+'.txt')
        ccf_collect= ascii.read(save_dir+'aan_'+target+'.txt')
        print(ccf_collect)
        #os.system('rm aan.txt')
    except IOError:
        continue 
    for ii in range(len(ccf_collect)):
        print(night,ccf_collect['ccf'][ii])
        hdu= fits.open(reduced_dir+night+'/'+str(ccf_collect['ccf'][ii]))
        rv= hdu[0].header['HIERARCH OHP DRS CCF RV']
        fwhm= hdu[0].header['HIERARCH OHP DRS CCF FWHM']
        try:
        	rhk_header.append(hdu[0].header['HIERARCH OHP DRS RHK'])
        except KeyError:
        	rhk_header.append(9999.99)
        try:
           sn1= hdu[0].header['HIERARCH OHP DRS CAL EXT SN1']
        except KeyError:
            Sn1.append(9999.99)
        BJD= hdu[0].header['HIERARCH OHP DRS BJD']
        hdu.close()
######################################################################################################e2dsA
        e2dsfits= replace_all(reduced_dir+night+'/'+str(ccf_collect['ccf'][ii]), d)
        hduii_lst = fits.open(e2dsfits)
        mode = hduii_lst[0].header['HIERARCH OHP INS FIBER']
        texp=  hduii_lst[0].header['HIERARCH OHP CCD DKTM']
        dprtype = hduii_lst[0].header['HIERARCH OHP DPR TYPE']
        berv= hduii_lst[0].header['HIERARCH OHP DRS BERV']
        try:
            rhk_header.append(hdu[0].header['HIERARCH OHP DRS RHK'])
       	except KeyError:
            rhk_header.append(9999.99)
        try:
            sn1= hdu[0].header['HIERARCH OHP DRS CAL EXT SN1']
        except KeyError:
            Sn1.append(9999.99)
        BJD= hdu[0].header['HIERARCH OHP DRS BJD']
        hdu.close()
######################################################################################################e2dsA
        e2dsfits= replace_all(reduced_dir+night+'/'+str(ccf_collect['ccf'][ii]), d)
        hduii_lst = fits.open(e2dsfits)
        mode = hduii_lst[0].header['HIERARCH OHP INS FIBER']
        texp=  hduii_lst[0].header['HIERARCH OHP CCD DKTM']
        dprtype = hduii_lst[0].header['HIERARCH OHP DPR TYPE']
        berv= hduii_lst[0].header['HIERARCH OHP DRS BERV']
        bv= hduii_lst[0].header['HIERARCH OHP TARG BV']
        print('bv',bv)
        ccdsigdet= hduii_lst[0].header['HIERARCH OHP DRS CCD SIGDET']
        if abs(rv)>200.:
           print('activity index is not computed')
           Mode.append(mode)
           bjd.append(BJD)
           Rv.append(rv)
           DPRtype.append(dprtype)
           Prot.append(9999.99)
           Sn1.append(sn1)
           HAL.append(9999.99)
           Sig_Hal.append(9999.99)     
           Prot.append(9999.99)
           R_hk.append(9999.99)
           Sig_r_hk.append(9999.99)
           S_mw.append(9999.99)
           Sig_s.append(9999.99)
           continue

        e2dsff,nx,ny=read_data(e2dsfits)
        e2dsff = e2dsff.astype('d')
        raw = e2dsff.copy()
             
        wave,param_ll,param_x = get_e2ds_ll(e2dsfits)
        wave = wave * (1+berv/2.99792458e5)
        wave = wave / (1+rv/2.99792458e5)
        bin=(wave[-1]-wave[0])/(len(wave)-1)
        noise = ccdsigdet*np.sqrt(ic_extmeanzone)
        flat_file= hduii_lst[0].header['HIERARCH OHP DRS CAL FLAT FILE']
        hduii_lst.close() 
######################################################################################################blaze_file
        blaze_file= flat_file.replace('flat','blaze')
        blaze,fx,fy= read_data(reduced_dir+night+'/'+blaze_file)
        blaze = blaze.astype('d')

   
######################################################################################################e2dsB       
        e2dsBfits= e2dsfits.replace('A','B')
        e2dsBff,Bnx,Bny= read_data(e2dsBfits)
        e2dsBff = e2dsBff.astype('d')

        if dprtype.split(",")[1] == 'WAVE':
            if os.path.exists(e2dsBfits)==0:
                 print('error: File '+e2dsBfits+' doesnot exist')
            else:
                 e2dsBffc,back=correct_mins2d(e2dsBff,9,453)
                 e2dsff=e2dsff-1.*back


        if dprtype.split(",")[1] == 'SKY':
            if os.path.exists(e2dsBfits)==0:
                 print('error: File '+e2dsBfits+' doesnot exist')
            else:
                 e2dsBffc,back= correct_backs2d(e2dsBff,5000,niter=7,coeff=5)
                 e2dsff=e2dsff-1.*back

        if dprtype.split(",")[1] == 'FP':
            if os.path.exists(e2dsBfits)==0:
                 print('error: File '+e2dsBfits+' doesnot exist')
            else:
                 e2dsBffc,back=correct_backs2d_FP(e2dsBff)
                 e2dsff=e2dsff-1.*back


                 

######################################################################################################       
    
        for i in range(ny):
            e2dsff[i] = e2dsff[i]/blaze[i]*sum(blaze[i])/nx

        ##########  Flux relatif Halpha #########

        # Halpha 6562.808
        xstart=6562.469
        xend=6563.147
        nH =np.logical_and( wave[36]> 6562.469, wave[36] < 6563.147)
        response = 1.
        rawH =np.compress(nH,raw[36])
        blazeH = np.compress(nH,blaze[36])
        dataH=np.compress(nH,e2dsff[36])
        somH=sum(dataH)/len(dataH)
        sig_H = sum((rawH+noise**2)*response**2/blazeH**2)*sum(blaze[36])**2/nx**2/len(dataH)**2


        # V band 3901
        xstart=6545.495
        xend=6556.245
        nV =np.logical_and(wave[36] > 6545.495, wave[36] < 6556.245)
        dataV=np.compress(nV,e2dsff[36])                   # attention proche du bord de blaze
        rawV =np.compress(nV,raw[36])
        blazeV = np.compress(nV,blaze[36])
        somV=sum(dataV)/len(dataV)
        sig_fV = sum((rawV+noise**2)/blazeV**2)*sum(blaze[36])**2/nx**2/len(dataV)**2
        #somV=0.

        # R band 4001
        xstart=6575.934
        xend=6584.684
        nR =np.logical_and(wave[36] > 6575.934, wave[36] < 6584.684)
        dataR=np.compress(nR,e2dsff[36])
        rawR =np.compress(nR,raw[36])
        blazeR =np.compress(nR,blaze[36])
        somR=sum(dataR)/len(dataR)
        sig_fR = sum((rawR+noise**2)/blazeR**2)*sum(blaze[36])**2/nx**2/len(dataR)**2

        Hal= float(somH/(somV+somR))
        sig_Hal_square = (sig_H)/(somV+somR)**2 + (somH)**2*(sig_fV+sig_fR)/(somV+somR)**4
        sig_Hal= np.sqrt(sig_Hal_square)

#######################################CaIIHK#######################################

    ## Ca II K 3933.664
        nK = np.logical_and(wave[0] > 3932.574, wave[0] < 3934.754)
        response = -np.absolute(wave[0]-3933.664)/1.09 + 1.
        response = np.compress(nK,response)
        dataK1 = np.compress(nK,e2dsff[0])
        rawK1 = np.compress(nK,raw[0])
        blazeK1 = np.compress(nK,blaze[0])
        CaIIK1 = sum(dataK1*response)/sum(response)
        sig_CaIIK1 = sum((rawK1+noise**2)*response**2/blazeK1**2)*sum(blaze[0])**2/nx**2/sum(response)**2

        nK = np.logical_and(wave[1] > 3932.574, wave[1] < 3934.754)
        response = -np.absolute(wave[1]-3933.664)/1.09 + 1.
        response = np.compress(nK,response)
        dataK2 = np.compress(nK,e2dsff[1])
        rawK2 = np.compress(nK,raw[1])
        blazeK2 = np.compress(nK,blaze[1])
        CaIIK2 = sum(dataK2*response)/sum(response)
        sig_CaIIK2 = sum((rawK2+noise**2)*response**2/blazeK2**2)*sum(blaze[1])**2/nx**2/sum(response)**2

        CaIIK = (CaIIK1+CaIIK2)/2.
        sig_CaIIK = (sig_CaIIK1+sig_CaIIK2)/4.

    # Ca II H 3968.470
        nH = np.logical_and(wave[1] > 3967.380, wave[1] < 3969.560)
        response = -np.absolute(wave[1]-3968.470)/1.09 + 1.
        response = np.compress(nH,response)
        dataH = np.compress(nH,e2dsff[1])
        rawH = np.compress(nH,raw[1])
        blazeH = np.compress(nH,blaze[1])
        CaIIH = sum(dataH*response)/sum(response)
        sig_CaIIH = sum((rawH+noise**2)*response**2/blazeH**2)*sum(blaze[1])**2/nx**2/sum(response)**2

    # V band 3901
        nV = np.logical_and(wave[0] > 3891.07, wave[0] < 3911.07)
        dataV = np.compress(nV,e2dsff[0])
        rawV = np.compress(nV,raw[0])
        blazeV = np.compress(nV,blaze[0])
        fV = sum(dataV)/len(dataV)
        sig_fV = sum((rawV+noise**2)/blazeV**2)*sum(blaze[0])**2/nx**2/len(dataV)**2


    # R band 4001
        nR = np.logical_and(wave[2] > 3991.07, wave[2] < 4011.07)
        dataR = np.compress(nR,e2dsff[2])
        rawR = np.compress(nR,raw[2])
        blazeR = np.compress(nR,blaze[2])
        fR = sum(dataR)/len(dataR)
        sig_fR = sum((rawR+noise**2)/blazeR**2)*sum(blaze[2])**2/nx**2/len(dataR)**2
	
	# S_SOPHIE 
        s_sophie = (CaIIK + CaIIH) / (fV + fR)
        sig_s_square = (sig_CaIIK+sig_CaIIH)/(fV+fR)**2 + (CaIIK+CaIIH)**2*(sig_fV+sig_fR)/(fV+fR)**4
        sig_s= np.sqrt(sig_s_square)

        # distinguer cas HE HR
        if mode=='             HR':
           a=0.753559947
           b=0.01656492054

        if mode =='             HE':
           a=0.7196755409
           b=0.02145828307

        s_mw= (s_sophie/a)-(b/a)

        if ((bv=='UNKNOWN')&(bv =='999.9')):
            r_hk= 9999.99
            sig_r_hk= 9999.99
            prot= 9999.99

        else:

            if (bv > 0.44) & (bv < 1.20):

               Ccf = 1.13*bv**3 - 3.91*bv**2 + 2.84*bv - 0.47
               if bv < 0.63:
                 Ccf = Ccf + 0.135*(0.63-bv) - 0.814*(0.63-bv)**2 + 6.03*(0.63-bv)**3
               Ccf = 10.**Ccf
               r_hk_raw= 1.34e-4*Ccf*s_mw
               r_phot = -4.898 + 1.918*bv**2 - 2.893*bv**3
               r_phot= 10.**r_phot
               sig_r_hk_raw = 1.34e-4*0.962*Ccf*sig_s
               sig_r_hk = sig_r_hk_raw/(r_hk_raw-r_phot)/np.log(10.)

               if r_hk_raw-r_phot>0.:
                  r_hk= np.log10(r_hk_raw-r_phot)

               else:
                  print(r_hk_raw-r_phot)
                  r_hk= 9999.99
                  sig_r_hk= 9999.99
                  prot= 9999.99

# Rotation period

               if bv < 1:
                    tau = 1.362 - 0.166*(1-bv) + 0.025*(1-bv)**2 - 5.323*(1-bv)**3
               else:
                    tau = 1.362 - 0.14*(1-bv)
               prot_raw = 0.324 - 0.400*(5+r_hk) - 0.283*(5+r_hk)**2 - 1.325*(5+r_hk)**3 + tau
               prot = 10**prot_raw

            else:

              r_hk= 9999.99
              sig_r_hk= 9999.99
              prot=9999.99

        Mode.append(mode)
        R_hk.append(r_hk)
        bjd.append(BJD)
        Rv.append(rv)
        DPRtype.append(dprtype)
        Sig_s.append(sig_s)
        S_mw.append(s_mw)
        Sig_r_hk.append(sig_r_hk)
        Prot.append(prot)
        HAL.append(Hal)
        Sig_Hal.append(sig_Hal)
        Sn1.append(sn1)
#                print(R_hk,ii,rhk_header[ii],rhk_header,bjd,Rv,Sn1,DPRtype,Mode,Sig_s,S_mw,Sig_r_hk,Prot,HAL,Sig_Hal) 

#         if (sn1 < 10) and (thorium == 0):
#             print('warning S/N is low: result might be unreliable!')

#         if (sn1 < 30) & (thorium == 1):
#             print('warning: S/N is low and simultaneous ThAr is used: result might be unrel)

print(len(bjd),len(Rv),len(Sn1),len(DPRtype),len(Mode),len(Sig_s),len(S_mw),len(Sig_r_hk),len(R_hk),len(Prot),len(HAL),len(Sig_Hal))
p=open('activity_'+target+'.rdb','w')
p.write('jdb\tvrad\tsn1\tdprtype\tmode\tsig_s\ts_mw\tsig_r_hk\tr_hk\tprot\tHal\tsig_Hal\trhk_header\n')
for i in range (len(R_hk)):
    p.write('%.8f\t%.8f\t%.8f\t%s\t%s\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n' %(bjd[i],Rv[i],Sn1[i],DPRtype[i],Mode[i],Sig_s[i],S_mw[i],Sig_r_hk[i],R_hk[i],Prot[i],HAL[i],Sig_Hal[i],rhk_header[i]))
p.close()


bjd = [x - 2400000 for x in bjd]
p=open('activity_test.rdb','w')
p.write('jdb\tvrad\tsig_s\ts_mw\tsig_r_hk\tr_hk\tHal\tsig_Hal\n')
for i in range (len(R_hk)):
    p.write('%.5f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n' %(bjd[i],Rv[i],Sig_s[i],S_mw[i],Sig_r_hk[i],R_hk[i],HAL[i],Sig_Hal[i]))
p.close()

