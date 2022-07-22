#################################################
#	calculate master calib (following SP1)	#
#################################################

#########################
#	imports		#
#########################

#from BD_mod import *
from scipy.optimize import curve_fit
import scipy.special as sp
from lmfit.models import SkewedGaussianModel
model = SkewedGaussianModel()

from numpy import *
import matplotlib.pyplot as plt

from comp_func import *
import pdb

plt.switch_backend('agg')

#################################
#	paths and variables	#
#################################

#namelist = ['Gl104', 'Gl109', 'Gl133', 'Gl15A', 'Gl169', 'Gl169.1A', 'Gl172', 'Gl176', 'Gl184', 'Gl2', 'Gl208', 'Gl21', 'Gl226', 'Gl251', 'Gl270', 'Gl275.1', 'Gl310', 'Gl328', 'Gl338A', 'Gl338B', 'Gl353', 'Gl361', 'Gl378', 'Gl38', 'Gl397.1A', 'Gl410', 'Gl424', 'Gl436', 'Gl47', 'Gl48', 'Gl480', 'Gl521', 'Gl540', 'Gl552', 'Gl617A', 'Gl625', 'Gl649', 'Gl655', 'Gl671', 'Gl687', 'Gl694', 'Gl70', 'Gl720A', 'Gl728', 'Gl731', 'Gl745A', 'Gl793', 'Gl806', 'Gl842.2', 'Gl863', 'Gl895', 'Gl908', 'Gl96']
namelist = ['Gl873','Gl406','Gl410','Gl411','Gl15A','Gl686','Gl514','Gl251','Gl48','Gl96','Gl388','Gl205','Gl687','Gl725A','Gl133','Gl270','Gl338A','Gl378','Gl793']

rem=['G239-25'] #removed stars - bookkeeping

namefiable =  ['HD185144','HD9407','HD89269A','Gl15A','Gl514','Gl686']

filepath = '/net/GSP/users/pcortes/SOPHIE/'
saverms = '/net/GSP/users/pcortes/SOPHIE/RMS_master.txt'
savermscomp = '/net/GSP/users/pcortes/SOPHIE/RMS_master_comp.txt'

#namefiable=namelist	#to run on all seven stars

sint=3.5 	#[3.1] #arange(2,4,0.5)#3.1 	#Seuil de rejection des etoiles (RMS)
s=6 		# Nombre de points de la mediane
pondcons=1	#range(1,4)#[2]	# Ponderation des constantes fiables
pondvar=1

def gaus(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

def skewed(x, a, mu, sigmag, alpha):
    #normal distribution
    normpdf = (1/(sigmag*np.sqrt(2*math.pi)))*np.exp(-(np.power((x-mu),2)/(2*np.power(sigmag,2))))
    normcdf = (0.5*(1+sp.erf((alpha*((x-mu)/sigmag))/(np.sqrt(2)))))
    return 2*a*normpdf*normcdf


#########################
#	NAIRA master	#
#########################

test_bjd=[]
test_vr=[]
test_name=[]

int1=1
tour=1
con_or=len(namefiable)

rms_con=[]
echec=0
while int1!=0 or echec==1:

	rms=[]
	rms_old=[]
	rms_in=[]
	int1=0
	int2=0

			#Creation de la correction a partir de la liste des constantes

	aj=zeros(len(namefiable),'d')
	for tps in range(4):
		k=0
		switch=0
		for cons in namefiable:
			print(cons)
			pond=pondvar
			if k<=con_or-1: pond=pondcons 
#			else : pond=int((sint-rms_con[k])*10.0)
			#Lecture
			#con_dat = rv_corr(cons)
			#thor_dat = thorium(cons)
			#bjd1 = con_dat[:,0]
			#vr1 = con_dat[:,1]
			#drift = thor_dat[:,4]
			#err = thor_dat[:,3]
			#pdb.set_trace()
			bjd1,vr1 = readrdb2(filepath+cons+'_NAIRA_corr_th_2_cti.rdb',1,2)	
			print(bjd1)
			err,xxx = readrdb2(filepath+cons+'_NAIRA_corr_th_2_cti.rdb',3,3)
			vrc1=vr1-average(vr1)
			#pdb.set_trace()
			bjd1=compress(abs(array(vrc1))<0.07,bjd1)   # Set up 0.02 WHY!!!???
			#pdb.set_trace()
			vr1=compress(abs(array(vrc1))<0.07,vr1)
			#pdb.set_trace()
			err=compress(abs(array(vrc1))<0.07,err)
			#drift=compress(abs(array(vrc1))<0.02,drift)	already corrected!
			vrc1=(vr1-average(vr1))*1000.0+aj[k]
			err=err*1000.
			k+=1
			if switch == 0: #Initialisation
				bjdo=bjd1.tolist()
				#print(bjdo)
				vrco=vrc1.tolist()
				erro=err.tolist()
				switch=1
			else:           #Ajout et tri
				for i in range(len(bjd1)):
					#print(cons)
					#pdb.set_trace()
					if bjd1[i]<bjdo[0]: 
						j=0
					else:
						temp=compress(bjd1[i]>array(bjdo),range(len(bjdo)))
						j=max(temp)+1
					for ec in range(pond): bjdo.insert(j,bjd1[i])
					for ec in range(pond): vrco.insert(j,vrc1[i])
					for ec in range(pond): erro.insert(j,err[i])

			for i in range(len(bjd1)):
				if int(bjd1[i])>57518 and int(bjd1[i])<57525:
					test_bjd.append(bjd1[i])
					test_vr.append(vrc1[i])
					test_name.append(cons)

#	elif method==3: #Mediane glissante
		vrmoy=[]
		date=[]
		span=[]
#			do3=int(bjdo[0])
#			c3=0
#			vo3=0
		for i in range(uint(s/2),uint(len(bjdo)-s/2)):
#				d3=int(average(array(bjdo[i-s/2:i+s/2])))
#s				v3=int(median(array(vrco[i-s/2:i+s/2])))
			vrmoy.append(median(array(vrco[i-s/2:i+s/2])))
			date.append(average(array(bjdo[i-s/2:i+s/2])))
			span.append(bjdo[i+s/2]-bjdo[i-s/2])
		vrmoy.append(median(array(vrco[-s/2:])))
		date.append(bjdo[-1])
		span.append(bjdo[-1]-bjdo[-s/2])
		i=0
		vrmoy_m3=[]
		date_m3=[]
		while int(date[i])==int(date[i+1]): i+=1
		vrmoy_m3.append(median(vrmoy[:i]))
		date_m3.append(average(date[:i]))
		for el in vrmoy[i+1:]: vrmoy_m3.append(el)
		for el in date[i+1:]: date_m3.append(el)
		vrmoy=vrmoy_m3
		date=date_m3
		date.pop(0)
		vrmoy.pop(0)
#		vrmoy=vrmoy[::-1]
#		date=date[::-1]
#Reajustement vertical des donnees a posteriori
		k=0
		for cons in namefiable:
			#Reajustement
			bjd1,vr1 = readrdb2(filepath+cons+'_NAIRA_corr_th_2_cti.rdb',1,2)	
			vrc1=vr1-average(vr1)
#			if cons=='HD221354.rdb' : print mean(vrc1),aj[k]
			bjd1=compress(abs(array(vrc1))<0.07,bjd1)  ## PREVIOUSLY 0.02 WHYYYYYYYYY
			vr1=compress(abs(array(vrc1))<0.07,vr1)
			vrc1=(vr1-average(vr1))*1000.0+aj[k]
			vrcon=interp(bjd1,date,vrmoy)
			a=average(array(vrcon)-array(vrc1))
#			if cons=='HD221354.rdb' : break
			#	bjd1,vrc1=dacor(bjd1,vrc1)
#			a1,chi=ajuv(bjd1,array(vrc1)/1000.0,date,array(vrmoy)/1000.0,-0.005,0.005,0.0001)
#			if abs(a-a1*1000.0)>0.1 : print "Attention !",tps,cons,a,a1*1000.0
#			print tps,name,a,min(chi)	
			aj[k]=float(aj[k])+a
			k+=1
	bjdmin=min(bjdo)
	bjdmax=max(bjdo)
	dafin=range(int(bjdmin),int(bjdmax))
#	dafin=range(int(min(date)),int(max(date)))

	vrfin=interp(dafin,date,vrmoy)
	if tour==1:
		dold=dafin
		vrold=vrfin

	newname=[]
	for name in namelist:
		bjd,vr = readrdb2(filepath+name+'_NAIRA_corr_th_2_cti.rdb',1,2)	
		if len(vr)<10: continue
		vr=(vr-median(vr))*1000.0		
		bjd=compress(abs(array(vr))<20,bjd)
		vr=compress(abs(array(vr))<20,vr)
		vrc1=zeros(len(vr),'d')
		vrcold=zeros(len(vr),'d')
		j=0	
		for i in range(len(dafin)):
			while int(bjd[j])==dafin[i]:
				vrc1[j]=vr[j]-vrfin[i]
				j+=1
				if j==len(bjd): break
			if j==len(bjd): break
		j=0	
		for i in range(len(dold)):
			while int(bjd[j])==dold[i]:
				vrcold[j]=vr[j]-vrold[i]
				j+=1
				if j==len(bjd): break
			if j==len(bjd): break
		rms.append(std(vrc1))
		rms_old.append(std(vrcold))
		rms_in.append(std(vr))
		if std(vrc1)<sint: 
			namefiable.append(name) # ajout de l'etoile dans les constantes
			int1+=1
			rms_con.append(std(vrc1))
		else: 
#			print name,len(bjd),stin[3],stold[3],sta1[3]
			newname.append(name)
	newcons=[]
	echec=0
	for cons in namefiable:
		bjd,vr = readrdb2(filepath+cons+'_NAIRA_corr_th_2_cti.rdb',1,2)	
		#snr, xxx = readrdb2(filepath+cons+'_NAIRA_corr_th_2_cti.rdb',8,8)
                vr=(vr-average(vr))*1000.0
		#pdb.set_trace()
		bjd=compress(abs(array(vr))<70,bjd)  # WHY 20?!!?!?!?!
		#pdb.set_trace()
		vr=compress(abs(array(vr))<70,vr)
#		bjd,vr=dacor(bjd,vr)
		vrc1=zeros(len(vr),'d')
		vrcold=zeros(len(vr),'d')
		j=0
		for i in range(len(dafin)):
			while int(bjd[j])==dafin[i]:
				vrc1[j]=vr[j]-vrfin[i]
				j+=1
				if j==len(bjd): break
			if j==len(bjd): break	
		rms.append(std(vrc1))
		rms_con.append(std(vrc1))
		rms_in.append(std(vr))
#	sint=max(array(rms_con))
#	print "Nouveau seuil :", sint
	dold=dafin
	vrold=vrfin
	namelist=newname
	tour+=1
#				print rms,rms_in,rms_old

dafin_NAIRA = dafin
vrfin_NAIRA = vrfin
print(vrfin)
print "Total number of constants (NAIRA): ", len(set(namefiable))
print set(namefiable)

#plt.figure(2)
#plt.clf()
plt.subplots(1,1,figsize=(8,5))

plt.plot(bjdo,vrco,'ko')
plt.plot(dafin,vrfin,'r',linewidth=3)
plt.ylim([-20,20])
plt.xlim([55600,59300])	
plt.xlabel('BJD - 2 400 000 [days]')
plt.ylabel('VR [m/s]')
#plt.show()
plt.savefig('master_NAIRA_comb.png', bbox_inches='tight')
#plt.clf()
plt.subplots(1,1,figsize=(8,5))

plt.plot(bjdo,vrco,'ko')
plt.plot(dafin,vrfin,'r',linewidth=3)
plt.ylim([-20,20])
plt.xlim([57930,59300])	
plt.xlabel('BJD - 2 400 000 [days]')
plt.ylabel('VR [m/s]')
#plt.show()
plt.savefig('master_NAIRA_comb_dat2.png', bbox_inches='tight')
#plt.clf()

#plt.clf()
plt.subplots(1,1,figsize=(8,5))

plt.plot(bjdo,vrco,'ko')
plt.plot(dafin,vrfin,'r',linewidth=3)
plt.ylim([-40,40])
plt.xlim([55500,59300])	
plt.xlabel('BJD - 2 400 000 [days]')
plt.ylabel('RV [m/s]')
plt.axvline(55872,color='#1F77B4', label='lamp change')
plt.axvline(56274,color='#FF7F0E',label='double scrambler')
plt.axvline(56690,color='#2CA02C',label='lamp change')
plt.axvline(56730,color='#9467BD',label='new calib. unit')
plt.axvline(56941,color='#8C564B', label='lamp current')
plt.axvline(57435,color='#E377C2', label='new thermal reg.')
plt.axvline(58419,color='#7F7F7F', label='pressure leak')
plt.axvline(58431,color='#7F7F7F')	#end pressure leak
plt.legend(loc='lower right')
#plt.show()
plt.savefig('master_NAIRA_comb_dates.png', bbox_inches='tight')
#plt.clf()



#Ecriture des fichiers

p=open("vrconstantes_NAIRA_comb.rdb",'w')
for i in range(len(vrfin)):
    p.write('%.4f\t%.4f\n' %(float(dafin[i]),vrfin[i]))
p.close()

fdat=open('vrconstantes_test.rdb','w')
for i in range(len(test_bjd)):
	fdat.write('%s\t%.4f\t%.4f\n' %(test_name[i],float(test_bjd[i]),test_vr[i]))
fdat.close()
#################################################
#	correct from master and save file	#
#################################################

#namelist = ['G239-25', 'Gl104', 'Gl109', 'Gl133', 'Gl15A', 'Gl169', 'Gl169.1A', 'Gl172', 'Gl176', 'Gl184', 'Gl2', 'Gl208', 'Gl21', 'Gl226', 'Gl251', 'Gl270', 'Gl275.1', 'Gl310', 'Gl328', 'Gl338A', 'Gl338B', 'Gl353', 'Gl361', 'Gl378', 'Gl38', 'Gl397.1A', 'Gl410', 'Gl411', 'Gl424', 'Gl436', 'Gl47', 'Gl48', 'Gl480', 'Gl514', 'Gl521', 'Gl540', 'Gl552', 'Gl617A', 'Gl625', 'Gl649', 'Gl655', 'Gl671', 'Gl686', 'Gl687', 'Gl694', 'Gl70', 'Gl720A', 'Gl728', 'Gl731', 'Gl745A', 'Gl793', 'Gl806', 'Gl842.2', 'Gl863', 'Gl895', 'Gl908', 'Gl96']
#namelist = ['Gl205']
for name in namelist:
	pc=open(filepath+name+'_NAIRA_combmaster_cor_newact.rdb','w')
	pc.write('jdb\tvrad\tsvrad\tfwhm\tcontrast\tbis_span\tha\tsig_ha\trhk\tsig_rhk\ts_mw\tsig_s\tna\tsig_na\tca\tsig_ca\tsnr\n')
	pc.write('---\t----\t-----\t----\t--------\t--------\t--\t------\t---\t-------\t----\t-----\t--\t------\t--\t------\t---\n')

	bjdt_NAIRA,vrt_NAIRA=readrdb2(filepath+name+'_NAIRA_corr_th_2_cti.rdb',1,2)
	errt_NAIRA,bisst_NAIRA=readrdb2(filepath+name+'_NAIRA_corr_th_2_cti.rdb',3,6)
	fwhmt_NAIRA,contt_NAIRA=readrdb2(filepath+name+'_NAIRA_corr_th_2_cti.rdb',4,5)
        snr_NAIRA, xxx = readrdb2(filepath+name+'_NAIRA_corr_th_2_cti.rdb',8,8)
	j=0
	k=0
	for i in range(len(bjdt_NAIRA)):
		if int(bjdt_NAIRA[i])>dafin_NAIRA[-1]:
			crop_len_NAIRA=len(bjdt_NAIRA)-i
			bjdt_NAIRA = bjdt_NAIRA[:i]
			vrt_NAIRA = vrt_NAIRA[:i]
			bisst_NAIRA = bisst_NAIRA[:i]
			errt_NAIRA = errt_NAIRA[:i]
			fwhmt_NAIRA = fwhmt_NAIRA[:i]
			contt_NAIRA = contt_NAIRA[:i]
	                snr_NAIRA = snr_NAIRA[:i]

			print crop_len_NAIRA," observations beyond rv master end cropped for ", name
			break

	vr_adapt_NAIRA=array(interp(bjdt_NAIRA,dafin_NAIRA,vrfin_NAIRA))
	vrc_NAIRA=vrt_NAIRA-vr_adapt_NAIRA/1000.
        #pdb.set_trace()
	#####
	# reading Halpha values to add to file
	# also read Na and He (named as Ca for DACE to read it)

	Halpha, err_Halpha=readrdb2(filepath+'activity_'+name+'.rdb',11,12)
	bjd_Ha, xxx=readrdb2(filepath+'activity_'+name+'.rdb',1,1)
	bjd_Ha = bjd_Ha - 2400000  # so it is consistent with the NAIRA file
	Na, err_Na= 999.99*np.ones(len(Halpha)), 999.99*np.ones(len(Halpha)) #readrdb2('/home/mhobson/Documents/NAIRA/update_2019-04-01/analysis/act/newback/'+name+'_Halpha_testback.rdb',14,15)
	He, err_He= 999.99*np.ones(len(Halpha)), 999.99*np.ones(len(Halpha)) #readrdb2('/home/mhobson/Documents/NAIRA/update_2019-04-01/analysis/act/newback/'+name+'_Halpha_testback.rdb',16,17)
	#print(bjd_Ha)
	Ha_N=[999.99]*len(bjdt_NAIRA)
	err_Ha_N=[999.99]*len(bjdt_NAIRA)
	Na_N=[999.99]*len(bjdt_NAIRA)
	err_Na_N=[999.99]*len(bjdt_NAIRA)
	He_N=[999.99]*len(bjdt_NAIRA)
	err_He_N=[999.99]*len(bjdt_NAIRA)
        snr_N = [999.99]*len(bjdt_NAIRA)
	
        for k in range(len(bjdt_NAIRA)):
		for j in range(len(bjd_Ha)):
			if round(bjd_Ha[j],4)==round(bjdt_NAIRA[k],4):
				#pdb.set_trace()
				Ha_N[k] = Halpha[j]
				err_Ha_N[k] = err_Halpha[j]
				Na_N[k] = Na[j]
				err_Na_N[k] = err_Na[j]
				He_N[k] = He[j]
				err_He_N[k] = err_He[j]
                                #snr_N[k] = snr_NAIRA[j]
	####

	#####
	# reading Smw, R'hk values to add to file

	logRhk, Smw=readrdb2(filepath+'activity_'+name+'.rdb',9,7)
	bjd_Rhk, err_S=readrdb2(filepath+'activity_'+name+'.rdb',1,6)
	bjd_Rhk = bjd_Rhk - 2400000 # so it is consistent with NAIRA file

	logRhk_N=[999.99]*len(bjdt_NAIRA)
	Smw_N=[999.99]*len(bjdt_NAIRA)
	err_S_N=[0.]*len(bjdt_NAIRA)
	err_Rhk_N=[0.0001]*len(bjdt_NAIRA)

	for k in range(len(bjdt_NAIRA)):
		for j in range(len(bjd_Rhk)):
			if round(bjd_Rhk[j],4)==round(bjdt_NAIRA[k],4):
				logRhk_N[k] = logRhk[j]
				Smw_N[k] = Smw[j]
				err_S_N[k] = err_S[j]
	####

	
	for k in range(len(bjdt_NAIRA)):
		pc.write('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.9f\t%.6f\t%.6f\t%.9f\t%.6f\t%.4f\t%.9f\t%.4f\t%.9f\t%.1f\n' %(bjdt_NAIRA[k],vrc_NAIRA[k],errt_NAIRA[k],fwhmt_NAIRA[k],contt_NAIRA[k],bisst_NAIRA[k],Ha_N[k],err_Ha_N[k],logRhk_N[k],err_Rhk_N[k],Smw_N[k],err_S_N[k],Na_N[k],err_Na_N[k],He_N[k],err_He_N[k],snr_NAIRA[k]))

	pc.close()

#################################
#	rms calculation		#
#################################

f= open(saverms, 'w')

labvec1 = ['star','NAIRA_th rms','NAIRA_th w rms','NAIRA_master w rms']
labvec2 = ['----','------------','--------------','------------------']

f.write('\t'.join(labvec1)+'\n')
f.write('\t'.join(labvec2)+'\n')

for star in namelist:
	err,vr1 = readrdb2(filepath+star+'_NAIRA_corr_th_2_cti.rdb',3,2)	
	naira_rms = (np.std(vr1-np.mean(vr1)))*1000.
	naira_w_rms = (w_rms(vr1, err))*1000.

	naira_rv_corr, naira_rv_err = readrdb2(filepath+star+'_NAIRA_combmaster_cor_newact.rdb',2,3)
	naira_cor_w_rms = (w_rms(naira_rv_corr, naira_rv_err))*1000.

	f.write('\t'.join([star, '%.3f'%naira_rms, '%.3f'%naira_w_rms, '%.3f'%naira_cor_w_rms])+'\n')

f.close()



