#################################################
#	Thorium correction for NAIRA rvs	#
#################################################

#########################
#	imports		#
#########################

import numpy as np
from comp_func import *

#################################
#	paths and variables	#
#################################

namelist = ['Gl873','Gl406','Gl410','Gl251','Gl388','Gl48','Gl96','Gl411','Gl687','Gl725A','Gl133','Gl251','Gl270','Gl338A','Gl378','Gl617A','Gl793','Gl205','HD185144', 'HD9407', 'HD89269A', 'Gl514', 'Gl686', 'Gl15A']

#namelist = ['G239-25', 'G244-47', 'Gl104', 'Gl105B', 'Gl107B', 'Gl109', 'Gl133', 'Gl134', 'Gl150.1A', 'Gl150.1B', 'Gl15A', 'Gl15B', 'Gl162', 'Gl169', 'Gl169.1A', 'Gl172', 'Gl176', 'Gl184', 'Gl2', 'Gl208', 'Gl21', 'Gl226', 'Gl251', 'Gl270', 'Gl275.1', 'Gl310', 'Gl328', 'Gl338A', 'Gl338B', 'Gl353', 'Gl361', 'Gl378', 'Gl38', 'Gl397.1A', 'Gl410', 'Gl411', 'Gl424', 'Gl436', 'Gl47', 'Gl48', 'Gl480', 'Gl49', 'Gl4A', 'Gl4B', 'Gl514', 'Gl521', 'Gl540', 'Gl552', 'Gl617A', 'Gl625', 'Gl63', 'Gl649', 'Gl655', 'Gl671', 'Gl686', 'Gl687', 'Gl694', 'Gl70', 'Gl720A', 'Gl728', 'Gl731', 'Gl745A', 'Gl793', 'Gl806', 'Gl842.2', 'Gl863', 'Gl87', 'Gl895', 'Gl908', 'Gl96', 'HD9407', 'HD221354','HD185144', 'HD89269A']

#namelist = ['Gl205']


save_path = '/net/GSP/users/pcortes/SOPHIE/'
obsl_path = '/net/GSP/users/pcortes/SOPHIE/'
dat_path =  '/net/GSP/users/pcortes/SOPHIE/'

#################################
#	apply correction	#
#################################

for star in namelist:
	obs = obsl_path+'obslist_'+star+'.txt'
	print(obs)
	dat = thorium_3(star, obs)	#uses secular acceleration correction when it exists
	#print(dat)
	bjd = dat[:,0]
	vrad = dat[:,1]
	err = dat[:,3]
	sn_26 = dat[:,5]
	sn_35 = dat[:,6]
        #print(bjd, vrad, err, sn_26, sn_35)
	nairafile = dat_path+star+'_cti_w7_s10_c5__NAIRA_v1.rdb'
	fwhm, contr = readrdb2(nairafile,6,7)
	biss, bjdo = readrdb2(nairafile,8,1)
	pc=open(save_path+star+'_NAIRA_corr_th_2_cti.rdb','w')
	pc.write('jdb\tvrad\tsvrad\tfwhm\tcontrast\tbis_span\tolddate\tsn_26\tsn_35\n')
	pc.write('---\t----\t-----\t----\t--------\t--------\t-------\t-----\t-----\n')	 
	for k in range(len(bjd)):
		for j in xrange(0, len(bjdo)):
			print(bjd[k], bjdo[j])
			if (round(bjd[k],5)==round(float(bjdo[j]),5)):
				pc.write('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n' %(bjd[k],vrad[k],err[k], fwhm[j],contr[j],biss[j],bjdo[j],sn_26[k],sn_35[k]))
			else: continue	
	pc.close()
