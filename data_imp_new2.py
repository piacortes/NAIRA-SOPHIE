#	new data import script - manual single-star version

#remember to delete previous run directories to free up space!

#doesn't keep s1d files as they're only used for Svalue and Halpha which isn't implemented on NAIRA for OHP, and I don't have the space for them!

#SP1 constants: HD185144, HD9407, HD221354, HD89269A

#list of active stars from Javiera
#ACTIVES WITH MORE THAN 10 SOPHIE+ POINTS:
#HD99491A
#HD101242
#HD115043 -> Anticorr. RV vs. bisector
#HD115404A
#HD115755
#HD133003 -> Corr RV with Halpha
#HD142093
#HD146362
#HD150706 -> SOPHIE-V target (Boisse et al. 2012). Corr. between RV, bisector, Ha and Sindex. KECK doesn't see the same.
#HD154345 -> logR'HK= -4.79 published as planet (Boisse et al. 2012) but the KECK guys say it's activity
#HD76780 -> might be some corn with logR'HK if you clean the data and with fwhm
#HD210460 -> was highly correlated with Halpha. Some Keck guys said (via email) they see a correlation with calcium (that we don't see).

#ACTIVES WITH LESS THAN 10 SOPHIE+ POINTS:
#HD73344
#HD75302

# UPDATES July 2021: using pfits instead of pcfitsio

#import fitsio
import os
import shutil
from astropy.io import fits as pf
import pdb
import argparse
#################################
#	files and directories	#
#################################

list_dir = '/net/GSP/users/pcortes/SOPHIE/'	#where the obslist files are stored
datapath = '/net/GSP/spirou/sophie/SOPHIE_data/reduced/'	#all data files
#listingpath='/net/GSP/spirou/sophie/SOPHIE_data/Neda_SOPHIE/listing_SOPHIE/'		#Chemin des listings
listingpath='/net/GSP/users/pcortes/SOPHIE/listing_SOPHIE/'


# Create the parser
parser = argparse.ArgumentParser(description='Target name')

# Add the arguments
parser.add_argument('target', type=str, help='the name of the star you want to create the obslists')

# Execute the parse_args() method
args = parser.parse_args()

starname = args.target

#########################
#	make obslist	#
#########################

#os.system('more '+listingpath+'Listing_ALL_201* | grep -w '+starname+' | grep -v sig | grep HR > '+list_dir+'obslist_'+starname+'.txt') #le grep -v sig permet d'eviter des bugs sur HD185144(sig, Dra)

os.system('grep -w '+starname+' '+listingpath+'Listing_ALL_20* | grep -v sig | grep HR > '+list_dir+'obslist_'+starname+'.txt') #le grep -v sig permet d'eviter des bugs sur HD185144(sig, Dra)

#grep -w Gl521 /home/altesse1NS/sophie/NEWLISTING3.1/listing_ALL_201* | grep -v sig | grep HR 
save_dir = list_dir+'data_'+starname

if os.path.exists(save_dir):
	os.system('rm -r '+save_dir)

os.mkdir(save_dir)

obslist = list_dir+'obslist_'+starname+'.txt'

f=open(obslist,"r")
cat=[]
for line in f:
    #print(line)
    #print(line[184:195])
    #if line[184:195]>55730:
    #cat.append(line.strip('\n'))
    cat.append(line.replace('\t',' '))
    #print(cat)
    #pdb.set_trace()
f.close()

#################################
#	define functions	#
#################################

def spec_cp(file_list, savedir, file_type):	#copies the spectra
	numb_files = 0.
	for line in file_list:
		try:
			shutil.copy(line, savedir)
			numb_files = numb_files + 1
		except Exception as e:
			print(e)
			line_new = line.replace(line[43:53],line[61:71])
			shutil.copy(line_new, savedir)
			numb_files = numb_files + 1
	print(str(numb_files)+' '+file_type+' files copied')


def ccf_cp(file_list, savedir, file_type):	#finds the correct ccf and copies it
	numb_files = 0.
	for line in file_list:
		line_M5 = line.replace('*','M5')
		line_K5 = line.replace('*','K5')
		line_G2 = line.replace('*','G2')
		if os.path.exists(line_M5):
			shutil.copy(line_M5, savedir)
			numb_files = numb_files + 1
		elif os.path.exists(line_K5):
			shutil.copy(line_K5, savedir)
			numb_files = numb_files + 1
		elif os.path.exists(line_G2):
			shutil.copy(line_G2, savedir)
			numb_files = numb_files + 1

	print(str(numb_files)+' '+file_type+' files copied')


def read_key(fitsfilename,keyname,hdu=-1):	#from hadmrFITS - can't figure out how to just import it!

	"Read the FITS file and return the value of the keyword."
  
	status=0
	p= fitsio.open(fitsfilename,0)
	if hdu!=-1:
		fitsio.movabs_hdu(p,hdu)
	n=fitsio.fits_get_hdrpos(p)[0]
	for i in range(n):
		rec=fitsio.fits_read_record(p,i+1)
		if rec[0:len(keyname)] == keyname:
			status=1
  
	if status:
		value=fitsio.fits_read_keyword(p,keyname)[0]
	else:
		value="'UNKNOWN'"
     
	fitsio.fits_close_file(p)
   
	if value[0]=="'":
		value=value[1:len(value)-1]
	else:  
		value=float(value)
	return value
	
def read_key2(fitsfilename, keyname):
        p = pf.open(fitsfilename)
        try:
            value = p[1].header[keyname]
        except Exception as e:
            value = p[0].header[keyname]
        return value
#################################################
#	turn list into lists of files to copy	#
#################################################

#l=open(obslist,'r')	#open the file

#initialise lists

listing_e2ds = []
listing_ccf = []
listing_blaze =[]
listing_wave =[]
listing_s1d = []

for line in cat:
	print(line)
	#pdb.set_trace()
	#name = line[89:119]	#keep only the part with the file name
	name = line[72:102]
	rv_drift = line[200:207]
	vrad0 = line[193:200]
	vrad = line[220:227]
	print(vrad, vrad0, rv_drift)
	#pdb.set_trace()
	date = name[7:17]
	#pdb.set_trace()
	print(name)
	print(date)
	name_e2ds = datapath  + str(date) + '/' + name.replace('.rdb:','/') + '_e2ds_A.fits'
	name_ccf = datapath + str(date) + '/' +  name.replace('.rdb:','/') + '_ccf_*_A.fits'
	name_calib = name_e2ds[0:54]
	#print(name_e2ds)
	#print(name_ccf)
	#print(name_calib)
	#pdb.set_trace()
	try:
		flat_file = read_key2(name_e2ds, 'HIERARCH OHP DRS CAL FLAT FILE')
		wave_file = read_key2(name_e2ds, 'HIERARCH OHP DRS CAL TH FILE')
		blaze_file = flat_file.replace('flat', 'blaze')
		listing_e2ds.append(name_e2ds)
		listing_ccf.append(name_ccf)
		listing_wave.append(name_calib+wave_file)
		listing_blaze.append(name_calib+blaze_file)
		#pdb.set_trace()
		name_s1d = datapath + str(date) + '/'+ name.replace('.rdb:','/') + '_s1d_A.fits'
		listing_s1d.append(name_s1d)
	except Exception as e1:
		try:
			date = line[-11:-1] #line[-10:]
			#print(date)
			#pdb.set_trace()
			name_e2ds = datapath  + str(date) + '/' + name.replace('.rdb:','/') + '_e2ds_A.fits'
			name_ccf = datapath + str(date) + '/' +  name.replace('.rdb:','/') + '_ccf_*_A.fits'
			name_calib = name_e2ds[0:54]
			flat_file = read_key2(name_e2ds, 'HIERARCH OHP DRS CAL FLAT FILE')
			wave_file = read_key2(name_e2ds, 'HIERARCH OHP DRS CAL TH FILE')
			blaze_file = flat_file.replace('flat', 'blaze')
			listing_e2ds.append(name_e2ds)
			listing_ccf.append(name_ccf)
			listing_wave.append(name_calib+wave_file)
			listing_blaze.append(name_calib+blaze_file)
			#pdb.set_trace()
			name_s1d = datapath + str(date) + '/'+ name.replace('.rdb:','/') + '_s1d_A.fits'
			listing_s1d.append(name_s1d)
			print(e1)
			print('File found')
		except Exception as e2:
        		print(e2)
		
#print name_e2ds	

#l.close()


#########################
#	copy files	#
#########################

print('obtaining observations for '+starname)

spec_cp(listing_e2ds, save_dir, 'e2ds')
ccf_cp(listing_ccf, save_dir, 'ccf')
spec_cp(listing_blaze, save_dir, 'blaze')
spec_cp(listing_wave, save_dir, 'wave')
#spec_cp(listing_s1d, save_dir, 's1d')

