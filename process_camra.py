import getopt, sys, os

import datetime

#sys.path.append('/home/users/cjwalden/my-packages')
import camra_utils as camra

try:
    opts, args = getopt.getopt(sys.argv[1:], "d:i:o:c:", ["date=","inpath=","outpath="])
except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized
        sys.exit(2)

data_date = datetime.datetime.now()
datestr = data_date.strftime('%Y%m%d')

#tracking_tag = 'CFARR_0002';
tracking_tag = 'AMOF_20220922221548';

campaign = 'woest';

yaml_project_file = '/Users/cw66/campaigns/woest/woest_project.yml'
yaml_instrument_file = '/Users/cw66/campaigns/woest/woest_instruments.yml'


for o, a in opts:
    if o == "-d":
        datestr = a;
    elif o == "-i":
        inpath = a;
    elif o == "-o":
        outpath = a;
    elif o == "-c":
        campaign = a;
    else:
        assert False, "unhandled option"

inpath = '/Users/cw66/Data/woest/radar-camra/'
#inpath = os.path.join('/gws/pw/j07/ncas_obs_vol2/cao/raw_data/ncas-radar-camra-1/data/campaign',campaign,'mom');
outpath = '/Users/cw66/Data/ncas-radar-camra-1/'
#outpath = os.path.join('/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-radar-camra-1',campaign);

print(tracking_tag);

camra.process_camra(datestr,inpath,outpath,yaml_project_file,yaml_instrument_file,tracking_tag);

