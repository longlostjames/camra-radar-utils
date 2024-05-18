import getopt, sys, os
import glob
import pyart

from datetime import datetime

#sys.path.append('/home/users/cjwalden/my-packages')
import camra_utils as camra

try:
    opts, args = getopt.getopt(sys.argv[1:], "d:i:c:o:", ["date=","inpath=","campaign=","outpath="])
except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized
        sys.exit(2)

data_date = datetime.now()
datestr = data_date.strftime('%Y%m%d')
inpath = '/gws/pw/j07/ncas_obs_vol2/cao/raw_data/ncas-radar-camra-1/data/'
outpath = '/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-radar-camra-1/'


for o, a in opts:
    if o == "-d":
        datestr = a;
    elif o == "-i":
        inpath = a;
    elif o == "-c":
        campaign = a.lower();
    elif o == "-o":
        outpath = a;
    else:
        assert False, "unhandled option"

dateyr = datestr[0:4]


rawpath = os.path.join(inpath,"campaign",campaign,"raw",datestr);
outpath = os.path.join(outpath,campaign,datestr);
os.chdir(rawpath);
files = [f for f in glob.glob('*raw.nc')]

print(files);

print(outpath);


if not os.path.exists(outpath):
    os.makedirs(outpath)

scan_types = ['ppi','rhi','fix','man']

for f in files:
    for elem in scan_types:
        scan_type = None
        if elem in f.lower():
            scan_type = elem
            break;

    print("{} {}".format(f, scan_type));


    Radar = camra.read_camra_raw(f);

    file_timestamp = datetime.strftime(Radar.metadata["time_coverage_start"],'%Y%m%d-%H%M%SZ');


    outfile = os.path.join(outpath,'ncas-radar-camra-1_{}_{}.nc'.format(file_timestamp,scan_type));

    pyart.io.write_cfradial(outfile, Radar, format='NETCDF4', time_reference=None);
