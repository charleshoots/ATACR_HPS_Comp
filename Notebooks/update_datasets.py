import os
import subprocess
from pathlib import Path
from time import gmtime
def open_terminal(directory):
    if not os.path.isdir(directory):
        raise ValueError(f"The directory {directory} does not exist.")
    command = f'osascript -e \'tell application "Terminal" to do script "cd \\"{directory}\\""\''
    subprocess.run(command, shell=True, check=True)

ProjectFolder = Path('~/Documents/GalapagosRecovery_Test')
ShareMount = Path('/Volumes/cruise/SR2410')
rsync124 = ShareMount / 'multibeam-em124' / 'rawdata'
rsync712 = ShareMount / 'multibeam-em712' / 'rawdata'
rsync_dir = rsync124
current_date = ''.join([str(gmtime().tm_year),str(gmtime().tm_mon).zfill(2),str(gmtime().tm_mday).zfill(2)])
dayfolder = ProjectFolder / ('Cocos_Ridge_' + current_date)
rsync_cmd = 'rsync -arv /Volumes/cruise/SR2410/multibeam-em124/rawdata/*{date}* {dest}'.format(date=current_date,dest=dayfolder)
cp_cmd_files = 'cp ~/Documents/Scripts/*.cmd ' + str(dayfolder)


if not dayfolder.is_dir():
  dayfolder.mkdir(parents=True,exist_ok=True)
os.system(rsync_cmd)
os.system('open ' + str(dayfolder))

os.system(cp_cmd_files)
os.system('./preprocess.cmd')

directory_path = str(dayfolder)
open_terminal(directory_path)