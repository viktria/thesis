#!/usr/bin/python3

import subprocess
import os
import time
import os.path

input_dir = "/bigdisk/users/kelvi/data/AF_full"
for input_file in os.listdir(input_dir):
#wait until not so many running jobs
    p1 = subprocess.Popen("squeue", stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["grep","-c","kelvi"], stdin=p1.stdout, stdout=subprocess.PIPE)
    for s in (str(p2.communicate())[3:-10]).split('\\n'):
        r_job = int(s)
    while r_job > 25 :
        subprocess.Popen(["echo", "waiting"])
        time.sleep(10)
        p3 = subprocess.Popen("squeue", stdout=subprocess.PIPE)
        p4 = subprocess.Popen(["grep","-c","kelvi"], stdin=p3.stdout, stdout=subprocess.PIPE)
        for s in (str(p4.communicate())[3:-10]).split('\\n'):
            r_job = int(s)
    else:
        input_file_path = os.path.join(input_dir, input_file)
        if os.path.isfile(input_file_path):
            print("Processing file:", input_file_path)
            subprocess.Popen(['sbatch', 
                            '-e','/bigdisk/users/kelvi/error/' + input_file +  '.err',
                            '-o','/bigdisk/users/kelvi/out_file/' + input_file +  '.out',
                            '/bigdisk/users/kelvi/scripts/AF_full.py', input_file_path])


