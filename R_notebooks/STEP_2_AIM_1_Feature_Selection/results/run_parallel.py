from subprocess import PIPE
import subprocess
import time


def submit_job_max_len(job_list, max_processes):
  sleep_time = 0.1
  processes = list()
  for command in job_list:
    print ('running {n} processes. Submitting {proc}.'.format(n=len(processes), proc=str(command)))
    processes.append(subprocess.Popen(command, shell=True, stdout=None, stdin=PIPE))
    while len(processes) >= max_processes:
      time.sleep(sleep_time)
      processes = [proc for proc in processes if proc.poll() is None]
  while len(processes) > 0:
    time.sleep(sleep_time)
    processes = [proc for proc in processes if proc.poll() is None]


# cmd = '/bin/bash run_what.sh {n}'
#job_list = ((cmd.format(n=i)).split() for i in range(100))

job_list     = ['Rscript /users/dgibbs/proj/gbm_sex/feature_search_M.R 1 0.5',
		'Rscript /users/dgibbs/proj/gbm_sex/feature_search_M.R 2 0.5',
		'Rscript /users/dgibbs/proj/gbm_sex/feature_search_M.R 3 0.5',
		'Rscript /users/dgibbs/proj/gbm_sex/feature_search_M.R 4 0.5',
		'Rscript /users/dgibbs/proj/gbm_sex/feature_search_M.R 5 0.66',
		'Rscript /users/dgibbs/proj/gbm_sex/feature_search_F.R 1 0.5',
		'Rscript /users/dgibbs/proj/gbm_sex/feature_search_F.R 2 0.5',
		'Rscript /users/dgibbs/proj/gbm_sex/feature_search_F.R 3 0.66',
		'Rscript /users/dgibbs/proj/gbm_sex/feature_search_F.R 4 0.5',
		'Rscript /users/dgibbs/proj/gbm_sex/feature_search_F.R 5 0.5']

submit_job_max_len(job_list, max_processes=14)
