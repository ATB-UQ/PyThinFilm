ssh raijin << 'END'
module load python

cd /home/564/mzs564/phaseTransition
python deposition_batch_daemon.py -i CBP-10nm-CBP-10nm-450K-slow.yml -c raijin
END

