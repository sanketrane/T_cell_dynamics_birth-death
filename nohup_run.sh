sh sh_files/run_neutral0.sh > log_files/log_n0.out 2> log_files/log_n0.err &

sh sh_files/run_neutral1.sh > log_files/log_n1.out 2> log_files/log_n1.err &

sh sh_files/run_neutral3.sh > log_files/log_n3.out 2> log_files/log_n3.err &

nohup sh sh_files/run_neutral7.sh > log_files/log_n7.out 2> log_files/log_n7.err 
