#!/bin/sh
#$ -S /bin/sh
#$ -cwd 
## #$ -l month
## #$ -l medium
#$ -l s_vmem=8G,mem_req=8G
#$ -l short
### #$ -l debug
#$ -pe def_slot 1
#$ -o ./log
#$ -e ./log

source ~/.bashrc

cd /Users/takeshik/study/scripts/CompMSPeaks

# cat ../lists/shokuhin.v2.taxid.mod.list |grep natural |grep Viridiplantae >../lists/shokuhin.v2.taxid.vplant.all.list

ruby table_compmspeaks.rb 600 700 ./shokuhin.v2.taxid.mod.list   p ./test_dir  >./test_output_table_compmspeaks.out

# 
