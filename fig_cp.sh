#########################################################################
# File Name: fig_cp.sh
# Author: Neo
# mail: liuniu@smail.nju.edu.cn
# Created Time: Sat Jun 30 22:56:59 2018
#########################################################################
#!/bin/bash

figs=("dec_nf_sf.eps" \
	  "decimation-offset-dist.eps" \
	  "error-SBL.eps" \
	  "inflation_plot.eps" \
	  "sess_num_dependent_sf.eps")

for fig in ${figs[@]}
do
	cp ../plots/${fig} ../notes/20180530/figures/
done
