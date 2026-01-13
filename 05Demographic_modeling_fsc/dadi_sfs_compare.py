#!/usr/bin/python
import matplotlib
matplotlib.use('Agg')

### for mac; only plot 2D sfs of data
import numpy
from numpy import array
import dadi
import os
import pylab
import glob

        
################ compare 2d SFS
pop_names = {
    "0": "BUG",
    "1": "ASH",
    "2": "COS",
    "3": "ERM",
    "4": "ALB",
    "5": "UTI"
}

obs_files = sorted(glob.glob("6Pop_jointDAFpop*.obs.dadi.fs"))
for obs_name in obs_files:
    ID = obs_name.replace("6Pop_jointDAFpop", "").replace(".obs.dadi.fs", "")
    sim_name = "modelFixT-031_maxL.gof_jointDAFpop%s.sim.dadi.fs" % ID
    fig_name = "jointDAFpop%s.svg" % ID
    
    obs = dadi.Spectrum.from_file(obs_name)
    sim = dadi.Spectrum.from_file(sim_name)
    
    pop_ids = ID.split("_")  # ["1", "0"]
    pop_labels = [pop_names[pop_ids[0]], pop_names[pop_ids[1]]]  # ["ASH", "BUG"]
    
    obs.pop_ids = pop_labels
    sim.pop_ids = pop_labels
    
    pylab.figure()
    dadi.Plotting.plot_2d_comp_multinom(sim, obs, vmin=1, resid_range=50)
    pylab.savefig(fig_name)
    pylab.close()
    
    print("Done: %s -> %s vs %s" % (fig_name, pop_labels[0], pop_labels[1]))





#dd = dadi.Misc.make_data_dict("COL_UTA.dadi")
#data = dadi.Spectrum.from_data_dict(dd, pop_ids=['UTA','COL'], projections=[40,39], polarized=False)
#data.to_file("COL_UTA.fs")


#import pylab
#pylab.figure()
#dadi.Plotting.plot_single_2d_sfs(fs, vmin=1)
#pylab.show()
##pylab.savefig('COL_UTA.png', dpi=150)


################ polt 1d SFS
#dd = dadi.Misc.make_data_dict("COL.dadi")
#data = dadi.Spectrum.from_data_dict(dd, pop_ids=['group1'], projections=[78], polarized=False)
#data.to_file("COL.fs")
#
#import pylab
#pylab.figure()
#dadi.Plotting.plot_1d_fs(data)
#pylab.savefig('COL.png', dpi=150)
