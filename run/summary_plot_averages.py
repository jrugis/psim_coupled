from __future__ import print_function

import matplotlib as mpl
mpl.use('Agg') # because there's no X11 display

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import sys
import struct
import ctypes
import time
import argparse


# flow_names should match the solution_values enum in global_defs.hpp with Qtot appended
FLOW_NAMES = ["Nal", "Kl", "Cll", "VOL", "Na", "K", "Cl", "HCO3", "H", "Va", "Vb", "Qtot"]


##################################################################
# ctypes lib for loading data
##################################################################

script_dir = os.path.abspath(os.path.dirname(__file__))
lib_path = os.path.join(script_dir, "_summary_plot_averages.so")
if not os.path.exists(lib_path):
    print("you must run 'make'")
    sys.exit(1)
_lib = ctypes.CDLL(lib_path)
_lib.load_summary_plot_data.restype = ctypes.c_int
_lib.load_summary_plot_data.argtypes = [
  ctypes.c_char_p,
  ctypes.c_int,
  ctypes.c_int,
  np.ctypeslib.ndpointer(dtype=np.float32),
  np.ctypeslib.ndpointer(dtype=np.float32),
  ctypes.c_int,
  np.ctypeslib.ndpointer(dtype=np.int32),
  ctypes.c_int,
  np.ctypeslib.ndpointer(dtype=np.int32),
]

##################################################################
# functions
##################################################################

def get_data_fluid_flow(fname, ntime):
  if not os.path.exists(fname):
    return None

  data = np.fromfile(fname, dtype=np.float32)  # load binary data
  data = np.reshape(data, (ntime, len(FLOW_NAMES)))  # reshape
  df = pd.DataFrame(data=data, columns=FLOW_NAMES)  # pandas dataframe with named columns

  return df

## get the time values associated with the saved data
def get_time_vals(fname):
  f = open(fname + ".dat", "r") # get the saved data stride
  for line in f:
    if "Tstride" in line:
      tstride = int(line.rstrip().split()[1])
      break
  f.close()
  tvals = [] # list of time values
  n = int(0)
  f = open(fname + ".out", "r")
  for line in f:
    if "<Acinus> t:" in line: # iterate through the saved time step values
      if(float(line.rstrip().split()[2]) == 0.0): continue  # skip start time 0.0
      n += 1
      if (n % tstride) == 0:
        tvals.append(float(line.rstrip().split()[2])) # add to list of time values
  f.close()
  return np.asarray(tvals) # return time values in a numpy array

# get cell node count
def get_node_count(fname):
  f = open(fname, "r")
  for line in f:
    if line.startswith("<CellMesh> vertices_count:"):
      v = line.rstrip().split()
  f.close()
  return int(v[2])

def load_results(dname, value, num_time, num_nodes):
    # load full results
    res = np.fromfile(f"{dname}_{value}.bin", dtype=np.float32)
    res = res.reshape(num_time, num_nodes)

    return res

# get average values at apical and basal cells for the given result file
def get_apical_basal_averages(res, apical_nodes, basal_nodes):
    # averages at apical/basal
    res_apical = res[:, apical_nodes].mean(axis=1)
    res_basal = res[:, basal_nodes].mean(axis=1)

    return res_apical, res_basal

##################################################################
# main program
##################################################################
def main():
    # parse command line arguments
    parser = argparse.ArgumentParser(description="Create summary plot of psim5 simulation")
    parser.add_argument("--no-ca", action="store_true", help="Don't plot ca")
    parser.add_argument("--no-cer", action="store_true", help="Don't plot cer")
    parser.add_argument("--no-ip3", action="store_true", help="Don't plot ip3")
    parser.add_argument("--no-ion", action="store_true", help="Don't plot ion results")
    #parser.add_argument("--no-ffr", action="store_false", help="(NOT WORKING) Don't plot fluid flow rate")
    parser.add_argument("--time-start", type=int, default=0, help="Starting time to display on graphs (default=0)")
    parser.add_argument("--cells", type=int, nargs='*', default=[1, 2, 3, 4, 5, 6, 7], help="Cells to plot (default is 1 2 3 4 5 6 7)")
    parser.add_argument("--nintra", type=int, default=8, help="Number of intracellular variables (default is 8)")
    parser.add_argument("--font-size", type=int, default=16, help="Font size for matplotlib (default is 16)")
    parser.add_argument("-o", "--output", default="summary_plot2.pdf", help="File name to save plot to (default is summary_plot2.pdf)")
    args = parser.parse_args()

    plt.rcParams.update({'font.size': args.font_size})
    ncells = len(args.cells)
    print("plotting cells:", args.cells)
    nintra = args.nintra
    dtypes = []
    if not args.no_ca:
        dtypes.append("ca")
    if not args.no_cer:
        dtypes.append("cer")
    if not args.no_ip3:
        dtypes.append("ip3")
    if not args.no_ion:
        dtypes.append("ion")
    #plot_ffr = False if args.no_ffr else True
    plot_ffr = False

    # get the x-axis time values
    x = get_time_vals("a1")

    print("create summary plot")

    nplots = len(dtypes) + 1 if plot_ffr else len(dtypes)
    if "ion" in dtypes:
        nplots += 4
    fig, plots = plt.subplots(nplots, ncells, sharex='col', squeeze=False)
    plt.subplots_adjust(wspace = 0.5)
    #fig.set_size_inches(ncells * 3.8, nplots * 2.5)
    fig.set_size_inches(ncells * 7.6, nplots * 5.0)
    fig.text(0.02, 0.96, os.getcwd(), fontsize=20)
    ylabels = ["" for _ in range(nplots)]

    # main loop over plots
    for celli, cell in enumerate(args.cells):
      dname = "a1c" + str(cell)
      print("plotting {}".format(dname))
      nodes = get_node_count(dname + ".out")
      if celli == 0 and plot_ffr:
        pass
      else:
        plots[-1, celli].set_xlabel(" time (s)")
      plots[0, celli].set_title("Cell " + str(cell))

      # load ca and ip3 data
      if "ca" in dtypes or "cer" in dtypes or "ip3" in dtypes:
          # first we load the list of apical and basal nodes
          apical_nodes = np.fromfile("apical_nodes_{}.bin".format(dname), dtype=np.int32)
          basal_nodes = np.fromfile("basal_nodes_{}.bin".format(dname), dtype=np.int32)
          print("  num apical nodes = {0}, num basal nodes = {1}".format(len(apical_nodes), len(basal_nodes)))

          # now load the result and calculate averages at apical and basal nodes
          if "ca" in dtypes:
              # for calcium
              ca_apical = np.empty(x.shape[0], np.float32)
              ca_basal = np.empty(x.shape[0], np.float32)
              filename = dname + "_ca.bin"
              status = _lib.load_summary_plot_data(filename.encode("utf-8"), x.shape[0], nodes, ca_apical, ca_basal,
                      len(apical_nodes), apical_nodes, len(basal_nodes), basal_nodes)
              if status < 0:
                  sys.exit("Error open result file: {}".format(result_file))
              if status > 0:
                  sys.exit("Error reading row in result file: {}".format(status))

          if "cer" in dtypes:
              # for cer
              cer_apical = np.empty(x.shape[0], np.float32)
              cer_basal = np.empty(x.shape[0], np.float32)
              filename = dname + "_cer.bin"
              status = _lib.load_summary_plot_data(filename.encode("utf-8"), x.shape[0], nodes, cer_apical, cer_basal,
                      len(apical_nodes), apical_nodes, len(basal_nodes), basal_nodes)
              if status < 0:
                  sys.exit("Error open result file: {}".format(result_file))
              if status > 0:
                  sys.exit("Error reading row in result file: {}".format(status))

          if "ip3" in dtypes:
              # for ip3
              ip_apical = np.empty(x.shape[0], np.float32)
              ip_basal = np.empty(x.shape[0], np.float32)
              filename = dname + "_ip3.bin"
              status = _lib.load_summary_plot_data(filename.encode("utf-8"), x.shape[0], nodes, ip_apical, ip_basal,
                      len(apical_nodes), apical_nodes, len(basal_nodes), basal_nodes)
              if status < 0:
                  sys.exit("Error open result file: {}".format(result_file))
              if status > 0:
                  sys.exit("Error reading row in result file: {}".format(status))

      # plot ca
      if "ca" in dtypes:
          ca_index = dtypes.index("ca")
          plots[ca_index, celli].plot(x, ca_apical, color='blue', label='apical')
          plots[ca_index, celli].plot(x, ca_basal, color='red', label='basal')
          plots[ca_index, celli].legend(loc='best')
          ylabels[ca_index] = "ca ($\mu$M)"

      # plot cer
      if "cer" in dtypes:
          cer_index = dtypes.index("cer")
          plots[cer_index, celli].plot(x, cer_apical, color='blue', label='apical')
          plots[cer_index, celli].plot(x, cer_basal, color='red', label='basal')
          plots[cer_index, celli].legend(loc='best')
          ylabels[cer_index] = "cer ($\mu$M)"

      # plot ip3
      if "ip3" in dtypes:
          ip_index = dtypes.index("ip3")
          plots[ip_index, celli].plot(x, ip_apical, color='blue', label='apical')
          plots[ip_index, celli].plot(x, ip_basal, color='red', label='basal')
          plots[ip_index, celli].legend(loc='best')
          ylabels[ip_index] = "ip3 ($\mu$M)"

      # load fluid flow results
      if plot_ffr or "ion" in dtypes:
          iondf = get_data_fluid_flow(f"{dname}_ion.bin", x.shape[0])
          # if fluid flow was disabled don't plot it
          if iondf is None:
              print("Warning: disabling fluid flow and volume plotting since no lumen results file exists")
              plot_ffr = False
              if "ion" in dtypes:
                  dtypes.remove("ion")

      # plot flow results
      if "ion" in dtypes:
          flow_start_index = dtypes.index("ion")

          # first plot volume
          plots[flow_start_index+0, celli].plot(x, iondf["VOL"], color='blue', label="volume")
          plots[flow_start_index+0, celli].legend(loc='best')
          ylabels[flow_start_index+0] = "volume ($\mu$m$^3$)"
          # next plot Nal, Kl, Cll
          plots[flow_start_index+1, celli].plot(x, iondf["Nal"], color='blue', label="Nal")
          plots[flow_start_index+1, celli].plot(x, iondf["Kl"], color='red', label="Kl")
          plots[flow_start_index+1, celli].plot(x, iondf["Cll"], color='green', label="Cll")
          plots[flow_start_index+1, celli].legend(loc='best')
          ylabels[flow_start_index+1] = "?"
          # next plot Na, K, Cl
          plots[flow_start_index+2, celli].plot(x, iondf["Na"], color='blue', label="Na")
          plots[flow_start_index+2, celli].plot(x, iondf["K"], color='red', label="K")
          plots[flow_start_index+2, celli].plot(x, iondf["Cl"], color='green', label="Cl")
          plots[flow_start_index+2, celli].legend(loc='best')
          ylabels[flow_start_index+2] = "?"
          # next plot Va, Vb
          plots[flow_start_index+3, celli].plot(x, iondf["Va"], color='blue', label="Va")
          plots[flow_start_index+3, celli].plot(x, iondf["Vb"], color='red', label="Vb")
          plots[flow_start_index+3, celli].legend(loc='best')
          ylabels[flow_start_index+3] = "?"
          # next plot FFR (Qtot??)
          plots[flow_start_index+4, celli].plot(x, iondf["Qtot"], color='blue', label="Qtot")
          plots[flow_start_index+4, celli].legend(loc='best')
          ylabels[flow_start_index+4] = "Qtot (QFFR?)"

    # setting x axis (time) lower limit
    for (m,n), subplot in np.ndenumerate(plots):
        subplot.set_xlim(args.time_start)

    # setting y axis labels
    for i in range(nplots):
        plots[i, 0].set_ylabel(ylabels[i])

    if plot_ffr:
      for i in range(1, ncells):
        plots[-1, i].axis('off')
        plots[-2, i].xaxis.set_tick_params(which='both', labelbottom=True)
      plots[-1, 0].set_title("Total fluid flow rate")
      plots[-1, 0].plot(x, ffr, color='blue')
      plots[-1, 0].set_xlabel(" time (s)")
      plots[-1, 0].set_ylabel(" Flow rate ($\mu m^3$/sec)")

    fig.savefig(args.output)


if __name__ == "__main__":
    main()
