#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# a file called sigma_history.csv needs to be read in, this file has
# AA,AB,BC,CC
#9.283210432623473,6.763211467579521,5.631940903351083,1.7123906905688424
#9.88083069938394,6.774568176135562,5.631940903351083,1.7839581716260162
#9.876952876357493,7.406732015011208,5.631940903351083,1.8651928252806624
#9.950843801080307,7.366710440159406,5.631940903351083,1.8651928252806624
# and I want to plot them all on the same graph with different colors and a legend

#out_dir = os.getcwd + '/output_analysis/pairsampler_results/'
# take the path to file as input argument
# python inpyt argument for the filename
import sys
if len(sys.argv) != 2:
    print("Usage: python plot_sigmas.py <filename>")
    sys.exit(1)
filename = sys.argv[1]
# if no filename is provided, use the default
# filename = 'output_analysis/pairsampler_results/sigma_history.csv'
# check if the file exists
if not os.path.isfile(filename):
    print(f"File {filename} does not exist.")
    sys.exit(1)

# read in the data
df = pd.read_csv(filename)
plt.plot(df['AA'], label='AA', color='blue')
plt.plot(df['AB'], label='AB', color='orange')
plt.plot(df['BC'], label='BC', color='green')
plt.plot(df['CC'], label='CC', color='red')

plt.show()