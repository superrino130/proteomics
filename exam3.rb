# import pylab
require 'matplotlib/pyplot'
plt = Matplotlib::Pyplot

# from pyteomics import tandem, pepxml, mzid, auxiliary as aux, pylab_aux as pa
require_relative 'rbteomics/tandem'
require_relative 'rbteomics/pepxml'
require_relative 'rbteomics/mzid'
# require_relative 'rbteomics/mass/mass'

# import pandas as pd
require 'pandas'
pd = Pandas
# import numpy as np
require 'numpy'
np = Numpy

require 'set'
require 'open-uri'
include Tandem

plt.figure()

tf = read('example.t.xml')

exit
