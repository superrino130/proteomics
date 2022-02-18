require 'matplotlib/pyplot'
plt = Matplotlib::Pyplot

# from pyteomics import tandem, pepxml, mzid, auxiliary as aux, pylab_aux as pa
# require_relative 'rbteomics/pepxml'
# require_relative 'rbteomics/mzid'
require_relative 'rbteomics/usi'
require_relative 'rbteomics/pylab_aux'

require 'pandas'
pd = Pandas
require 'numpy'
np = Numpy

require 'set'
require 'open-uri'

plt.figure()

# spectrum = Usi.proxi(
#   'mzspec:PXD004732:01650b_BC2-TUM_first_pool_53_01_01-3xHCD-1h-R2:scan:41840',
#   'backend' => 'massive')
# peptide = 'WNQLQAFWGTGK'
spectrum = Usi.proxi(
  'mzspec:PXD005175:CRC_iTRAQ_06:scan:11803:VEYTLGEESEAPGQR/3',
  'backend' => 'jpost')
peptide = 'VEYTLGEESEAPGQR'

Pylab_aux.annotate_spectrum(spectrum, peptide, 'precursor_charge' => 2, 'backend' => 'spectrum_utils',
  'ion_types' => 'aby', 'title' => peptide)
plt.show()