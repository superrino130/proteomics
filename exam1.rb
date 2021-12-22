# import os
# from urllib.request import urlretrieve
# import gzip
require 'zlib'
# import matplotlib.pyplot as plt
# import numpy as np
require 'numpy'
np = Numpy
# from pyteomics import fasta, parser, mass, achrom, electrochem, auxiliary
require_relative 'rbteomics/fasta'
require_relative 'rbteomics/auxiliary/utils'
require_relative 'rbteomics/parser'
require 'set'

if FileTest.exist?('yeast.fasta.gz')
  # PASS
end

puts 'Cleaving the proteins with trypsin...'
unique_peptides = Set.new
Zlib::GzipReader.open('yeast.fasta.gz') do |gz|
  gz.read.split('>').each_with_index do |x, i|
    f = Bio::FastaFormat.new(x)
    new_peptides = cleave(f.seq, 'trypsin')
    unique_peptides.merge(new_peptides)
  end
end
puts "Done, #{unique_peptides.size} sequences obtained!"
