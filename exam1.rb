# import os
# from urllib.request import urlretrieve
# import gzip
require 'zlib'
# import matplotlib.pyplot as plt
# import numpy as np
require 'numpy'
np = Numpy
# from pyteomics import fasta, parser, mass, achrom, electrochem, auxiliary
require './rbteomics/fasta'
# require './rbteomics/auxiliary/'
require 'set'

if FileTest.exist?('yeast.fasta.gz')
  puts
  # PASS
end

puts 'Cleaving the proteins with trypsin...'
unique_peptides = Set.new
Zlib::GzipReader.open('yeast.fasta.gz') do |gz|
  # fasta.FASTA(gzfile).each do |description, sequence|

  # end
  
  # gz.read.split(">").each_with_index do |x, i|
  #   next if i == 0
  #   s = ""
  #   x.split("\n")[1..-1].each do |y|
  #     s << y unless y.nil?
  #   end
  #   unique_peptides.merge(s.split(/([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))/).to_set)
  # end
end
