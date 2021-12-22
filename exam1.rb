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
require 'open-uri'

if FileTest.exist?('yeast.fasta.gz').!
  puts 'Downloading the FASTA file for Saccharomyces cerevisiae...'
  URI.open("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000002311/UP000002311_559292.fasta.gz") {|f|
    IO.copy_stream(f, 'yeast.fasta.gz')
  }
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

peptides = unique_peptides.map{ {'sequence': _1} }

puts 'Parsing peptide sequences...'

