# import os
# from urllib.request import urlretrieve
# import gzip
require 'zlib'
# import matplotlib.pyplot as plt
# import numpy as np
require 'matplotlib/pyplot'
plt = Matplotlib::Pyplot
require 'numpy'
np = Numpy
# from pyteomics import fasta, parser, mass, achrom, electrochem, auxiliary
require_relative 'rbteomics/fasta'
require_relative 'rbteomics/auxiliary/utils'
require_relative 'rbteomics/mass/mass'
require_relative 'rbteomics/parser'
require_relative 'rbteomics/electrochem'
require_relative 'rbteomics/achrom'
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

peptides = unique_peptides.map{ {'sequence' => _1} }

puts 'Parsing peptide sequences...'
peptides.each do |peptide|
  peptide['parsed_sequence'] = parse(peptide['sequence'], show_unmodified_termini=true)
  peptide['length'] = length(peptide['parsed_sequence'])
end
puts 'Done!'

peptides = peptides.select{ _1['length'] <= 100 }

puts "peptides.size = #{peptides.size}"

puts 'Calculating the mass, charge and m/z...'
peptides.each do |peptide|
  peptide['charge'] = charge(peptide['parsed_sequence'], pH=2.0).round
  peptide['mass'] = calculate_mass(peptide['parsed_sequence'])
  peptide['m/z'] = calculate_mass(peptide['parsed_sequence'], {'charge' => peptide['charge']})
end
puts 'Done!'

puts 'Calculating the retention time...'
peptides.each do |peptide|
  peptide['RT_RP'] = calculate_RT(
    peptide['parsed_sequence'],
    RCs_zubarev)
peptide['RT_normal'] = calculate_RT(
    peptide['parsed_sequence'],
    RCs_yoshida_lc)
end

plt.figure()
plt.hist(peptides.map{ _1['m/z'] },
    bins = 2000,
    range=[0,4000])
plt.xlabel('m/z, Th')
plt.ylabel('# of peptides within 2 Th bin')

plt.show()
