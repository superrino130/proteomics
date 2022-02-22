# from pyteomics import mgf, pepxml, mass
require_relative 'rbteomics/mgf'
require_relative 'rbteomics/pepxml'
require_relative 'rbteomics/mass/mass'

# import os
# from urllib.request import urlopen, Request
# import pylab
require 'matplotlib/pyplot'
plt = Matplotlib::Pyplot
require 'set'
require 'open-uri'

['mgf', 'pep.xml'].each do |fname|
  if FileTest.exist?('example.' + fname).!
    if fname == 'mgf'
      url = 'https://pyteomics.readthedocs.io/en/latest/_downloads/5f9796f834db02d33534508ba9266f0d/example.mgf'
    elsif fname == 'pep.xml'
      url = 'https://pyteomics.readthedocs.io/en/latest/_downloads/9ec3a75036fb2adf2acd1f33ab7fbb9a/example.pep.xml'
    end
    URI.open(url) {|f|
      puts 'Downloading ' + 'example.' + fname + '...'
      IO.copy_stream(f, 'example.' + fname)
    }
  end
end

def fragments(peptide, types: ['b', 'y'], maxcharge: 1)
  Fiber.new do
    1.upto(peptide.size - 2) do |i|
      types.each do |ion_type|
        1.upto(maxcharge) do |charge|
          if 'abc'.include?(ion_type[0])
            Fiber.yield fast_mass(peptide[0...i], ion_type: ion_type, charge: charge)
          else
            Fiber.yield fast_mass(peptide[i..], ion_type: ion_type, charge: charge)
          end
        end
      end
    end  
  end
end

spectrum = Mgf::Read.call('example.mgf')
psms = Pepxml::Read.call('example.pep.xml')



Mgf::Read.call('example.mgf') do |spectra|
  # File.open('example.pep.xml') do |psms|
    spectrum = _next(spectra)
    p spectrum

  # end
end

exit
# with mgf.read('example.mgf') as spectra, pepxml.read('example.pep.xml') as psms:
#     spectrum = next(spectra)
#     psm = next(psms)




exit


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
  peptide['parsed_sequence'] = parse(peptide['sequence'], show_unmodified_termini: true)
  peptide['length'] = length(peptide['parsed_sequence'])
end
puts 'Done!'

peptides = peptides.select{ _1['length'] <= 100 }

puts "peptides.size = #{peptides.size}"

puts 'Calculating the mass, charge and m/z...'
peptides.each do |peptide|
  peptide['charge'] = charge(peptide['parsed_sequence'], pH=2.0).round
  peptide['mass'] = calculate_mass(peptide['parsed_sequence'])
  peptide['m/z'] = calculate_mass(peptide['parsed_sequence'], **{'charge' => peptide['charge']})
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
    bins: 2000,
    range: [0, 4000])
plt.xlabel('m/z, Th - Ruby')
plt.ylabel('# of peptides within 2 Th bin')

plt.show()

plt.figure()
plt.hist(peptides.map{ _1['charge'] },
    bins: 20,
    range: [0, 10])
plt.xlabel('charge, e - Ruby')
plt.ylabel('# of peptides')

plt.show()

x = peptides.map{ _1['RT_RP'] }
y = peptides.map{ _1['RT_normal'] }
heatmap, xbins, ybins = np.histogram2d(x, y, bins=100)
heatmap[heatmap == 0] = np.nan
a, b, r, stderr = linear_regression(x, y: y)

plt.figure()
plt.imshow(heatmap)
plt.xlabel('RT on RP, min - Ruby')
plt.ylabel('RT on normal phase, min')
plt.title("All tryptic peptides, RT correlation = #{r}")

plt.show()

x = peptides.map{ _1['m/z'] }
y = peptides.map{ _1['RT_RP'] }
heatmap, xbins, ybins = np.histogram2d(x, y,
    bins: [150, 2000],
    range: [[0, 4000], [0, 150]])
heatmap[heatmap == 0] = np.nan
a, b, r, stderr = linear_regression(x, y: y)

plt.figure()
plt.imshow(heatmap,
    aspect: 'auto',
    origin: 'lower')
plt.xlabel('m/z, Th - Ruby')
plt.ylabel('RT on RP, min')
plt.title("All tryptic peptides, correlation = #{r}")

plt.show()

close_mass_peptides = peptides.select{ 700.0 <= _1['m/z'] && _1['m/z'] <= 701.0 }
x = close_mass_peptides.map{ _1['RT_RP'] }
y = close_mass_peptides.map{ _1['RT_normal'] }
a, b, r, stderr = linear_regression(x, y: y)

plt.figure()
plt.scatter(x, y)
plt.xlabel('RT on RP, min - Ruby')
plt.ylabel('RT on normal phase, min')
plt.title("Tryptic peptides with m/z=700-701 Th\nRT correlation = #{r}")

plt.show()