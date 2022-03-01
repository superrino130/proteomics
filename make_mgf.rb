
class Charge < Array
  def initialize(*args, **kwargs)
    begin
      self << args[0].to_i
      return
    rescue => exception
      if args[0].instance_of?(String)
        begin
          num, sign = /^(\d+)(\+|-)$/.match(args[0])
          self << super(sign + num, args[1..-1], kwargs)
          return
        rescue => exception
          # pass
        end
      end
      raise PyteomicsError.new(*e.args)
    end
  end

  def __str__()
    self.abs.to_s + (self < 0 ? '-' : '+')
  end
end

class Ion < String
  @_pattern = /([abcxyz]\d+(\-H2O|\-NH3)?)([\+|-]\d+)/

  def initialize(...)
    __init__(...)
  end

  def __init__(*args, **kwargs)
    if args || args[0].instance_of?(String)
      begin
        @ion_type, @neutral_loss, charge = args[0].match(@_pattern)
      rescue => exception
        raise PyteomicsError.new("Malformed ion string, must match the regex #{_pattern}")
      end
    end
  end
end

class ChargeList < Array
  def initialize(...)
    __init__(...)
  end

  def __init__(*args, **kwargs)
    if ['', 0, nil, [], {}].include?(args).! || args[0].instance_of?(String)
      delim = /(?:,\s*)|(?:\s*and\s*)/
      args[0].split(delim.match(args[0]).to_s).each do |e|
        e = Charge.new(e)[0].to_i
        self << e
      end
    else
      begin
        super(arges[0].uniq.sort, *args[1..-1], **kwargs)
      rescue => exception
        super(*args, **kwargs)
      end
      self[0] = Charge.new(self[0])[0].to_i
    end
  end

  def __str__()
    if self.size > 1
      self[0..-2].join(', ') + " and #{self[-1]}"
    elsif self
      self[0].to_s
    else
      super.__str__(self)
    end
  end
end

def _parse_charge(s, list_only: false)
  if list_only.!
    begin
      return Charge.new(s)[0].to_i
    rescue => exception
      # pass
    end
  end
  return ChargeList.new(s)
end

def _parse_ion(ion_text)
  begin
    return Ion(ion_text)
  rescue => exception
    warn "Could not parse ion string: #{ion_text} #{exception.args[0]})"
  end
end


module Mgf
  class MGFBase < Array
    def initialize(path)
      @path = path
      self.clear
    end

    def self.parse_precursor_charge(charge_text, list_only: false)
      _parse_charge(charge_text, list_only: list_only)
    end
  
    def self.parse_peak_charge(charge_text, list_only: false)
      _parse_charge(charge_text, list_only: false)
    end
  
    def self.parse_peak_ion(ion_text)
      _parse_ion(ion_text)
    end
  end

  class MGF < MGFBase
    def initialize(path)
      super
    end

    def parser
      h = {}
      d = []
      di = -1
      File.open(@path) do |f|
        f.readlines.each do |line|
          if line.chomp.empty?
            next
          elsif line.include?('BEGIN IONS')
            di += 1
            d[di] = {'intensity array' => [], 'm/z array' => [], 'params' => {}}
          elsif line.include?('END IONS')
            # PASS
          elsif di < 0
            k, v = line.chomp.split('=')
            h[k.downcase] = v if k.nil?.!
          else
            if line.include?('=')
              k, v = line.chomp.split('=')
              d[di]['params'][k.downcase] = v
            else
              ma, ia = line.split.map(&:to_f)
              d[di]['intensity array'] << ia
              d[di]['m/z array'] << ma
            end
          end
        end
      end
      d.each do |elems|
        elems.keys.each do |key|
          if key == 'params'
            elems[key] = (elems[key].to_a + h.to_a).sort_by{ _1[0] }.to_h
          end
        end
        elems['charge array'] = [] if elems.include?('charge array').!
        self << elems
      end
      fix
    end

    def fix
      self.each do |elems|
        elems['params'].keys.each do |key|
          case key
          when 'pepmass'
            pm = elems['params'][key].split.map(&:to_f)
            # if pm.size > 1
              ch = MGFBase.parse_precursor_charge(elems['params']['charge'], list_only: true)
            # else
            #   ch = [MGFBase.parse_precursor_charge(elems['params']['charge'], list_only: false)]
            # end
            if ch.size == pm.size
              elems['params'][key] = pm
              elems['params']['charge'] = ch
            elsif ch.size > pm.size
              elems['params'][key] = pm + [nil] * (ch.size - pm.size)
              elems['params']['charge'] = ch[0, pm.size]
            end
          when 'rtinseconds'
            elems['params'][key] = elems['params'][key].to_f
          end
        end
      end
    end
  end

  class IndexedMGF < MGFBase

  end
end

include Mgf

path = 'tests/test.mgf'
# path = 'example.mgf'
d = Mgf::MGF.new(path).parser
p d

# p d == [
#   {
#     "intensity array"=>[73.0, 44.0, 67.0, 291.0, 54.0, 49.0], 
#     "m/z array"=>[846.6, 846.8, 847.6, 1640.1, 1640.6, 1895.5], 
#     "params"=>{
#       "charge"=>[2], 
#       "com"=>"Based on http://www.matrixscience.com/help/data_file_help.html", 
#       "it_mods"=>"Oxidation (M)", 
#       "itol"=>"1", 
#       "itolu"=>"Da", 
#       "mass"=>"Monoisotopic", 
#       "mods"=>"Carbamidomethyl (C)", 
#       "pepmass"=>[983.6, nil], 
#       "title"=>"Spectrum 1", 
#       "useremail"=>"leu@altered-state.edu", 
#       "username"=>"Lou Scene"
#     }, 
#     "charge array"=>[]
#   }, 
#   {
#     "intensity array"=>[237.0, 128.0, 108.0, 1007.0, 974.0, 79.0], 
#     "m/z array"=>[345.1, 370.2, 460.2, 1673.3, 1674.0, 1675.3], 
#     "params"=>{
#       "charge"=>[2, 3], 
#       "com"=>"Based on http://www.matrixscience.com/help/data_file_help.html", 
#       "it_mods"=>"Oxidation (M)", 
#       "itol"=>"1", 
#       "itolu"=>"Da", 
#       "mass"=>"Monoisotopic", 
#       "mods"=>"Carbamidomethyl (C)", 
#       "pepmass"=>[1084.9, 1234.0], 
#       "rtinseconds"=>25.0, 
#       "scans"=>"3", 
#       "title"=>"Spectrum 2", 
#       "useremail"=>"leu@altered-state.edu", 
#       "username"=>"Lou Scene"
#     }, 
#     "charge array"=>[]
# }]

require 'matplotlib/pyplot'
plt = Matplotlib::Pyplot

spectrum = d[0]

plt.figure()
plt.title('Theoretical and experimental spectra for ')
        # + psm['search_hit'][0]['peptide'])
plt.xlabel('m/z, Th')
plt.ylabel('Intensity, rel. units')

plt.bar(spectrum['m/z array'], spectrum['intensity array'], width: 0.1, linewidth: 2,
  edgecolor: 'black')

# theor_spectrum = list(fragments(psm['search_hit'][0]['peptide'],
#     maxcharge=psm['assumed_charge']))
# plt.bar(theor_spectrum,
#         [spectrum['intensity array'].max()]*len(theor_spectrum),
#         width=0.1, edgecolor='red', alpha=0.7)
plt.show()