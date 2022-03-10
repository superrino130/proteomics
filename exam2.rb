require_relative 'rbteomics/mass/mass'
require_relative 'rbteomics/xml'
require_relative 'rbteomics/_schema_defaults'

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

module Pepxml
  class PepXML < Array
    def initialize(source, ...)
      @source = source

      @file_format ||= 'pepXML'
      @_root_element ||= 'msms_pipeline_analysis'
      @_default_schema ||= Pepxml_schema_defaults
      @_default_version ||= '1.15'
      @_default_iter_tag ||= 'spectrum_query'
      @_indexed_tags ||= Set.new(['spectrum_query'])
      @_indexed_tag_keys ||= {'spectrum_query' => 'spectrum'}
      @_default_id_attr ||= 'spectrum'
      @_structures_to_flatten ||= Set.new(['search_score_summary', 'modification_info'])
      @_convert_items ||= {
        'float' => Set.new(['calc_neutral_pep_mass', 'massdiff', 'probability', 'variable', 'static']),
        'int' => Set.new(['start_scan', 'end_scan', 'index', 'num_matched_peptides']),
        'bool' => Set.new(['is_rejected']),
        'floatarray' => Set.new(['all_ntt_prob'])
      }
  
      self.clear
      read
    end

    def convert_items(item, value)
      case item
      when 'retention_time_sec', 'precursor_neutral_mass', 'Morpheus Score', 'PSM q-value', 'calc_neutral_pep_mass', 'massdiff'
        value.to_f
      when 'assumed_charge', 'start_scan', 'end_scan', 'index', 'hit_rank', 'num_tol_term', 'num_tot_proteins', 'num_matched_ions', 'tot_num_ions'
        value.to_i
      when 'is_rejected'
        value == '0' ? false : true
      else
        value
      end
    end

    def read
      ss_flg = false
      c1_cnt = -1
      c2_cnt = -1
      c3_cnt = -1

      IO.foreach(@source) do |line|
        if line.include?('<' + @_default_iter_tag)
          c1_cnt += 1
          c2_cnt = -1
          c3_cnt = -1
          if c1_cnt > 0
            return
          end
          self << {}
          line.split.each do |param|
            if param.include?('=')
              key = param.split('=')[0]
              value = param.match(Regexp.new(key + '=\"(.*?)\"'))[1]
              self[c1_cnt][key] = convert_items(key, value)
            end
          end
        end
        if c1_cnt >= 0
          if line.include?('<search_hit')
            self[c1_cnt]['search_hit'] ||= []
            c2_cnt += 1
            c3_cnt = -1
            self[c1_cnt]['search_hit'] << {}
            key = ''
            line.split[1..].each do |param|
              if param.include?('="')
                key = param.split('=')[0]
                if param.match?(Regexp.new(key + '=\"(.*?)\"'))
                  value = param.match(Regexp.new(key + '=\"(.*?)\"'))[1]
                else
                  value = param.sub(key + '="', '')
                end
                if key == 'hit_rank'
                  c3_cnt += 1
                  self[c1_cnt]['search_hit'][c2_cnt]['proteins'] ||= []
                  self[c1_cnt]['search_hit'][c2_cnt]['proteins'] << {}
                end
                if ['protein', 'protein_descr', 'peptide_next_aa', 'peptide_prev_aa', 'num_tol_term'].include?(key)
                  self[c1_cnt]['search_hit'][c2_cnt]['proteins'][c3_cnt][key] = convert_items(key, value)
                else
                  self[c1_cnt]['search_hit'][c2_cnt][key] = convert_items(key, value)
                end
              else
                self[c1_cnt]['search_hit'][c2_cnt]['proteins'][c3_cnt][key] += ' ' + param
                self[c1_cnt]['search_hit'][c2_cnt]['proteins'][c3_cnt][key].chop! if self[c1_cnt]['search_hit'][c2_cnt]['proteins'][c3_cnt][key][-1] == '>'
                self[c1_cnt]['search_hit'][c2_cnt]['proteins'][c3_cnt][key].chop! if self[c1_cnt]['search_hit'][c2_cnt]['proteins'][c3_cnt][key][-1] == '"'
              end
            end
          elsif line.include?('<search_score')
            self[c1_cnt]['search_hit'][c2_cnt]['search_score'] ||= {}
            key = line.match(Regexp.new('name=\"(.*?)\"'))[1]
            value = line.match(Regexp.new('value=\"(.*?)\"'))[1]
            self[c1_cnt]['search_hit'][c2_cnt]['search_score'][key] = convert_items(key, value)
          end
        end
      end
    end
  end
end

def fragments(peptide, types: ['b', 'y'], maxcharge: 1)
  arr = []
    1.upto(peptide.size - 2) do |i|
      types.each do |ion_type|
        1.upto(maxcharge) do |charge|
          if 'abc'.include?(ion_type[0])
            arr << fast_mass(peptide[0...i], ion_type: ion_type, charge: charge)
          else
            arr << fast_mass(peptide[i..], ion_type: ion_type, charge: charge)
          end
        end
      end
    end
  arr
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
      IO.foreach(@path, chomp: true) do |line|
        if line.empty?
          next
        elsif line.include?('BEGIN IONS')
          di += 1
          d[di] = {'intensity array' => [], 'm/z array' => [], 'params' => {}}
        elsif line.include?('END IONS')
          # PASS
        elsif di < 0
          k, v = line.split('=')
          h[k.downcase] = v if k.nil?.!
        else
          if line.include?('=')
            k, v = line.split('=')
            d[di]['params'][k.downcase] = v
          else
            ma, ia = line.split.map(&:to_f)
            d[di]['intensity array'] << ia
            d[di]['m/z array'] << ma
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
end

include Mgf
include Pepxml

path = 'example.mgf'
spectrum = Mgf::MGF.new(path).parser[0]

path = 'example.pep.xml'
psm = Pepxml::PepXML.new(path)[0]

plt.figure()
plt.title('Theoretical and experimental spectra for ' + psm['search_hit'][0]['peptide'])
plt.xlabel('m/z, Th by Ruby')
plt.ylabel('Intensity, rel. units')

plt.bar(spectrum['m/z array'], spectrum['intensity array'], width: 0.1, linewidth: 2, edgecolor: 'black')

theor_spectrum = fragments(psm['search_hit'][0]['peptide'], types: ['b', 'y'], maxcharge: psm['assumed_charge']).select{ _1 < 700 }
plt.bar(theor_spectrum, [spectrum['intensity array'].max] * theor_spectrum.size, width: 0.1, edgecolor: 'red', alpha: 0.7)
plt.show()