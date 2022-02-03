# from __future__ import division
# import math
# from .. import parser
require_relative '../parser'
# from ..auxiliary import PyteomicsError, _nist_mass, BasicComposition
require_relative '../auxiliary/constants'
require_relative '../auxiliary/structures'
# from itertools import chain, product, combinations_with_replacement
# from collections import defaultdict
# try:
#     from urllib import urlopen
# except ImportError:
#     from urllib.request import urlopen
# from datetime import datetime
# import re
# import operator
require 'set'
# from lxml import etree
# from ..xml import _local_name
require 'time'
# from ..xml import xpath

NIST_mass = Nist_mass.dup
STD_aa_comp = {}
STD_ion_comp = {}

# @_isotope_string = /^([A-Z][a-z+]*)(?:\[(\d+)\])?$/
# @_atom = /([A-Z][a-z+]*)(?:\[(\d+)\])?([+-]?\d+)?/
# @_formula = /^(#{@_atom})*$/
ISOTOPE_string = /^([A-Z][a-z+]*)(?:\[(\d+)\])?$/
ATOM_string = /([A-Z][a-z+]*)(?:\[(\d+)\])?([+-]?\d+)?/
FORMULA_string = /^(([A-Z][a-z+]*)(?:\[(\d+)\])?([+-]?\d+)?)*$/

def _make_isotope_string(element_name, isotope_num)
  if isotope_num == 0
    element_name
  else
    "#{element_name}[#{isotope_num}]"
  end
end

def _parse_isotope_string(label)
  m = ISOTOPE_string.match(label)
  element_name, num = m[1], m[2]
  isotope_num = ['', 0, nil, false, [], {}].include?(num).! ? num.to_i : 0
  [element_name, isotope_num]
end
  
class Composition < BasicComposition
  attr_reader :_kw_sources

  def _from_parsed_sequence(parsed_sequence, aa_comp)
    self.clear
    comp = Hash.new(0)
    parsed_sequence.each do |aa|
      if aa_comp.include?(aa)
        aa_comp[aa].each do |elem, cnt|
          comp[elem] += cnt
        end
      else
        begin
          mod, aa = _split_label(aa)
          (aa_comp[mod].to_a + aa_comp[aa].to_a).each do | elem, cnt|
            comp[elem] += cnt
          end
        rescue => exception
          raise PyteomicsError.new("No information for #{aa} in `aa_comp`")
        end
      end
    end
    _from_composition(comp)
  end

  def _from_split_sequence(split_sequence, aa_comp)
    clear()
    comp = Hash.new(0)
    split_sequence.each do |group|
      i = 0
      while i < group.size
        jj = 0
        group.size.next.downto(0) do |j|
          jj = j
          begin
            label = group[i...j].join('')
            aa_comp[label].each do |elem, cnt|
              comp[elem] += cnt
            end
          rescue => exception
            next
          else
            i = j
            break
          end
        end
        raise PyteomicsError.new("Invalid group starting from position #{i + 1}: #{group}") if jj == 0
      end
    end
    _from_composition(comp)
  end

  def _from_sequence(sequence, aa_comp)
    parsed_sequence = parse(
        sequence,
        show_unmodified_termini: true,
        'labels' => aa_comp
      )
    _from_parsed_sequence(parsed_sequence, aa_comp)
  end

  def _from_formula(formula, mass_data)
    if FORMULA_string.match(formula).!
      raise PyteomicsError.new('Invalid formula: ' + formula)
    end
    formula.scan(ATOM_string).each do |elem, isotope, number|
      if mass_data.include?(elem).!
        raise PyteomicsError.new('Unknown chemical element: ' + elem)
      end
      # self[_make_isotope_string(elem, int(isotope) if isotope else 0)] += int(number) if number else 1
      self[_make_isotope_string(elem, ['', 0, nil, false, []].include?(isotope) ? 0 : isotope.to_i)] += number.nil? ? 1 : number.to_i
    end
  end

  def _from_composition(comp)
    comp.each do |isotope_string, num_atoms|
      element_name, isotope_num = _parse_isotope_string(isotope_string)
      self[_make_isotope_string(element_name, isotope_num)] = num_atoms
    end
  end
  
  def initialize(...)
    @_kw_sources = Set.new(['formula', 'sequence', 'parsed_sequence', 'split_sequence', 'composition'])

    __init__(...)
  end

  def __init__(*args, **kwargs)
    self.default = 0
    aa_comp = kwargs['aa_comp'] || STD_aa_comp
    mass_data = kwargs['mass_data'] || NIST_mass
    kw_given = @_kw_sources & kwargs.keys
    if kw_given.size > 1
      raise PyteomicsError.new("Only one of #{@_kw_sources.join(', ')} can be specified!\nGiven: #{kw_given.join(', ')}")
    elsif kw_given.empty?.!
      kwa = kw_given.first
      kw_given.delete(kwa)
      method(('_from_' + kwa).to_sym).call(kwargs[kwa], kwa == 'formula' ? mass_data : aa_comp)
    elsif ['', 0, nil, false, [], {}].include?(args).!
      if args[0].is_a?(Hash)
          _from_composition(args[0])
      elsif args[0].is_a?(String)
        begin
          _from_sequence(args[0], aa_comp)
        rescue => exception
          begin
            _from_formula(args[0], mass_data)
          rescue => exception
            raise PyteomicsError.new("Could not create a Composition object from string: '#{args[0]}': not a valid sequence or formula")
          end
        end
      else
        begin
          _from_sequence(tostring(args[0], show_unmodified_termini: true), aa_comp)
        rescue => exception
          raise PyteomicsError.new("Could not create a Composition object from '#{args[0]}'. A Composition object must be specified by sequence, parsed or split sequence, formula or dict.")
        end
      end
    else
      _from_composition(kwargs)
    end

    ion_comp = kwargs['ion_comp'] || STD_ion_comp
    if kwargs.include?('ion_type')
      ion_comp[kwargs['ion_type']].each do |k, v|
        self[k] = v
      end
    end

    charge = self['H+']
    if kwargs.include?('charge')
      if ['', 0, nil, false, [], {}].include?(charge).!
        raise PyteomicsError.new("Charge is specified both by the number of protons and 'charge' in kwargs")
      end
      charge = kwargs['charge']
      self['H+'] = charge
    end
  end

  def mass(**kwargs)
    composition = self.dup
    mass_data = kwargs['mass_data'] || NIST_mass

    mass = 0.0
    average = kwargs['average'] || false
    composition.each do |isotope_string, amount|
      element_name, isotope_num = _parse_isotope_string(isotope_string)
      if ['', 0, nil, false, [], {}].include?(isotope_num) && average
        # if isotope_num.! && average
        mass_data[element_name].each do |isotope, data|
          mass += (amount * data[0] * data[1]) if ['', 0, nil, false, [], {}].include?(isotope).!
        end
      else
        mass += (amount * mass_data[element_name][isotope_num][0])
      end
    end

    charge = kwargs['charge'] || composition['H+']
    if ['', 0, nil, false, [], {}].include?(charge).!
      mass += mass_data['H+'][0][0] * charge if ['', 0, nil, false, [], {}].include?(composition['H+'])
      mass /= charge
    end
    mass
  end
end

STD_aa_comp.merge!({
  'A' => Composition.new(**{'H' => 5, 'C' => 3, 'O' => 1, 'N' => 1}),
  'C' => Composition.new(**{'H' => 5, 'C' => 3, 'S' => 1, 'O' => 1, 'N' => 1}),
  'D' => Composition.new(**{'H' => 5, 'C' => 4, 'O' => 3, 'N' => 1}),
  'E' => Composition.new(**{'H' => 7, 'C' => 5, 'O' => 3, 'N' => 1}),
  'F' => Composition.new(**{'H' => 9, 'C' => 9, 'O' => 1, 'N' => 1}),
  'G' => Composition.new(**{'H' => 3, 'C' => 2, 'O' => 1, 'N' => 1}),
  'H' => Composition.new(**{'H' => 7, 'C' => 6, 'N' => 3, 'O' => 1}),
  'I' => Composition.new(**{'H' => 11, 'C' => 6, 'O' => 1, 'N' => 1}),
  'J' => Composition.new(**{'H' => 11, 'C' => 6, 'O' => 1, 'N' => 1}),
  'K' => Composition.new(**{'H' => 12, 'C' => 6, 'N' => 2, 'O' => 1}),
  'L' => Composition.new(**{'H' => 11, 'C' => 6, 'O' => 1, 'N' => 1}),
  'M' => Composition.new(**{'H' => 9, 'C' => 5, 'S' => 1, 'O' => 1, 'N' => 1}),
  'N' => Composition.new(**{'H' => 6, 'C' => 4, 'O' => 2, 'N' => 2}),
  'P' => Composition.new(**{'H' => 7, 'C' => 5, 'O' => 1, 'N' => 1}),
  'Q' => Composition.new(**{'H' => 8, 'C' => 5, 'O' => 2, 'N' => 2}),
  'R' => Composition.new(**{'H' => 12, 'C' => 6, 'N' => 4, 'O' => 1}),
  'S' => Composition.new(**{'H' => 5, 'C' => 3, 'O' => 2, 'N' => 1}),
  'T' => Composition.new(**{'H' => 7, 'C' => 4, 'O' => 2, 'N' => 1}),
  'V' => Composition.new(**{'H' => 9, 'C' => 5, 'O' => 1, 'N' => 1}),
  'W' => Composition.new(**{'C' => 11, 'H' => 10, 'N' => 2, 'O' => 1}),
  'Y' => Composition.new(**{'H' => 9, 'C' => 9, 'O' => 2, 'N' => 1}),
  'U' => Composition.new(**{'H' => 5, 'C' => 3, 'O' => 1, 'N' => 1, 'Se' => 1}),
  'O' => Composition.new(**{'H' => 19, 'C' => 12, 'O' => 2, 'N' => 3}),
  'H-' => Composition.new(**{'H' => 1}),
  '-OH' => Composition.new(**{'O' => 1, 'H' => 1}),
})
STD_ion_comp.merge!({
  'M' => Composition.new(**{'formula' => ''}),
  'M-H2O' => Composition.new(**{'formula' => 'H-2O-1'}),
  'M-NH3' => Composition.new(**{'formula' => 'N-1H-3'}),
  'a' => Composition.new(**{'formula' => 'H-2O-1' + 'C-1O-1'}),
  'a-H2O' => Composition.new(**{'formula' => 'H-2O-1' + 'C-1O-1' + 'H-2O-1'}),
  'a-NH3' => Composition.new(**{'formula' => 'H-2O-1' + 'C-1O-1' + 'N-1H-3'}),
  'b' => Composition.new(**{'formula' => 'H-2O-1'}),
  'b-H2O' => Composition.new(**{'formula' => 'H-2O-1' + 'H-2O-1'}),
  'b-NH3' => Composition.new(**{'formula' => 'H-2O-1' + 'N-1H-3'}),
  'c' => Composition.new(**{'formula' => 'H-2O-1' + 'NH3'}),
  'c-1' => Composition.new(**{'formula' => 'H-2O-1' + 'NH3' + 'H-1'}),
  'c-dot' => Composition.new(**{'formula' => 'H-2O-1' + 'NH3' + 'H1'}),
  'c+1' => Composition.new(**{'formula' => 'H-2O-1' + 'NH3' + 'H1'}),
  'c+2' => Composition.new(**{'formula' => 'H-2O-1' + 'NH3' + 'H2'}),
  'c-H2O' => Composition.new(**{'formula' => 'H-2O-1' + 'NH3' + 'H-2O-1'}),
  'c-NH3' => Composition.new(**{'formula' => 'H-2O-1'}),
  'x' => Composition.new(**{'formula' => 'H-2O-1' + 'CO2'}),
  'x-H2O' => Composition.new(**{'formula' => 'H-2O-1' + 'CO2' + 'H-2O-1'}),
  'x-NH3' => Composition.new(**{'formula' => 'H-2O-1' + 'CO2' + 'N-1H-3'}),
  'y' => Composition.new(**{'formula' => ''}),
  'y-H2O' => Composition.new(**{'formula' => 'H-2O-1'}),
  'y-NH3' => Composition.new(**{'formula' => 'N-1H-3'}),
  'z' => Composition.new(**{'formula' => 'H-2O-1' + 'ON-1H-1'}),
  'z-dot' => Composition.new(**{'formula' => 'H-2O-1' + 'ON-1'}),
  'z+1' => Composition.new(**{'formula' => 'H-2O-1' + 'ON-1H1'}),
  'z+2' => Composition.new(**{'formula' => 'H-2O-1' + 'ON-1H2'}),
  'z+3' => Composition.new(**{'formula' => 'H-2O-1' + 'ON-1H3'}),
  'z-H2O' => Composition.new(**{'formula' => 'H-2O-1' + 'ON-1H-1' + 'H-2O-1'}),
  'z-NH3' => Composition.new(**{'formula' => 'H-2O-1' + 'ON-1H-1' + 'N-1H-3'}),
})

def calculate_mass(*args, **kwargs)
  composition = kwargs.include?('composition') ? Composition.new(kwargs['composition']) : Composition.new(*args, **kwargs)
  composition.mass(**kwargs)
end

def most_probable_isotopic_composition(*args, **kwargs)
  composition = kwargs.include?('composition') ? kwargs['composition'].to_h : Composition.new(*args, **kwargs)

  composition.each do |isotope_string, _|
    element_name, isotope_num = _parse_isotope_string(isotope_string)
    if isotope_num != 0
      i = composition[isotope_string] || 0
      composition.delete(isotope_string)
      composition[element_name] += i
    end
  end

  mass_data = kwargs['mass_data'] || NIST_mass
  elements_with_isotopes = kwargs['elements_with_isotopes']
  isotopic_composition = Composition.new

  composition.each do |element_name, _|
    if elements_with_isotopes.! || elements_with_isotopes.include?(element_name)
      # first_iso, second_iso = sorted([(i[0], i[1][1]) for i in mass_data[element_name].items() if i[0]], key=lambda x: -x[1])[:2]
      first_iso, second_iso = mass_data[element_name].select{ |k, v| ['', 0, nil, false].include?(k).! }.map{ |k, v| [k, v[1]] }.sort_by{ |x| -x[1] }[0...2]

      first_iso_str = _make_isotope_string(element_name, first_iso[0])
      isotopic_composition[first_iso_str] = composition[element_name].ceil.to_i * first_iso[1]

      second_iso_str = _make_isotope_string(element_name, second_iso[0])
      isotopic_composition[second_iso_str] = composition[element_name] - isotopic_composition[first_iso_str]
    else
      isotopic_composition[element_name] = composition[element_name]
    end
  end
  [isotopic_composition, isotopic_composition_abundance('composition' => isotopic_composition, 'mass_data' => mass_data)]
end

def isotopic_composition_abundance(*args, **kwargs)
  composition = kwargs.include?('composition') ? Composition.new(kwargs['composition']) : Composition.new(*args, **kwargs)

  isotopic_composition = Hash.new{ |h, k| h[k] = {} }

  composition.each do |element, _|
    element_name, isotope_num = _parse_isotope_string(element)

    if isotopic_composition.include?('element_name') && (isotope_num == 0 || isotopic_composition[element_name].include?(0))
      raise PyteomicsError.new("Please specify the isotopic states of all atoms of #{element_name} or do not specify them at all.")
    else
      isotopic_composition[element_name][isotope_num] = composition[element]
    end
  end
  mass_data = kwargs['mass_data'] || NIST_mass
  num1, num2, denom = 1, 1, 1
  isotopic_composition.each do |element_name, isotope_dict|
    num1 *= (1..isotope_dict.values.sum).inject(1, :*)
    isotope_dict.each do |isotope_num, isotope_content|
      denom *= (1..isotope_content).inject(1, :*)
      num2 *= mass_data[element_name][isotope_num][1] ** isotope_content if isotope_num != 0
    end
  end
  num2 * (num1 / denom)
end

# Fiber
def isotopologues(*args, **kwargs)
    iso_threshold = kwargs.delete('isotope_threshold') || 5e-4
    overall_threshold = kwargs.delete('overall_threshold') || 0.0
    mass_data = kwargs['mass_data'] || NIST_mass
    elements_with_isotopes = kwargs['elements_with_isotopes']
    report_abundance = kwargs['report_abundance'] || false
    composition = kwargs.include?('composition') ? Composition.new(kwargs['composition']) : Composition.new(*args, **kwargs)
    other_kw = kwargs.dup
    Composition.new._kw_sources.each do |k|
      other_kw.delete(k)
    end
    dict_elem_isotopes = {}
    composition.each do |element, _|
      if elements_with_isotopes.nil? || elements_with_isotopes.include?(element)
        element_name, isotope_num = _parse_isotope_string(element)
        isotopes = mass_data[element_name].select{ |k, v| k != 0 && v[1] >= iso_threshold }
        list_isotopes = isotopes.map{ |k, v| _make_isotope_string(element_name, k) }
        dict_elem_isotopes[element] = list_isotopes
      else
        dict_elem_isotopes[element] = [element]
      end
    end
    all_isotoplogues = []
    dict_elem_isotopes.each do |element, list_isotopes|
      n = composition[element]
      list_comb_element_n = []
      list_isotopes.repeated_combination(n).each do |elementXn|
        list_comb_element_n. << elementXn
      end
      all_isotoplogues << list_comb_element_n
    end
    all_isotoplogues[0].product(*all_isotoplogues[1..-1]).each do |isotopologue|
      ic = Composition.new(isotopologue.flatten.join(''), **other_kw)
      if report_abundance || overall_threshold > 0.0
        abundance = isotopic_composition_abundance('composition' => ic, **other_kw)
        if abundance > overall_threshold
          if report_abundance
            yield [ic, abundance]
          else
            yield ic
          end
        end
      else
        yield ic
      end
    end
end

STD_aa_mass = {
  'G' => 57.02146,
  'A' => 71.03711,
  'S' => 87.03203,
  'P' => 97.05276,
  'V' => 99.06841,
  'T' => 101.04768,
  'C' => 103.00919,
  'L' => 113.08406,
  'I' => 113.08406,
  'J' => 113.08406,
  'N' => 114.04293,
  'D' => 115.02694,
  'Q' => 128.05858,
  'K' => 128.09496,
  'E' => 129.04259,
  'M' => 131.04049,
  'H' => 137.05891,
  'F' => 147.06841,
  'U' => 150.95364,
  'R' => 156.10111,
  'Y' => 163.06333,
  'W' => 186.07931,
  'O' => 237.14773,
}

def fast_mass(sequence, ion_type: nil, charge: nil, **kwargs)
  aa_mass = kwargs['aa_mass'] || STD_aa_mass
  begin
    mass = sequence.split('').map{ aa_mass[_1] }.sum
  rescue => exception
    raise PyteomicsError.new('No mass data for residue: ' + e.args[0].to_s)
  end

  mass_data = kwargs['mass_data'] || NIST_mass
  mass += mass_data['H'][0][0] * 2 + mass_data['O'][0][0]

  if ion_type
    begin
      icomp = (kwargs['ion_comp'] || STD_ion_comp)[ion_type]
    rescue => exception
      raise PyteomicsError.new("Unknown ion type: #{ion_type}")
    end

    mass += icomp.map{ |element, num| mass_data[element][0][0] * num }.sum
  end
  mass = (mass + mass_data['H+'][0][0] * charge) / charge if charge.nil?.! && charge != 0
  mass
end

def fast_mass2(sequence, ion_type: nil, charge: nil, **kwargs)
  aa_mass = kwargs['aa_mass'] || STD_aa_mass
  mass_data = kwargs['mass_data'] || NIST_mass
  aa_mass['H-'] = mass_data['H'][0][0] if aa_mass.include?('H-').!
  aa_mass['-OH'] = mass_data['H'][0][0] + mass_data['O'][0][0] if aa_mass.include?('-OH').!
  begin
    comp = amino_acid_composition(sequence, show_unmodified_termini: true, allow_unknown_modifications: true, 'labels' => aa_mass)
  rescue => exception
    raise PyteomicsError.new("Mass not specified for label(s): #{(Set.new(parse(sequence)) - aa_mass).to_a.join(', ')}")
  end

  begin
    mass = 0
    comp.each do |aa, num|
      if aa_mass.include?(aa)
        mass += aa_mass[aa] * num
      else
        mod, x = _split_label(aa)
        mass += (aa_mass[mod] + aa_mass[x]) * num
      end
    end
  rescue => exception
    raise PyteomicsError.new("Unspecified mass for modification: '#{e.args[0]}'")
  end

  if ion_type
    begin
      icomp = (kwargs['ion_comp'] || STD_ion_comp)[ion_type]
    rescue => exception
      raise PyteomicsError.new("Unknown ion type: #{ion_type}")      
    end

    mass += icomp.map{ |element, num| mass_data[element][0][0] * num }.sum
  end

  mass = (mass + mass_data['H+'][0][0] * charge) / charge if charge.nil?.! && charge != 0
  mass
end

class Unimod
  def initialize(...)
    __init__(...)
  end

  def __init__(source: 'http://www.unimod.org/xml/unimod.xml')
    def process_mod(mod)
      d = mod.attrib
      new_d = {}
      ['date_time_modified', 'date_time_posted'].each do |key|
        new_d[key] = Time.parse(d.delete(key)).strftime("%Y-%m-%d %H:%M:%S")
      end
      comp = Composition.new()
      self._xpath('delta', mod).each do |delta|
        ['avge_mass', 'mono_mass'].each do |key|
          new_d[key] = delta.attrib.delete(key).to_f
        end
        self._xpath('element', delta).each do |elem|
          e_d = elem.attrib
          amount = e_d.delete('number').to_i
          label = e_d.delete('symbol')
          isotope, symbol = /^(\d*)(\D+)$/.match(label)
          if isotope.nil?
            isotope = 0
          else
            isotope = isotope.to_i
          end
          comp += Composition.new(**{'formula' => _make_isotope_string(symbol, isotope), 'mass_data' => self._massdata}) * amount
        end
      end
      new_d['composition'] = comp
      new_d['record_id'] = d.delete('record_id').to_i
      new_d['approved'] = d.delete('approved') == '1'
      new_d.merge!(d)
      spec = []
      self._xpath('specificity', mod).each do |sp|
        sp_d = sp.attrib
        sp_new_d = {}
        sp_new_d['hidden'] = sp_d.delete('hidden') == '1'
        sp_new_d['spec_group'] = sp_d.delete('spec_group').to_i
        sp_new_d.merge!(sp_d)
        notes = []
        self._xpath('*', sp).each do |note|
          notes << note.text.strip if note.text || note.text.strip
        end
        sp_new_d['note'] = notes.join("\n") if notes
        spec << sp_new_d
      end
      new_d['specificity'] = spec

      alt_names = []
      self._xpath('alt_name', mod).each do |alt_name|
        alt_names << alt_name.text
      end
      new_d['alt_names'] = alt_names if alt_names

      refs = []
      self._xpath('xref', mod).each do |ref|
        ref_d = {}
        ref.iterchildren().each do |sub|
          ref_d[_local_name(sub)] = sub.text
        end
        ['text', 'source', 'url'].each do |key|
          ref_d[key] = nil if ref_d.include?(key).!
        end
        refs << ref_d
      end
      new_d['refs'] = refs
      return new_d
    end

    if source.instance_of?(String)
      self._tree = etree.parse(urlopen(source))
    else
      self._tree = etree.parse(source)
    end
    self._massdata = self._mass_data()
    self._mods = []
    self._id = {}
    enumerate(self._xpath('/unimod/modifications/mod')).each do |i, mod|
      mod_dict = process_mod(mod)
      self._mods << mod_dict
      self._id[mod_dict['record_id']] = i
    end
  end

  def _xpath(path, element: nil)
    return xpath(self._tree, path, 'umod') if element.nil?
    xpath(element, path, 'umod')
  end

  def _mass_data()
    massdata = Hash.new({})
    elements = self._xpath('/unimod/elements/elem').map{ _1.attrib }
    avg = {}
    elements.each do |elem|
      i, label = /^(\d*)(\D+)$/.match(elem['title'])
      if i.nil?
        iso = 0
      else
        iso = i.to_i
      end
      massdata[label][iso] = [elem['mono_mass'].to_f, iso == 0 ? 1.0 : 0.0]
      avg[label] = elem['avge_mass'].to_f if iso.nil?.!
    end
    massdata.each do |elem, isotopes|
      isotopes[isotopes[0][0].round.to_i] = isotopes[0]
      if isotopes.size == 3
        m1, m2 = isotopes.sort[1..-1].map{ _1[1][0] }
        m_avg = avg[elem]
        a = (m2 - m_avg) / (m2 - m1)
        b = (m_avg - m1) / (m2 - m1)
        isotopes.sort[1..-1].zip([a, b]).each do |state, abundance|
          isotopes[state] = [isotopes[state][0], abundance]
        end
      end
    end
    massdata
  end

  # @property
  def mods()
    self._mods
  end

  # @property
  def mass_data()
    self._massdata
  end

  def by_title(title, strict: true)
    f = {true: operator.eq, false: operator.contains}
    func = f[strict]
    result = self._mods.select{ func(_1['title'], title) }
    return result[0] if result.size == 1
    result
  end

  def by_name(name, strict: true)
    f = {true: operator.eq, talse: operator.contains}
    func = f[strict]
    result = self._mods.select{ func(_1['full_name'], name) }
    return result[0] if result.size == 1
    result
  end

  def by_id(i)
    i = i.to_i if i.instance_of?(String)
    self._mods[self._id[i]]
  end

  # __getitem__ = by_id
end