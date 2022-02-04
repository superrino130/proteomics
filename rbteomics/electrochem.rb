# from __future__ import division
# from . import parser
require_relative 'parser'
# from .auxiliary import PyteomicsError
require_relative 'auxiliary/structures'
# from collections import Iterable, Counter

PK_lehninger = {
  'E' =>   [[4.25,  -1]],
  'R' =>   [[12.48,  1]],
  'Y' =>   [[10.07, -1]],
  'D' =>   [[3.65,  -1]],
  'H' =>   [[6.00,  +1]],
  'K' =>   [[10.53, +1]],
  'C' =>   [[8.18,  -1]],
  'H-' =>  [[9.69,  +1]],
  '-OH' => [[2.34,  -1]],
}

PK_sillero = {
  'E' =>   [[4.5,  -1]],
  'R' =>   [[12.0, +1]],
  'Y' =>   [[10.0, -1]],
  'D' =>   [[4.0,  -1]],
  'H' =>   [[6.4,  +1]],
  'K' =>   [[10.4, +1]],
  'C' =>   [[9.0,  -1]],
  'H-' =>  [[8.2,  +1]],
  '-OH' => [[3.2,  -1]],
}

PK_dawson = {
  'E' =>   [[4.3,  -1]],
  'R' =>   [[12.0, +1]],
  'Y' =>   [[10.1, -1]],
  'D' =>   [[3.9,  -1]],
  'H' =>   [[6.0,  +1]],
  'K' =>   [[10.5, +1]],
  'C' =>   [[8.3,  -1]],
  'H-' =>  [[8.2,  +1]],
  '-OH' => [[3.2,  -1]],
}

PK_rodwell = {
  'E' =>   [[4.25, -1]],
  'R' =>   [[11.5, +1]],
  'Y' =>   [[10.7, -1]],
  'D' =>   [[3.86, -1]],
  'H' =>   [[6.0,  +1]],
  'K' =>   [[11.5, +1]],
  'C' =>   [[8.33, -1]],
  'H-' =>  [[8.0,  +1]],
  '-OH' => [[3.1,  -1]],
}

PK_bjellqvist = {
  'E' =>   [[4.45, -1]],
  'R' =>   [[12.0, +1]],
  'Y' =>   [[10.0, -1]],
  'D' =>   [[4.05, -1]],
  'H' =>   [[5.98, +1]],
  'K' =>   [[10.0, +1]],
  'C' =>   [[9.0,  -1]],
  'H-' =>  [[7.5,  +1]],
  '-OH' => [[3.55, -1]],
}

PK_nterm_bjellqvist = {
  'H-' => {
      'A' => [[7.59, +1]],
      'M' => [[7.0,  +1]],
      'S' => [[6.93, +1]],
      'P' => [[8.36, +1]],
      'T' => [[6.82, +1]],
      'V' => [[7.44, +1]],
      'E' => [[7.7,  +1]]
  }
}

PK_cterm_bjellqvist = {
  '-OH' => {
      'D' => [[4.55, -1]],
      'E' => [[4.75, -1]]
  }
}

Hydropathicity_KD = {
  "A" => 1.800,
  "R" => -4.500,
  "N" => -3.500,
  "D" => -3.500,
  "C" => 2.500,
  "Q" => -3.500,
  "E" => -3.500,
  "G" => -0.400,
  "H" => -3.200,
  "I" => 4.500,
  "L" => 3.800,
  "K" => -3.900,
  "M" => 1.900,
  "F" => 2.800,
  "P" => -1.600,
  "S" => -0.800,
  "T" => -0.700,
  "W" => -0.900,
  "Y" => -1.300,
  "V" => 4.200,
}

def charge(sequence, pH, **kwargs)
  peptide_dict, pK = _prepare_charge_dict(sequence, **kwargs)

  pH_list = pH.instance_of?(Array) ? pH : [pH]

  charge_list = _charge_for_dict(peptide_dict, pH_list, pK)
  pH.instance_of?(Array).! ? charge_list[0] : charge_list
end

def _prepare_charge_dict(sequence, **kwargs)
  nterm = cterm = n_aa = c_aa = nil
  pK = kwargs['pK'] || PK_lehninger.dup
  pK_nterm = kwargs['pK_nterm'] || {}
  pK_cterm = kwargs['pK_cterm'] || {}

  if sequence.instance_of?(Hash)
    peptide_dict = sequence.dup
    sequence.each do |k, v|
      if k[-1] == '-'
        if v > 1 || nterm
          raise PyteomicsError.new("More that one N-terminal group in #{sequence}")
        end
        nterm = k
      end
      if k[0] == '-'
        if v > 1 || cterm
          raise PyteomicsError.new("More that one C-terminal group in #{sequence}")
        end
        cterm = k
      end
      if k[0...5] == 'nterm'
        if v > 1 || n_aa
          raise PyteomicsError.new("More that one N-terminal residue in #{sequence}")
        end
        n_aa = k[5..-1]
        peptide_dict[n_aa] = (peptide_dict[n_aa] || 0) + 1
      end
      if k[0...5] == 'cterm'
        if v > 1 || c_aa
          raise PyteomicsError.new("More that one C-terminal residue in #{sequence}")
        end
        c_aa = k[5..-1]
        peptide_dict[c_aa] = (peptide_dict[c_aa] || 0) + 1
      end
    end
    if nterm.nil? || cterm.nil?
      raise PyteomicsError.new('Peptide must have two explicit terminal groups')
    end
    if (n_aa.nil? || c_aa.nil?) && (pK_nterm.empty?.! || pK_cterm.empty?.!)
      raise PyteomicsError.new("Two terminal residues must be present in peptide (designated as 'ntermX' and 'ctermX', where 'X' is the one-letter residue label). Use 'term_aa=true' when calling 'parser.amino_acid_composition'.")
    end
  elsif [String, Array].include?(sequence.class)
    if sequence.instance_of?(String)
      if sequence.match?(/\A\p{Lu}+\z/)
        parsed_sequence = [STD_nterm] + sequence.split('') + [STD_cterm]
      else
        parsed_sequence = parse(sequence, show_unmodified_termini: true)
      end
    elsif sequence.instance_of?(Array)
      if sequence[0][-1] != '-' || sequence[-1][0] != '-'
        raise PyteomicsError.new('Parsed sequences must contain terminal groups at 0-th and last positions.')
      end
      parsed_sequence = sequence
    end
    n_aa = parsed_sequence[1]
    c_aa = parsed_sequence[-2]
    nterm = parsed_sequence[0]
    cterm = parsed_sequence[-1]
    peptide_dict = parsed_sequence.tally
  else
    raise PyteomicsError.new("Unsupported type of sequence: #{sequence}")
  end

  if pK_nterm.include?(nterm)
    if pK_nterm[nterm].include?(n_aa)
      pK[nterm] = pK_nterm[nterm][n_aa]
    end
  end
  if pK_cterm.include?(cterm)
    if pK_cterm[cterm].include?(c_aa)
      pK[cterm] = pK_cterm[cterm][c_aa]
    end
  end
  [peptide_dict, pK]
end

def _charge_for_dict(peptide_dict, pH_list, pK)
  charge_list = []
  pH_list.each do |pH_value|
    charge = 0
    peptide_dict.each_key do |aa|
      (pK[aa] || []).each do |ionizable_group|
        charge += peptide_dict[aa] * ionizable_group[1] * (
          1.0 / (1.0 + 10 ** (ionizable_group[1] * (pH_value - ionizable_group[0]))))
      end
    end
    charge_list << charge
  end

  charge_list
end

def pI(sequence, pI_range: [0.0, 14.0], precision_pI: 0.01, **kwargs)
  pK = kwargs['pK'] || PK_lehninger.dup
  pK_nterm = {}
  pK_cterm = {}
  if sequence.instance_of?(String) || sequence.instance_of?(Array)
    pK_nterm = kwargs['pK_nterm'] || {}
    pK_cterm = kwargs['pK_cterm'] || {}
  elsif sequence.instance_of?(Hash) && (kwargs.include?('pK_nterm') || kwargs.include?('pK_cterm'))
    raise PyteomicsError.new("Can not use terminal features for #{sequence}")
  end
  peptide_dict, pK = _prepare_charge_dict(sequence, **{'pK' => pK, 'pK_cterm' => pK_cterm, 'pK_nterm' => pK_nterm})
  left_x, right_x = pI_range
  left_y = _charge_for_dict(peptide_dict, [left_x], pK)[0]
  right_y = _charge_for_dict(peptide_dict, [right_x], pK)[0]
  while (right_x - left_x) > precision_pI
    if left_y * right_y > 0
      return left_y.abs < right_y.abs ? left_x : right_x
    end
    middle_x = (left_x + right_x) / 2.0
    middle_y = _charge_for_dict(peptide_dict, [middle_x], pK)[0]
    if middle_y * left_y < 0
      right_x = middle_x
      right_y = middle_y
    else
      left_x = middle_x
      left_y = middle_y
    end
  end
  (left_x + right_x) / 2.0
end

def gravy(sequence, hydropathicity: Hydropathicity_KD)
  begin
    return sequence.map{hydropathicity[_1]}.sum / sequence.size
  rescue => exception
    raise PyteomicsError.new("Hydropathicity for amino acid {} not provided.".format(e.args[0]))    
  end
end
