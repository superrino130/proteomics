# from __future__ import division
# from . import parser
# from .auxiliary import PyteomicsError
require_relative 'auxiliary/structures'
# from collections import Iterable, Counter

def charge(sequence, pH, **kwargs)
  peptide_dict, pK = _prepare_charge_dict(sequence, **kwargs)

  pH_list = pH.instance_of?(Array) ? pH : [pH]

  charge_list = _charge_for_dict(peptide_dict, pH_list, pK)
  pH.instance_of?(Array).! ? charge_list[0] : charge_list
end

def _prepare_charge_dict(sequence, **kwargs)
  nterm = cterm = n_aa = c_aa = nil
  pK = (kwargs['pK'] || @pK_lehninger).dup
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
    if (n_aa is None or c_aa is None) and (pK_nterm or pK_cterm)
      raise PyteomicsError.new('Two terminal residues must be present in peptide (designated as "ntermX" and "ctermX", where "X" is the one-letter residue label). Use ``term_aa=True`` when calling `parser.amino_acid_composition`.')
    end
  elsif [String, Array].include?(sequence.class)
    if sequence.instance_of?(String)
      if sequence.match?(/\A\p{Lu}+\z/)
        parsed_sequence = [STD_nterm] + sequence.to_a + [STD_cterm]
      else
        parsed_sequence = parse(sequence, show_unmodified_termini=true)
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
          1. / (1. + 10 ** (ionizable_group[1] * (pH_value - ionizable_group[0]))))
      end
    end
    charge_list << charge
  end

  charge_list
end




@pK_lehninger = {
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
