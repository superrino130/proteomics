# import numpy as np
require 'numpy'
# np = Numpy
# from .auxiliary import linear_regression, PyteomicsError
require_relative 'auxiliary/structures'
require_relative 'auxiliary/math'
# from . import parser
require_relative 'parser'

def get_RCs(sequences, rts, lcp: -0.21, term_aa: false, **kwargs)
  labels = kwargs['labels']

  peptide_dicts = sequences.map{ peptide.instance_of?(Hash).! ? amino_acid_composition(_1, false, term_aa, allow_unknown_modifications=true, labels=labels) : _1 }
  
  compositon_array = []
  peptide_dicts.each_key do |pdict|
    loglen = Numpy.log(length(pdict))
  end
end

def get_RCs_vary_lcp(sequences, rts, term_aa: false, lcp_range: [-1.0, 1.0], **kwargs)
  labels = kwargs['labels']

  best_r = -1.1
  best_RC_dict = {}
  lcp_accuracy = kwargs['lcp_accuracy'] || 0.1

  min_lcp = lcp_range[0]
  max_lcp = lcp_range[1]
  step = (max_lcp - min_lcp) / 10.0

  peptide_dicts = sequences.map{ _1.instance_of?(Hash).! ? amino_acid_composition(_1, false, term_aa, allow_unknown_modifications=true, labels=labels) : _1 }
  while step > lcp_accuracy
    lcp_grid = Numpy.arange(min_lcp, max_lcp, (max_lcp - min_lcp) / 10.0)
    lcp_grid.each do |lcp|
      rc_dict = get_RCs(peptide_dicts, rts, lcp, term_aa, labels=labels)
      regression_coeffs = peptide_dicts.keys.map{ linear_regression(rts, [calculate_RT(_1, RC_dict)]) }
      if regression_coeffs[2] > best_r
        best_r = regression_coeffs[2]
        best_RC_dict = rc_dict.to_h
      end
    end
    min_lcp = best_rc_dict['lcp'] - step
    max_lcp = best_rc_dict['lcp'] + step
    step = (max_lcp - min_lcp) / 10.0
  end
  best_rc_dict
end

def calculate_RT(peptide, rc_dict, raise_no_mod: true)
  amino_acids = rc_dict['aa'].keys.select{ (_1[0...5] == 'nterm' || _1[0...5] == 'cterm').! }

  term_aa = false
  rc_dict['aa'].each_key do |aa|
    if aa[0...5] == 'nterm' || aa[0...5] == 'cterm'
      term_aa = true
      break
    end
  end

  if peptide.instance_of?(Hash)
    peptide_dict = peptide.dup
  else
    peptide_dict = amino_acid_composition(peptide, show_unmodified_termini: false, term_aa: term_aa, allow_unknown_modifications: true, 'labels' => amino_acids)
  end
  rt = 0.0
  peptide_dict.each_key do |aa|
    if rc_dict['aa'].include?(aa).!
      if aa.size == 1
        raise PyteomicsError.new("No RC for residue '#{aa}'")
      end
      if ['', 0, nil, false, [], {}].include?(raise_no_mod).! && rc_dict['aa'].include?(aa[-1])
        rt += peptide_dict[aa] * rc_dict['aa'][aa[-1]]
      else
        raise PyteomicsError.new(
          "Residue '#{aa}' not found in RC_dict. Set raise_no_mod=False to ignore this error and use the RC for '#{aa[-1]} instead.")
      end
    else
      rt += peptide_dict[aa] * rc_dict['aa'][aa]
    end
  end
  
  length_correction_term = (
    1.0 + (rc_dict['lcp'] || 0) * Numpy.log(length(peptide_dict)))
  rt *= length_correction_term

  rt += rc_dict['const'] || 0
end


RCs_guo_ph2_0 = {
  'aa' => {
    'K' => -2.1,
    'G' => -0.2,
    'L' =>  8.1,
    'A' =>  2.0,
    'C' =>  2.6,
    'E' =>  1.1,
    'D' =>  0.2,
    'F' =>  8.1,
    'I' =>  7.4,
    'H' => -2.1,
    'M' =>  5.5,
    'N' => -0.6,
    'Q' =>  0.0,
    'P' =>  2.0,
    'S' => -0.2,
    'R' => -0.6,
    'T' =>  0.6,
    'W' =>  8.8,
    'V' =>  5.0,
    'Y' =>  4.5,
    'H-' => 0.0,
    '-OH' =>0.0
  },
  'lcp' => 0.0,
  'const' => 0.0
}

RCs_guo_ph7_0 = {
  'aa' => {
    'K' => -0.2,
    'G' => -0.2,
    'L' =>  9.0,
    'A' =>  2.2,
    'C' =>  2.6,
    'E' => -1.3,
    'D' => -2.6,
    'F' =>  9.0,
    'I' =>  8.3,
    'H' =>  2.2,
    'M' =>  6.0,
    'N' => -0.8,
    'Q' =>  0.0,
    'P' =>  2.2,
    'S' => -0.5,
    'R' =>  0.9,
    'T' =>  0.3,
    'W' =>  9.5,
    'V' =>  5.7,
    'Y' =>  4.6,
    'H-' => 0.0,
    '-OH' =>0.0
  },
  'lcp' => 0.0,
  'const' => 0.0
}

RCs_meek_ph2_1 = {
  'aa' => {
    'K' => -3.2,
    'G' => -0.5,
    'L' => 10.0,
    'A' => -0.1,
    'C' => -2.2,
    'E' => -7.5,
    'D' => -2.8,
    'F' => 13.9,
    'I' => 11.8,
    'H' =>  0.8,
    'M' =>  7.1,
    'N' => -1.6,
    'Q' => -2.5,
    'P' =>  8.0,
    'S' => -3.7,
    'R' => -4.5,
    'T' =>  1.5,
    'W' => 18.1,
    'V' =>  3.3,
    'Y' =>  8.2,
    'H-' => 0.0,
    '-OH' =>0.0
  },
  'lcp' => 0.0,
  'const' => 0.0
}

RCs_meek_ph7_4 = {
  'aa' => {
    'K' =>  0.1,
    'G' =>  0.0,
    'L' =>  8.8,
    'A' =>  0.5,
    'C' => -6.8,
    'E' =>-16.9,
    'D' => -8.2,
    'F' => 13.2,
    'I' => 13.9,
    'H' => -3.5,
    'M' =>  4.8,
    'N' =>  0.8,
    'Q' => -4.8,
    'P' =>  6.1,
    'S' =>  1.2,
    'R' =>  0.8,
    'T' =>  2.7,
    'W' => 14.9,
    'V' =>  2.7,
    'Y' =>  6.1,
    'H-' => 0.0,
    '-OH' =>0.0
  },
  'lcp' => 0.0,
  'const' => 0.0
}

RCs_browne_tfa = {
  'aa' => {
    'K' => -3.7,
    'G' => -1.2,
    'L' => 20.0,
    'A' =>  7.3,
    'C' => -9.2,
    'E' => -7.1,
    'D' => -2.9,
    'F' => 19.2,
    'I' =>  6.6,
    'H' => -2.1,
    'M' =>  5.6,
    'N' => -5.7,
    'Q' => -0.3,
    'P' =>  5.1,
    'S' => -4.1,
    'pS' =>-6.5,
    'R' => -3.6,
    'T' =>  0.8,
    'pT' =>-1.6,
    'W' => 16.3,
    'V' =>  3.5,
    'Y' =>  5.9,
    'pY' => 3.5,
    'H-' => 0.0,
    '-OH' =>0.0
  },
  'lcp' => 0.0,
  'const' => 0.0
}

RCs_browne_hfba = {
  'aa' => {
    'K' => -2.5,
    'G' => -2.3,
    'L' => 15.0,
    'A' =>  3.9,
    'C' =>-14.3,
    'E' => -7.5,
    'D' => -2.8,
    'F' => 14.7,
    'I' => 11.0,
    'H' =>  2.0,
    'M' =>  4.1,
    'N' => -2.8,
    'Q' =>  1.8,
    'P' =>  5.6,
    'S' => -3.5,
    'pS' =>-7.6,
    'R' =>  3.2,
    'T' =>  1.1,
    'pT' =>-3.0,
    'W' => 17.8,
    'V' =>  2.1,
    'Y' =>  3.8,
    'pY' =>-0.3,
    'H-' => 0.0,
    '-OH' =>0.0
  },
  'lcp' => 0.0,
  'const' => 0.0
}

RCs_palmblad = {
  'aa' => {
    'K' => -0.66,
    'G' => -0.29,
    'L' =>  2.28,
    'A' =>  0.41,
    'C' => -1.32,
    'E' => -0.26,
    'D' =>  0.04,
    'F' =>  2.68,
    'I' =>  2.70,
    'H' =>  0.57,
    'M' =>  0.98,
    'N' => -0.54,
    'Q' =>  1.02,
    'P' =>  0.97,
    'S' => -0.71,
    'R' => -0.76,
    'T' =>  0.37,
    'W' =>  4.68,
    'V' =>  2.44,
    'Y' =>  2.78,
    'H-' => 0.0,
    '-OH' =>0.0
  },
  'lcp' => 0.0,
  'const' => 0.0
}

RCs_yoshida = {
  'aa' =>{
    'K' =>  2.77,
    'G' => -0.16,
    'L' => -2.31,
    'A' =>  0.28,
    'C' =>  0.80,
    'camC' =>  0.80,
    'E' =>  1.58,
    'D' =>  2.45,
    'F' => -2.94,
    'I' => -1.34,
    'H' =>  3.44,
    'M' => -0.14,
    'N' =>  3.25,
    'Q' =>  2.35,
    'P' =>  0.77,
    'S' =>  2.53,
    'R' =>  3.90,
    'T' =>  1.73,
    'W' => -1.80,
    'V' => -2.19,
    'Y' => -0.11,
    'H-' => 0.0,
    '-OH' =>0.0
  },
  'lcp' => 0.0,
  'const' => 0.0
}

RCs_yoshida_lc = {
  'aa' => {
    'A' => 1.29,
    'C' => 0.94,
 'camC' => 0.94,
    'D' => 3.89,
    'E' => 4.40,
    'F' => -4.18,
    'G' => 1.29,
    'H' => 7.57,
    'I' => -2.65,
    'K' => 7.33,
    'L' => -3.93,
    'M' => -1.48,
    'N' => 6.65,
    'P' => 1.03,
    'Q' => 6.68,
    'R' => 7.08,
    'S' => 5.09,
    'T' => 3.46,
    'V' => -2.52,
    'W' => -1.87,
    'Y' => -0.46,
    'H-' => 0.0,
    '-OH' => 0.0
  },
  'const' => 0.0,
  'lcp' => -0.2
}

RCs_zubarev = {
  'aa' => {
    'A' => 6.73,
    'E' => 5.66,
    'C' => 3.25,
    'D' => 5.64,
    'G' => 2.35,
    'F' => 27.43,
    'I' => 20.50,
    'H' => -0.66,
    'K' => -4.47,
    'M' => 17.39,
    'L' => 23.38,
    'N' => 2.57,
    'Q' => 2.93,
    'P' => 5.66,
    'S' => 3.58,
    'R' => -2.55,
    'T' => 4.88,
    'Y' => 13.22,
    'W' => 31.27,
    'V' => 13.05,
   'camC' => 3.25,
    'C' => 3.25,
   'oxM' => -7.61,
   '-OH' => 0.0,
   'H-' => 0.0
   },
  'const' => 0.53,
  'lcp' => -0.21
}

RCs_gilar_atlantis_ph3_0 = {
  'aa' => {
    'K' => 15.90,
    'R' => 13.64,
    'H' => 12.94,
    'E' => 2.97,
    'P' => 4.77,
    'Q' => 5.43,
    'D' => 3.20,
    'C*' => 4.87,
    'C' => 4.87,
    'N' => 3.91,
    'A' => 3.34,
    'G' => 3.33,
    'S' => 3.04,
    'T' => 2.71,
    'V' => 1.75,
    'I' => 0.65,
    'M' => 1.13,
    'L' => 0.13,
    'F' => -1.17,
    'Y' => -0.22,
    'W' => -2.47
  },
  'lcp' => 0.0,
  'const' => 21.33
}

RCs_gilar_atlantis_ph4_5 = {
  'aa' => {
    'K' => 15.49,
    'R' => 13.33,
    'H' => 12.19,
    'E' => 6.93,
    'P' => 5.89,
    'Q' => 5.68,
    'D' => 5.31,
    'C*' => 5.23,
    'C' => 5.23,
    'N' => 4.07,
    'A' => 3.6,
    'G' => 3.46,
    'S' => 2.62,
    'T' => 2.33,
    'V' => 1.42,
    'I' => 0.84,
    'M' => 0.34,
    'L' => 0.29,
    'F' => -1.21,
    'Y' => -1.62,
    'W' => -2.08
  },
  'lcp' => 0.0,
  'const' => 23.95
}

RCs_gilar_atlantis_ph10_0 = {
  'aa' => {
    'K' => 25.23,
    'R' => 23.38,
    'H' => 5.94,
    'E' => 0.59,
    'P' => 4.00,
    'Q' => 3.53,
    'D' => -0.84,
    'C*' => 3.52,
    'C' => 3.52,
    'N' => 3.26,
    'A' => 3.64,
    'G' => 3.02,
    'S' => 2.28,
    'T' => 1.74,
    'V' => 1.05,
    'I' => 1.51,
    'M' => -0.61,
    'L' => 0.25,
    'F' => -0.17,
    'Y' => -0.79,
    'W' => 0.23
  },
  'lcp' => 0.0,
  'const' => 13.78
}

RCs_gilar_beh = {
  'aa' => {
    'K' => 9.49,
    'R' => 8.56,
    'H' => 8.40,
    'E' => 5.95,
    'P' => 4.73,
    'Q' => 4.65,
    'D' => 4.97,
    'C' => 3.47,
    'C*' => 3.47,
    'N' => 3.50,
    'A' => 2.90,
    'G' => 2.63,
    'S' => 2.14,
    'T' => 2.19,
    'V' => 1.71,
    'I' => 1.30,
    'M' => 1.40,
    'L' => 0.73,
    'F' => -0.09,
    'Y' => -0.40,
    'W' => 0.11
  },
  'lcp' => 0.0,
  'const' => 18.41
}

RCs_gilar_beh_amide = {
  'aa' => {
    'K' => 7.19,
    'R' => 6.68,
    'H' => 6.16,
    'E' => 6.11,
    'P' => 3.18,
    'Q' => 5.19,
    'D' => 6.02,
    'C*' => 3.71,
    'C' => 3.71,
    'N' => 4.16,
    'A' => 2.64,
    'G' => 3.12,
    'S' => 3.17,
    'T' => 3.41,
    'V' => 0.83,
    'I' => -0.69,
    'M' => -0.12,
    'L' => -1.24,
    'F' => -1.93,
    'Y' => 0.46,
    'W' => -2.11
  },
  'lcp' => 0.0,
  'const' => 24.26
}

RCs_gilar_rp = {
  'aa' => {
    'K' => -1.015,
    'R' => -0.681,
    'H' => -1.937,
    'E' => 1.475,
    'P' => 3.496,
    'Q' => 1.228,
    'D' => 1.326,
    'C' => 1.832,
    'C*' => 1.832,
    'N' => 0.299,
    'A' => 2.322,
    'G' => 1.172,
    'S' => 1.165,
    'T' => 1.894,
    'V' => 5.695,
    'I' => 8.343,
    'M' => 5.128,
    'L' => 9.069,
    'F' => 10.877,
    'Y' => 5.603,
    'W' => 12.183
  },
  'lcp' => 0.0,
  'const' => -3.696
}

RCs_krokhin_100A_fa = {
  'aa' => {
    'K' => -5.08,
    'G' => -0.07,
    'L' =>  9.89,
    'A' =>  1.63,
    'C' =>  0.7,
  'camC' =>  0.7,
    'E' =>  1.75,
    'D' =>  0.95,
    'F' =>  11.92,
    'I' =>  9.06,
    'H' => -5.05,
    'M' =>  6.96,
    'N' => -0.59,
    'Q' =>  0.2,
    'P' =>  1.98,
    'S' => 0.27,
    'R' => -3.55,
    'T' =>  1.37,
    'W' =>  13.67,
    'V' =>  5.72,
    'Y' =>  5.97
  },
  'lcf' => 0.0,
  'const' => 0.0
}

RCs_krokhin_100A_tfa = {
  'aa' => {
    'K' => -3.53,
    'G' => -0.35,
    'L' =>  9.44,
    'A' =>  1.11,
    'C' =>  0.04,
  'camC' =>  0.04,
    'E' =>  1.08,
    'D' =>  -0.22,
    'F' =>  11.34,
    'I' =>  7.86,
    'H' => -3.04,
    'M' =>  6.57,
    'N' => -1.44,
    'Q' =>  -0.53,
    'P' =>  1.62,
    'S' => -0.33,
    'R' => -2.58,
    'T' =>  0.48,
    'W' =>  13.12,
    'V' =>  4.86,
    'Y' =>  5.4
  },
  'lcf' => 0.0,
  'const' => 0.0
}