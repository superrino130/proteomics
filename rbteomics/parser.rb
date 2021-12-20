# import re
# from collections import deque
# import itertools as it
# from .auxiliary import PyteomicsError, memoize, BasicComposition, cvstr, cvquery

STD_amino_acids = ['Q', 'W', 'E', 'R', 'T', 'Y', 'I', 'P', 'A', 'S',
  'D', 'F', 'G', 'H', 'K', 'L', 'C', 'V', 'N', 'M']
"""modX labels for the 20 standard amino acids."""

STD_nterm = 'H-'
"""modX label for the unmodified N-terminus."""

STD_cterm = '-OH'
"""modX label for the unmodified C-terminus."""

STD_labels = std_amino_acids + [std_nterm, std_cterm]
"""modX labels for the standard amino acids and unmodified termini."""


def is_term_mod(label)
  _nterm_mod = /[^-]+-$/
  _cterm_mod = /-[^-]+$/
  _ntrem_mod.match(label) || !!_cterm_mod.match(label)
end

def match_modX(label)
  _modX_single = /^([^A-Z-]*)([A-Z])$/
  _modX_single.match(label)
end

def is_modX(label)
  !!match_modX(label)
end

def length(sequence, **kwargs)
  return 0 if ['', [], {}].include?(sequence).!
  if sequence.instance_of?(String) || sequence.instance_of?(Array)
    if sequence.instance_of?(String)
      parsed_sequence = parse(...)
    else
      parsed_sequence = sequence
    end
    num_term_groups = 0
    num_term_groups += 1 if is_term_mod(parsed_sequence[0])
    num_term_groups += 1 if is_term_mod(parsed_sequence[-1])
    return parsed_sequence.size - num_term_groups
  elsif sequence.instance_of?(Hash)
    return sequence.select{ |k, v| is_term_mod(k) }.sum
  end

  raise PyteomicsError.new('Unsupported type of sequence.')
end

def _split_label(label)
  begin
    mod, x = match_modX(label)
  rescue => exception
    raise PyteomicsError.new("Cannot split a non-modX label: #{label}")    
  end
  if mod.!
    [x,]
  else
    [mod, x]
  end
end

# _modX_sequence = re.compile(r'^([^-]+-)?((?:[^A-Z-]*[A-Z])+)(-[^-]+)?$')
# _modX_group = re.compile(r'[^A-Z-]*[A-Z]')
# _modX_split = re.compile(r'([^A-Z-]*)([A-Z])')
# _modX_single = re.compile(r'^([^A-Z-]*)([A-Z])$')

def parse(sequence, show_unmodified_termini=false, split=false, allow_unknown_modifications=false, **kwargs)
  sequence = sequence.to_s
  _modX_sequence = /^([^-]+-)?((?:[^A-Z-]*[A-Z])+)(-[^-]+)?$/

  begin
    n, body, c = _modX_sequence.match(sequence)
  rescue => exception
    raise PyteomicsError.new('Not a valid modX sequence: ' + sequence)
  end

  labels = kwargs['labels']
  if c.nil? && !!n
    labels = STD_labels if labels.nil?
    if n != STD_nterm && labels.include?(n).!
      c = '-' + body
      body = n[0..-2]
      n = nil
    end
  end

  _modX_split = /([^A-Z-]*)([A-Z])/
  _modX_group = /[^A-Z-]*[A-Z]/
  if ['', nil, 0, [], {}].include?(split).!
    parsed_sequence = body.scan(_modX_split).map{ _1[0] ? _1 : [_1[1],] }
  else
    parsed_sequence = body.scan(_modX_group)
  end
  nterm, cterm = (['', nil, [], {}, 0].include?(n).! ? n : STD_nterm), (['', nil, [], {}, 0].include?(c).! ? c : STD_cterm)

  if !!labels
    labels = Set.new(labels)
    [n, c].zip([STD_nterm, STD_cterm]).each do |term, std_term|
      if ['', nil, 0, [], {}].include?(term).! && labels.include?(term).! && allow_unknown_modifications.!
        raise PyteomicsError.new("Unknown label: #{term}")
      end
    end
    parsed_sequence.each do |group|
      if ['', nil, 0, [], {}].include?(split).!
        mod, x = group.size == 2 ? group : ['', group[0]]
      else
        mod, X = _modX_split.match(group)
      end
      if (mod.! && labels.include?(x).!) || ((labels.include?(mod + x) || labels.include?(x)).! && (labels.include?(mod) || allow_unknown_modifications))
        raise PyteomicsError.new("Unknown label: #{group}")
      end
    end
  end

  if ['', nil, 0, [], {}].include?(show_unmodified_termini).! || nterm != STD_nterm
    if ['', nil, 0, [], {}].include?(split).!
      parsed_sequence[0] = [nterm, ] + parsed_sequence[0]
    else
      parsed_sequence.insert(0, nterm)
    end
  end
  if ['', nil, 0, [], {}].include?(show_unmodified_termini).! || cterm != STD_cterm
    if ['', nil, 0, [], {}].include?(split).!
      parsed_sequence[-1] = parsed_sequence[-1] + [cterm,]
    else
      parsed_sequence << cterm
    end
  end
  parsed_sequence
end

def valid(*args, **kwargs)
  begin
    parse(...)
  rescue => exception
    return false
  end
  true
end

def fast_valid(sequence, labels=Set.new(STD_labels))
  Set.new(sequence).subset?(labels)
end

def tostring(parsed_sequence, show_unmodified_termini=true)
  parsed_sequence = parsed_sequence.to_a
  labels = []
  nterm = parsed_sequence[0]
  cterm = parsed_sequence[-1]

  if nterm.instance_of?(String)
    if nterm != STD_nterm || show_unmodified_termini
      labels << nterm
    end
    labels.concat(parsed_sequence[1..-2])
    if parsed_sequence.size > 1 && (cterm != STD_cterm || show_unmodified_termini)
      labels << cterm
    end
  else
    if parsed_sequence.size == 1
      g = nterm
      if nterm[0] == STD_nterm && show_unmodified_termini.!
        g = g[1..-1]
      end
      if nterm[-1] == STD_cterm && show_unmodified_termini.!
        g = g[0..-2]
      end
      return g.join('')
    end
    if nterm[0] != STD_nterm || show_unmodified_termini
      labels << nterm.join('')
    else
      labels << nterm[1..-1].join('')
    end
    parsed_sequence[1..-1].map{ labels << _1.join('') }
    if parsed_sequence.size > 1
      if cterm[-1] != STD_cterm || show_unmodified_termini
        labels << cterm.join('')
      else
        labels << cterm[0..-2]
      end
    end
  end
  labels.join('')
end

def amino_acid_composition(sequence, show_unmodified_termini=false, term_aa=false, allow_unknown_modifications=false, **kwargs)
  labels = kwargs['labels']

  if sequence.instance_of?(String)
    parsed_sequence = parse(sequence, show_unmodified_termini, allow_unknown_modifications=allow_unknown_modifications, labels=labels)
  elsif sequence.instance_of?(Array)
    if ['', nil, 0, [], {}].include?(sequence).! && sequence[0].instance_of?(Array)
      parsed_sequence = parse(sequence.to_s || true), show_unmodified_termini, allow_unknown_modifications=allow_unknown_modifications, labels=labels)
    else
      parsed_sequence = sequence
    end
  else
    raise PyteomicsError.new("Unsupported type of a sequence. Must be str or list, not #{sequence.class}"))
  end

  aa_dict = BasicComposition.new()

  if ['', nil, 0, [], {}].include?(term_aa).!
    nterm_aa_position = is_term_mod(parsed_sequence[0]) ? 1 : 0
    cterm_aa_position = (
      is_term_mod(parsed_sequence[-1]) ? parsed_sequence.size - 2 : parsed_sequence.size - 1
    )
    if parsed_sequence.size > 1
      aa_dict['cterm' + parsed_sequence.delete(cterm_aa_position)] = 1
    end
    aa_dict['nterm' + parsed_sequence.delete(nterm_aa_position)] = 1
  end

  parsed_sequence.each do |aa|
    aa_dict[aa] += 1
  end

  aa_dict
end