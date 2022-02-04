# import re
# from collections import deque
# import itertools as it
# from .auxiliary import PyteomicsError, memoize, BasicComposition, cvstr, cvquery
require 'set'
require_relative 'auxiliary/utils'
require_relative 'auxiliary/structures'


STD_amino_acids = ['Q', 'W', 'E', 'R', 'T', 'Y', 'I', 'P', 'A', 'S',
  'D', 'F', 'G', 'H', 'K', 'L', 'C', 'V', 'N', 'M']

STD_nterm = 'H-'

STD_cterm = '-OH'

STD_labels = STD_amino_acids + [STD_nterm, STD_cterm]

def is_term_mod(label)
  _nterm_mod = /[^-]+-$/
  _cterm_mod = /-[^-]+$/
  _nterm_mod.match(label) || _cterm_mod.match(label)
end

def match_modX(label)
  _modX_single = /^([^A-Z-]*)([A-Z])$/
  _modX_single.match(label)

  m = _modX_single.match(label)
  mod = m[1]
  x = m[0].sub(mod,'')
  [mod, x]
end

def is_modX(label)
  match_modX(label)[0].nil?.!
end

def length(sequence, **kwargs)
  return 0 if ['', nil, 0, false, [], {}].include?(sequence)
  if sequence.instance_of?(String) || sequence.instance_of?(Array)
    if sequence.instance_of?(String)
      parsed_sequence = parse(sequence, **kwargs)
    else
      parsed_sequence = sequence
    end
    num_term_groups = 0
    num_term_groups += 1 if is_term_mod(parsed_sequence[0])
    num_term_groups += 1 if is_term_mod(parsed_sequence[-1])
    return parsed_sequence.size - num_term_groups
  # elsif sequence.instance_of?(Hash)
  elsif sequence.is_a?(Hash)
    return sequence.sum{ |k, v| is_term_mod(k).! ? v : 0 }
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
    [x]
  else
    [mod, x]
  end
end

# _modX_sequence = re.compile(r'^([^-]+-)?((?:[^A-Z-]*[A-Z])+)(-[^-]+)?$')
# _modX_group = re.compile(r'[^A-Z-]*[A-Z]')
# _modX_split = re.compile(r'([^A-Z-]*)([A-Z])')
# _modX_single = re.compile(r'^([^A-Z-]*)([A-Z])$')

def parse(sequence, show_unmodified_termini: false, split: false, allow_unknown_modifications: false, **kwargs)
  sequence = sequence.to_s
  return if sequence.empty?
  _modX_sequence = Regexp.compile('^([^-]+-)?((?:[^A-Z-]*[A-Z])+)(-[^-]+)?$')

  begin
    m = _modX_sequence.match(sequence)
    n, body, c = m[1], m[2], m[3] if m.nil?.!
  rescue => exception
    raise PyteomicsError.new('Not a valid modX sequence: ' + sequence)
  end

  labels = kwargs['labels']
  if c.nil? && n.nil?.!
    labels = STD_labels if labels.nil?
    if n != STD_nterm && labels.include?(n).!
      c = '-' + body
      body = n[0..-2]
      n = nil
    end
  end
  
  _modX_split = Regexp.compile('([^A-Z-]*)([A-Z])')
  _modX_group = Regexp.compile('[^A-Z-]*[A-Z]')
  if split
    parsed_sequence = body.scan(_modX_split).map{ ['', 0, nil, false, [], {}].include?(_1[0]).! ? _1 : [_1[1]] }
  else
    parsed_sequence = body.scan(_modX_group)
  end
  nterm, cterm = (['', nil, 0, false, [], {}].include?(n).! ? n : STD_nterm), (['', nil, 0, false, [], {}].include?(c).! ? c : STD_cterm)
  if labels.nil?.!
    if labels.instance_of?(String)
      labels = Set.new(labels.split(''))
    elsif labels.instance_of?(Hash)
      labels = labels.keys.to_set
    else
      labels = Set.new(labels)
    end
    [n, c].zip([STD_nterm, STD_cterm]).each do |term, std_term|
      if ['', 0, nil, false, [], {}].include?(term).! && labels.include?(term).! && allow_unknown_modifications.!
        return PyteomicsError.new("Unknown label: #{term}")
      end
    end
    parsed_sequence.each do |group|
      if split
        mod, x = group.size == 2 ? group : ['', group[0]]
      else
        m = _modX_split.match(group)
        mod = m[1]
        x = m[0].sub(mod,'')
      end
      if (mod == '' && labels.include?(x).!) ||
         (labels.include?(mod + x) || (labels.include?(x) && (labels.include?(mod) || allow_unknown_modifications))).!
        return PyteomicsError.new("Unknown label: #{group}")
      end
    end
  end

  if show_unmodified_termini || nterm != STD_nterm
    if split
      parsed_sequence[0] = [nterm] + parsed_sequence[0].compact
    else
      parsed_sequence.insert(0, nterm)
    end
  end
  if show_unmodified_termini || cterm != STD_cterm
    if split
      parsed_sequence[-1] = parsed_sequence[-1].compact + [cterm]
    else
      parsed_sequence << cterm
    end
  end
  parsed_sequence
end

def valid(...)
  m = ''
  begin
    m = parse(...)
  rescue => exception
    return false
  end
  if m.instance_of?(PyteomicsError)
    false
  else
    true
  end
end

def fast_valid(sequence, labels: Set.new(STD_labels))
  if sequence.instance_of?(String)
    Set.new(sequence.split('')).subset?(labels)
  else
    Set.new(sequence).subset?(labels)
  end
end

def tostring(parsed_sequence, show_unmodified_termini: true)
  if parsed_sequence.is_a?(Hash)
    parsed_sequence = parsed_sequence.keys
  elsif parsed_sequence.is_a?(String)
    parsed_sequence = parsed_sequence.split('')
  end
  labels = []
  nterm = parsed_sequence[0]
  cterm = parsed_sequence[-1]

  if nterm.instance_of?(String)
    if nterm != STD_nterm || show_unmodified_termini
      labels << nterm
    end
    labels.concat(parsed_sequence[1...-1])
    if parsed_sequence.size > 1 && (cterm != STD_cterm || show_unmodified_termini)
      labels << cterm
    end
  else
    if parsed_sequence.size == 1
      g = nterm
      if nterm[0] == STD_nterm && show_unmodified_termini.!
        g = g[1..]
      end
      if nterm[-1] == STD_cterm && show_unmodified_termini.!
        g = g[0...-1]
      end
      return g.join('')
    end
    if nterm[0] != STD_nterm || show_unmodified_termini
      labels << nterm.join('')
    else
      labels << nterm[1..].join('')
    end
    # parsed_sequence[1...-1].map{ labels << _1.join('') }
    parsed_sequence[1...-1].each do |e|
      if e.instance_of?(String)
        labels << e
      else
        labels << e.join('')
      end
    end
    if parsed_sequence.size > 1
      if cterm[-1] != STD_cterm || show_unmodified_termini
        if cterm.instance_of?(String)
          labels << cterm
        else
          labels << cterm.join('')
        end
      else
        labels << cterm[0...-1].join('')
      end
    end
  end
  labels.join('')
end

def amino_acid_composition(sequence, show_unmodified_termini: false, term_aa: false, allow_unknown_modifications: false, **kwargs)
  labels = kwargs['labels']

  if sequence.instance_of?(String)
    parsed_sequence = parse(sequence, show_unmodified_termini: show_unmodified_termini, allow_unknown_modifications: allow_unknown_modifications, 'labels' => labels)
  elsif sequence.instance_of?(Array)
    if ['', nil, 0, [], {}].include?(sequence).! && sequence[0].instance_of?(Array)
      parsed_sequence = parse(sequence.to_s || true, show_unmodified_termini: show_unmodified_termini, allow_unknown_modifications: allow_unknown_modifications, 'labels' => labels)
    else
      parsed_sequence = sequence.dup
    end
  else
    raise PyteomicsError.new("Unsupported type of a sequence. Must be str or list, not #{sequence.class}")
  end

  aa_dict = BasicComposition.new()

  if ['', nil, false, 0, [], {}].include?(term_aa).!
    nterm_aa_position = is_term_mod(parsed_sequence[0]) ? 1 : 0
    cterm_aa_position = (
      is_term_mod(parsed_sequence[-1]) ? parsed_sequence.size - 2 : parsed_sequence.size - 1
    )
    if parsed_sequence.size > 1
      aa_dict['cterm' + parsed_sequence.delete_at(cterm_aa_position)] = 1
    end
    aa_dict['nterm' + parsed_sequence.delete_at(nterm_aa_position)] = 1
  end

  parsed_sequence.each do |aa|
    aa_dict[aa] += 1
  end

  aa_dict
end

class Deque
  attr_accessor :que
  def initialize(args, maxlen: nil)
    @que = args
    @maxlen = maxlen
  end

  def push(args)
    @que.push args
    while @maxlen.nil?.! && @que.size > @maxlen
      @que.shift
    end
  end

  def unshift(args)
    @que.unshift args
    while @maxlen.nil?.! && @que.size > @maxlen
      @que.pop
    end
  end
end

#memoize()
def cleave(...)
  Set.new(_cleave(...))
end

def _cleave(sequence, rule, missed_cleavages: 0, min_length: nil, semi: false, except: nil)
  if EXPASY_rules.include?(rule)
    rule = EXPASY_rules[rule]
  elsif PSIMS_rules.include?(rule)
    rule = PSIMS_rules[rule]
  elsif PSIMS_index.include?(rule)
    rule = PSIMS_index[rule]
  end
  except = EXPASY_rules[except] || except
  peptides = []
  ml = missed_cleavages + 2
  trange = (0...ml)
  cleavage_sites = Deque.new([0], maxlen: ml)
  if min_length.nil?
    min_length = 1
  end
  cl = 1
  if except.nil?.!
    exceptions = Set.new
    i = 0
    while m = sequence.match(except, i)
      break if m.empty?
      i = m.end(0)
      exceptions << i
    end
    exceptions << nil
  end

  @md = Struct.new(:text, :begin, :end)
  def find_iter(rgx, text)
    idx = 0
    while (m = rgx.match(text, idx))
      yield @md.new(m[0], m.begin(0), m.end(0))
      idx = m.end(0)
      break if m[0].empty?
    end
  end
  rule = Regexp.new(rule) if rule.instance_of?(String)
  a = enum_for(:find_iter, rule, sequence).map(&:end) + [nil]

  a.each do |i|
    next if !!except && exceptions.include?(i)
    cleavage_sites.push(i)
    cl += 1 if cl < ml
    trange.to_a[0...cl - 1].each do |j|
      seq = sequence[cleavage_sites.que[j]...cleavage_sites.que[-1]]
      if ['', nil, false, 0, [], {}].include?(seq).! && seq.size >= min_length
        peptides << seq
        if ['', nil, false, 0, [], {}].include?(semi).!
          (min_length...seq.size - 1).map{ peptides << seq[0..._1] }
          (1...seq.size - min_length + 1).map{ peptides << seq[_1..-1] }
        end
      end
    end
  end
  peptides
end

def num_sites(sequence, rule, **kwargs)
  _cleave(sequence, rule, **kwargs).size - 1
end

EXPASY_rules = {
  'arg-c' =>         /R/,
  'asp-n' =>         /\w(?=D)/,
  'bnps-skatole'  => /W/,
  'caspase 1' =>     /(?<=[FWYL]\w[HAT])D(?=[^PEDQKR])/,
  'caspase 2' =>     /(?<=DVA)D(?=[^PEDQKR])/,
  'caspase 3' =>     /(?<=DMQ)D(?=[^PEDQKR])/,
  'caspase 4' =>     /(?<=LEV)D(?=[^PEDQKR])/,
  'caspase 5' =>     /(?<=[LW]EH)D/,
  'caspase 6' =>     /(?<=VE[HI])D(?=[^PEDQKR])/,
  'caspase 7' =>     /(?<=DEV)D(?=[^PEDQKR])/,
  'caspase 8' =>     /(?<=[IL]ET)D(?=[^PEDQKR])/,
  'caspase 9' =>     /(?<=LEH)D/,
  'caspase 10' =>    /(?<=IEA)D/,
  'chymotrypsin high specificity'  => /([FY](?=[^P]))|(W(?=[^MP]))/,
  'chymotrypsin low specificity' =>
      /([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))/,
  'clostripain' =>   /R/,
  'cnbr' =>          /M/,
  'enterokinase' =>  /(?<=[DE]{3})K/,
  'factor xa' =>     /(?<=[AFGILTVM][DE]G)R/,
  'formic acid' =>   /D/,
  'glutamyl endopeptidase' => /E/,
  'granzyme b' =>    /(?<=IEP)D/,
  'hydroxylamine' => /N(?=G)/,
  'iodosobenzoic acid' => /W/,
  'lysc' =>          /K/,
  'ntcb' =>          /\w(?=C)/,
  'pepsin ph1.3' =>  /((?<=[^HKR][^P])[^R](?=[FL][^P]))|((?<=[^HKR][^P])[FL](?=\w[^P]))/,
  'pepsin ph2.0' =>  /((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|((?<=[^HKR][^P])[FLWY](?=\w[^P]))/,
  'proline endopeptidase' => /(?<=[HKR])P(?=[^P])/,
  'proteinase k' =>  /[AEFILTVWY]/,
  'staphylococcal peptidase i' => /(?<=[^E])E/,
  'thermolysin' =>   /[^DE](?=[AFILMV])/,
  'thrombin' =>      /((?<=G)R(?=G))|((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))/,
  'trypsin' =>       /([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))/,
  'trypsin_exception' => /((?<=[CD])K(?=D))|((?<=C)K(?=[HY]))|((?<=C)R(?=K))|((?<=R)R(?=[HR]))/,
}

PSIMS_rules = {
  Cvstr.new('2-iodobenzoate', accession: 'MS:1001918') => /(?<=W)/,
  Cvstr.new('Arg-C', accession: 'MS:1001303') => /(?<=R)(?!P)/,
  Cvstr.new('Asp-N', accession: 'MS:1001304') => /(?=[BD])/,
  Cvstr.new('Asp-N ambic', accession: 'MS:1001305') => /(?=[DE])/,
  Cvstr.new('CNBr', accession: 'MS:1001307') => /(?<=M)/,
  Cvstr.new('Chymotrypsin', accession: 'MS:1001306') => /(?<=[FYWL])(?!P)/,
  Cvstr.new('Formic acid', accession: 'MS:1001308') => /((?<=D))|((?=D))/,
  Cvstr.new('Lys-C', accession: 'MS:1001309') => /(?<=K)(?!P)/,
  Cvstr.new('Lys-C/P', accession: 'MS:1001310') => /(?<=K)/,
  Cvstr.new('PepsinA', accession: 'MS:1001311') => /(?<=[FL])/,
  Cvstr.new('TrypChymo', accession: 'MS:1001312') => /(?<=[FYWLKR])(?!P)/,
  Cvstr.new('Trypsin', accession: 'MS:1001251') => /(?<=[KR])(?!P)/,
  Cvstr.new('Trypsin/P', accession: 'MS:1001313') => /(?<=[KR])/,
  Cvstr.new('V8-DE', accession: 'MS:1001314') => /(?<=[BDEZ])(?!P)/,
  Cvstr.new('V8-E', accession: 'MS:1001315') => /(?<=[EZ])(?!P)/,
  Cvstr.new('glutamyl endopeptidase', accession: 'MS:1001917') => /(?<=[^E]E)/,
  Cvstr.new('leukocyte elastase', accession: 'MS:1001915') => /(?<=[ALIV])(?!P)/,
  Cvstr.new('proline endopeptidase', accession: 'MS:1001916') => /(?<=[HKR]P)(?!P)/,
}

PSIMS_index = Cvquery.new.__call__(PSIMS_rules)

def isoforms(sequence, **kwargs)
  def main(group)
    if group[-1][0] == '-'
      i = -2
    else
      i = -1
    end
    [group.size + i, group[i]]
  end

  def apply_mod(label, mod)
    group = label.dup.to_a
    m = main(group)[0]
    c = true
    if m == 0 && is_term_mod(mod).!
      group.insert(0, mod)
    elsif mod[0] == '-' && (group[-1] == STD_cterm || (group[-1][0] == '-' && @override))
      group[-1] = mod
    elsif mod[-1] == '-' && (group[0] == STD_nterm || (group[0][-1] == '-' && @override))
      group[0] = mod
    elsif is_term_mod(mod).!
      if ['', nil, false, 0, [], {}].include?(m).! && group[m - 1][-1] != '-'
        if ['', nil, false, 0, [], {}].include?(@override).!
          group[m - 1] = mode
        else
          c = false
        end
      else
        group.insert(m, mod)
      end
    else
      c = false
    end
    if c
      return group
    end
  end

  variable_mods = kwargs['variable_mods'] || {}
  varmods_term, varmods_non_term = [], []
  variable_mods.sort.each do |m, r|
    if is_term_mod(m)
      varmods_term << [m, r]
    else
      varmods_non_term << [m, r]
    end
  end
  fixed_mods = kwargs['fixed_mods'] || {}
  parse_kw = {}
  if kwargs.include?('labels')
    parse_kw['labels'] = kwargs['labels'].to_a + fixed_mods.keys
  end
  parsed = parse(sequence, show_unmodified_termini: true, split: true, **parse_kw)
  @override = kwargs['override'] || false
  @show_unmodified_termini = kwargs['show_unmodified_termini'] || false
  @max_mods = kwargs['max_mods']
  format_ = kwargs['format'] || 'str'

  fixed_mods.each do |cmod, res|
    parsed.each_with_index do |group, i|
      if res == true || res.include?(main(group)[1])
        parsed[i] = apply_mod(group, cmod) || parsed[i]
      end
    end
  end

  @states = [[parsed[0]]]
  m0 = main(parsed[0])[1]
  varmods_non_term.each do |m, r|
    if r == true || r.include?(m0) || r.include?('nterm' + m0) || parsed.size == 1 && r.include?('cterm' + m0)
      applied = apply_mod(parsed[0], m)
      if applied.nil?.!
        @states[0] << applied
      end
    end
  end
  more_states = []
  varmods_term.each do |m, r|
    if r == true || r.include?(m0)
      if m[-1] == '-' || parsed.size == 1
        @states[0].each do |group|
          applied = apply_mod(group, m)
          if applied.nil?.!
            more_states << applied
          end
        end
      end
    end
  end
  @states[0].concat(more_states)

  parsed[1...-1].each do |group|
    gstates = [group]
    varmods_non_term.each do |m, r|
      if r == true || r.include?(group[-1])
        applied = apply_mod(group, m)
        if applied.nil?.!
          gstates << applied
        end
      end
    end
    @states << gstates
  end

  if parsed.size > 1
    @states << [parsed[-1]]
    m1 = main(parsed[-1])[1]
    varmods_non_term.each do |m, r|
      if r == true || r.include?(m1) || r.include?('cterm' + m1) || parsed.size == 1 && r.include?('nterm' + m1)
        applied = apply_mod(parsed[-1], m)
        if applied.nil?.!
          @states[-1] << applied
        end
      end
    end
    more_states = []
    varmods_term.each do |m, r|
      if r == true || r.include?(m1)
        if m[0] == '-' || parsed.size == 1
          @states[-1].each do |group|
            applied = apply_mod(group, m)
            if applied.nil?.!
              more_states << applied
            end
          end
        end
      end
    end
    @states[-1].concat(more_states)
  end

  # sites = states.select{ _1[0].size > 1 }
  # sites = states.flatten(1).select{ _1[0].size > 1 }
  @sites = []
  @states.each_with_index do |x, i|
    @sites << [i, x] if x.size > 1
  end
  if @max_mods.nil? || @max_mods > @sites.size
    @possible_states = @states[0].product(*@states[1..])
  else
    @new_state_lists = []
    def state_lists()
      @max_mods.next.times do |m|
        @sites.combination(m).each do |comb|
          skel = @states.map{ [_1[0]] }.dup
          comb.each do |i, e|
            skel[i] = e[1..]
          end
          @new_state_lists << skel
        end
      end
    end
    state_lists()
    @possible_states = @new_state_lists.map{ |e| e[0].product(*e[1..]) }.inject(&:concat).flatten(1)
  end

  if format_ == 'split'
    @new_states = []
    def strip_std_terms
      @possible_states.each do |ps|
        if @show_unmodified_termini.!
          if ps[0][0] == STD_nterm
            ps[0] = ps[0][1..].dup
          end
          if ps[-1][-1] == STD_cterm
            ps[-1] = ps[-1][0...-1].dup
          end
        end
        @new_states << ps
      end
    end
    strip_std_terms()
    return @new_states
  elsif format_ == 'str'
    return @possible_states.map{ tostring(_1, show_unmodified_termini: @show_unmodified_termini) }
  else
    raise PyteomicsError.new("Unsupported value of 'format': #{format_}")
  end
end

def coverage(protein, peptides)
  require 'numpy'
  protein = protein.gsub(/[^A-Z]/, '')
  mask = Numpy.zeros(pritein.size, dtype: Numpy.int8)
  peptides.each do |peptide|
    s = peptide.gsub(/[^A-Z]/, '')
    indices = protein.scan(%r[?=#{peptide.gsub(/[^A-Z]/, '')}]).map{ text.index(_1) }
    indices.each do |i|
      mask[i...i + peptide.size] = 1
    end
  end
  mask.sum(dtype: float) / mask.size
end