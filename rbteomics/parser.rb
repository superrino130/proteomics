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
  _nterm_mod.match(label) || !!_cterm_mod.match(label)
end

def match_modX(label)
  _modX_single = /^([^A-Z-]*)([A-Z])$/
  _modX_single.match(label)
end

def is_modX(label)
  !!match_modX(label)
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
  elsif sequence.instance_of?(Hash)
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

def parse(sequence, show_unmodified_termini: false, splitflg: false, allow_unknown_modifications: false, **kwargs)
  sequence = sequence.to_s
  return if sequence.empty?
  _modX_sequence = /^([^-]+-)?((?:[^A-Z-]*[A-Z])+)(-[^-]+)?$/

  begin
    m = _modX_sequence.match(sequence)
    n, body, c = m[1], m[2], m[3]
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
  if splitflg
    parsed_sequence = body.scan(_modX_split).map{ ['', 0, nil, false, [], {}].include?(_1[0]).! ? _1 : [_1[1]] }
  else
    parsed_sequence = body.scan(_modX_group)
  end
  nterm, cterm = (['', nil, 0, false, [], {}].include?(n).! ? n : STD_nterm), (['', nil, 0, false, [], {}].include?(c).! ? c : STD_cterm)
  if labels.nil?.!
    if labels.instance_of?(String)
      labels = Set.new(labels.split(''))
    else
      labels = Set.new(labels)
    end
    [n, c].zip([STD_nterm, STD_cterm]).each do |term, std_term|
      if ['', nil, 0, [], {}].include?(term).! && labels.include?(term).! && allow_unknown_modifications.!
        raise PyteomicsError.new("Unknown label: #{term}")
      end
    end
    parsed_sequence.each do |group|
      if splitflg
        mod, x = group.size == 2 ? group : ['', group[0]]
      else
        m = _modX_split.match(group)
        mod = m[1]
        x = m[0]
      end
      if (mod != '' && labels.include?(x).!) || ((labels.include?(mod + x) || labels.include?(x)).! && (labels.include?(mod) || allow_unknown_modifications))
        raise PyteomicsError.new("Unknown label: #{group}")
      end
    end
  end

  if show_unmodified_termini || nterm != STD_nterm
    if splitflg
      parsed_sequence[0] = [nterm] + parsed_sequence[0].compact
    else
      parsed_sequence.insert(0, nterm)
    end
  end
  if show_unmodified_termini || cterm != STD_cterm
    if splitflg
      parsed_sequence[-1] = parsed_sequence[-1].compact + [cterm]
    else
      parsed_sequence << cterm
    end
  end
  parsed_sequence
end

def valid(...)
  begin
    parse(...)
  rescue => exception
    return false
  end
  true
end

def fast_valid(sequence, labels: Set.new(STD_labels))
  if sequence.instance_of?(String)
    Set.new(sequence.split('')).subset?(labels)
  else
    Set.new(sequence).subset?(labels)
  end
end

def tostring(parsed_sequence, show_unmodified_termini: true)
  parsed_sequence = parsed_sequence.to_a
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
    parsed_sequence[1...-1].map{ labels << _1.join('') }
    if parsed_sequence.size > 1
      if cterm[-1] != STD_cterm || show_unmodified_termini
        labels << cterm.join('')
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
    parsed_sequence = parse(sequence, show_unmodified_termini: show_unmodified_termini, allow_unknown_modifications: allow_unknown_modifications, **{'labels' => labels})
  elsif sequence.instance_of?(Array)
    if ['', nil, 0, [], {}].include?(sequence).! && sequence[0].instance_of?(Array)
      parsed_sequence = parse(sequence.to_s || true, show_unmodified_termini: show_unmodified_termini, allow_unknown_modifications: allow_unknown_modifications, **{'labels' => labels})
    else
      parsed_sequence = sequence
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
      aa_dict.defaultdict['cterm' + parsed_sequence.delete_at(cterm_aa_position)] = 1
    end
    aa_dict.defaultdict['nterm' + parsed_sequence.delete_at(nterm_aa_position)] = 1
  end

  parsed_sequence.each do |aa|
    aa_dict.defaultdict[aa] += 1
  end

  aa_dict.defaultdict
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
  icleave(...).resume.map{ _2 }.to_set
end

def icleave(sequence, rule, missed_cleavages: 0, min_length: nil, max_length: nil, semi: false, except: nil, regex: false)
    rule2 = rule.inspect.gsub(/\//, '')
    rule2 = rule2.gsub(/\"/, '')
    if regex.!
      if EXPASY_rules.include?(rule2)
        rule = EXPASY_rules[rule2]
      elsif PSIMS_rules.include?(rule2)
        rule = PSIMS_rules[rule2]
      elsif PSIMS_index.include?(rule2)
        rule = PSIMS_index[rule2]
      elsif rule2.match(/[a-z]/)
        warn "Interpreting the rule as a regular expression: #{rule2}. Did you mistype the rule? Specify `regex=True` to silence this warning."
      end
    end
    except = EXPASY_rules[except] || except
    peptides = []
    ml = missed_cleavages + 2
    trange = (0...ml)
    cleavage_sites = Deque.new([0], maxlen: ml)
    if min_length.nil?
      min_length = 1
    end
    if max_length.nil?
      max_length = sequence.size
    end
    cl = 1
    if except.nil?.!
      excepts = sequence.scan(except).map{ sequence.index(_1) }
    end
    rule3 = %r[#{'^'+rule2}]

    a = []
    sequence.scan(rule).flatten.compact.uniq.each do |k|
      sequence.size.times do |i|
        a << i + k.size if sequence[i..].match(rule3)
      end
    end
    Fiber.new do
      a.uniq.each do |i|
      next if except.nil?.! && excepts.include?(i)
      cleavage_sites.push(i)
      cl += 1 if cl < ml
      trange.to_a[0...cl - 1].each do |j|
        seq = sequence[cleavage_sites.que[j]...cleavage_sites.que[-1]]
        lenseq = seq.size
        if i.nil?.!
          start = i - lenseq
        else
          start = sequence.size - lenseq
        end
        if ['', 0, nil, false, [], {}].include?(seq).! && min_length <= lenseq && lenseq <= max_length
          Fiber.yield [start, seq]
          if semi
            (min_length...[lenseq, max_length].min).each do |k|
              Fiber.yield [start, seq[0...k]]
            end
            ([1, lenseq - max_length].max...lenseq - min_length + 1).each do |k|
              Fiber.yield [start + k, seq[k..]]
            end
          end
        end
      end
    end
  end
end

def xcleave(...)
  icleave(...).resume
end

def num_sites(sequence, rule, **kwargs)
  icleave(sequence, rule, **kwargs).resume.size - 1
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

PSIMS_index = @cvquery.__call__(PSIMS_rules)

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
    group = label.to_a
    m = main(group)[0]
    c = true
    if m == 0 && is_term_mod(mod)
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
      group
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
    parse_kw['labels'] = kwargs['labels'].to_a + fixed_mods.to_a
  end
  parsed = parse(sequence, show_unmodified_termini: true, splitflg: true, **parse_kw)
  @override = kwargs['override'] || false
  show_unmodified_termini = kwargs['show_unmodified_termini'] || false
  @max_mods = kwargs['@max_mods']
  format_ = kwargs['format'] || 'str'

  fixed_mods.each do |cmod, res|
    parsed.each_with_index do |group, i|
      if !!res || res.include?(main(group)[1])
        parsed[i] = apply_mod(group, cmod) || parsed[i]
      end
    end
  end

  states = [[parsed[0]]]
  m0 = main(parsed[0])[1]
  varmods_non_term.each do |m, r|
    if !!r || r.include?(m0) || r.include?('nterm' + m0) || parsed.size == 1 && r.include?('cterm' + m0)
      applid = apply_mod(parsed[0], m)
      if !!applid
        states[0] << applid
      end
    end
  end
  more_states = []
  varmods_term.each do |m, r|
    if !!r || r.include?(m0)
      if m[-1] == '-' || parsed.size == 1
        states[0].each do |group|
          applid = apply_mod(group, m)
          if !!applid
            more_states << applid
          end
        end
      end
    end
  end
  states[0].concat(more_states)

  parsed[1...-1].each do |group|
    gstates = [group]
    varmods_non_term.each do |m, r|
      if !!r || r.include?(group[-1])
        applid = apply_mod(group, m)
        if !!applid
          gstates << applid
        end
      end
    end
    states << gstates
  end

  if parsed.size > 1
    states << [parsed[-1]]
    m1 = main(parsed[-1])[1]
    varmods_non_term.each do |m, r|
      if !!r || r.include?(m1) || r.include?('cterm' + m1) || parsed.size == 1 && r.include?('nterm' + m1)
        applid = apply_mod(parsed[-1], m)
        if !!applid
          states[-1] << applid
        end
      end
    end
    more_states = []
    varmods_term.each do |m, r|
      if !!r || r.include?(m1)
        if m[0] == '-' || parsed.size == 1
          states[-1].each do |group|
            applid = apply_mod(group, m)
            if !!applid
              more_states << applid
            end
          end
        end
      end
    end
    states[-1].concat(more_states)
  end

  sites = states.select{ _1[1].size > 1 }
  if @max_mods.nil? || @max_mods > sites.size
    @possible_states = states.inject(:product).map(&:flatten)
  else
    def state_lists
      Fiber.new do
        (0..@max_mods).to_a.each do |m|
          sites.combination(m).each do |comb|
            skel = states.map{ [_1[0]] }
            comb.each do |i, e|
              skel[i] = e[1..-1]
            end
            yield skel
          end
        end
      end
    end
    @possible_states = state_lists.resume.map{ _1.inject(:product).map(&:flatten) }.inject(&:concat).flatten
  end

  if format_ == 'split'
    def strip_std_terms
      Fiber.new do
        @possible_states.each do |ps|
          ps = ps.to_a
          if show_unmodified_termini.!
            if ps[0][0] == STD_nterm
              ps[0] = ps[0][1..-1]
            end
            if ps[-1][-1] == STD_cterm
              ps[-1] = ps[-1][-1]
            end
          end
          yield ps
        end
      end
    end
    strip_std_terms.resume
  elsif format_ == 'str'
    @possible_states.map{ tostring(_1, show_unmodified_termini: show_unmodified_termini) }
  else
    raise PyteomicsError("Unsupported value of 'format': #{format_}")
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