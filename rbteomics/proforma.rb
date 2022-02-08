require_relative 'mass/mass'
require_relative 'auxiliary/file_helpers'
require_relative 'auxiliary/structures'
require_relative '_schema_defaults'
require 'set'

# from pyteomics.mass.mass import calculate_mass
# import re
# import warnings
# from collections import deque
# from functools import partial

# try:
#     from enum import Enum
# except ImportError:
#     # Python 2 doesn't have a builtin Enum type
#     Enum = object

# from pyteomics.mass import Composition, std_aa_mass, Unimod, nist_mass
# from pyteomics.auxiliary import PyteomicsError, BasicComposition
# from pyteomics.auxiliary.utils import add_metaclass

# try:
#     from psims.controlled_vocabulary.controlled_vocabulary import (load_psimod, load_xlmod, load_gno, obo_cache)
# except ImportError:
#     def _needs_psims(name):
#         raise ImportError("Loading %s requires the `psims` library. To access it, please install `psims`" % name)

#     load_psimod = partial(_needs_psims, 'PSIMOD')
#     load_xlmod = partial(_needs_psims, 'XLMOD')
#     load_gno = partial(_needs_psims, 'GNO')
#     obo_cache = None

module Proforma
  module_function

  Water_mass = calculate_mass('H2O')
  STD_aa_mass = STD_aa_mass.dup
  STD_aa_mass['X'] = 0

  element_symbols = Set.new(Nist_mass)
  element_symbols.delete('e*')
  element_symbols << 'e'

  class ProFormaError < PyteomicsError
    def initialize(...)
      __init__(...)
    end

    def __init__(message, **kwargs)
      index = kwargs['index'] || nil
      parser_state = kwargs['parser_state'] || nil
      super(message, index, parser_state)
      @message = message
      @index = index
      @parser_state = parser_state
    end
  end

  class PrefixSavingMeta
    def initialize(...)
      __new__(...)
    end

    def __new__(mcs, name, parents, attrs)
      new_type = type.__new__(mcs, name, parents, attrs)
      prefix = attrs["prefix_name"]
      if prefix
        new_type.prefix_map[prefix.downcase] = new_type
      end
      short = attrs["short_prefix"]
      if short
        new_type.prefix_map[short.downcase] = new_type
      end
      new_type
    end

    def find_by_tag(tag_name)
      if tag_name.nil?
        raise ValueError, "tag_name cannot be None!"
      end
      tag_name = tag_name.lower()
      return @prefix_map[tag_name]
    end
  end

  tagtypeenum = Struct.new(
    :unimod,
    :psimod,
    :massmod,
    :generic,
    :info,
    :gnome,
    :xlmod,
    :formula,
    :glycan,
    :localization_marker,
    :position_label,
    :group_placeholder
  )
  TagTypeEnum = tagtypeenum.new(
    unimod = 0,
    psimod = 1,
    massmod = 2,
    generic = 3,
    info = 4,
    gnome = 5,
    xlmod = 6,
    formula = 7,
    glycan = 8,
    localization_marker = 9,
    position_label = 10,
    group_placeholder = 999
  )

  Sentinel = Object.new

  #@add_metaclass(PrefixSavingMeta)
  class TagBase
    attr_reader :type, :value, :extra, :group_id

    def initialize(...)
      __init__(...)
    end

    def __init__(type, value, extra: nil, group_id: nil)
      @prefix_name = nil
      @short_prefix = nil
      @prefix_map = {}

      @type = type
      @value = value
      @extra = extra
      @group_id = group_id
    end

    def __str__
      part = self._format_main()
      if @extra
        rest = @extra.map { |e| e.to_s }
        label = ([part] + rest).join('|')
      else
        label = part
      end
      if @group_id
        label = "#{label}#{@group_id}"
      end
      return label.to_s
    end

    def __repr__
      "#{self.class}(#{@value.inspect}, #{@extra.inspect}, #{@group_id.inspect})"
    end

    def __eq__(other)
      if other.nil?
        return false
      end
      if other.instance_of?(String)
        return self.to_s == other
      end
      (@type == other.class) && (@value == other.value) && (@extra == other.extra) && (@group_id == other.group_id)
    end
  
    def __ne__(other)
      self != other
    end

    def find_tag_type(tag_type)
      out = []
      if @type == tag_type
        out << self
      end
      if [0, '', nil, false, [], {}].include?(@extra)
        return out
      end
      @extra.each do |e|
        if e.type == tag_type
          out << e
        end
      end
      out
    end

    def self.parse(buffer)
      process_tag_tokens(buffer)
    end
  end

  class GroupLabelBase < TagBase
    def __str__
      part = _format_main
      if @extra.empty?.!
        rest = @extra.map{ |e| e.to_s }
        label = ([part] + rest).join('|')
      else
        label = part
      end
      return label.to_s
    end
  end
  
  class PositionLabelTag < GroupLabelBase
    def initialize(...)
      __init__(...)
    end

    def __init__(value: nil, extra: nil, group_id: nil)
      assert group_id.nil?.!
      value = group_id
      super(TagTypeEnum.position_label, value, extra, group_id)
    end

    def _format_main
      "#{self.group_id}"
    end
  end

  class LocalizationMarker < GroupLabelBase
    def initialize(...)
      __init__(...)
    end

    def __init__(value, extra: nil, group_id: nil)
      assert group_id.nil?.!
      super(TagTypeEnum.localization_marker, value.to_f, extra, group_id)
    end

    def _format_main
      # "#{@group_id}(#{@value:.4g})"
      "#{@group_id}(#{@value})"
    end
  end

  class InformationTag < TagBase
    def initialize(...)
      __init__(...)
    end

    def __init__(value, extra: nil, group_id: nil)
      @prefix_name = 'INFO'
      super(TagTypeEnum.info, value.to_s, extra, group_id)
    end

    def _format_main
      @value.to_s
    end
  end

  class MassModification < TagBase
    attr_reader :_significant_figures

    def initialize(...)
      __init__(...)
    end

    def __init__(value, extra: nil, group_id: nil)
      if value.instance_of?(String)
        sigfigs = value.split('.')[-1].gsub(/0+$/, '').size
      else
        sigfigs = 4
      end
      @_significant_figures = sigfigs
      super(TagTypeEnum.massmod, value.to_f, extra, group_id)
    end

    def _format_main
      if @value >= 0:
        sprintf("+%#.#{@_significant_figures}", @value).gsub(/0+$/, '').gsub(/.+$/, '')
      else
        sprintf("%#.#{@_significant_figures}", @value).gsub(/0+$/, '').gsub(/.+$/, '')
      end
    end

    #@property
    def mass
      @value
    end
  end

  class ModificationResolver
    def initialize(...)
      __init__(...)
    end

    def __init__(name, **kwargs)
      @name = name
      @_database = nil
    end

    def load_database
      raize NtoImplementedError
    end

    #@property
    def database
      if @_database.!
        @_database = load_database()
      end
      @_database
    end

    def resolve(**kwargs)
      name = kwargs['name'] || nil
      id = kwargs['id'] || nil 
      raise NotImplementedError
    end

    def __call__(**kwargs)
      name = kwargs['name'] || nil
      id = kwargs['id'] || nil
      resolve(**kwargs)
    end
  end
  class UnimodResolver < ModificationResolver
  
  
  
  end
  class PSIModResolver < ModificationResolver;end
  class XLMODResolver < ModificationResolver;end
  class GNOResolver < ModificationResolver;end
  class GenericResolver < ModificationResolver;end
  class ModificationBase < TagBase;end
  class FormulaModification < ModificationBase;end
  class GlycanModification < ModificationBase;end
  class UnimodModification < ModificationBase;end
  class PSIModModification < ModificationBase;end
  class GNOmeModification < ModificationBase;end
  class XLMODModification < ModificationBase;end
  class GenericModification < ModificationBase;end
  def split_tags(tokens);end
  def find_prefix(tokens);end
  def process_marker(tokens);end
  def process_tag_tokens(tokens);end
  class ModificationRule
      # Not started

  end
  class StableIsotope
      # Not started

  end
  class IntersectionEnum < Enumerator
      # Not started
  end

  class TaggedInterval
    attr_reader :start, :_end, :tags, :ambiguous

    def initialize(...)
      __init__(...)
    end

    def __init__(start, _end: nil, tags: nil, ambiguous: false)
      @start = start
      @_end = _end
      @tags = tags
      @ambiguous = ambiguous
    end

    def __eq__(other)
      if other.nil?
          return false
      end
      @start == other.start && @_end == other._end && @tags == other.tags
    end

    def __ne__(other)
      self != other
    end

    def __str__
      return "(#{@start}-#{@_end})#{@tags.inspect}"
    end

    def __repr__
      "#{self.class}(#{@start}, #{@_end}, #{@tags})"
    end

    def as_slice
      @start...@_end
    end

    def copy
      __class__(@start, @_end, @tags)
    end

    def _check_slice(qstart, qend, warn_ambiguous)
      valid = qstart <= @start && qend >= @_end
      _case = valid ? IntersectionEnum.new.full_contained_in_interval : IntersectionEnum.new.no_overlap
      if valid.!
        valid = qstart <= @start && qend > @start
        if valid
          _case = IntersectionEnum.new.start_overlap
          if warn_ambiguous
            warn "Slice bisecting interval #{self}"
          end
        end
      end
  
      if valid.!
        valid = qstart < @_end and qend > @_end
        if valid
          _case = IntersectionEnum.new.end_overlap
          if warn_ambiguous
            warn "Slice bisecting interval #{self}"
          end
        end
      end
  
      if valid
        valid = qstart >= @start && qend < @_end
        if valid
          _case = IntersectionEnum.new.full_contains_interval
          if warn_ambiguous
            warn "Slice bisecting interval #{self}"
          end
        end
      end
      return [valid, _case]
    end

    def _update_coordinates_sliced(start: nil, _end: nil, warn_ambiguous: true)
      if _end.nil?
        qend = @_end + 1
      else
        qend = _end
      end
      if start.nil?
        qstart = @start - 1
      else
        qstart = @start
      end

      valid, intersection_type = _check_slice(qstart, qend, warn_ambiguous)
      if @ambiguous && [IntersectionEnum.full_contained_in_interval, IntersectionEnum.no_overlap].include?(intersection_type).!
        raise ValueError("Cannot bisect an ambiguous interval")
      end
      return nil if valid.!
      _new = self.dup
      if start.nil?.!
        diff = @start - start
        if diff < 0
          diff = 0
        end
        @start = diff
      end
      if _end.nil?.!
        width = [@_end, _end].min - @start
      else
        width = @_end - [start, @start].max
      end
      _new._end = _new.start + width
      return _new
    end
  end

  class ChargeState
      # Not started
  end
  class TokenBuffer
    def initialize(...)
      __init__(...)
    end

    def __init__(initial: nil)
      @buffer = initial.nil? ? [] : initial.to_a
      @boundaries = []
    end

    def append(c)
      @buffer << c
    end

    def reset
      if @buffer.empty?.!
        @buffer = []
      end
      if @boundaries.empty?.!
        @boundaries = []
      end
    end

    def __bool__
      [0, '', nil, false, [], {}].include?(@buffer).!
    end

    def __iter__
      @buffer.to_enum
    end

    def __getitem__(i)
      @buffer[i]
    end

    def __len__
      @buffer.size
    end

    def tokenize
      i = 0
      pieces = []
      (@boundaries + [self.size]).each do |k|
        piece = @buffer[i...k]
        i = k
        pieces << piece
      end
      pieces
    end

    def _transform(value)
      value
    end

    def process
      if @boundaries
        value = tokenize.map{ |v| _transform(v) }
      else
        value = _transform(@buffer)
      end
      reset()
      value
    end

    def bound
      k = self.size
      @boundaries << k
      k
    end

    def __call__
      process
    end
  end

  class NumberParser < TokenBuffer
    def _transform(value)
      value.join('').to_i
    end
  end
  class StringParser < TokenBuffer
    # Not started
  end

  class TagParser < TokenBuffer
    def initialize(...)
      __init__(...)
    end

    def __init__(initial: nil, group_ids: nil)
      super(initial: initial)
      if group_ids.nil?.!
        @group_ids = group_ids.to_set
      else
        @group_ids = Set.new
      end
    end

    def _transform(value)
      tag = process_tag_tokens(value)
      if tag.group_id
        @group_ids << tag.group_id
      end
      return tag
    end

    def process
      value = super
      if value.instance_of?(Array).!
        value = [value]        
      end
      value
    end
  end

  parserstateenum = Struct.new(
      :before_sequence,
      :tag_before_sequence,
      :global_tag,
      :fixed_spec,
      :labile_tag,
      :sequence,
      :tag_in_sequence,
      :interval_tag,
      :tag_after_sequence,
      :stable_isotope,
      :post_tag_before,
      :unlocalized_count,
      :post_global,
      :post_global_aa,
      :post_interval_tag,
      :post_tag_after,
      :charge_state_start,
      :charge_state_number,
      :charge_state_adduct_start,
      :charge_state_adduct_end,
      :inter_chain_cross_link_start,
      :chimeric_start,
      :interval_initial,
      :done
  )
  ParserStateEnum = parserstateenum.new(
    before_sequence = 0,
    tag_before_sequence = 1,
    global_tag = 2,
    fixed_spec = 3,
    labile_tag = 4,
    sequence = 5,
    tag_in_sequence = 6,
    interval_tag = 7,
    tag_after_sequence = 8,
    stable_isotope = 9,
    post_tag_before = 10,
    unlocalized_count = 11,
    post_global = 12,
    post_global_aa = 13,
    post_interval_tag = 14,
    post_tag_after = 15,
    charge_state_start = 16,
    charge_state_number = 17,
    charge_state_adduct_start = 18,
    charge_state_adduct_end = 19,
    inter_chain_cross_link_start = 20,
    chimeric_start = 21,
    interval_initial = 22,
    done = 999
  )

  BEFORE = ParserStateEnum.before_sequence
  TAG_BEFORE = ParserStateEnum.tag_before_sequence
  FIXED = ParserStateEnum.fixed_spec
  GLOBAL = ParserStateEnum.global_tag
  ISOTOPE = ParserStateEnum.stable_isotope
  LABILE = ParserStateEnum.labile_tag
  SEQ = ParserStateEnum.sequence
  TAG = ParserStateEnum.tag_in_sequence
  INTERVAL_TAG = ParserStateEnum.interval_tag
  INTERVAL_INIT = ParserStateEnum.interval_initial
  TAG_AFTER = ParserStateEnum.tag_after_sequence
  POST_TAG_BEFORE = ParserStateEnum.post_tag_before
  POST_TAG_AFTER = ParserStateEnum.post_tag_after
  UNLOCALIZED_COUNT = ParserStateEnum.unlocalized_count
  POST_GLOBAL = ParserStateEnum.post_global
  POST_GLOBAL_AA = ParserStateEnum.post_global_aa
  POST_INTERVAL_TAG = ParserStateEnum.post_interval_tag
  CHARGE_START = ParserStateEnum.charge_state_start
  CHARGE_NUMBER = ParserStateEnum.charge_state_number
  ADDUCT_START = ParserStateEnum.charge_state_adduct_start
  ADDUCT_END = ParserStateEnum.charge_state_adduct_end
  DONE = ParserStateEnum.done
  
  VALID_AA = 'QWERTYIPASDFGHKLCVNMXUOJZB'.split('').to_set
  
  def parse(sequence)
    labile_modifications = []
    fixed_modifications = []
    unlocalized_modifications = []
    intervals = []
    isotopes = []

    n_term = nil
    c_term = nil

    i = 0
    n = sequence.size

    positions = []
    state = BEFORE
    depth = 0

    current_aa = nil
    current_tag = TagParser.new()
    current_interval = nil
    current_unlocalized_count = NumberParser.new()
    current_aa_targets = TokenBuffer.new()

    charge_buffer = nil
    adduct_buffer = nil

    while i < n
      c = sequence[i]
      i += 1
      if state == BEFORE
        if c == '['
          state = TAG_BEFORE
          depth = 1
        elsif c == '{'
          state = LABILE
          depth = 1
        elsif c == '<'
          state = FIXED
        elsif VALID_AA.include?(c)
          current_aa = c
          state = SEQ
        else
          raise ProFormaError.new("Error In State #{state}, unexpected #{c} found at index #{i}")
        end
      elsif state == SEQ || state == INTERVAL_INIT
        if state == INTERVAL_INIT
          state = SEQ
          if c == '?'
            if current_interval.nil?.!
              current_interval.ambiguous = true
            end
            next
          end
        end
        if VALID_AA.include?(c)
          if current_aa.nil?.!
            positions << [[current_aa, current_tag  ? current_tag() : nil]]
          end
          current_aa = c
        elsif c == '['
          state = TAG
          if current_tag
            current_tag.bound()
          end
          depth = 1
        elsif c == '('
          if current_interval.nil?.!
            raise ProFormaError.new("Error In State #{state}, nested range found at index #{i}. Nested ranges are not yet supported by ProForma.")
          end
          current_interval = TaggedInterval.new(positions.size + 1)
          state = INTERVAL_INIT
        elsif c == ')'
          positions << [[current_aa, current_tag ? current_tag() : None]]
          current_aa = nil
          if current_interval.nil?
            raise ProFormaError.new("Error In State #{state}, unexpected #{c} found at index #{i}")
          else
            current_interval._end = positions.size
            if i < n && sequence[i] == '['
              i += 1
              depth = 1
              state = INTERVAL_TAG
            else
              intervals << current_interval
              current_interval = nil
            end
          end
        elsif c == '-'
          state = TAG_AFTER
          if i >= n || sequence[i] != '['
            raise ProFormaError.new("Missing Closing Tag")
            i += 1
            depth = 1
          elsif c == '/'
            state = CHARGE_START
            charge_buffer = NumberParser.new()
          elsif c == '+'
            raise ProFormaError.new("Error In State #{state}, #{c} found at index #{i}. Chimeric representation not supported")
          else
            raise ProFormaError.new("Error In State #{state}, unexpected #{c} found at index #{i}")
          end
        elsif state == TAG || state == TAG_BEFORE || state == TAG_AFTER || state == GLOBAL || state == INTERVAL_TAG
          if c == '['
            depth += 1
            current_tag << c
          elsif c == ']'
            depth -= 1
            if depth <= 0
              depth = 0
              if state == TAG
                state = SEQ
              elsif state == TAG_BEFORE
                state = POST_TAG_BEFORE
              elsif state == TAG_AFTER
                c_term = current_tag()
                state = POST_TAG_AFTER
              elsif state == GLOBAL
                state = POST_GLOBAL
              elsif state == INTERVAL_TAG
                state = POST_INTERVAL_TAG
                depth = 0
              end
            else
              current_tag << c
            end
          else
            current_tag << c
          end
        elsif state == FIXED
          if c == '['
            state = GLOBAL
          else
            state = ISOTOPE
            current_tag.reset()
            current_tag << c
          end
        elsif state == ISOTOPE
          if c != '>'
            current_tag << c
          else
            isotopes << StableIsotope.new(current_tag.join(''))
            current_tag.reset()
            state = BEFORE
          end
        elsif state == LABILE
          if c == '{'
            depth += 1
          elsif c == '}'
            depth -= 1
            if depth <= 0
              depth = 0
              labile_modifications << current_tag()[0]
              state = BEFORE
            end
          else
            current_tag.append(c)
          end
        elsif state == POST_INTERVAL_TAG
          if c == '['
            current_tag.bound()
            state = INTERVAL_TAG
          elsif VALID_AA.include?(c)
            current_aa = c
            current_interval.tags = current_tag()
            intervals << current_interval
            current_interval = nil
            state = SEQ
          elsif c == '-'
            state = TAG_AFTER
            if i >= n || sequence[i] != '['
              raise ProFormaError.new("Missing Closing Tag")
            end
            i += 1
            depth = 1
          elsif c == '/'
            state = CHARGE_START
            charge_buffer = NumberParser()
          elsif c == '+'
            raise ProFormaError.new("Error In State #{state}, #{c} found at index #{i}. Chimeric representation not supported")
          else
            raise ProFormaError.new("Error In State #{state}, unexpected #{c} found at index #{i}")
          end
        elsif state == POST_TAG_BEFORE
          if c == '?'
            unlocalized_modifications << current_tag()[0]
            state = BEFORE
          elsif c == '-'
            n_term = current_tag()
            state = BEFORE
          elsif c == '^'
            state = UNLOCALIZED_COUNT
          else
            raise ProFormaError.new("Error In State #{state}, unexpected #{c} found at index #{i}")
          end
        elsif state == UNLOCALIZED_COUNT
          if c.isdigit()
            current_unlocalized_count << c
          elsif c == '['
              state = TAG_BEFORE
              depth = 1
              tag = current_tag()[0]
              multiplicity = current_unlocalized_count()
              multiplicity.times do |i|
                unlocalized_modifications << tag
              end
          elsif c == '?'
              state = BEFORE
              tag = current_tag()[0]
              multiplicity = current_unlocalized_count()
              multiplicity.times do |i|
                unlocalized_modifications << tag
              end
          else
              raise ProFormaError.new("Error In State #{state}, unexpected #{c} found at index #{i}")
          end
        elsif state == POST_GLOBAL
          if c == '@'
            state = POST_GLOBAL_AA
          else
            raise ProFormaError.new("Error In State #{state}, fixed modification detected without target amino acids found at index #{i}")
          end
        elsif state == POST_GLOBAL_AA
          if VALID_AA.include?(c)
            current_aa_targets << c
          elsif c == ','
            # pass
          elsif c == '>'
            fixed_modifications << ModificationRule.new(current_tag()[0], current_aa_targets())
              state = BEFORE
          else
            raise ProFormaError.new("Error In State #{state}, unclosed fixed modification rule")
          end
        elsif state == POST_TAG_AFTER
          if c == '/'
              state = CHARGE_START
              charge_buffer = NumberParser.new()
          elsif c == '+'
            raise ProFormaError.new("Error In State #{state}, #{c} found at index #{i}. Chimeric representation not supported")
          end
        elsif state == CHARGE_START
          if '+-'.include?(c)
              charge_buffer << c
              state = CHARGE_NUMBER
          elsif c.isdigit()
              charge_buffer << c
              state = CHARGE_NUMBER
          elsif c == '/'
              state = ParserStateEnum.inter_chain_cross_link_start
              raise ProFormaError.new("Inter-chain cross-linked peptides are not yet supported")
          else
              raise ProFormaError.new("Error In State #{state}, unexpected #{c} found at index #{i}")
          end
        elsif state == CHARGE_NUMBER
          if c.isdigit()
              charge_buffer << c
          elsif c == "["
              state = ADDUCT_START
              adduct_buffer = StringParser()
          else
              raise ProFormaError.new("Error In State #{state}, unexpected #{c} found at index #{i}")
          end
        elsif state == ADDUCT_START
          if c.isdigit() || "+-".include?(c) || element_symbols.include?(c)
            adduct_buffer << c
          elsif c == ','
            adduct_buffer.bound()
          elsif c == ']'
            state = ADDUCT_END
          end
        elsif state == ADDUCT_END
          if c == '+'
              raise ProFormaError.new("Error In State #{state}, #{c} found at index #{i}. Chimeric representation not supported")
          end
        else
          raise ProFormaError.new("Error In State #{state}, unexpected #{c} found at index #{i}")
        end
      end
    end
    if charge_buffer
      charge_number = charge_buffer()
      if adduct_buffer
        adducts = adduct_buffer()
      else
        adducts = nil
      end
      charge_state = ChargeState.new(charge_number, adducts)
    else
      charge_state = nil
    end
    if current_aa
      positions << [current_aa, current_tag ? current_tag() : nil]
    end
    if [ISOTOPE, TAG, TAG_AFTER, TAG_BEFORE, LABILE].include?(state)
      raise ProFormaError.new("Error In State #{state}, unclosed group reached end of string!")
    end
    return [positions, {
      'n_term' => n_term,
      'c_term' => c_term,
      'unlocalized_modifications' => unlocalized_modifications,
      'labile_modifications' => labile_modifications,
      'fixed_modifications' => fixed_modifications,
      'intervals' => intervals,
      'isotopes' => isotopes,
      'group_ids' => lambda { |x| current_tag.group_ids[x] },
      'charge_state' => charge_state,
    }]
  end

  def to_proforma(sequence, n_term: nil, c_term: nil, unlocalized_modifications: nil,
    labile_modifications: nil, fixed_modifications: nil, intervals: nil,
    isotopes: nil, charge_state: nil, group_ids: nil)
  end

  class ProFormaProperty
    def initialize(...)
      __init__(...)
    end

    def __init__(name)
      @name = name
    end

    def __get__(obj, cls)
      obj.properties[@name]
    end

    def __set__(obj, value)
      obj.properties[@name] = value
    end

    def __repr__
      "#{self.class}(#{@name.inspect})"
    end
  end

  class ProForma
    def initialize(...)
      __init__(...)
    end

    def __init__(sequence, properties)
      @sequence = sequence
      @properties = properties

      @isotopes = ProFormaProperty.new('isotopes')
      @charge_state = ProFormaProperty.new('charge_state')
  
      @intervals = ProFormaProperty.new('intervals')
      @fixed_modifications = ProFormaProperty.new('fixed_modifications')
      @labile_modifications = ProFormaProperty.new('labile_modifications')
      @unlocalized_modifications = ProFormaProperty.new('unlocalized_modifications')
  
      @n_term = ProFormaProperty.new('n_term')
      @c_term = ProFormaProperty.new('c_term')
  
      @group_ids = _ProFormaProperty('group_ids')  
    end

    def __str__
      to_proforma(@sequence, **@properties)
    end

    def __repr__
      "#{self.class}(#{@sequence}, #{@properties})"
    end

    def __getitem__(i)
      if i.instance_of?(Range)
        props = @properties.dup
        ivs = []
        props['intervals'].each do |iv|
          iv = iv._update_coordinates_sliced(i.start, i.stop)
          next if iv.nil?
          ivs << iv
        end
        props['intervals'] = ivs
        __class__(@sequence[i], props)
      else
        @sequence[i]
      end
    end

    def __eq__(other)
      if other.instance_of?(String)
        self.to_s == other
      elsif other.nil?
        false
      else
        @sequence == other.sequence && @properties == other.properties
      end
    end

    def __ne__(other)
      self != other
    end

    def self.parse(string)
      self.new(*parse(string))
    end

    #@property
    def mass
      mass = 0.0

      fixed_modifications = @properties['fixed_modifications']
      fixed_rules = {}
      fixed_modifications.each do |rule|
        rule.targets.each do |aa|
          fixed_rules[aa] = rule.modification_tag.mass
        end
      end

      @sequence.each do |position|
        aa = position[0]
        begin
          mass += std_aa_mass[aa]          
        rescue => exception
          if exception.instance_of?(KeyError)
            warn "#{aa} does not have an exact mass"
          end
        end
        if fixed_rules.include?(aa)
          mass += fixed_rules[aa]
        end
        tags = position[1]
        if tags.empty?.!
          tags.each do |tag|
            begin
              mass += tag.mass
            rescue => exception
              next if [AttributeError, KeyError].include?(exception.class)
            end
          end
        end
      end
      @properties['labile_modifications'].each do |mod|
        mass += mod.mass
      end
      @properties['unlocalized_modifications'].each do |mod|
        mass += mod.mass
      end
      if @properties['n_term']
        @properties['n_term'].each do |mod|
          begin
            mass += mod.mass            
          rescue => exception
            next if [AttributeError, KeyError].include?(exception.class)
          end
        end        
      end
      mass += calculate_mass('formula' => "H")
      if @properties['c_term']
        @properties['c_term'].each do |mod|
          begin
            mass += mod.mass            
          rescue => exception
            next if [AttributeError, KeyError].include?(exception.class)            
          end
        end
      end
      mass += calculate_mass('formula' => "OH")
      @properties['intervals'].each do |iv|
        begin
          mass += iv.tag.mass
        rescue => exception
          next if [AttributeError, KeyError].include?(exception.class)            
        end
      end
      mass
    end

    def find_tags_by_id(tag_id, include_position: true)
      if tag_id.start_with('#').!
        tag_id = '#' + tag_id
      end
      matches = []
      @sequence.to_enum.each do |i, j|
        _token, tags = j
        if tags.empty?.!
          tags.each do |tag|
            if tag.group_id == tag_id
              if include_position
                  matches << [i, tag]
              else
                  matches << tag
              end
            end
          end
        end
      end
      @properties['intervals'].each do |iv|
        if iv.tag.group_id == tag_id
          matches << include_position ? [iv, iv.tag] : iv.tag
        end
      end
      @properties['unlocalized_modifications'].each do |ulmod|
        if ulmod.group_id == tag_id
          matches << include_position ? ['unlocalized_modifications', ulmod] : ulmod
        end
      end
      @properties['labile_modifications'].each do |lamod|
        if lamod.group_id == tag_id
          matches << include_position ? ['labile_modifications', lamod] : lamod
        end
      end
      matches
    end
  end
end
