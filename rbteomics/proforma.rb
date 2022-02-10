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

  module PrefixSavingMeta

    def __new__(prefix_map)
      new_class = Class.new
      if self.respond_to?("prefix_name")
        prefix_map["prefix_name"] = new_class
      end
      if self.respond_to?("short_prefix")
        prefix_map["short_prefix"] = new_class
      end
      prefix_map
    end

    def find_by_tag(tag_name)
      if tag_name.nil?
        raise ValueError, "tag_name cannot be None!"
      end
      tag_name = tag_name.downcase
      prefix_map = get_prefix_map()
      return prefix_map[tag_name]
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
    include PrefixSavingMeta
    attr_reader :type, :value, :extra, :group_id

    @@prefix_name = nil
    @@short_prefix = nil
    @@prefix_map = {}

    def initialize(...)
      __init__(...)
    end

    def __init__(type, value, extra: nil, group_id: nil)
      @type = type
      @value = value
      @extra = extra
      @group_id = group_id
      @@prefix_map = __new__(@@prefix_map)
    end

    def get_prefix_map
      @@prefix_map
    end

    def to_s
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

    def inspect
      "#{self.class}(#{@value.inspect}, #{@extra.inspect}, #{@group_id.inspect})"
    end

    def ==(other)
      if other.nil?
        return false
      end
      if other.instance_of?(String)
        return self.to_s == other
      end
      (@type == other.class) && (@value == other.value) && (@extra == other.extra) && (@group_id == other.group_id)
    end
  
    def !=(other)
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
      Proforma.process_tag_tokens(buffer)
    end
  end

  class GroupLabelBase < TagBase
    def to_s
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
      "#{@group_id}"
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
    @@prefix_name = 'INFO'

    def initialize(...)
      __init__(...)
    end

    def __init__(value, extra: nil, group_id: nil)
      super(TagTypeEnum.info, value.to_s, extra, group_id)
    end

    def _format_main
      @value.to_s
    end
  end

  class MassModification < TagBase
    attr_reader :_significant_figures

    @@prefix_name = "Obs"
    
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
      super(TagTypeEnum.massmod, value.to_f, extra: extra, group_id: group_id)
    end

    def _format_main
      if @value >= 0
        sprintf("+%#.#{@_significant_figures}f", @value).gsub(/0+$/, '').gsub(/.+$/, '')
      else
        sprintf("%#.#{@_significant_figures}f", @value).gsub(/0+$/, '').gsub(/.+$/, '')
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
    def initialize(...)
      __init__(...)
    end

    def __init__(**kwargs)
      super("unimod", **kwargs)
      @_database = kwargs["database"]
      @strict = kwargs["strict"] || true
    end

    def load_database
      Unimod.new
    end

    def resolve(**kwargs)
      name = kwargs['name'] || nil
      id = kwargs['id'] || nil

      strict = kwargs["strict"] || @strict
      exhaustive = kwargs["exhaustive"] || true
      if name.nil?.!
        defn = @database.by_title(name, strict: strict)
        if defn.!
          defn = @database.by_name(name, strict: strict)
        end
        if defn.! && exhaustive && strict
          defn = @database.by_title(name, strict: false)
          if defn.!
            defn = @database.by_name(name, strict: False)
          end
        end
        if defn && defn.instance_of?(Array)
          warn "Multiple matches found for #{name.inspect} in Unimod, taking the first, #{defn[0]['record_id']}."
          defn = defn[0]
        end
        if defn.!
          raise KeyError(name)
        end
      elsif id.nil?.!
        defn = database.by_id(id)
        if defn.!
          raise KeyError(id)
        end
      else
        raise ValueError("Must provide one of `name` or `id`")
      end
      return {
          'composition' => defn['composition'],
          'name' => defn['title'],
          'id' => defn['record_id'],
          'mass' => defn['mono_mass'],
          'provider' => @name
      }
    end
  end

  class PSIModResolver < ModificationResolver
    def initialize(...)
      __init__(...)
    end

    def __init__(**kwargs)
      super('psimod', **kwargs)
      @_database = kwargs["database"]
    end

    def load_database
      load_psimod()
    end

    def resolve(**kwargs)
      name = kwargs['name'] || nil
      id = kwargs['id'] || nil 
      if name.nil?.!
        defn = @database[name]
      elsif id.nil?.!
        defn = @database[sprintf("MOD:%#.#5d", id)]
      else
        raise ValueError("Must provide one of `name` or `id`")
      end
      mass = defn.DiffMono.to_f
      composition = Composition.new(defn.DiffFormula.strip.replace(" ", ''))
      return {
          'mass' => mass,
          'composition' => composition,
          'name' => defn.name,
          'id' => defn.id,
          'provider' => @name
      }
    end
  end

  class XLMODResolver < ModificationResolver
    def initialize(...)
      __init__(...)
    end

    def __init__(**kwargs)
      super('xlmod', **kwargs)
      @_database = kwargs["database"]
    end

    def load_database
      load_xlmod()
    end

    def resolve(**kwargs)
      name = kwargs['name'] || nil
      id = kwargs['id'] || nil 
      if name.nil?.!
        defn = @database[name]
      elsif id.nil?.!
        defn = @database[sprintf("XLMOD:%#.#5d", id)]
      else
        raise ValueError("Must provide one of `name` or `id`")
      end
      mass = defn['monoIsotopicMass'].to_f
      if defn.include?('deadEndFormula')
        composition = Composition.new(defn['deadEndFormula'].replace(" ", '').replace("D", "H[2]"))
      elsif defn.include?('bridgeFormula')
        composition = Composition.new(defn['bridgeFormula'].replace(" ", '').replace("D", "H[2]"))
      end
      return {
          'mass' => mass,
          'composition' => composition,
          'name' => defn.name,
          'id' => defn.id,
          'provider' => @name
      }
    end
  end

  class GNOResolver < ModificationResolver
    @@mass_pattern = Regexp.compile("(\d+(:?\.\d+)) Da")

    def initialize(...)
      __init__(...)
    end

    def __init__(**kwargs)
      super('gnome', **kwargs)
      @_database = kwargs["database"]
    end

    def load_database
      load_gno()
    end

    def get_mass_from_glycan_composition(term)
      val = term['GNO:00000202']
      monosaccharides = BasicComposition.new()
      composition = Composition.new()
      if val != ''
        tokens = val.scan(/([A-Za-z0-9]+)\((\d+)\)/)
        mass = 0.0
        tokens.each do |symbol, count|
          count = count.to_i
          begin
            mono_mass, mono_comp = GlycanModification.valid_monosaccharides[symbol]
            mass += mono_mass * count
            composition += mono_comp * count
            monosaccharides[symbol] += count        
          rescue => exception
            next if exception.instance_of?(KeyError)
          end
        end
        return [mass, monosaccharides, composition]
      end
      return [nil, nil, nil]
    end
  end

  class GenericResolver < ModificationResolver
    def initialize(...)
      __init__(...)
    end

    def __init__(resolvers, **kwargs)
      super('generic', **kwargs)
      @resolvers = resolvers.to_a
    end

    def load_database
      nil
    end

    def resolve(**kwargs)
      name = kwargs['name'] || nil
      id = kwargs['id'] || nil
      defn = nil
      @resolvers.each do |resolver|
        begin
          defn = resolver('name' => name, 'id' => id, **kwargs)          
        rescue => exception
          next if exception.instance_of?(KeyError)          
        end
      end
      if defn.nil?
        if name.nil?
          raise KeyError(id)
        elsif id.nil?
          raise KeyError(name)
        else
          raise ValueError("Must provide one of `name` or `id`")
        end
      end
      defn
    end
  end

  class ModificationBase < TagBase
    attr_reader :_definition

    @@_tag_type = nil

    def initialize(...)
      __init__(...)
    end

    def __init__(value, extra: nil, group_id: nil)
      super(@@_tag_type, value, extra: extra, group_id: group_id)
      @_definition = nil
    end
  
    #@property
    def definition
      if @_definition.nil?
        @_definition = resolve()
      end
      @_definition
    end

    #@property
    def mass
      @definition['mass']
    end

    #@property
    def composition
      @definition['composition']
    end

    #@property
    def id
      @definition['id']
    end

    #@property
    def name
      @definition['name']
    end

    #@property
    def provider
      @definition['provider']
    end

    def _populate_from_definition(definition)
      @_definition = definition
    end

    def _format_main
      "#{self.prefix_name}:#{self.value}"
    end

    def _parse_identifier
      tokens = value.split(":", 2)
      if tokens.size > 1
        value = tokens[1]
      else
        value = @value
      end
      if value.isdigit()
        id = value.to_i
        name = nil
      else
        name = value
        id = nil
      end
      return [name, id]
    end

    def resolve
      keys = _parse_identifier()
      resolver(*keys)
    end
  end

  class FormulaModification < ModificationBase
    @@prefix_name = "Formula"

    def initialize(...)
      __init__(...)
    end

    def __init__(value,  extra: nil, group_id: nil)
      @isotope_pattern = Regexp.compile('\[(?P<isotope>\d+)(?P<element>[A-Z][a-z]*)(?P<quantity>[\-+]?\d+)\]')
      @_tag_type = TagTypeEnum.formula
      super
    end

    def _normalize_isotope_notation(match)
      parts = match.groupdict()
      "#{parts['element']}[#{parts['isotope']}]#{parts['quantity']}"
    end

    def resolve
      normalized = @value.replace(' ', '')
      if normalized.include?('[')
        normalized = isotope_pattern.sub(@_normalize_isotope_notation, normalized)
      end
      composition = Composition.new('formula' => normalized)
      return {
          "mass" => composition.mass(),
          "composition" => composition,
          "name" => @value
      }
    end
  end

  class GlycanModification < ModificationBase
    @@prefix_name = "Glycan"


    def initialize(...)
      __init__(...)
    end

    def __init__(value,  extra: nil, group_id: nil)
      @_tag_type = TagTypeEnum.glycan
  
      @valid_monosaccharides = {
          "Hex" => [162.0528, Composition.new("C6H10O5")],
          "HexNAc" => [203.0793, Composition.new("C8H13N1O5")],
          "HexS" => [242.009, Composition.new("C6H10O8S1")],
          "HexP" => [242.0191, Composition.new("C6H11O8P1")],
          "HexNAcS" => [283.0361, Composition.new("C8H13N1O8S1")],
          "dHex" => [146.0579, Composition.new("C6H10O4")],
          "NeuAc" => [291.0954, Composition.new("C11H17N1O8")],
          "NeuGc" => [307.0903, Composition.new("C11H17N1O9")],
          "Pen" => [132.0422, Composition.new("C5H8O4")],
          "Pent" => [132.0422, Composition.new("C5H8O4")],
          "Fuc" => [146.0579, Composition.new("C6H10O4")]
      }
  
      @tokenizer = Regexp.compile("([A-Za-z]+)\s*(\d*)\s*")
      @monomer_tokenizer = Regexp.compile(valid_monosaccharides.keys.sort_by{ |x| -x.size }.join('|'))
      super
    end

    #@property
    def monosaccharides
      @definition['monosaccharides']
    end

    def resolve
      composite = BasicComposition.new()
      @value.scan(@tokenizer).each do |tok, cnt|
        if [0, '', nil, false, [], {}].include?(cnt).!
          cnt = cnt.to_i
        else
          cnt = 1
        end
        if @valid_monosaccharides.include?(tok).!
          parts = tok.scan(@monomer_tokenizer)
          t = 0
          parts.each do |ps|
            break if @valid_monosaccharides.include?(ps).!
            t += ps.size
          end
          if t != tok.size
            raise ValueError("#{tok.inspect} is not a valid monosaccharide name")
          else
            parts[0...-1].each do |ps|
              composite[ps] += 1
            end
            composite[parts[-1]] += cnt
          end
        else
          composite[tok] += cnt
        end
      end
      mass = 0
      chemcomp = Composition.new()
      composite.each do |key, cnt|
        m, c = @valid_monosaccharides[key]
        mass += m * cnt
        chemcomp += c * cnt
      end
      return {
          "mass" => mass,
          "composition" => chemcomp,
          "name" => @value,
          "monosaccharides" => composite
      }
    end
  end

  class UnimodModification < ModificationBase
    @@resolver = UnimodResolver.new()

    @@prefix_name = "UNIMOD"
    @@short_prefix = "U"
    @@_tag_type = TagTypeEnum.unimod

    def self.resolver
      @@resolver
    end
  end

  class PSIModModification < ModificationBase
    @@resolver = PSIModResolver.new()

    @@prefix_name = "MOD"
    @@short_prefix = 'M'
    @@_tag_type = TagTypeEnum.psimod

    def self.resolver
      @@resolver
    end
  end

  class GNOmeModification < ModificationBase
    @@resolver = GNOResolver.new()

    @@prefix_name = "GNO"
    @@short_prefix = 'G'
    @@_tag_type = TagTypeEnum.gnome

    def self.resolver
      @@resolver
    end

    #@property
    def monosaccharides
      @definition['monosaccharides']
    end
  end

  class XLMODModification < ModificationBase
    @@resolver = XLMODResolver.new()

    @@prefix_name = "XLMOD"
    # short_prefix = 'XL'
    @@_tag_type = TagTypeEnum.xlmod

    def self.resolver
      @@resolver
    end
  end

  class GenericModification < ModificationBase
    @@_tag_type = TagTypeEnum.generic
    # @@resolver = GenericResolver.new([
    #       partial(UnimodModification.resolver, exhaustive: false),
    #       PSIModModification.resolver,
    #       XLMODModification.resolver,
    #       GNOmeModification.resolver,
    #       partial(UnimodModification.resolver, strict: false)
    # ])

    def self.resolver
      @@resolver
    end

    def _format_main
      @value
    end

    def resolve
      keys = _parse_identifier()
      defn = nil
      begin
        defn = UnimodModification.resolver(*keys)
      rescue => exception
        if exception.instance_of?(KeyError)
          # PASS
        end
      end
      if defn.nil?.!
        return defn
      end
      raise KeyError(keys)
    end
  end

  def split_tags(tokens)
    starts = [0]
    ends = []
    tokens.each_with_index do |c, i|
      if c == '|'
        ends << i
        starts << i + 1
      elsif i != 0 && c == '#'
        ends << i
        starts << i
      end
    end
    ends << tokens.size
    out = []
    starts.each_with_index do |start, i|
      eend = ends[i]
      tag = tokens[start...eend]
      next if tag.size == 0
      out << tag
    end
    out
  end

  def find_prefix(tokens)
    tokens.each_with_index do |c, i|
      if c == ':'
        return [tokens[0...i].join(''), tokens[i + 1..].join('')]
      end
    end
    [nil, tokens.join('')]
  end

  def process_marker(tokens)
    if tokens[1...3] == 'XL'
      return PositionLabelTag.new(nil, group_id: tokens.join(''))
    else
      group_id = None
      value = None
      tokens.each_with_index do |c, i|
        if c == '('
          group_id = tokens[0...i].join('')
          if tokens[-1] != ')'
            raise Exception("Localization marker with score missing closing parenthesis")
          end
          value = tokens[i + 1...-1].join('').to_f
          return LocalizationMarker.new(value, group_id: group_id)
        else
          group_id = ''.join(tokens)
          return PositionLabelTag.new(group_id: group_id)
        end
      end
    end
  end

  def process_tag_tokens(tokens)
    parts = split_tags(tokens)
    main_tag = parts[0]
    if ['+', '-'].include?(main_tag[0])
      main_tag = main_tag.join('')
      main_tag = MassModification.new(main_tag)
    elsif main_tag[0] == '#'
      main_tag = process_marker(main_tag)
    else
      prefix, value = find_prefix(main_tag)
      if prefix.nil?
        main_tag = GenericModification.new(value)
      else
        tag_type = TagBase.new(self.class, value).find_by_tag(prefix)
        main_tag = tag_type.new(value)
      end
    end
    if parts.size > 1
      extras = []
      parts[1..].each do |part|
        prefix, value = find_prefix(part)
        if prefix.nil?
          if value[0] == "#"
            marker = process_marker(value)
            if marker.instance_of?(PositionLabelTag)
              main_tag.group_id = value
            else
              main_tag.group_id = marker.group_id
              extras << marker
            end
          else
            extras << GenericModification.new(value.join(''))
          end
        else
          tag_type = TagBase.new.find_by_tag(prefix)
          extras << tag_type(value)
        end
      end
      main_tag.extra = extras
    end
    main_tag
  end

  class ModificationRule
    attr_reader :modification_tag, :targets

    def initialize(...)
      __init__(...)
    end

    def __init__(modification_tag, targets: nil)
      @modification_tag = modification_tag
      @targets = targets
    end

    def ==(other)
      return false if other.nil?
      @modification_tag == other.modification_tag && @targets == other.targets
    end

    def !=(other)
      self != other
    end

    def to_s
      "<[#{@modification_tag}]@#{@targets.join(',')}>"
    end

    def inspect
      "#{self.class}(#{@modification_tag.inspect}, #{@targets})"
    end
  end

  class StableIsotope
    attr_reader :isotope

    def initialize(...)
      __init__(...)
    end

    def __init__(isotope)
      @isotope = isotope
    end

    def ==(other)
      return false if other.nil?
      @isotope == other.isotope
    end

    def !=(other)
      self != other
    end

    def to_s
      "<#{@isotope}>"
    end

    def inspect
      "#{self.class}(#{@isotope})"
    end
  end

  intersectionenum = Struct.new(
    :no_overlap,
    :full_contains_interval,
    :full_contained_in_interval,
    :start_overlap,
    :end_overlap
  )

  IntersectionEnum = intersectionenum.new(
    no_overlap = 0,
    full_contains_interval = 1,
    full_contained_in_interval = 2,
    start_overlap = 3,
    end_overlap = 4
  )
  
  class TaggedInterval
    attr_accessor :start, :eend, :tags, :ambiguous

    def initialize(...)
      __init__(...)
    end

    def __init__(start, eend: nil, tags: nil, ambiguous: false)
      @start = start
      @eend = eend
      @tags = tags
      @ambiguous = ambiguous
    end

    def __eq__(other)
      if other.nil?
          return false
      end
      @start == other.start && @eend == other.eend && @tags == other.tags
    end

    def __ne__(other)
      self != other
    end

    def __str__
      return "(#{@start}-#{@eend})#{@tags.inspect}"
    end

    def __repr__
      "#{self.class}(#{@start}, #{@eend}, #{@tags})"
    end

    def as_slice
      @start...@eend
    end

    def copy
      __class__(@start, @eend, @tags)
    end

    def _check_slice(qstart, qend, warn_ambiguous)
      valid = qstart <= @start && qend >= @eend
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
        valid = qstart < @eend and qend > @eend
        if valid
          _case = IntersectionEnum.new.end_overlap
          if warn_ambiguous
            warn "Slice bisecting interval #{self}"
          end
        end
      end
  
      if valid
        valid = qstart >= @start && qend < @eend
        if valid
          _case = IntersectionEnum.new.full_contains_interval
          if warn_ambiguous
            warn "Slice bisecting interval #{self}"
          end
        end
      end
      return [valid, _case]
    end

    def _update_coordinates_sliced(start: nil, eend: nil, warn_ambiguous: true)
      if eend.nil?
        qend = @eend + 1
      else
        qend = eend
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
      if eend.nil?.!
        width = [@eend, eend].min - @start
      else
        width = @eend - [start, @start].max
      end
      _new.eend = _new.start + width
      return _new
    end
  end

  class ChargeState
    attr_reader :charge, :adducts

    def initialize(...)
      __init__(...)
    end

    def __init__(charge, adducts: nil)
      if adducts.nil?
        adducts = []
      end
      @charge = charge
      @adducts = adducts
    end

    def to_s
      tokens = [@charge.to_s]
      if adducts.empty?.!
        tokens << "["
        tokens << @adducts.map{ |adduct| adduct.to_s }.join(',')
        tokens << "]"
      end
      tokens.join('')
    end

    def inspect
      "#{self.class}(#{@charge}, #{@adducts})"
    end
  end

  class TokenBuffer
    def initialize(...)
      __init__(...)
    end

    def __init__(initial: nil)
      @buffer = initial.nil? ? [] : initial.to_a
      @boundaries = []
    end

    def <<(c)
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

    def bool
      [0, '', nil, false, [], {}].include?(@buffer).!
    end

    def __iter__
      @buffer.to_enum
    end

    def [](i)
      @buffer[i]
    end

    def size
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
      if [0, '', nil, false, [], {}].include?(@boundaries).!
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
    def _transform(value)
      value.join('')
    end
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
      tag = Proforma.process_tag_tokens(value)
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
      :charge_state_adducteend,
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
    charge_state_adducteend = 19,
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
  ADDUCTeend = ParserStateEnum.charge_state_adducteend
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
            positions << [current_aa, current_tag.bool ? current_tag.process : nil]
          end
          current_aa = c
        elsif c == '['
          state = TAG
          if current_tag.bool
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
          positions << [current_aa, current_tag.bool ? current_tag.process : nil]
          current_aa = nil
          if current_interval.nil?
            raise ProFormaError.new("Error In State #{state}, unexpected #{c} found at index #{i}")
          else
            current_interval.eend = positions.size
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
          end
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
                c_term = current_tag.process
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
            isotopes << StableIsotope.new(current_tag.process.join(''))
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
              labile_modifications << current_tag.process[0]
              state = BEFORE
            end
          else
            current_tag << c
          end
      elsif state == POST_INTERVAL_TAG
          if c == '['
            current_tag.bound()
            state = INTERVAL_TAG
          elsif VALID_AA.include?(c)
            current_aa = c
            current_interval.tags = current_tag.process
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
            unlocalized_modifications << current_tag.process[0]
            state = BEFORE
          elsif c == '-'
            n_term = current_tag.process
            state = BEFORE
          elsif c == '^'
            state = UNLOCALIZED_COUNT
          else
            raise ProFormaError.new("Error In State #{state}, unexpected #{c} found at index #{i}")
          end
      elsif state == UNLOCALIZED_COUNT
          if c.match(/[0-9]/)
            current_unlocalized_count << c
          elsif c == '['
              state = TAG_BEFORE
              depth = 1
              tag = current_tag.process[0]
              multiplicity = current_unlocalized_count()
              multiplicity.times do |i|
                unlocalized_modifications << tag
              end
          elsif c == '?'
              state = BEFORE
              tag = current_tag.process[0]
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
            fixed_modifications << ModificationRule.new(current_tag.process[0], targets: current_aa_targets.process)
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
          elsif c.match(/[0-9]/)
              charge_buffer << c
              state = CHARGE_NUMBER
          elsif c == '/'
              state = ParserStateEnum.inter_chain_cross_link_start
              raise ProFormaError.new("Inter-chain cross-linked peptides are not yet supported")
          else
              raise ProFormaError.new("Error In State #{state}, unexpected #{c} found at index #{i}")
          end
      elsif state == CHARGE_NUMBER
          if c.match(/[0-9]/)
              charge_buffer << c
          elsif c == "["
              state = ADDUCT_START
              adduct_buffer = StringParser()
          else
              raise ProFormaError.new("Error In State #{state}, unexpected #{c} found at index #{i}")
          end
      elsif state == ADDUCT_START
          if c.match(/[0-9]/) || "+-".include?(c) || element_symbols.include?(c)
            adduct_buffer << c
          elsif c == ','
            adduct_buffer.bound()
          elsif c == ']'
            state = ADDUCTeend
          end
      elsif state == ADDUCTeend
          if c == '+'
            raise ProFormaError.new("Error In State #{state}, #{c} found at index #{i}. Chimeric representation not supported")
          end
      else
        raise ProFormaError.new("Error In State #{state}, unexpected #{c} found at index #{i}")
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
      positions << [current_aa, current_tag.bool ? current_tag.process : nil]
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
      'charge_state' => charge_state
    }]
  end

  def to_proforma(sequence, n_term: nil, c_term: nil, unlocalized_modifications: nil,
    labile_modifications: nil, fixed_modifications: nil, intervals: nil,
    isotopes: nil, charge_state: nil, group_ids: nil)
    primary = []
    sequence.each do |aa, tags|
      if tags.nil?
        primary << aa.to_s
      else
        primary << aa.to_s.concat(tags.map{ |t| "[#{t.to_s}]" }.join(''))
      end
    end
    if intervals.nil?.!
      intervals.sort_by{ |x| x.start }.each do |iv|
        if iv.ambiguous
          primary[iv.start] = '(?' + primary[iv.start]
        else
          primary[iv.start] = '(' + primary[iv.start]
        end

        terminator = '{0!s})'.format(primary[iv.end - 1])
        if iv.tags
          terminator += iv.tags.map{ |t| "[#{t.to_s}]"}.join('')
        end
        primary[iv.eend - 1] = terminator
      end
    end
    if n_term
      primary.unshift n_term.map{ |t| "[#{t.to_s}]" }.join('').concat('-')
    end
    if c_term
      primary << '-'.concat(c_term.map{ |t| "[#{t.to_s}" }.join(''))
    end
    if charge_state
      primary << "/#{charge_state.to_s}"
    end
    if labile_modifications
      primary.unshift labile_modifications.map{ |m| "{{#{m.to_s}}}" }
    end
    if unlocalized_modifications
      primary.unshift "?"
      primary.unshift unlocalized_modifications.map{ |m| "#{m.to_s}" }
    end
    if isotopes
      primary.unshift isotopes.map{ |m| "#{m.to_s}" }
    end
    if fixed_modifications
      primary.unshift fixed_modifications.map{ |m| "#{m.to_s}" }
    end
    primary.join('')
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

    def inspect
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
    end

    @@isotopes = ProFormaProperty.new('isotopes')
    @@charge_state = ProFormaProperty.new('charge_state')

    @@intervals = ProFormaProperty.new('intervals')
    @@fixed_modifications = ProFormaProperty.new('fixed_modifications')
    @@labile_modifications = ProFormaProperty.new('labile_modifications')
    @@unlocalized_modifications = ProFormaProperty.new('unlocalized_modifications')

    @@n_term = ProFormaProperty.new('n_term')
    @@c_term = ProFormaProperty.new('c_term')

    @@group_ids = ProFormaProperty.new('group_ids')  

    def to_s
      to_proforma(@sequence, **@properties)
    end

    def inspect
      "#{self.class}(#{@sequence}, #{@properties})"
    end

    def [](i)
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

    def ==(other)
      if other.instance_of?(String)
        self.to_s == other
      elsif other.nil?
        false
      else
        @sequence == other.sequence && @properties == other.properties
      end
    end

    def !=(other)
      self != other
    end

    def self.parse(string)
      Proforma.parse(string)
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
