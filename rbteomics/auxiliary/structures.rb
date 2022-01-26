# import re
# from collections import defaultdict, Counter
# import warnings

# try:
#     basestring
#     PY2 = True
# except NameError:
#     basestring = (str, bytes)
#     PY2 = False
BaseString ||= String
PY2 ||= false

# _UNIT_CV_INTERN_TABLE = dict()
UNIT_CV_INTERN_TABLE ||= {}

def clear_unit_cv_table()
  UNIT_CV_INTERN_TABLE.clear
end

def _intern_unit_or_cv(unit_or_cv)
  return nil if unit_or_cv.nil?
  begin
    return UNIT_CV_INTERN_TABLE[unit_or_cv]
  rescue => exception
    UNIT_CV_INTERN_TABLE[unit_or_cv] = unit_or_cv
    return UNIT_CV_INTERN_TABLE[unit_or_cv]
  end
end

class PyteomicsError < Exception
  def initialize(...)
    __init__(...)
  end

  def __init__(msg, *values)
    @message = msg
    @values = values
  end

  def to_s
    if [0, nil, false].include?(@values) || @values.empty?
      return "Pyteomics error, message: #{@message.inspect}"
    else
      return "Pyteomics error, message: #{@message.inspect}, #{@values.inspect})"
    end
  end
end

class Charge < Integer
  def initialize(...)
    __new__(...)
  end

  def __new__(*args, **kwargs)
    begin
      return super
    rescue => exception
      if args[0].instance_of?(String)
        begin
          num, sign = /^(\d+)(\+|-)$/.match(args[0])
          return super(sign + num, args[1..-1], kwargs)
        rescue => exception
          # pass
        end
      end
      raise PyteomicsError.new(*e.args)
    end
  end

  def __str__()
    self.abs.to_s + (self < 0 ? '-' : '+')
  end
end

class Ion < String
  @_pattern = /([abcxyz]\d+(\-H2O|\-NH3)?)([\+|-]\d+)/

  def initialize(...)
    __init__(...)
  end

  def __init__(*args, **kwargs)
    if args || args[0].instance_of?(String)
      begin
        @ion_type, @neutral_loss, charge = args[0].match(@_pattern)
      rescue => exception
        raise PyteomicsError.new("Malformed ion string, must match the regex #{_pattern}")
      end
    end
  end
end

class ChargeList < Array
  def initialize(...)
    __init__(...)
  end

  def __init__(*args, **kwargs)
    if ['', 0, nil, [], {}].include?(args).! || args[0].instance_of?(String)
      delim = /(?:,\s*)|(?:\s*and\s*)/
      self[0].nil? ? self[0] = Charge.new(delim.match(arges[0])) : self[0] += Charge.new(delim.match(args[0]))
    else
      begin
        super(arges[0].uniq.sort, *args[1..-1], **kwargs)
      rescue => exception
        super(*args, **kwargs)
      end
      self[0] = Charge.new(self[0])
    end
  end

  def __str__()
    if self.size > 1
      self[0..-2].join(', ') + " and #{self[-1]}"
    elsif self
      self[0].to_s
    else
      super.__str__(self)
    end
  end
end

def _parse_charge(s, list_only: false)
  if list_only.!
    begin
      return Charge.new(s)
    rescue => exception
      # pass
    end
  end
  return ChargeList.new(s)
end

def _parse_ion(ion_text)
  begin
    return Ion(ion_text)
  rescue => exception
    warn "Could not parse ion string: #{ion_text} #{exception.args[0]})"
  end
end

class BasicComposition < Hash
  def initialize(...)
    __init__(...)
  end

  def __init__(*args, **kwargs)
    self.default = 0

    if args
      args.each do |x|
        self[x] += 1
      end
    end
    if kwargs
      kwargs.each do |k, v|
        self[k] = v
      end
    end
  end

  def to_s
    # return '{}({})'.format(type(self).__name__, dict.__repr__(self))
    "#{self.class}(#{self.inspect})"
  end

  def __str__()
    to_s
  end
  def __repr__()
    to_s
  end

  def __repr_pretty__(p, cycle)
    if [0, '', nil, []].include?(cycle).!
      p.text("#{type().__name__} object with a cyclic reference'")
    end
    p.text(__str__())
  end

  def __add__(other)
    result = self.dup
    other.each do |elem, cnt|
      result[elem] += cnt
    end
    result
  end

  def +(...)
    __add__(...)
  end

  def __iadd__(other)
    other.each do |elem, cnt|
      self[elem] += cnt
    end
    @defaultdict
  end

  def __radd__(other)
    raise "PASS"
  end

  def __sub__(other)
    result = self.dup
    other.each do |elem, cnt|
      result[elem] -= cnt
    end
    result
  end

  def -(...)
    __sub__(...)
  end

  def __isub__(other)
    other.each do |elem, cnt|
      @defaultdict[elem] -= cnt
    end
    @defaultdict
  end

  def __rsub__(other)
    raise "PASS"
  end

  def __mul__(other)
    if other.instance_of?(Integer).!
      raise PyteomicsError.new('Cannot multiply Composition by non-integer', other)
    end
    self.each do |k, v|
      self[k] *= other
    end
  end

  def *(...)
    __mul__(...)
  end

  def __imul__(other)
    if other.instance_of?(Integer).!
      raise PyteomicsError.new('Cannot multiply Composition by non-integer', other)
    end
    @defaultdict.each do |k, v|
      @defaultdict[k] *= other
    end
    @defaultdict
  end

  def __rmul__(other)
    raise "PASS"
  end

  def __eq__(other)
    return false if other.instance_of?(Hash).!
    self_items = @defaultdict.select{ [0, '', nil, []].include?(_1[1]).! }
    other_items = other.select{ [0, '', nil, []].include?(_1[1]).! }
    self_items == other_items
  end

  def __missing__(key)
    0
  end

  def __setitem__(key, value)
    if value.instance_of?(Float)
      value = value.round
    elsif value.instance_of?(Integer).!
      raise PyteomicsError.new("Only integers allowed as values in Composition, got #{value.class}.")
    end
    if value != 0 # reject 0's
      self.merge!({key => value})
    elsif self.include?(key)
      self.delete(key)
    end
  end

  def []=(...)
    __setitem__(...)
  end

  def copy()
    @defaultdict.class.concat(@defaultdict)
  end

  def __reduce__()
    class_, args, state, list_iterator, dict_iterator = super.__reduce__()
    args = []
    return [class_, args, state, list_iterator, dict_iterator]
  end  
end

class MappingOverAttributeProxy
  def initialize(...)
    __init__(...)
  end

  def __init__(obj)
    @obj = obj
  end

  def __getitem__(key)
    @obj[key]
  end

  def __setitem__(key, value)
    @obj[key] = value
  end

  def __contains__(key)
    @obj.include?(key)
  end

  def __repr__()
    "#{@obj.class}(#{@obj.inspect})"
  end
end

class Unitint < Integer
  def initialize(value, unit_info)
    __new__(value, unit_info)
  end

  def __new__(value, unit_info)
    @inst = value.to_i
    @unit_info = unit_info
  end

  def __reduce__
    [self.__class__, [@inst, @unit_info]]
  end

  def _repr_pretty_(p, cycle)
    base = super.__repr__
    if @unit_info
      string = "#{base} #{@unit_info}"
    else
      string = base
    end
    p.text(string)
  end
end

class UnitFloat < Float
  def initialize(...)
    __new__(...)
  end

  def __new__(value, unit_info: nil)
    @value = value.to_f
    @unit_info = unit_info
    return self
  end

  #@property
  def __dict__
    _MappingOverAttributeProxy(self)
  end

  def __reduce__
    [self.class, [self.to_f, @unit_info]]
  end

  def _repr_pretty_(p, cycle)
    base = super.__repr__
    if @unit_info
      string = "#{base} #{@unit_info}"
    else
      string = base
    end
    p.text(string)
  end
end

class Unitstr < String
  def initialize(...)
    if PY2.!
      __slots__ = ["unit_info"]
    end
    __new__(...)
  end

  def __new__(value, unit_info: nil)
    # if PY2 && value.instance_of?(U)
    inst = super.new(cls, value)
    inst.unit_info = unit_info
    inst
  end

  @property
  def __dict__
    _MappingOverAttributeProxy.new
  end

  def __reduce__
    [self.class, [self.to_s, self.unit_info]]
  end

  def _repr_pretty_(p, cycle)
    base = super.__repr__
    if [0, '', nil, [], {}].include?(self.unit_info).!
      string = "#{base} #{self.unit_info}"
    else
      string = base
    end
    p.text(string)
  end
end

class Cvstr < String
  attr_accessor :accession, :unit_accession
  @@_cache = {}
  def initialize(value, ...)
    super(value)
    __new__(value, ...)
  end

  def __new__(value, accession: nil, unit_accession: nil)
    if @@_cache[value]
      @accession = accession
      @unit_accession = unit_accession
      return self
    end
    @accession = _intern_unit_or_cv(accession)
    @unit_accession = _intern_unit_or_cv(unit_accession)
    @@_cache[value] = self
    self
  end

  @property
  def __dict__
    _MappingOverAttributeProxy.new
  end

  def __reduce__
    [self.class, [self.to_s, self.accession, self.unit_accession]]
  end
end

class CVQueryEngine
  def _accession(key)
    key.respond_to?(:accession) ? key.accession : nil
  end

  def _query_dict(data, accession)
    data.each do |key, value|
      if self.accession(key) == accession
        if value.instance_of?(String).! || value != ''
          return value
        else
          return key
        end
      elsif value.instance_of?(Hash)
        inner = self._query_dict(value, accession)
        return inner if !!inner
      elsif value.instance_of?(Array)
        inner = self._query_sequence(value, accession)
        return inner if !!inner
      elsif self.accession(value) == accession
        return value
      end
    end
  end

  def _query_sequence(data, accession)
    data.each do |value|
      if value.instance_of?(Hash)
        inner = self._query_dict(value, accession)
        return inner if !!inner
      elsif value.instance_of?(Array)
        inner = self._query_sequence(value, accession)
        return inner if !!inner
      elsif _accession(value) == accession
        return value
      end
    end
  end

  def query(data, accession)
    if accession.nil?
      raise TypeError("'accession' cannot be None")
    end
    return self._query_dict(data, accession)
  end

  def _is_empty(value)
    return value == '' if value.instance_of?(String)
    false
  end

  def _walk_dict(data, index)
    data.each do |key, value|
      accession = _accession(key)
      if accession
        if self._is_empty(value).!
          index[accession] = value
        else
          index[accession] = key
        end
      elsif value.instance_of?(Hash)
        _walk_dict(value, index)
      elsif value.instance_of?(Array)
        _walk_sequence(value, index)
      end
      accession = _accession(value)
      index[accession] = value if accession
    end
    index
  end

  def _walk_sequence(data, index)
    data.each do |value|
      if value.instance_of?(Hash)
        _walk_dict(value, index)
      elsif value.instance_of?(Array)
        _walk_sequence(value, index)
      else
        accession = _accession(value)
        index[accession] = value if accession
      end
    end
  end

  def index(data)
    index = _walk_dict(data, {})
  end

  def __call__(data, accession: nil)
    if accession.nil?
      self.index(data)
    else
      self.query(data, accession)
    end
  end
end

Cvquery ||= CVQueryEngine