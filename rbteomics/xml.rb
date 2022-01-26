# import re
# import socket
# from traceback import format_exc
# import warnings
# from collections import OrderedDict, namedtuple
# from itertools import islice
# from lxml import etree
# import numpy as np
require 'numpy'

# from .auxiliary import FileReader, PyteomicsError, basestring, _file_obj, HierarchicalOffsetIndex
# from .auxiliary import unitint, unitfloat, unitstr, cvstr
# from .auxiliary import _keepstate_method as _keepstate
# from .auxiliary import BinaryDataArrayTransformer
# from .auxiliary import TaskMappingMixin, IndexedReaderMixin, IndexSavingMixin
require_relative 'auxiliary/structures'
require_relative 'auxiliary/file_helpers'
BaseString ||= String


# try:  # Python 2.7
#     from urllib2 import urlopen, URLError
# except ImportError:  # Python 3.x
#     from urllib.request import urlopen, URLError
require 'set'
require 'open-uri'

def _local_name(element)
  tag = element.tag
  if ['', 0, nil, false, [], {}].include?(tag).! && tag[0] == '{'
    return tag.rpartition('}')[2]
  end
  tag
end

def xsd_parser(schema_url)
  # Not started
end

module XMLValueConverter
  def self.duration_str_to_float(s)
    if s.star_with?('P').!
      begin
        return Unitfloat.new(s, 'duration')
      rescue => exception
        return Unitstr.new(s, 'duration')
      end
    end
    _duration_parser = Regexp.compile('(?P<sign>-?)P(?:(?P<years>\d+\.?\d*)Y)?(?:(?P<months>\d+\.?\d*)M)?(?:(?P<days>\d+\.?\d*)D)?(?:T(?:(?P<hours>\d+\.?\d*)H)?(?:(?P<minutes>\d+\.?\d*)M)?(?:(?P<seconds>\d+\.?\d*)S)?)?')
    match = self._duration_parser.search(s)
    if match.nil?.!
      matchdict = match.groupdict()
      hours = matchdict['hours'].to_f || 0
      minutes = matchdict['minutes'].to_f || 0
      seconds = matchdict['seconds'].to_f || 0
      minutes += hours * 60.0
      minutes += seconds / 60.0
      return Unitfloat.new(minutes, 'minute')
    else
      return Unitstr.new(s, 'duration')
    end
  end

  def self.str_to_bool(s)
    return true if ['true', '1', 'y'].include?(s.downcase)
    return false if ['false', '0', 'n'].include?(s.downcase)
    raise PyteomicsError.new('Cannot convert string to bool: ' + s)
  end

  def str_to_num(s, numtype)
    ['', nil].include?(s).! ? numtype(s) : nil
  end

  def self.to(t)
    def convert_from(s)
      return str_to_num(s, t)
    end
    convert_from
  end

  def self.converters
    return {
      'ints': self.to(Unitint),
      'floats': self.to(Unitfloat),
      'bools': self.str_to_bool,
      'intlists': lambda { |x| Numpy.fromstring(x.replace('\n', ' '), dtype: int, sep: ' ') },
      'floatlists': lambda { |x| Numpy.fromstring(x.replace('\n', ' '), sep: ' ') },
      'charlists': list,
      'duration': duration_str_to_float
    }
  end
end

SXMLParam ||= Struct.new(:xmlparam, [:name, :value, :type])

class XMLParam
  def initialize(sxmlparam: nil)
    @sxmlparam = sxmlparam
  end
end

class XML < FileReader
  # Not started
end

Pattern_path = Regexp.compile('([\w/*]*)(.*)')

def get_rel_path(element, names)
  # Not started
end

def xpath(tree, path, ns: nil)
  # Not started
end

def _make_version_info
  # Not started
end

class ByteCountingXMLScanner < File_obj
  # Not started
end

class TagSpecificXMLByteIndex
  def initialize(...)
    @_default_indexed_tags = []
    @_default_keys = {}
    @_scanner_class = ByteCountingXMLScanner
    __init__(...)
  end

  def __init__(source, indexed_tags: nil, keys: nil)
    keys = @_default_keys.dup if keys.nil?
    indexed_tags = @_default_indexed_tags if indexed_tags.nil?
    @indexed_tags = indexed_tags
    @indexed_tag_keys = keys
    @source = source
    @offsets = HierarchicalOffsetIndex.new
    build_index()
  end

  def __getstate__
    # Not started
  end

  def __setstate__(state)
    # Not started
  end

  def __getitem__(key)
    @offsets[key]
  end

  def build_index
    scanner = @_scanner_class.new(@source, @indexed_tags)
    @offsets = scanner.build_byte_index(@indexed_tag_keys)
    @offsets
  end

  def items
    @offsets
  end

  def keys
    @offsets.keys
  end

  def __iter__
    @keys.to_enum
  end

  def __len__
    items().sum{ |_, group| group.size }
  end

  def self.build(source, indexed_tags: nil, keys: nil)
    indexer = cls(source, indexed_tags, keys)
    indexer.offsets
  end
end

def ensure_bytes_single(string)
  # Not started
end

def ensure_bytes(strings)
  # Not started
end

def _flatten_map(hierarchical_map)
  all_records = []
  hierarchical_map.each do |key, records|
    all_records.concat(records)
  end
  all_records.sort!{ |x| x[1] }
  OrderedDict(all_records)
end

class IndexedXML < XML
  prepend IndexedReaderMixin

  def initialize(...)
    @_indexed_tags = Set.new
    @_indexed_tag_keys = {}
    @_use_index = true
    __init__(...)
  end

  def __init__(source, *args, **kwargs)
    read_schema = kwargs['read_schema'] || false
    iterative = kwargs['iterative'] || true
    build_id_cache = kwargs['build_id_cache'] || false
    use_index = kwargs['use_index'] || nil
    tags = kwargs['indexed_tags']
    tag_index_keys = kwargs['indexed_tag_keys']

    @_indexed_tags = tags if tags.nil?.!
    @_indexed_tag_keys = tag_index_keys if tag_index_keys.nil?.!

    @_use_index = use_index if use_index.nil?.!
    if ['', 0, nil, false, [], {}].include?(use_index).!
      build_id_cache = false
      if ['', 0, nil, false, [], {}].include?(@_default_iter_path).! && @_default_iter_path != @_default_iter_tag
        warn('_default_iter_path differs from _default_iter_tag and index is enabled. _default_iter_tag will be used in the index, mind the consequences.')
      end
    end
    super(source, read_schema, iterative, build_id_cache, *args, **kwargs)

    @_offset_index = nil
    _build_index()
  end

  # @property
  def default_index
    @_offset_index[@_default_iter_tag]
  end

  def __reduce_ex__(protocol)
    reconstructor, args, state = XML.new.__reduce_ex__(protocol)
    args += [False]
    [reconstructor, args, state]
  end

  def __getstate__
    state = super
    state['_indexed_tags'] = @_indexed_tags
    state['_indexed_tag_keys'] = @_indexed_tag_keys
    state['_use_index'] = @_use_index
    state['_offset_index'] = @_offset_index
    state
  end

  def __setstate__(state)
    super
    @_indexed_tags = state['_indexed_tags']
    @_indexed_tag_keys = state['_indexed_tag_keys']
    @_use_index = state['_use_index']
    @_offset_index = state['_offset_index']
  end

  # @_keepstate
  def _build_index
    if ['', 0, nil, false, [], {}].include?(@_indexed_tags) || ['', 0, nil, false, [], {}].include?(@_use_index)
      return
    @_offset_index = TagSpecificXMLByteIndex.build(
      @_source, @_indexed_tags, @_indexed_tag_keys)
  end

  #@_keepstate
  def _find_by_id_reset(elem_id, id_key: nil)
    _find_by_id_no_reset(elem_id, id_key: id_key)
  end

  #@_keepstate
  def get_by_id(elem_id, **kwargs)
    id_key = kwargs['id_key'] || nil
    element_type = kwargs['element_type'] || nil
    begin
      index = @_offset_index
      if element_type.nil?
        offset, element_type = index.find_no_type(elem_id)
      else
        offset = index.find(elem_id, element_type)
      end
      @_source.seek(offset)
      if id_key.nil?
        id_key = @_indexed_tag_keys[element_type]
      end
      elem = _find_by_id_no_reset(elem_id, id_key: id_key)        
    rescue => exception
      elem = _find_by_id_reset(elem_id, id_key: id_key)      
    end
    data = _get_info_smart(elem, **kwargs)
  end

  def __contains__(key)
    @_offset_index[@_default_iter_tag].include?(key)
  end

  def __len__
    @_offset_index[@_default_iter_tag].size
  end

  def iterfind(path, **kwargs)
    if @_indexed_tags.include?(path) && @_use_index
      return IndexedIterfind.new(path, **kwargs)
    end
    Iterfind.new(path, **kwargs)
  end
end

class MultiProcessingXML < IndexedXML
  include TaskMappingMixin
  def _task_map_iterator
    @_offset_index[self._default_iter_tag].to_enum
  end
end

class IndexSavingXML < IndexedXML
  prepend IndexSavingMixin
  
  def initialize(...)
    @_index_class = HierarchicalOffsetIndex.new
  end

  def _read_byte_offsets
    File.open(@_byte_offset_filename, 'r') do |f|
      index = @_index_class.load(f)
      if index.schema_version.nil?
        raise TypeError("Legacy Offset Index!")
      end
      @_offset_index = index
    end
  end
end

class ArrayConversionMixin < BinaryDataArrayTransformer
  def initialize(...)
    @_dtype_dict = {}
    @_array_keys = ['m/z array', 'intensity array']
    __init__(...)
  end

  def __init__(*args, **kwargs)
    @_dtype_dict = {'None' => nil}
    dtype = kwargs.delete('dtype') || nil
    if dtype.is_a(Hash)
      @_dtype_dict.merge!(dtype)
    elsif dtype.nil?.!
      @_dtype_dict = @_array_keys.map{ |k| [k, dtype] }.to_h
      @_dtype_dict['None'] = dtype
    end
    super
  end

  def __getstate__
    state = super
    state['_dtype_dict'] = @_dtype_dict
    state
  end

  def __setstate__(state)
    super
    @_dtype_dict = state['_dtype_dict']
  end

  def _convert_array(k, array)
    dtype = @_dtype_dict[k]
    if dtype.nil?.!
      return array.astype(dtype)
    end
    array
  end

  def _finalize_record_conversion(array, record)
    key = record.key
    return _convert_array(key, array)
  end
end

class Iterfind
  def initialize(...)
    @md = Struct.new(:start, :stop, :step)
    __init__(...)
  end

  def __init__(parser, tag_name, **kwargs)
    @parser = parser
    @tag_name = tag_name
    @config = kwargs
    @_iterator = nil
  end

  def __repr__
    if @config != ''
      @config = ", " + repr(@config)
    else
      @config = ''
    end
    template = "#{self.class}(#{__repr__(@tag_name)}#{@config})"
  end

  def __iter__
    self
  end

  def _make_iterator
    @parser._iterfind_impl(@tag_name, **@config)
  end

  def __next__
    if @_iterator.nil?
      @_iterator = _make_iterator()
    end
    return _next(@_iterator)
  end

  def _next
    __next__()
  end

  #@property
  def is_indexed
    false
  end

  def reset
    @_iterator = nil
    @parser.reset()
  end

  def __enter__
    self
  end

  def __exit__(*args, **kwargs)
    reset()
  end

  def map(*args,**kwargs)
    raise NotImplementedError.new("This query isn't indexed, it cannot be mapped with multiprocessing")
  end

  def _get_by_index(idx)
    self.reset()
    value = _next(islice(idx, idx + 1))
  end

  def _get_by_slice(slc)
    self.reset()
    value = islice(slc.start, slc.stop, slc.step).to_a
  end

  def __getitem__(i)
    if i.is_a?(Struct)
        return _get_by_slice(i)
    end
    return _get_by_index(i)
  end
end

class IndexedIterfind < Iterfind
  prepend TaskMappingMixin
  
  def initialize(...)
    __init__(...)
  end

  def __init__(parser, tag_name, **kwargs)
    super
  end

  def _task_map_iterator
    _index.to_enum
  end

  def _index
    @parser.index(@tag_name)
  end

  def _get_reader_for_worker_spec
    @parser
  end

  def _yield_from_index
    _task_map_iterator.each do |key, _|
      yield @parser.get_by_id(key, @config)
    end
  end

  def _make_iterator
    if is_indexed
      return _yield_from_index()
    end
    warn "Non-indexed iterator created from #{self}"
    super
  end

  #@property
  def is_indexed
    if @parser.respond_to?(:_index)
      if @parser._index.nil?.!
        index = @parser._index
        if index.is_a?(HierarchicalOffsetIndex)
          return index.include?(@tag_name) && ['', 0, nil, false, [], {}].include?(index[@tag_name]).!
        end
      end
    end
    false
  end

  def _get_by_index(idx)
    key = @index.form_index(idx)
    @parser.get_by_id(key)
  end

  def _get_by_slice(slc)
    keys = @index.form_slice(slc)
    @parser.get_by_ids(keys)
  end

  def __len__
    @index.size
  end
end