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

class XMLValueConverter
  # Not started
end

XMLParam = Struct.new(:xmlparam, [:name, :value, :type])

class XMLParam(xmlparam: ["name", "value", "type"])
  def initialize(**kwargs)
    @xmlparam = kwargs['XMLParam']
    @name = kwargs['name']
    @value = kwargs['value']
    @type = kwargs['type']
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
    items().sum{ |_, group| group.size })
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
  # Not started
end

class IndexedXML < XML
  prepend IndexedReaderMixin

  def initialize(...)
    @_indexed_tags = Set.new
    @_indexed_tag_keys = {}
    @_use_index = true
    __init__(...)
  end

  def __init__(source, read_schema: false, iterative: true, build_id_cache: false,
    use_index: nil, *args, **kwargs)
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
    # Not started
  end

  def __getstate__
    # Not started
  end

  def __setstate__(state)
    # Not started
  end

  # @_keepstate
  def _build_index
    if ['', 0, nil, false, [], {}].include?(@_indexed_tags) || ['', 0, nil, false, [], {}].include?(@_use_index)
      return
    @_offset_index = TagSpecificXMLByteIndex.build(
      @_source, @_indexed_tags, @_indexed_tag_keys)
  end
end

class MultiProcessingXML < IndexedXML
  include TaskMappingMixin
  # Not started
end

class IndexSavingXML < IndexedXML
  prepend IndexSavingMixin
  # Not started
end

class ArrayConversionMixin < BinaryDataArrayTransformer)
  # Not started
end

class Iterfind
  # Not started
end

class IndexedIterfind < Iterfind
  prepend TaskMappingMixin
  # Not started
end