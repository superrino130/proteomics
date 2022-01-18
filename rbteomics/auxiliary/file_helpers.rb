# import sys
# import codecs
# import re
# from functools import wraps
# from contextlib import contextmanager
# from collections import OrderedDict, defaultdict
# import json
# import multiprocessing as mp
# import threading
# import warnings
# import os
# from abc import ABCMeta

Basestring = String

require 'pandas'
require 'numpy'

# try:
#     import dill
# except ImportError:
#     dill = None
#     try:
#         import cPickle as pickle
#     except ImportError:
#         import pickle
#     serializer = pickle
# else:
#     serializer = dill

# try:
#     from queue import Empty
# except ImportError:
#     from Queue import Empty

# try:
#     from collections.abc import Sequence
# except ImportError:
#     from collections import Sequence

# from .structures import PyteomicsError
# from .utils import add_metaclass
require './rbteomics/auxiliary/structures'
require './rbteomics/auxiliary/utils'
require "pathname"

def _keepstate(func)
#   @wraps(func)
  def wrapped(*args, **kwargs)
    positions = [args.map{ (_1.respond_to?(:seek) ? (_1.seek && (_1.respond_to?(:tell) ? _1.arg : nill.class)) : nil) }]
    args.zip(positions).each do |arg, pos|
      arg.seek(0) if !!pos
    end
    res = func(*args, **kwargs)
    args.zip(positions).each do |arg, pos|
      if !!pos
        begin
          arg.seek(pos)
        rescue => exception
          # PASS
        end
      end
    end
    res
  end
  wrapped
end

def _keepstate_method(func)
#   @wraps(func)
  def wrapped(*args, **kwargs)
    position = self.tell
    self.seek(0)
    begin
      return func(*args, **kwargs)
    ensure
      self.seek(position)
    end
  end
  wrapped
end

# class _file_obj
class File_obj
  def initialize(...)
    __init__(...)
  end

  def __init__(f, mode, encoding: nil)
    @_file_spec = nil
    @mode = mode
    if f.nil?
      @file = {r: 'gets', a: 'puts', w: 'warn'}[mode[0]]
      @_file_spec = nil
    elsif f.instance_of?(Basestring)
      @file = File.open(f, mode)
      @file.set_encoding(encoding) if encoding.nil?.!
      @_file_spec = f
    else
      @_file_spec = f
      @file = f
    end
    @encoding = @file.respond_to?(:encoding) ? @file.encoding : encoding
    @close_file = @file != f
  end

  def __enter__()
    self
  end

  def __reduce_ex__(protocol)
    [self.class, [@file_spec, @mode, @encoding]]
  end

  def __exit__(*args, **kwargs)
    if @close_file.! || @_file_spec.nil?
      return # do nothing
    end
    exit = @file.respond_to?(:__exit__) ? @file.__exit__ : nil
    if !!exit
      return exit(*args, **kwargs)
    else
      exit = @file.respond_to?(:close) ? @file.close : nil
      exit() if !!exit
    end
  end

  def __getattr__(attr)
    @file.attr if @file.respond_to?(:attr)
  end

  def __iter__()
    iter(@file)
  end
end

module NoOpBaseReader
  def __init__(*args, **kwargs)
    # PASS
  end
end

class IteratorContextManager
  include NoOpBaseReader
  def initialize(*args, **kwargs)
    __init__(args, kwargs)
  end

  def __init__(*args, **kwargs)
    @_func = kwargs.delete('parser_func')
    @_args = args
    @_kwargs = kwargs
    reset() if self.class == IteratorContextManager
    super
  end

  def __getstate__()
    state = {}
    state['_iterator_args'] = @_args
    state['_iterator_kwargs'] = @_kwargs
    state
  end

  def __setstate__(state)
    @_args = state['_iterator_args']
    @_kwargs = state['_iterator_kwargs']
  end

  def reset
    begin
      @_reader = _func(*self._args, **self._kwargs)
    rescue => exception
      __exit__(*sys.exc_info())
      raise "Error in IteratorContextManager"
    end
  end

  def __enter__
    self
  end

  def __exit__(*args, **kwargs)
    # PASS
  end

  def __iter__
    self
  end

  def __next__
    _next(@_reader)
  end

#   _next = __next__()
end

# @add_metaclass(ABCMeta)
class FileReader < IteratorContextManager
  def initialize(...)
    __init__(...)
  end

  def __init__(source, **kwargs)
    func = kwargs['parser_func']
    # super(*kwargs[:args], parser_func: func, **kwargs[:kwargs])
    # super(*kwargs[:args], parser_func: func, **kwargs[:kwargs])
    super(kwargs[:args], parser_func: func, **kwargs[:kwargs])
    @_pass_file = kwargs['pass_file']
    @_source_init = source
    @_mode = kwargs['mode']
    @_encoding = kwargs['encoding']
    reset()
  end

  def reset
    if self.respond_to?(:_source)
      @_source.__exit__(nil, nil, nil)
    end
    @_source = File_obj.new(@_source_init, @_mode, @_encoding)
    begin
      if @_pass_file
        @_reader = _func(@_source, @_args, @_kwargs)
      else
        @_reader = _func(@_args, @kwargs)
      end
    rescue => exception
      @_source.__exit__(exception.message, exception.backtrace)
      raise
    end
  end

  def __exit__(*args, **kwargs)
    @_souce.__exit__(args, kwargs)
  end

  def __getattr__(attr)
    if attr == '_source'
      raise AttributeError.new
    end
    return @_source.respond_to?(:attr)
  end
end

def remove_bom(bstr)
  bstr.gsub(/\xFE\xFF|\xFF\xFE/n,"").gsub(/\xFE\xFF|\xFF\xFE/n,"").gsub!(/\x00$/, '')
end

module IndexedReaderMixin
  include NoOpBaseReader

  attr_reader :offset_index
  
  def index
    @offset_index
  end

  def default_index
    @offset_index
  end

  def __len__
    @offset_index.size
  end

  def __contaions__(key)
    @offset_index.include?(key)
  end

  def _item_form_offsets(other)
    raise NotImplementedError
  end

  def get_py_id(elem_id)
    index = @default_index
    if index.nil?
      raise PyteomicsError.new('Access by ID requires building an offset index.')
    end
    offsets = index[elem_id]
    _item_form_offsets(offsets)
  end

  def get_by_ids(ids)
    ids.map{ get_by_id(_1) }
  end

  def get_by_index(i)
    begin
      key = @default_index.from_index(i, false)
    rescue => exception
      raise PyteomicsError.new('Positional access requires building an offset index.')
    end
    get_by_id(key)
  end

  def get_by_indexes(indexes)
    indexes.map{ get_by_index(_1) }
  end

  def get_by_index_slice(s)
    begin
      keys = @default_index.from_slice(s, false)
    rescue => exception
      raise PyteomicsError.new('Positional access requires building an offset index.')
    end
    get_by_ids(keys)
  end

  def get_by_key_slice(s)
    keys = @default_index.between(s.start, s.stop)
    if s.step && s.step != 0
      keys = keys.select.with_index { |num, i| i % s.keys == 0 }
    end
    get_by_ids(keys)
  end

  def __getitem__(key)
    return get_by_id(key) if key.instance_of?(Basestring)
    return get_by_index(key) if key.instance_of?(Integer)
    if key.instance_of?(Sequence)
      return [] if key.!
      return get_by_indexes(key) if key[0].instance_of?(Integer)
      return get_by_ids(key) if key[0].instance_of?(Basestring)
    end
    # if isinstance(key, slice):
    if key.instance_of?(Array)
      return get_by_index_slice(key) if key[0].instance_of?(Integer)
      return get_by_key_slice(key) if key[0].instance_of?(Basestring)
      return self if key.nil?
    end
    raise PyteomicsError.new("Unsupported query key: #{key}")
  end
end

class RTLocator
  def initialize(reader)
    __init__(reader)
  end

  def __init__(reader)
    @_reader = reader
  end

  def _get_scan_by_time(time)
    if @_reader.default_index.!
      raise .new("This method requires the index. Please pass `use_index=True` during initialization")
    end

    scan_ids = [@_reader.default_index]
    lo = 0
    hi = scan_ids.size

    best_match = nil
    best_error = Float::INFINITY
    best_time = nil
    best_id = nil

    if time == Float::INFINITY
      scan = @_reader.get_by_id(scan_ids[-1])
      return [scan_ids[-1], scan, @_reader._get_time(scan)]
    end

    while hi != lo
      mid = (hi + lo) / 2
      sid = scan_ids[mid]
      scan = @_reader.get_by_id(sid)
      scan_time = @_reader._get_time(scan)
      err = (scan_time - time).abs
      if err < best_error
        best_error = err
        best_match = scan
        best_time = scan_time
        best_id = sid
      end
      if scan_time == time
        return [sid, scan, scan_time]
      elsif (hi - lo) == 1
        return [best_id, best_match, best_time]
      elsif scan_time > time
        hi = mid
      else
        lo = mid
      end
    end
  end

  def __getitem__(key)
    return _get_scan_by_time(key)[1] if [Integer, Float].include?(key.class)
    return key.map{ _get_scan_by_time(t)[1] } if key.instance_of?(Sequence)
    # if isinstance(key, slice):
    raise "there is no slice class"
  end
end

class TimeOrderedIndexedReaderMixin
  include IndexedReaderMixin

  def time()
      @_time
  end

  def initialize(*args, **kwargs)
    __init__(args, kwargs)
  end

  def __init__(*args, **kwargs)
    super
  end

  #staticmethod
  def _get_time(scan)
    raise NotImplementedError
  end
end

class IndexedTextReader < FileReader
  include IndexedReaderMixin

  def initialize(source, **kwargs)
    @delimiter = nil
    @label = nil
    @block_size = 1000000
    @label_group = 1
    __init__(source, kwargs)
  end

  def __init__(source, **kwargs)
    super(source, mode: 'rb', encoding: nil, **kwargs)
    @encoding = encoding
    ['delimiter', 'label', 'block_size', 'label_group'].each do |attr|
      if kwargs.include?(attr)
        eval("@#{attr} = #{kwargs.delete(attr)}")
      end
    end
    @_offset_index = nil
    if kwargs.delete('_skip_index').!
        @_offset_index = build_byte_index()
    end
  end

  def __getstate__
    state = super
    state['offset_index'] = @_offset_index
    state
  end

  def __setstate__(state)
    super
    @_offset_index = state['offset_index']
  end

  def _chunk_iterator
    # 未実装
  end

  def _generate_offsets
    Fiber.new do
      i = 0
      pattern = Regexp.compile(remove_bom(@label.encode(@encoding)))
      _chunk_iterator().each do |chunk|
        match = pattern.search(chunk)
        if match
          label = match[@label_group]
          yield [i, label.decode(@encoding), match]
        end
        i += chunk.size
      end
      yield [i, nil, nil]
    end
  end

  def build_byte_index
    index = OffsetIndex.new
    g = _generate_offsets()
    last_offset = 0
    last_label = nil
    g.each do |offset, label, keyline|
      if !!last_label
        index[last_label] = [last_offset, offset]
      end
      last_label = label
      last_offset = offset
    end
    raise "last_label.nil?" if last_label.nil?
    index
  end

  def _read_lines_from_offsets(start, _end)
    _source.seek(start)
    lines = _source.read(_end - start).decode(@encoding).split('\n')
    lines
  end
end

module IndexSavingMixin
  include NoOpBaseReader
  _index_class = NotImplementedError

  #property
  def _byte_offset_filename
    begin
      path = _source.name      
    rescue => exception
      return nil
    end
    name = Pathname(path).sub_ext('')
    ext = Pathname(path).extname
    byte_offset_filename = "#{name}-#{ext[1..-1]}-byte-offsets.json"
  end

  def _check_has_byte_offset_file
    path = @_byte_offset_filename
    return false if path.nil?
    return Pathname(path).exist?
  end

  #classmethod
  def self.prebuild_byte_offset_file(path)
    inst = self.new(path)
    inst.write_byte_offsets()
  end

  def write_byte_offsets
    File.open(@_byte_offset_filename, 'w') do |f|
      f.write(@_offset_index)
    end
  end

  @_keepstate_method
  def _build_index
    return if @_use_index.!
    begin
      _read_byte_offsets()
    rescue => exception
      super
    end
  end

  def _read_byte_offsets
    File.open(@_byte_offset_filename) do |f|
      index = _index_class.load(f)
      @_offset_index = index
    end
  end
end

def _file_reader(_mode: 'r')
  def decorator(_func)
    # @wraps(_func)
    def helper(*args, **kwargs)
      if args
        return FileReader(args[0], mode: _mode, parser_func: _func, pass_file: true, args: args[1..-1], kwargs: kwargs,
          encoding: kwargs.delete('encoding') || nil)
      end
      source = kwargs.delete('source') || nil
      FileReader(source, mode: _mode, parser_func: _func, pass_file: true, args: [], kwargs: kwargs, encoding: kwargs.delete('encoding') || nil)
    end
    helper
  end
  decorator
end

