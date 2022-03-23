# import codecs
# from functools import wraps
# from contextlib import contextmanager
# from collections import OrderedDict, defaultdict
require 'json'
# import multiprocessing as mp
# import threading
# from abc import ABCMeta

BaseString ||= String

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
require_relative 'structures'
require_relative 'utils'
# require_relative '../xml'
require "pathname"

module File_helpers
  module_function
  def _keepstate(func)
    #   @wraps(func)
    orig = method(func)
    define_singleton_method(func) do |*args, **kwargs|
      positions = args.map{ |arg| (arg.respond_to?(:seek) ? arg.seek : nil) && (arg.respond_to?(:tell) ? arg.tell : nill.class) }
      args.zip(positions).each do |arg, pos|
        arg.seek(0) if pos.nil?.!
      end
      res = orig.call(*args, **kwargs)
      args.zip(positions).each do |arg, pos|
        if pos.nil?.!
          begin
            arg.seek(pos)
          rescue => exception
            # PASS
          end
        end
      end
      return res
    end
  end
    
  def _keepstate_method(func)
    orig = method(func)
    define_singleton_method(func) do |*args, **kwargs|
      # position = self.tell # tell : binary mode
      self.seek(0) if self.respond_to?(:seek)
      begin
        orig.call(*args, **kwargs)
      ensure
        # self.seek(position)
      end
    end
  end
  
  class File_obj
    def initialize(...)
      __init__(...)
    end
  
    def __init__(f, mode, encoding: nil)
      @_file_spec = nil
      @mode = mode
      if f.nil?
        @file = {'r' => 'gets', 'a' => 'puts', 'w' => 'warn'}[mode[0]]
        @_file_spec = nil
      elsif f.instance_of?(BaseString)
        # @file = File.open(f, mode)
        @file = File.open(f)
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
        exit() if exit.nil?.!
      end
    end
  
    def __getattr__(attr)
      @file.attr if @file.respond_to?(:attr)
    end
  
    def __iter__()
      iter(@file)
    end
  
    def lines
      @file
    end
  end
  
  module NoOpBaseReader
    def __init__(...)
      # PASS
    end
  end
  
  module IteratorContextManager
    include NoOpBaseReader
    module_function
    def __init__(*args, **kwargs)
      @_func ||= kwargs.delete('parser_func')
      @_args ||= args
      @_kwargs ||= kwargs
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
        @_reader = @_func.call(*@_args, **@_kwargs)
      rescue => exception
        # __exit__(*sys.exc_info())
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
  module FileReader
    include IteratorContextManager
    module_function
    def __init__(source, **kwargs)
      func = kwargs['parser_func']
      super(*kwargs['args'], 'parser_func' => func, **kwargs['kwargs'])
      @_pass_file ||= kwargs['pass_file']
      @_source_init ||= source
      @_mode ||= kwargs['mode']
      @_encoding ||= kwargs['encoding']
      reset
    end
  
    def reset
      if self.respond_to?(:_source)
        @_source.__exit__(nil, nil, nil)
      end
      @_source = File_obj.new(@_source_init, @_mode, encoding: @_encoding).lines.readlines
      begin
        if @_pass_file.!
          # @_reader = @_func.call(@_source, *@_args, **@_kwargs)
          @_reader = @_func.call(@_source, **@_kwargs)
        else
          @_reader = @_func.call(*@_args, **@_kwargs)
        end
      rescue => exception
        # @_source.__exit__(exception.message, exception.backtrace)
        raise
      end
    end
  
    def __exit__(*args, **kwargs)
      @_source.__exit__(args, kwargs)
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
  
    module_function
    
    attr_reader :_offset_index
    
    def index
      @_offset_index
    end
  
    def default_index
      @_offset_index
    end
  
    def __len__
      @_offset_index.size
    end
  
    def __contaions__(key)
      @_offset_index.include?(key)
    end
  
    def _item_form_offsets(other)
      raise NotImplementedError.new
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
      return get_by_id(key) if key.instance_of?(BaseString)
      return get_by_index(key) if key.instance_of?(Integer)
      if key.instance_of?(Sequence)
        return [] if key.!
        return get_by_indexes(key) if key[0].instance_of?(Integer)
        return get_by_ids(key) if key[0].instance_of?(BaseString)
      end
      # if isinstance(key, slice):
      if key.instance_of?(Array)
        return get_by_index_slice(key) if key[0].instance_of?(Integer)
        return get_by_key_slice(key) if key[0].instance_of?(BaseString)
        return self if key.nil?
      end
      raise PyteomicsError.new("Unsupported query key: #{key}")
    end
  end
  
  class RTLocator
    def initialize(...)
      __init__(...)
    end
  
    def __init__(reader)
      @_reader = reader
    end
  
    def _get_scan_by_time(time)
      if @_reader.default_index.!
        raise PyteomicsError.new("This method requires the index. Please pass `use_index=True` during initialization")
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
      if key.start.nil?
        start_index = @_reader.default_index.from_index(0)
      else
        start_index = _get_scan_by_time(key.start)[0]
      end
      if key.stop.nil?
        stop_index = @_reader.default_index.from_index(-1)
      else
        stop_index = _get_scan_by_time(key.stop)[0]
      end
      raise "there is no slice class"
    end
  end
  
  module TimeOrderedIndexedReaderMixin
    include IndexedReaderMixin
  
    module_function
  
    #@property
    def time
        @_time
    end
  
    def __init__(*args, **kwargs)
      super
      @_time = RTLocator.new(self)
    end
  
    #staticmethod
    def self._get_time(scan)
      raise NotImplementedError.new
    end
  end
  
  module IndexedTextReader
    include IndexedReaderMixin
    include FileReader
  
    module_function
  
    def __init__(source, **kwargs)
      encoding = kwargs.delete('encoding') || 'utf-8'
      super(source, 'mode' => 'rb', 'encoding' => nil, **kwargs)
      @delimiter = nil
      @label = nil
      @block_size = 1000000
      @label_group = 1
    
      @encoding = encoding
      ['delimiter', 'label', 'block_size', 'label_group'].each do |attr|
        if kwargs.include?(attr)
          eval("@#{attr} = #{kwargs.delete(attr)}")
        end
      end
      @_offset_index = nil
      if (kwargs.delete('_skip_index') || false).!
          @_offset_index = self.build_byte_index()
      end
    end
  
    def __getstate__
      state = @c1.__getstate__ || @c2.__getstate__
      state['offset_index'] = @_offset_index
      state
    end
  
    def __setstate__(state)
      @c1.__setstate__(state)
      @c2.__setstate__(state)
      @_offset_index = state['offset_index']
    end
  
    def _chunk_iterator
      # Not started
    end
  
    def _generate_offsets
      Fiber.new do
        i = 0
        pattern = Regexp.compile(remove_bom(@label.encode(@encoding)))
        _chunk_iterator().each do |chunk|
          match = pattern.search(chunk)
          if match
            label = match[@label_group]
            Fiber.yield [i, label.decode(@encoding), match]
          end
          i += chunk.size
        end
        Fiber.yield [i, nil, nil]
      end
    end
  
    def build_byte_index
      index = OffsetIndex.new
      g = _generate_offsets()
      last_offset = 0
      last_label = nil
      g.each do |offset, label, keyline|
        if last_label.nil?.!
          index[last_label] = [last_offset, offset]
        end
        last_label = label
        last_offset = offset
      end
      raise "last_label.nil?" if last_label.nil?
      index
    end
  
    def _read_lines_from_offsets(start, _end)
      @_source.seek(start)
      lines = @_source.read(_end - start).decode(@encoding).split('\n')
      lines
    end
  end
  
  module IndexSavingMixin
    include File_helpers
    include NoOpBaseReader
  
    module_function
  
    def __init__(...)
      @_index_class ||= NotImplementedError
      super
    end
  
    #property
    def _byte_offset_filename
      begin
        # path = _source.name      
        path = @_source_init  
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
  
    #@_keepstate_method
    def in_build_index
      return if @_use_index.!
      begin
        _read_byte_offsets()
      rescue => exception
        super
      end
    end
    def _build_index
      _keepstate_method(:in_build_index)
      in_build_index
    end
  
    def _read_byte_offsets
      File.open(@_byte_offset_filename) do |f|
        index = @@_index_class.load(f)
        @_offset_index = index
      end
    end
  end
  
  def _file_reader(_mode: 'r')
    def decorator(_func)
      # @wraps(_func)
      def helper(*args, **kwargs)
        if args
          return FileReader.new(args[0], 'mode' => @_mode, 'parser_func' => @_func, 'pass_file' => true, 'args' => args[1..], 'kwargs' => kwargs,
            'encoding' => kwargs.delete('encoding') || nil)
        end
        source = kwargs.delete('source') || nil
        FileReader.new(source, 'mode' => @_mode, 'parser_func' => @_func, 'pass_file' => true, 'args' => [], 'kwargs' => kwargs, 'encoding' => kwargs.delete('encoding') || nil)
      end
      helper
    end
    decorator
  end
  
  def _file_writer(_mode: 'a')
    def docerator(_func)
  
    end
  
  end
  
  module WritableIndex
    Schema_version = [1, 0, 0]
    Schema_version_tag_key = "@pyteomics_schema_version"
  
    module_function
  
    def _serializable_container
      {'index' => self.constants}
    end
  
    def save(fp)
      container = _serializable_container
      container[Schema_version_tag_key] = Schema_version
      JSON.dump(container, fp)
    end
  
    def self.load(cls, fp)
      container = json.load(fp)
      version_tag = container[@_schema_version_tag_key]
      if version_tag.nil?
        inst = cls
        inst.Schema_version = nil
        return inst
      end
      version_tag = [version_tag]
      index = container['index']
      if version_tag < cls.Schema_version
        inst = cls(index)
        inst.Schema_version = version_tag
        return inst
      end
      return cls(index)
    end
  end
  
  class OffsetIndex
    include WritableIndex
  
    def initialize(...)
      __init__(...)
    end
    
    def __init__(*args, **kwargs)
      @hash = {}
      if args.empty?.!
        if args.size == 1 && args[0].is_a?(String)
          args[0].split('').each do |s|
            @hash[s] = 1
          end
        else
          args.each do |x|
            if x.is_a?(Array)
              x.each do |y|
                if y.is_a?(Array)
                  @hash[y[0]] = y[1]
                else
                  @hash[y] = 1
                end
              end
            else
              @hash[x] = 1
            end
          end
        end
      elsif kwargs.empty?.!
        @hash.merge!(kwargs)
      end
      @_index_sequence = nil
    end
  
    def _invalidate
      @_index_sequence = nil
    end
  
    #@property
    def index_sequence
      if @_index_sequence.nil?
        @_index_sequence = @hash.to_a
      end
      @_index_sequence
    end
  
    def __setitem__(key, value)
      _invalidate()
      return super
    end
  
    def pop(*args, **kwargs)
      _invalidate()
      return super
    end
  
    def find(key, *args, **kwargs)
      return @hash[key]
    end
  
    def from_index(index, include_value: false)
      items = @hash.sort
      if include_value
        items[index]
      else
        items[index][0]
      end
    end
  
    def from_slice(spec, include_value: false)
      items = @hash.sort
      items[spec].map{ |k, v| include_value ? [k, v] : k }
    end
  
    def between(start, stop, include_value: false)
      keys = @hash.keys.sort
      if start.nil?.!
        begin
          start_index = keys.index(start)
        rescue => exception
          raise KeyError start
        end
      else
        start_index = 0
      end
      if stop.nil?.!
        begin
          stop_index = keys.index(stop)
        rescue => exception
          raise KeyError stop
        end
      else
        stop_index = keys.size - 1
      end
      if start.nil? || stop.nil?
        #pass  # won't switch indices
      else
        start_index, stop_index = [start_index, stop_index].min, [start_index, stop_index].max
      end
  
      if include_value
        return keys[start_index...stop_index + 1].map{ |k| [k, @hash[k]] }
      end
      return keys[start_index...stop_index + 1]
    end
  
    def __repr__
      "#{self.class}(#{self.to_a})"
    end
  
    def _integrity_check
      indices = self.values
      sorted_indices = self.values.sort
      indices == sorted_indices
    end
  
    def sort
      sorted_pairs = self.sort_by{ |_, v| v }
      self.clear
      self._invalidate()
      sorted_pairs.each do |key, value|
        self[key] = value
      end
      return self
    end
  end
  
  module IndexSavingTextReader
    include IndexedTextReader
    include IndexSavingMixin
    @@_index_class = OffsetIndex
  end
  
  module HierarchicalOffsetIndex
    include WritableIndex
    module_function
    @@_inner_type = OffsetIndex
    
    def __init__(base: nil)
      @mapping = Hash.new #(self._inner_type)
      (base || {}).each do |key, value|
        @mapping[key] = @_inner_type.new(value)
      end
    end
  
    def _integrity_check
      self.each do |key, value|
        if value._integrity_check().!
          return false
        end
      end
      true
    end
  
    def sort
      self.each do |key, value|
        value.sort!
      end
    end
  
    def __getitem__(key)
      @mapping[key]
    end
  
    def __setitem__(key, value)
      @mapping[key] = value
    end
  
    def __iter__
      @mapping.to_enum
    end
  
    def __len__
      self.sum{ |key, group| group.size }
    end
  
    def __contains__(key)
      @mapping.include?(key)
    end
  
    def find(key, element_type: nil)
      if element_type.nil?
        self.keys.each do |element_type|
          begin
            return find(key, element_type)
          rescue => exception
            next
          end
        end
        raise KeyError.new(key)
      else
        self[element_type][key]
      end
    end
  
    def find_no_type(key)
      self.keys.each do |element_type|
        begin
          return [find(key, element_type), element_type]
        rescue => exception
          next
        end
      end
      raise KeyError.new(key)
    end
  
    def update(*args, **kwargs)
      @mapping.update(*args, **kwargs)
    end
  
    def pop(key, default: nil)
      @mapping.delete(key) || default
    end
  
    def keys
      @mapping.keys
    end
  
    def values
      @mapping.values
    end
  
    def items
      @mapping
    end
  
    def _serializable_container
      encoded_index = {}
      container = {
        'keys': keys
      }
      self.each do |key, offset|
        if offset.is_a?(Hash)
          encoded_index[key] = offset.to_a
        elsif offset.is_a?(String)
          encoded_index[key] = offset.split('')
        else
          encoded_index[key] = offset
        end
      end
      container['index'] = encoded_index
      container
    end
  end
  
  def _make_chain(reader, readername, full_output: false)
    @full_output = full_output
    def concat_results(*args, **kwargs)
      results = args.map{ |arg| reader(arg, **kwargs) }
      if args.map{ |a| a.instance_of?(Pandas.DataFrame)}.all?
        return Pandas.concat(results)
      end
      return Numpy.concatenate(results)
    end
  
    def _iter(files, kwargs)
      files.each do |f|
        File.open(f) do |r|
          r.each do |item|
            yield item
          end
        end
      end
    end
  
    def chain(*files, **kwargs)
      _iter(files, kwargs)
    end
  
    def from_iterable(files, **kwargs)
      _iter(files, kwargs)
    end
  
    #@contextmanager
    def _chain(...)
      yield chain(...)
    end
  
    #@contextmanager
    def _from_iterable(...)
      yield from_iterable(...)
    end
  
    dispatch = lambda do |x = nil, *args, **kwargs|
      if x.nil?
        @dispatch_from_iterable.call(*args, **kwargs)
      elsif x == 'from_iterable' || x == :from_iterable
        @dispatch_from_iterable.call(*args, **kwargs)
      end
    end
  
    @dispatch_from_iterable = lambda do |*args, **kwargs|
      if kwargs['full_output'] || @full_output
        return concat_results(*args, **kwargs)
      end
      return _chain(*args, **kwargs)
    end
  
    # dispatch.__doc__ = "Chain :py:func:`#{readername}` for several files. Positional arguments should be file names or file objects. Keyword arguments are passed to the :py:func:`#{readername}` function."
    # dispatch_from_iterable.__doc__ = "Chain :py:func:`#{readername}` for several files. Keyword arguments are passed to the :py:func:`#{readername}` function./nParameters/n----------/nfiles : iterable/n    Iterable of file names or file objects.\n"
    dispatch
  end
  
  def _check_use_index(source, use_index, default)
    begin
      if use_index.nil?.!
        use_index = ['', 0, nil, false, [], {}].include?(use_index).!
      end
  
      begin
        f = File.open(source)
      rescue => exception
        if source.is_a?(BaseString)
          return (use_index.nil?.! ? use_index : default)
        end        
      end
  
      begin
        File.open(source) do |f|
          @seekable = f.methods.include?(:seek)
        end
      rescue => exception
        @seekable = nil
      end
  
      if f.nil?.!
        binary = f.read.encoding == Encoding::ASCII_8BIT
      else
        binary = nil
      end
  
      if @seekable == false
        if binary
          raise PyteomicsError.new("Cannot work with non-seekable file in binary mode: #{source}.")
        end
        if use_index
          warn "Cannot use indexing as #{source} is not seekable. Setting `use_index` to False."
          use_index = false
        end
      elsif binary.nil?.!
        if use_index.nil?.! && binary != use_index
          warn "use_index is #{use_index}, but the file mode is #{source.mode}. Setting `use_index` to #{binary}"
        end
        use_index = binary
      elsif use_index.nil?
        warn "Could not check mode on #{source}. Specify `use_index` explicitly to avoid errors."
      end
  
      if use_index.nil?.!
        return use_index
      end
  
      return default
    rescue => exception
      if exception.is_a?(PyteomicsError)
        raise
      else
        if use_index.nil?
          warn "Could not check mode on #{source}. Reason: #{exception.message}. Specify `use_index` explicitly to avoid errors."
          return default
        end
        warn exception.message
        return use_index
      end
    end
  end
  
  class FileReadingProcess
    # Not started
  end
  
  begin
    require 'etc'
    NPROC = Etc.nprocessers  
  rescue => exception
    NPROC = 4
  end
  QUEUE_TIMEOUT = 4
  QUEUE_SIZE = 1e7.to_i
  
  module TaskMappingMixin
    include NoOpBaseReader
    def __init__(*args, **kwargs)
      @_queue_size ||= kwargs.delete('queue_size') || QUEUE_SIZE
      @_queue_timeout ||= kwargs.delete('timeout') || QUEUE_TIMEOUT
      @_nproc ||= kwargs.delete('processes') || NPROC
      super
    end
  
    def _get_reader_for_worker_spec
      self
    end
  
    def _build_worker_spec(target, args, kwargs)
      serialized = []
      [[_get_reader_for_worker_spec(), 'reader'], [target, 'target'], [args, 'args'], [kwargs, 'kwargs']].each do |obj, objname|
        begin
          serialized << serializer.dumps(obj)
        rescue => exception
          msg = "Could not serialize #{objname} #{obj} with #{serializer.__name__}."
          if serializer != dill
              msg += ' Try installing `dill`.'
          end
          raise PyteomicsError.new(msg)
        end
      end
      serialized
    end
  
    def _spawn_workers(specifications, in_queue, out_queue, processes)
      reader_spec, target_spec, args_spec, kwargs_spec = specifications
      workers = []
      processes.times do
        worker = FileReadingProcess.new(
          reader_spec, target_spec, in_queue, out_queue, args_spec, kwargs_spec)
        workers << worker
      end
      workers
    end
  
    def _spawn_feeder_thread(in_queue, iterator, processes)
      def feeder()
        iterator.each do |key|
          in_queue << key
        end
        processes.times do
          in_queue << nil
        end
      end
  
      feeder_thread = threading.Thread(target=feeder)
      feeder_thread.daemon = true
      feeder_thread.start()
      feeder_thread
    end
  
    def map(**_kwargs)
      target = _kwargs['target'] || nil
      processes = _kwargs['processes'] || -1
      args = _kwargs['args'] || nil
      kwargs = _kwargs['kwargs'] || nil
  
      if @_offset_index.nil?
        raise PyteomicsError.new('The reader needs an index for map() calls. Create the reader with `use_index=True`.')
      end
  
      if processes < 1
        processes = @_nproc
      end
      iterator = _task_map_iterator()
  
      if args.nil?
        args = []
      else
        if args.is_a?(Hash)
          args = args.to_a
        elsif args.is_a?(String)
          args = args.split('')
        else
          args = args.to_a
        end
      end
      if kwargs.nil?
        kwargs = {}
      else
        kwargs = kwargs.to_h
      end
      kwargs.merge!(_kwargs)
  
      serialized = _build_worker_spec(target, args, kwargs)
  
      # in_queue = mp.Queue(self._queue_size)
      # out_queue = mp.Queue(self._queue_size)
      in_queue = []
      out_queue = []
  
      workers = _spawn_workers(serialized, in_queue, out_queue, processes)
      feeder_thread = _spawn_feeder_thread(in_queue, iterator, processes)
      # for worker in workers:
      #   worker.start()
  
      def iterate()
        while true
          begin
            result = out_queue.get(True, @_queue_timeout)
            yield result
          rescue => exception
            if workers.map{ |w| w.is_done }.all?
              break
            else
              next
            end
        
          end
        end
        feeder_thread.join()
        workers.each do |worker|
          worker.join()
        end
      end
      return iterate()
    end
  
    def _task_map_iterator
      @_offset_index.keys.to_enum
    end
  end
  
  class ChainBase
    def initialize(...)
      __init__(...)
    end
  
    def __init__(*sources, **kwargs)
      @sources = sources
      @kwargs = kwargs
      @_iterator = nil
    end
  
    def self.from_iterable(cls, sources, **kwargs)
      cls.new(*sources, **kwargs)
    end
  
    def self._make_chain(sequence_maker, *args)
      if sequence_maker.instance_of?(Class)
        tp = Class.new(self) {
          def sequence_maker
            sequence_maker
          end
          define_method(sequence_maker.to_s, instance_method(:sequence_maker))
        }
      else
        tp = Class.new(self)
        tp.define_singleton_method(sequence_maker.to_s, args[0])
      end
      return tp
    end
  
    def sequence_maker(file)
      raise NotImplementedError.new
    end
  
    def _create_sequence(file)
      sequence_maker(file, **@kwargs)
    end
  
    def _iterate_over_series
      @sources.each do |f|
        _create_sequence(f) do |r|
          r.each do |item|
            yield item
          end
        end
      end
    end
  
    def __enter__
      @_iterator = _iterate_over_series().to_enum
      self
    end
  
    def __exit__(*args, **kwargs)
      @_iterator = nil
    end
  
    def __iter__
      self
    end
  
    def __next__
      if @_iterator.nil?
        @_iterator = _iterate_over_series()
      end
      _next(@_iterator)
    end
  
    def _next
      __next__()
    end
  
    def map(**_kwargs)
      target = _kwargs['target'] || nil
      processes = _kwargs['processes'] || -1
      queue_timeout = _kwargs['queue_timeout'] || _QUEUE_TIMEOUT
      args = _kwargs['args'] || nil
      kwargs = _kwargs['kwargs'] || nil
      @sources.each do |f|
        _create_sequence(f) do |r|
          r.map(target, processes, queue_timeout, args, kwargs, **_kwargs).each do |result|
            yield result
          end
        end
      end
    end
  end
  
  class TableJoiner < ChainBase
    def concatenate(results)
      if results.map{ |a| a.is_a?(Pandas.DataFrame) }.all?
        return Pandas.concat(results)
      end
      if results[0].is_a?(Numpy.ndarray)
        return Numpy.concatenate(results)
      else
        return Numpy.array(results.map{ |a| a.map{ |b| b } })
      end
    end
  
    def _iterate_over_series
      results = @sources.map{ |f| _create_sequence(f) }
      return concatenate(results)
    end
  end
end
