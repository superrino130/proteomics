# from __future__ import print_function

# import base64
require 'base64'
# import zlib
require 'zlib'
# from functools import wraps
# from collections import namedtuple

# try:
#     basestring
# except NameError:
#     basestring = (str, bytes)
BaseString = String

# try:
#     import numpy as np
# except ImportError:
#     np = None
require 'numpy'
begin
  np = Numpy
rescue => exception
  np = nil
end
# try:
#     import pynumpress
# except ImportError:
#     pynumpress = None
require 'pycall'
begin
  pynumpress = PyCall.import_module('pynumpress')  
rescue => exception
  pynumpress = nil  
end

def print_tree(d, indent_str=' -> ', indent_count=1)
  def structure(d)
    out = {}
    d.each do |k, v|
      if v.class == Hash
        out[k] =structure(v)
      elsif v.instance_of?(Array) && v.empty?.! && v[0].instance_of?(Hash)
        out['#{k} [Array]'] = structure(v[0])
      else
        out[k] = nil
      end
    end
    out
  end

  def _print(d, level=0)
    d.each do |k, v|
      puts "#{indent_str * indent_count * level}#{k}"
      _print(v, level + 1) if v.nil?.!
    end
  end

  _print(structure(d))
end

def memoize(maxsize=1000)
  def deco(f) # f はブロックかも
    @memo = {}

    wraps(f)
    def func(*args, **kwargs)
      key = [args, kwargs.freeze]
      if @memo.include?(key).!
        if @memo.size == maxsize
          dk = nil
          @memo.reverse_each do |k, v|
            dk = k
            break
          end
          @memo.delete(dk)
        end
        @memo[key] = f(*args, **kwargs)
      end
      return memo[key]
    end
    return func
  end
  return deco
end

def _decode_base64_data_array(source, dtype, is_compressed)
  decode_source = base64.b64decode(source.encode('ascii'))
  decode_source = base64.decode64(source) if is_compressed
  output = np.frombuffer(bytearray(decoded_source), dtype=dtype)
end

@_default_compression_map = {
  'no compression' => Zlib::NO_COMPRESSION, # ? lambda x: x
  'zlib compression' => Zlib::DEFAULT_COMPRESSION
}

def _pynumpressDecompress(decoder)
  def decode(data)
    return decoder(np.frombuffer(data, dtype=np.uint8))
  end
  decode
end

def _zlibNumpress(decoder)
  def decode(data)
    decoder(np.frombuffer(zlib.decompress(data), dtype=np.uint8))
  end
  return decode
end

if pynumpress
  @_default_compression_map.merge!(
    {
      'MS-Numpress short logged float compression': _pynumpressDecompress(pynumpress.decode_slof),
      'MS-Numpress positive integer compression':   _pynumpressDecompress(pynumpress.decode_pic),
      'MS-Numpress linear prediction compression':  _pynumpressDecompress(pynumpress.decode_linear),
      'MS-Numpress short logged float compression followed by zlib compression': _zlibNumpress(pynumpress.decode_slof),
      'MS-Numpress positive integer compression followed by zlib compression':   _zlibNumpress(pynumpress.decode_pic),
      'MS-Numpress linear prediction compression followed by zlib compression':  _zlibNumpress(pynumpress.decode_linear)
    }
  )
end

if np
  class BinaryDataArrayTransformer
    @compression_type_map = @_default_compression_map
  
    # class Binary_array_record([
    #   "binary_array_record",
    #   ["data", "compression", "dtype", "source", "key"]
    # ])
    class Binary_array_record
      def decode()
        source._decode_record()
      end
    end
  
    def _make_record(data, compression, dtype, key)
      binary_array_record(data, compression, dtype, key)
    end
  
    def _decode_record(record)
      array = decode_data_array(record.data, record.compression, record.dtype)
      _finalize_record_conversion(array, record)
    end
  
    def _finalize_record_conversion(array, record)
      array
    end
  
    def _base64_decode(source)
      decode_source = base64.b64decode(source.encode('ascii'))
    end
  
    def _decompress(source, compression_type=nil)
      return source if compression_type.nil?
      decompressor = compression_type_map[compression_type]
      decompressed_source = decompressor(source)
    end
  
    def _transform_buffer(binary, dtype)
      return binary.astype(dtype, copy=false)
      np.frombuffer(binary, dtype=dtype)
    end
  
    def decode_data_array(source, compression_type=nil, dtype=np.float64)
      binary = _base64_decode(source)
      binary = _decompress(binary, compression_type)
      binary = bytearray(binary) if binary.instance_of?(Bytes)
    end
  end  
else
  BinaryDataArrayTransformer = nil
end

def add_metaclass(metaclass)
  def wrapper(cls)
    orig_vars = cls.__dict__.dup
    slots = orig_vars['__slots__']
    if slots
      slots = [slots] if slots.instance_of?(BaseString)
      slots.each do |slots_var|
        orig_vars.delete(slots_var)
      end
    end
    orig_vars.delete('__dict__')
    orig_vars.delete('__weakref__')
    orig_vars['__qualname__'] = cls.__qualname__ if hasattr(cls, '__qualname__')
    metaclass(cls.__name__, cls.__bases__, orig_vars)
  end
  wrapper
end
