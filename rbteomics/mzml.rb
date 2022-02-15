require 'rexml/document'
require_relative 'xml'
require_relative 'auxiliary/file_helpers'
require_relative 'auxiliary/target_decoy'
require_relative 'auxiliary/utils'
require_relative '_schema_defaults'
require 'set'
require 'numpy'

module Mzml
  NON_STANDARD_DATA_ARRAY = 'non-standard data array'

  STANDARD_ARRAYS = Set.new([
    'm/z array',
    'intensity array',
    'charge array',
    'signal to noise array',
    'time array',
    'wavelength array',
    'flow rate array',
    'pressure array',
    'temperature array',
    'mean charge array',
    'resolution array',
    'baseline array',
    'noise array',
    'sampled noise m/z array',
    'sampled noise intensity array',
    'sampled noise baseline array',
    'ion mobility array',
    'deconvoluted ion mobility drift time array',
    'deconvoluted inverse reduced ion mobility array',
    'deconvoluted ion mobility array',
    'raw ion mobility drift time array',
    'raw inverse reduced ion mobility array',
    'raw ion mobility array',
    'mean inverse reduced ion mobility array',
    'mean ion mobility array',
    'mean ion mobility drift time array',
    'mass array',
    'scanning quadrupole position lower bound m/z array',
    'scanning quadrupole position upper bound m/z array',
  ])

  class MzML < Xml::XML
    # (xml.ArrayConversionMixin, aux.TimeOrderedIndexedReaderMixin, xml.MultiProcessingXML, xml.IndexSavingXML)

    @@file_format = 'mzML'
    @@_root_element = 'mzML'
    @_default_schema = Mzml_schema_defaults
    @@_default_version = '1.1.0'
    @@_default_iter_tag = 'spectrum'
    @@_structures_to_flatten = Set.new(['binaryDataArrayList', 'referenceableParamGroupRef'])
    @@_indexed_tags = Set.new(['spectrum', 'chromatogram'])

    def initialize(...)
      __init__(...)
    end

    def __init__(*args, **kwargs)
      @decode_binary = kwargs.delete('decode_binary') || true
      @acm = Xml::ArrayConversionMixin.new(*args, **kwargs)
      @tim = TimeOrderedIndexedReaderMixin.new(*args, **kwargs)
      @mpx = Xml::MultiProcessingXML.new(*args, **kwargs)
      @isx = Xml::IndexSavingXML.new(*args, **kwargs)
    end
  end

  
  
  
  
  
  # def _decode_peaks(info, peaks_data)
  #   compressed = info['compressionType'] == 'zlib'
  #   dt = info['precision'] == '32' ? Numpy.float32 : Numpy.float64
  #   dtype = Numpy.dtype([PyCall::Tuple.(['m/z array', dt]), PyCall::Tuple.(['intensity array', dt])]).newbyteorder('>')
  #   data = _decode_base64_data_array(peaks_data, dtype, compressed)
  # end

  # class IteratorQueue
  #   def initialize(...)
  #     __init__(...)
  #   end

  #   def __init__(iterator)
  #     @queue = PriorityQueue.new
  #     @iterator = iterator
  #     @last_index = -1
  #     @producer = consume(iterator)
  #   end

  #   def insert_item(scan)
  #     @queue << [scan['num'].to_i, scan]
  #   end

  #   def __iter__
  #     @producer
  #   end

  #   def consume(iterator)
  #     iterator.each do |scan|
  #       scan.delete("scan")
  #       if scan['msLevel'] != 1
  #         insert_item(scan)
  #       else
  #         insert_item(scan)
  #         barrier = scan['num'].to_i
  #         while true
  #           idx, item = @queue.pop
  #           if idx >= barrier
  #             insert_item(item)
  #             break
  #           end
  #           yield item
  #         end
  #       end
  #     end
  #     while @queue.empty?.!
  #       idx, item = @queue.pop
  #       yield item
  #     end
  #   end
  # end

  # class MzXML < Xml::XML
  #   # (xml.ArrayConversionMixin, aux.TimeOrderedIndexedReaderMixin, xml.MultiProcessingXML, xml.IndexSavingXML)
  #   @@_root_element = 'mzXML'
  #   @@_default_iter_tag = 'scan'
  #   @@_indexed_tags = Set.new(['scan'])
  #   @@_indexed_tag_keys = {'scan' => 'num'}
  #   @@_default_version = nil
  #   @@_default_schema = Mzxml_schema_defaults
  #   @_default_id_attr = 'num'

  #   def initialize(...)
  #     __init__(...)
  #   end

  #   def __init__(*args, **kwargs)
  #     @decode_binary = kwargs.delete('decode_binary') || true
  #     @acm = Xml::ArrayConversionMixin.new(*args, **kwargs)
  #     @tim = TimeOrderedIndexedReaderMixin.new(*args, **kwargs)
  #     @mpx = Xml::MultiProcessingXML.new(*args, **kwargs)
  #     @isx = Xml::IndexSavingXML.new(*args, **kwargs)
  #   end

  #   def __getstate__
  #     [@acm, @tim, @mpx, @isx].each do |obj|
  #       if obj.respond_to?(:__getstate__)
  #         state = obj.__getstate__
  #         break
  #       end
  #     end
  #     state['decode_binary'] = @decode_binary
  #     state
  #   end

  #   def __setstate__(state)
  #     [@acm, @tim, @mpx, @isx].each do |obj|
  #       if obj.respond_to?(:__setstate__)
  #         state = obj.__setstate__(state)
  #         break
  #       end
  #     end
  #     @decode_binary = state['decode_binary']
  #   end

  #   def _get_info_smart(element, **kw)
  #     name = Xml._local_name(element)

  #     kwargs = dict(kw)
  #     rec = kwargs.delete('recursive') || nil
  #     if 'mzXML'.include?(name)
  #       info = _get_info(element, 'recursive' => rec.nil?.! ? rec : false, **kwargs)
  #     else
  #       info = _get_info(element, 'recursive' => rec.nil?.! ? rec : true, **kwargs)
  #     end
  #     if info.include?('num') && info.is_a?(Hash)
  #       info['id'] = info['num']
  #     end
  #     if info.include?('peaks') && info.is_a?(Hash)
  #       _decode_peaks(info)
  #     end
  #     info
  #   end

  #   def _determine_compression(info)
  #     if info['compressionType'] == 'zlib'
  #       return 'zlib compression'
  #     end
  #     return "no compression"
  #   end

  #   def _determine_dtype(info)
  #     dt = info['precision'] == '32' ? Numpy.float32 : Numpy.float64
  #     endianess = ['network', 'big'].include?(info['byteOrder']) ? ">" : "<"
  #     dtype = Numpy.dtype([PyCall::Tuple.(['m/z array', dt]), PyCall::Tuple.(['intensity array', dt])]).newbyteorder(endianess)
  #   end

  #   def _finalize_record_conversion(array, record)
  #     key = record.key
  #     _convert_array(key, array[key])
  #   end

  #   def _decode_peaks(info)
  #     if info['peaks'].is_a?(Hash).! && info['peaks'].is_a?(Array).!
  #       compression_type = _determine_compression(info)
  #       dtype = _determine_dtype(info)
  #       binary = info.delete('peaks')
  #       if @decode_binary.!
  #         @_array_keys.each do |k|
  #           record = _make_record(binary, compression_type, dtype, k)
  #           info[k] = record
  #         end
  #       else
  #         peak_data = decode_data_array(binary, compression_type, dtype)
  #         @_array_keys.each do |k|
  #           info[k] = _convert_array(k, peak_data[k])
  #         end
  #       end
  #     else
  #       if @decode_binary.!
  #         arrays = info.delete('peaks')[0]
  #         @_array_keys.each do |k|
  #           info[k] = arrays[k]
  #         end
  #       else
  #         peak_data = info.delete('peaks')[0]
  #         @_array_keys.each do |k|
  #           info[k] = _convert_array(k, peak_data[k] || Numpy.array([]))
  #         end
  #       end
  #     end
  #   end

  #   def iterfind(path, **kwargs)
  #     if path == 'scan'
  #       [@acm, @tim, @mpx, @isx].each do |obj|
  #         if obj.respond_to?(:iterfind)
  #           generator = obj.iterfind(path, **kwargs)
  #           break
  #         end
  #       end
  #       IteratorQueue.new(generator).each do |item|
  #         yield item
  #       end
  #     else
  #       [@acm, @tim, @mpx, @isx].each do |obj|
  #         if obj.respond_to?(:iterfind)
  #           obj.iterfind(path, **kwargs).each do |item|
  #             yield item
  #           end
  #           break
  #         end
  #       end
  #     end
  #   end

  #   def next(...)
  #     iterfind(...)
  #   end

  #   def _get_time(scan)
  #     scan['retentionTime']
  #   end
  # end

  # def read(source, read_schema: false, iterative: true, use_index: false, dtype: nil, huge_tree: false)
  #   MzXML.new(source, 'read_schema' => read_schema, 'iterative' => iterative,
  #     'use_index' => use_index, 'dtype' => dtype, 'huge_tree' => huge_tree)
  # end

  # def iterfind(source, path, **kwargs)
  #   MzXML.new(source, **kwargs).iterfind(path, **kwargs)
  # end

  # @@version_info = Xml._make_version_info(MzXML)

  # # chain = aux._make_chain(read, 'read')
  # @@chain = ChainBase._make_chain(MzXML)
end
