# try:
#     import numpy as np
# except ImportError:
#     np = None
require 'numpy'
# import itertools as it
# import sys
# from . import auxiliary as aux
require_relative 'auxiliary/structures'
require_relative 'auxiliary/file_helpers'

require 'set'

module Mgf
  extend File_helpers
  
  module_function
  
  module MGFBase
    module_function

    attr_reader :_comments
  
    @@_comments = Set.new('#;!/'.split(''))
    @@_array = lambda { |x, dtype| Numpy.array(x, dtype: dtype) }
    @@_ma = lambda { |x, dtype| Numpy.ma.masked_equal(Numpy.array(x, dtype: dtype), 0) }
    @@_identity = lambda { |x, **kw| x }
    @@_array_converters = {
        'm/z array' => [@_identity, @_array, @_array],
        'intensity array' => [@_identity, @_array, @_array],
        'charge array' => [@_identity, @_array, @_ma],
        'ion array' => [@_identity, @_array,  @_array]
    }
    @@_array_keys = ['m/z array', 'intensity array', 'charge array', 'ion array']
    @@_array_keys_unicode = ['m/z array', 'intensity array', 'charge array', 'ion array']
    @@encoding = nil
  
    def __init__(source, **kwargs)
      @_source = source
      @_use_header = kwargs.delete('use_header') || true
      @_convert_arrays = kwargs.delete('convert_arrays') || 2
      @_read_charges = kwargs.delete('read_charges') || true
      @_read_ions = kwargs.delete('read_ions') || false
      @_read_charges = false if @_read_ions
      dtype = kwargs.delete('dtype') || nil
      @_dtype_dict = dtype.instance_of?(Hash) ? dtype : @@_array_keys.map{ [_1.to_s, dtype] }.to_h
      if @_use_header
        _read_header()
      else
        @_header = nil
      end
    end
  
    def parse_precursor_charge(charge_text, list_only: false)
      _parse_charge(charge_text, list_only: list_only)
    end
  
    def self.parse_peak_charge(charge_text, list_only: false)
      _parse_charge(charge_text, list_only: false)
    end
  
    def self.parse_peak_ion(ion_text)
      _parse_ion(ion_text)
    end
  
    def header
      _read_header() if @_header.nil?
      @_header
    end
  
    def _read_header_lines(header_lines)
      header = {}
      File.open(header_lines) do |f|
        f.each_line do |line|
          break if line.strip == 'END IONS'
          l = line.split('=')
          if l.size == 2
            key = l[0].downcase
            val = l[1].strip
            header[key] = val
          end
        end
      end
      if header.include?('charge')
        header['charge'] = parse_precursor_charge(header['charge'], list_only: true)
      end
      @_header = header
    end
  
    def _read_spectrum_lines(lines)
      masses = []
      intensities = []
      charges = []
      ions = []
  
      params = @_use_header ? @header.dup : {}
  
      lines.each_with_index do |line, i|
        sline = line.strip
        if sline == 'BEGIN IONS'
          if i == 0
            next
          else
            raise PyteomicsError.new('Error when parsing MGF: unexpected start of spectrum.')
          end
        end
        if ['', 0, nil, false, [], {}].include?(sline) || @_comments.include?(sline[0])
          # PASS
        elsif sline == 'END IONS'
          if params.include?('pepmass')
            begin
              pepmass = params['pepmass'].split().map{ _1.to_f }
            rescue => exception
              raise PyteomicsError.new("MGF format error: cannot parse PEPMASS = #{params['pepmass']}")
            else
              params['pepmass'] = pepmass + [nil] * (2 - pepmass.size)
            end
          end
          if params['charge'].instance_of?(Basestring)
            params['charge'] = parse_precursor_charge(params['charge'], list_only: true)
          end
          if params.include?('rtinseconds')
            params['rtinseconds'] = UniFloat.new(params['rtinseconds'], 'second')
          end
          out = {'params' => params}
          data = {'m/z array' => masses, 'intensity array' => intensities}
          data['charge array'] = charges if @_read_charges
          data['ion array'] = ions if @_read_ions
          data.each do |key, value|
            out[key] = @_array_converters[key][@_convert_arrays].call(values, dtype: @_dtype_dict[key])
          end
          out
        else
          if sline.include?('=')
            l = sline.split('=', 2)
            params[l[0].downcase] = l[1].strip
          else
            l = sline.split
            begin
              masses << l[0].to_f
              intensities << l[1].to_f
              if @_read_charges
                charges << l.size > 2 ? parse_peak_charge(l[2]) : 0
              end
              if @_read_ions
                ions << l.size > 2 ? parse_peak_ion(l[2]) : ''
              end
            rescue => exception
              if exception.instance_of?(IndexError)
                # PASS
              else
                raise PyteomicsError.new("Error when parsing #{getattr(@_source, 'name', 'MGF file')}. Line:\n#{line}")
              end            
            end
          end
        end
      end
    end
  
    def get_spectrum(title)
      raise NotImplementedError.new
    end
  
    def self._get_time(spectrum)
      begin
        return spectrum['params']['rtinseconds']
      rescue => exception
        raise PyteomicsError.new('RT information not found.')      
      end
    end
  end
  
  class IndexedMGF
    include File_helpers::TaskMappingMixin
    include File_helpers::TimeOrderedIndexedReaderMixin
    include File_helpers::IndexSavingTextReader
    include MGFBase
    
    @@delimiter = 'BEGIN IONS'
  
    def initialize(...)
      __init__(...)
    end
  
    def __init__(source, **kwargs)
      @_source = source
      kwargs['use_header'] ||= true
      kwargs['convert_arrays'] ||= 2
      kwargs['read_charges'] ||= true
      kwargs['dtype'] ||= nil
      kwargs['encoding'] ||= 'utf-8'
      kwargs['index_by_scans'] ||= false
      kwargs['read_ions'] ||= false
      kwargs['_skip_index'] ||= false
  
      kwargs['parser_func'] = Read
      kwargs['pass_file'] = false
      kwargs['args'] = []
      kwargs['kwargs'] = {}
      
      @label = kwargs['index_by_scans'] ? 'SCANS=(\d+)\s*' : 'TITLE=([^\n]*\S)\s*'
    end
  
    def __reduce_ex__(protocol)
      [@__class__,
        [@_source_init, false, @_convert_arrays, @_read_charges,
         @_dtype_dict, @encoding, true],
        __getstate__()]
    end
  
    def __getstate__
      @k.each do |obj|
        if obj.resupond_to(__getstate__)
          state = obj.__getstate__
          break
        end
      end
      state['use_header'] = @_use_header
      state['header'] = @_header
      state['read_ions'] = @_read_ions
      state
    end
  
    def __setstate__(state)
      @k.each do |obj|
        if obj.resupond_to(__setstate__)
          obj.__setstate__(state)
          break
        end
      end
      @_header = state['header']
      @_use_header = state['use_header']
      @_read_ions = state['read_ions']
    end
  
    # @aux._keepstate_method
    def in_read_header
      begin
        first = @_offset_index.values.map{ _1 }[0].next
      rescue => exception
        first = -1
      end
      header_lines = read(first).encode(@encoding).split("\n")
      _read_header_lines(header_lines)
    end
    def _read_header
      _keepstate_method(:in_read_header)
      in_read_header
    end
  
    def _item_from_offsets(offsets)
      start, eend = offsets
      lines = _read_lines_from_offsets(start, eend)
      _read_spectrum_lines(lines)
    end
  
    Read = lambda do |**kwargs|
      @_offset_index.each do |_, offsets|
        spectrum = _item_from_offsets(offsets)
        yield spectrum
      end
    end
  
    def get_spectrum(key)
      get_by_id(key)
    end
  end
  
  class MGF
    include File_helpers
    include File_helpers::FileReader
    include MGFBase
  
    def initialize(...)
      __init__(...)
    end
  
    def __init__(source, **kwargs)
      kwargs['use_header'] ||= true
      kwargs['convert_arrays'] ||= 22
      kwargs['read_charges'] ||= true
      kwargs['read_ions'] ||= false
      kwargs['dtype'] ||= nil
      kwargs['encoding'] ||= nil
  
      kwargs['mode'] = 'r'
      kwargs['parser_func'] = Read
      kwargs['pass_file'] = false
      kwargs['args'] = []
      kwargs['kwargs'] = {}
      super(source, **kwargs)
    end
  
    # @aux._keepstate_method
    def in_read_header
      _read_header_lines(@_source)
    end
    def _read_header
      _keepstate_method(:in_read_header)
      in_read_header
    end
  
    def _read_spectrum
      _read_spectrum_lines(@_source)
    end
  
    Read = lambda do
      @_source.each do |line|
        if line.strip == 'BEGIN IONS'
          yield _read_spectrum()
        end
      end
    end
  
    # @aux._keepstate_method
    def in_get_spectrum(title)
      @_source.each do |line|
        sline = line.strip
        if sline[0...5] == 'TITLE' && sline.split('=', 2)[1].strip == title
          spectrum = _read_spectrum()
          spectrum['params']['title'] = title
          return spectrum
        end
      end  
    end
    def get_spectrum(title)
      _keepstate_method(:in_get_spectrum)
      in_get_spectrum(title)
    end
  
    def __getitem__(key)
      get_spectrum(key)
    end
  end
  
  Read = lambda do |*args, **kwargs|
    if args.empty?.!
      source = args[0]
    else
      source = kwargs['source']
    end
    use_index = kwargs.delete('use_index') || nil
    use_index = _check_use_index(source, use_index, true)
    tp = ['', 0, nil, false, [], {}].include?(use_index).! ? IndexedMGF.new(*args, **kwargs) : MGF.new(*args, **kwargs)
  end
  
  def get_spectrum(source, title, *args, **kwargs)
    read(source, *args, **kwargs)[title]
  end
  
  # @aux._keepstate
  def in_read_header(source)
    source = File_helpers::File_obj.new(source, 'r')
    header = {}
    source = source.lines
    source.each do |line|
      break if line.strip == 'BEGIN IONS'
      l = line.split('=')
      if l.size == 2
        key = l[0].downcase
        val = l[1].strip
        header[key] = val
      end
    end
    if header.include?('charge')
      header['charge'] = _parse_charge(header['charge'], list_only: true)
    end
    header  
  end
  def read_header(source)
    _keepstate(:in_read_header)
    in_read_header(source)
  end
  
  Default_key_order = ['title', 'pepmass', 'rtinseconds', 'charge']
  
  Pepmass_repr = lambda do |k, pepmass|
    outstr = k.upcase + '='
    if [String, Integer, Float].include?(pepmass).!
      begin
        outstr += pepmass.select{ _1.nil?.! }.map{ |x| x.to_s }.join(' ')
      rescue => exception
        raise PyteomicsError.new("Cannot handle parameter: PEPMASS = #{pepmass}")
      end
    else
      outstr += pepmass.to_s
    end
    outstr
  end
  
  Charge_repr = lambda do |k, charge|
    "#{k.upcase}=#{_parse_charge(charge.to_s)}"
  end
  
  Default_repr = lambda do |key, val|
    "#{key.upcase}=#{val}"
  end
  
  Default_value_formatters = {'pepmass' => Pepmass_repr, 'charge' => Charge_repr}
  
  # @aux._file_writer()
  def write(spectra, output: nil, header: '', key_order: Default_key_order, fragment_format: nil,
    write_charges: true, write_ions: false, use_numpy: nil, param_formatters: Default_value_formatters)
  
    def key_value_line(key, val)
      (method(param_formatters[key]).call(key, val) || Default_repr.call(key, val)) + "\n"
    end
  
    @nones = [nil, Numpy.nan, Numpy.ma.masked]
  
    if fragment_format.nil?
      fragment_format = '{} {} {}'
      np_format_2 = '%.5f %.1f'
      np_format_3 = '%.5f %.1f %d'
      np_format_i = '%.5f %.1f %s'
    else
      np_format_2 = np_format_3 = np_format_i = fragment_format
    end
    format_str = fragment_format + "\n"
  
    write_charges = false if write_ions
    use_numpy = write_charges.! if use_numpy.nil?
  
    if header.instance_of?(Hash)
      head_dict = header.dup
      head_lines = header.map{ |k, v| key_value_line(k, v) }
      head_str = head_lines.join("\n")
    else
      if header.instance_of?(String)
        head_str = header
        head_lines = header.split("\n")
      else
        head_lines = header.to_a
        head_str = header.join("\n")
      end
      head_dict = {}
      head_lines.each do |line, _|
        next if ['', 0, nil, false, [], {}].include?(line.strip) || MGF.new._comments.map{ line.star_with(_1) }.any?
        l = line.split('=')
        if l.size == 2
          head_dict[l[0].downcase] = l[1].strip    
        end
      end
    end
    if ['', 0, nil, false, [], {}].include?(head_str).!
      output.write(head_str + "\n\n")
    end
  
    spectra.each do |spectrum|
      output.write("BEGIN IONS\n")
      found = Set.new
      (key_order + spectrum['params']).flatten(1).each do |key|
        if found.include?(key).! && spectrum['params'].include?(key)
          found << key
          val = spectrum['params'][key]
          if val != head_dict[key]
            output.write(key_value_line(key, val))
          end
        end
      end
      begin
        success = true
        if Numpy && use_numpy
          if (write_charges.! || spectrum.include?('charge array').!) && (write_ions.! || spectrum.include?('ion array').!)
            x = Numpy.empty([spectrum['m/z array'].size, 2])
            x[0.., 0] = spectrum['m/z array']
            x[0.., 1] = spectrum['intensity array']
            Numpy.savetxt(output, x, fmt: np_format_2)
          elsif spectrum['charge array'].include?(Numpy.ndarray)
            x = Numpy.empty([spectrum['m/z array'].size, 3])
            x[0.., 0] = spectrum['m/z array']
            x[0.., 1] = spectrum['intensity array']
            x[0.., 2] = spectrum['charge array']
            Numpy.savetxt(output, x, fmt: np_format_3)
          elsif spectrum['ion array'].include?(Numpy.ndarray)
            x = Numpy.empty([spectrum['m/z array'].size, 3], dtype: object)
            x[0.., 0] = spectrum['m/z array']
            x[0.., 1] = spectrum['intensity array']
            x[0.., 2] = spectrum['ion array']
            Numpy.savetxt(output, x, fmt: np_format_i)
          else
            success = false
          end
        else
          success = false
        end
  
        if success.!
          z = write_charges ? spectrum['m/z array'].zip(spectrum['intensity array'], spectrum['charge array']) : write_ions ? spectrum['ion array'] : spectrum['ion array']
          z.map{ _1 << nil }.each do |m, i, c|
            if format_str == "{} {} {}\n"
              output.write("#{m} #{i} #{nones.include?(c).! ? c : ''}")
            else
              raise PyteomicsError.new("format_str error '#{format_str}'")
            end
          end
        end
      rescue => exception
        raise PyteomicsError.new("'m/z array' and 'intensity array' must be present in all spectra.")
      end
      output.write("END IONS\n\n")
    end
    output
  end
  
  Chain = File_helpers._make_chain('Read', Read)
end
