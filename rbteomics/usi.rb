# import json
# import warnings
# import threading
# import multiprocessing

# from collections import namedtuple, defaultdict

# try:
#     from multiprocessing.dummy import Pool as ThreadPool
# except ImportError:
#     ThreadPool = None

# try:
#     from urllib2 import Request, urlopen
# except ImportError:
#     from urllib.request import Request, urlopen

# try:
#     import numpy as np

#     def coerce_array(array_data):
#         return np.array([float(v) for v in array_data])

# except ImportError:

#     def coerce_array(array_data):
#         return [float(v) for v in array_data]

# from .auxiliary import PyteomicsError




require 'numpy'
require_relative 'auxiliary/structures'
require_relative 'auxiliary/file_helpers'
require_relative 'auxiliary/utils'

require 'set'
require 'open-uri'
require 'delegate'
require 'rexml/document'

module Usi
  module_function
  USI = Struct.new(:USI, :protocol, :dataset, :datafile, :scan_identifier_type, :scan_identifier, :interpretation) do
    include Usi
    def to_s
      self.select{ |x| x.nil?.! }.join(':')
    end

    def parse
      self.protocol, self.dataset, self.datafile, self.scan_identifier_type, self.scan_identifier, self.interpretation = _usi_parser(self.USI)
      self
    end
  end

  # def coerce_array(array_data)
  #   Numpy.array(array_data.map{ |v| v.to_f })
  # end

  def cast_numeric(value)
    if value.instance_of?(Integer) || value.instance_of?(Float)
      return value
    else
      begin
        return Float(value)
      rescue => exception
        # PASS
      end
      begin
        return Integer(value)
      rescue => exception
        raise exception
      end
    end
  end

  def _usi_parser(usi)
    tokens = usi.split(':', 6)
    protocol = tokens[0]
    dataset = tokens[1]
    datafile = tokens[2]
    scan_identifier_type = tokens[3]
    scan_identifier = tokens[4]
    interpretation = tokens[5] || nil
    [protocol, dataset, datafile, scan_identifier_type, scan_identifier, interpretation]
  end

  class PROXIBackend
    def initialize(...)
      __init__(...)
    end

    def __init__(name, url_template, **kwargs)
      kwargs['version'] = '0.1' if kwargs.include?('version').!
      @name = name
      @url_template = url_template
      @options = kwargs
    end

    def inspect
      "#{self.class}(#{@options})"
    end

    def _request(usi)
      url = @url_template.sub('{usi!s}', usi)
      @options.each do |k, v|
        url = url.sub("{#{k}}", v.to_s) if url.include?(k)
      end
      # req = Request(url)
      URI.open(url) do |response|
        if response.status[0] != '200'
          raise "PROXI Service Response Code #{response.status}"
        end
        data = response.read().decode('utf-8')
        data = json.loads(data)
        return data
      end
    end

    def get(usi)
      data = _request(usi)
      result = _coerce(data)
    end

    def _coerce(data)
      if data.instance_of?(Array)
        data = data[0]
      end
      result = {}
      result['attributes'] = data.delete('attributes') || []
      result['attributes'].each do |attrib|
        if attrib.include?('value') && attrib['value'].instance_of?(String) && attrib['value'][0].isdigit()
          begin
            attrib['value'] = cast_numeric(attrib['value'])
          rescue => exception
            # PASS
          end
        end
      end
      result['m/z array'] = coerce_array(data.delete('mzs') || [])
      result['intensity array'] = coerce_array(data.delete('intensities') || [])
      data.each do |key, value|
        if result.include?(key)
          raise "Attempting to set explicit value for #{key.inspect}"
        end
        result[key] = value
      end
      result
    end

    def __call__(usi)
      get(usi)
    end
  end

  class PeptideAtlasBackend < PROXIBackend
    def initialize(...)
      @_url_template ||= "http://www.peptideatlas.org/api/proxi/v{version}/spectra?resultType=full&usi={usi!s}"
      __init__(...)
    end

    def __init__(**kwargs)
      super('PeptideAtlas', @_url_template, **kwargs)
    end
  end

  class MassIVEBackend < PROXIBackend
    def initialize(...)
      @_url_template ||= "http://massive.ucsd.edu/ProteoSAFe/proxi/v{version}/spectra?resultType=full&usi={usi}"
      __init__(...)
    end

    def __init__(**kwargs)
      super('MassIVE', @_url_template, **kwargs)
    end
  end

  class PRIDEBackend < PROXIBackend
    def initialize(...)
      @_url_template ||= "http://wwwdev.ebi.ac.uk/pride/proxi/archive/v{version}/spectra?resultType=full&usi={usi}"
      __init__(...)
    end

    def __init__(**kwargs)
        super('PRIDE', @_url_template, **kwargs)
    end
  end

  class JPOSTBackend < PROXIBackend
    def initialize(...)
      @_url_template ||= 'https://repository.jpostdb.org/proxi/spectra?resultType=full&usi={usi}'
      __init__(...)
    end

    def __init__(**kwargs)
      super('jPOST', @@_url_template, **kwargs)
      kwargs.delete('version')
    end
  end

  class ProteomeExchangeBackend < PROXIBackend
    def initialize(...)
      @_url_template ||= 'http://proteomecentral.proteomexchange.org/api/proxi/v{version}/spectra?resultType=full&usi={usi!s}'
      __init__(...)
    end

    def __init__(**kwargs)
        super('ProteomeExchange', @_url_template, **kwargs)
    end
  end

  class PROXIAggregator
    def initialize(...)
      @_coalesce_resolution_methods = ["first"]
      __init__(...)
    end
    # Not started
  end

  @@proxies = {
    'peptide_atlas' => PeptideAtlasBackend,
    'massive' => MassIVEBackend,
    'pride' => PRIDEBackend,
    'jpost' => JPOSTBackend,
    'proteome_exchange' => ProteomeExchangeBackend
  }

  @@default_backend = 'peptide_atlas'

  AGGREGATOR_KEY = 'aggregator'
  AGGREGATOR = PROXIAggregator

  def proxi(usi, **kwargs)
    backend = kwargs['backend'] || @@default_backend

    if backend.instance_of?(String)
      if backend == AGGREGATOR_KEY
        AGGREGATOR.new(usi)
      elsif @@proxies.include?(backend)
        @@proxies[backend].new(**kwargs).get(usi)
      else
        raise PyteomicsError.new("Unknown PROXI backend name: #{backend}.")
      end
    elsif backend.instance_of?(Class) && ((backend.ancestors.include?(PROXIBackend)) || (backend.ancestors.include?(PROXIAggregator)))
      backend = backend.new(**kwargs).get(usi)
    elsif callable?(backend)
      backend(usi)
    else
      raise TypeError "Unrecognized backend type: #{backend.class}"
    end
  end
end
