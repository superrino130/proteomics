# import warnings
# from . import auxiliary as aux
# from . import xml, _schema_defaults
require_relative 'xml'
require_relative '_schema_defaults'
require_relative 'auxiliary/file_helpers'
require_relative 'auxiliary/target_decoy'
require 'set'

class MzIdentML
  # (xml.MultiProcessingXML, xml.IndexSavingXML)
  # Not started
end

def read(source, **kwargs)
  # Not started
end

def iterfind(source, path, **kwargs)
  # Not started
end

# Version_info = _make_version_info(MzIdentML)

def get_by_id(source, elem_id, **kwargs)
  # Not started
end

Chain = ChainBase._make_chain(MzIdentML)

@is_decoy = lambda do |psm, prefix: nil|
  # Not started
end

def dataframe(*args, **kwargs)
  # Not started
end

def filter_df(*args, **kwargs)
  # Not started
end

fdr = _make_fdr(@is_decoy, nil)
@_key = lambda { |x| x['SpectrumIdentificationItem'].min{ |sii| sii['mascot:expectation value'] } }
@qvalues = _make_qvalues(Chain, @is_decoy, nil, @_key)
filter = lambda do |x = nil|
  if x.nil?
    _make_filter(Chain, @is_decoy, nil, @_key, @qvalues)
  elsif x == 'chain' || x == :chain
    y = _make_filter(Chain, @is_decoy, nil, @_key, @qvalues)
    _make_chain(y, 'filter', true)
  end
end
# filter = _make_filter(Chain, @is_decoy, nil, @_key, @qvalues)
# filter.chain = _make_chain(filter, 'filter', true)
