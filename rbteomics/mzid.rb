# import warnings
# from . import auxiliary as aux
# from . import xml, _schema_defaults
require_relative 'xml'
require_relative '_schema_defaults'
require_relative 'auxiliary/file_helpers'
require_relative 'auxiliary/target_decoy'
require 'set'

module MZID
  class MzIdentML
    # (xml.MultiProcessingXML, xml.IndexSavingXML)
    # Not started
  
    def initialize(...)
      __init__(...)
    end
  
    def __init__(*args, **kwargs)
      @file_format = 'mzIdentML'
      @_root_element = 'MzIdentML'
      @_default_schema = Mzid_schema_defaults
      @_default_version = '1.1.0'
      @_default_iter_tag = 'SpectrumIdentificationResult'
      @_structures_to_flatten = ['Fragmentation'].to_set
      @_indexed_tags = ['SpectrumIdentificationResult', 'SpectrumIdentificationItem',
                       'SearchDatabase', 'SourceFile', 'SpectraData', 'Sample',
                       'DBSequence',  'Peptide', 'PeptideEvidence',
                       'Measure', 'TranslationTable', 'MassTable', 'Enzyme',
                       'Organization', 'AnalysisSoftware', 'BibliographicReference', 'Person', 'Provider',
                       'SpectrumIdentificationList', 'SpectrumIdentificationProtocol', 'SpectrumIdentification',
                       'ProteinDetectionList', 'ProteinDetectionProtocol', 'ProteinDetection',
                       'ProteinDetectionHypothesis', 'ProteinAmbiguityGroup',
      ].to_set
  
      @_element_handlers = XML::XML._element_handlers.dup
      @_element_handlers.merge!({
          "Modification": XML::XML::Promote_empty_parameter_to_name,
          "SpectrumIDFormat": XML::XML::Promote_empty_parameter_to_name,
          "FileFormat": XML::XML::Promote_empty_parameter_to_name,
          "Role": XML::XML::Promote_empty_parameter_to_name
      })
  
  
      @mpx = XML::MultiProcessingXML.new(*args, **kwargs)
      @isx = XML::IndexSavingXML.new(*args, **kwargs)
    end

    # def size
      
    # end
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
  
end
