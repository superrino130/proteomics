# import warnings
# from . import auxiliary as aux
# from . import xml, _schema_defaults
require_relative 'Xml'
require_relative '_schema_defaults'
require_relative 'auxiliary/file_helpers'
require_relative 'auxiliary/target_decoy'
require 'set'
require 'delegate'
require 'pandas'

module Mzid
  module_function
  
  class MzIdentML
    # (xml.MultiProcessingXML, xml.IndexSavingXML)
    @@file_format = 'mzIdentML'
    @@_root_element = 'MzIdentML'
    @@_default_schema = Mzid_schema_defaults
    @@_default_version = '1.1.0'
    @@_default_iter_tag = 'SpectrumIdentificationResult'
    @@_structures_to_flatten = ['Fragmentation'].to_set
    @@_indexed_tags = ['SpectrumIdentificationResult', 'SpectrumIdentificationItem',
                     'SearchDatabase', 'SourceFile', 'SpectraData', 'Sample',
                     'DBSequence',  'Peptide', 'PeptideEvidence',
                     'Measure', 'TranslationTable', 'MassTable', 'Enzyme',
                     'Organization', 'AnalysisSoftware', 'BibliographicReference', 'Person', 'Provider',
                     'SpectrumIdentificationList', 'SpectrumIdentificationProtocol', 'SpectrumIdentification',
                     'ProteinDetectionList', 'ProteinDetectionProtocol', 'ProteinDetection',
                     'ProteinDetectionHypothesis', 'ProteinAmbiguityGroup',
    ].to_set

    @@_element_handlers = Xml::XML._element_handlers.dup
    @@_element_handlers.merge!({
        "Modification": Xml::XML::Promote_empty_parameter_to_name,
        "SpectrumIDFormat": Xml::XML::Promote_empty_parameter_to_name,
        "FileFormat": Xml::XML::Promote_empty_parameter_to_name,
        "Role": Xml::XML::Promote_empty_parameter_to_name
    })

    def initialize(...)
      __init__(...)
    end
  
    def __init__(*args, **kwargs)
      kwargs['retrieve_refs'] = true if kwargs.include?('retrieve_refs').!
  
      @mpx = SimpleDelegator.new(Xml::MultiProcessingXML).new(*args, **kwargs)
      @isx = SimpleDelegator.new(Xml::IndexSavingXML).new(*args, **kwargs)
      @k = [@mpx, @isx]
    end

    def _get_info_smart(element, **kwargs)
      name = Xml._local_name(element)
      kwargs = kwargs.to_h
      rec = kwargs.delete("recursive") || nil

      _get_info(element, 'recursive' => (rec if rec.nil?.! ? rec : name != @_root_element), **kwargs)
    end

    def _retrieve_refs(info, **kwargs)
      info.to_h.each do |k, v|
        if k.end_with?('_ref')
          begin
            by_id = get_by_id(v, retrieve_refs: true)            
          rescue => exception
            warn 'Ignoring unresolved reference: ' + v            
          else
            info.merge!(by_id)
            info.delete(k)
            info.delete('id')
          end
        end
      end
    end

    def method_missing(method, ...)
      @k.each do |x|
        if x.respond_to?(method)
          return x.method(method).call(...)
        else
          if x.method_missing(method, ...) != 'dummy'
            return x.method_missing(method, ...)
          end
        end
      end
    end
  end
  
  def read(source, **kwargs)
    kwargs = kwargs.dup
    kwargs['retrieve_refs'] = true if kwargs.include?('retrieve_refs').!
    kwargs['build_id_cache'] = kwargs['build_id_cache'] || kwargs['retrieve_refs']
    MzIdentML.new(source, **kwargs)
  end
  
  def iterfind(source, path, **kwargs)
    kwargs = kwargs.dup
    kwargs['build_id_cache'] = kwargs['build_id_cache'] || kwargs['retrieve_refs']
    MzIdentML.new(source, **kwargs).iterfind(path, **kwargs)
  end
  
  # Version_info = _make_version_info(MzIdentML)
  
  def get_by_id(source, elem_id, **kwargs)
    MzIdentML.new(source, **kwargs).get_by_id(elem_id, **kwargs)
  end
  
  Chain = ChainBase._make_chain(MzIdentML)
  
  @is_decoy = lambda do |psm, prefix: nil|
    psm['SpectrumIdentificationItem'].all{ |sii| sii['PeptideEvidenceRef'].all{ |pe| ['', 0, nil, false, [], {}].include?(pe['isDecoy']).! } }
  end
  
  def dataframe(*args, **kwargs)
    data = []

    sep = kwargs.delete('sep') || nil
    Chain.new(*args, **kwargs) do |f|
      f.each do |item|
        info = {}
        item.each do |k, v|
          if [String, Integer, Float].include?(v.class)
            info[k] = v
          end
        end
        sii = item['SpectrumIdentificationItem'][0] || nil
        if sii.nil?.!
          info.merge!(sii.select{ |k, v| [String, Integer, Float].include?(v) }.to_h)
          evref = sii['PeptideEvidenceRef']
          if evref.nil?.! && evref.empty?.!
            prot_descr, accessions, isd, starts, ends, lengths = [], [], [], [], [], []
            evref.each do |d|
              prot_descr << d['protein description']
              accessions << d['accession']
              isd << d['isDecoy']
              starts << d['start']
              ends << d['end']
              lengths << d['length']
            end
            isd = isd.all?{ ['', 0, nil, false, [], {}].include?(_1).! }
            if sep.nil?.!
              if prot_descr.all?{ |prd| prd.instance_of?(String) }
                prot_descr = sep.join(prot_descr)
              end
              if accessions.all?{ |acc| acc.instance_of?(String) }
                accessions = sep.join(accessions)
              end
              if prot_descr.all?{ |prd| prd.nil? }
                prot_descr = nil
              end
              if accessions.all?{ |acc| acc.nil? }
                accessions = nil
              end

              info.merge!(evref[0].select{ |k, v| [String, Integer, Float, Array].include?(v.class) }.to_h)
              info['protein description'] = prot_descr
              info['accession'] = accessions
              info['isDecoy'] = isd
              info['start'] = starts
              info['end'] = ends
              info['length'] = lengths
            end
            data << info
          end
        end
      end
    end
    df = Pandas.DataFrame(data)
  end
  
  def filter_df(*args, **kwargs)
    kwargs['key'] = 'mascot:expectation value' if kwargs.include?('key').!
    kwargs['is_decoy'] = 'isDecoy' if kwargs.include?('is_decoy').!
    if args.all?{ |arg| arg.instance_of?(Pandas.DataFrame) }
      df = Pandas.concat(args)
    else
      df = dataframe(*args, **kwargs)
    end
    filter(df, **kwargs)
  end
  
  fdr = Target_decoy::Make_fdr.call(@is_decoy, nil)
  @_key = lambda { |x| x['SpectrumIdentificationItem'].min{ |sii| sii['mascot:expectation value'] } }
  @qvalues = Target_decoy::Make_qvalues.call(Chain, @is_decoy, nil, @_key)
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
