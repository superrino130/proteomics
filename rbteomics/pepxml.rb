# from lxml import etree
require 'rexml/document'
# from . import xml, auxiliary as aux, _schema_defaults
require_relative 'xml'
require_relative 'auxiliary/file_helpers'
require_relative 'auxiliary/target_decoy'
require_relative '_schema_defaults'
require 'set'
class PepXML
# class PepXML < IndexedXML
  # prepend IndexSavingMixin
  # include TaskMappingMixin
  #class PepXML(xml.MultiProcessingXML, xml.IndexSavingXML):
  #class MultiProcessingXML(IndexedXML, TaskMappingMixin):
  #class IndexSavingXML(IndexSavingMixin, IndexedXML):

  def initialize
    @_index_class = HierarchicalOffsetIndex.new
    @file_format = 'pepXML'
    @_root_element = 'msms_pipeline_analysis'
    @_default_schema = _schema_defaults._pepxml_schema_defaults
    @_default_version = '1.15'
    @_default_iter_tag = 'spectrum_query'
    @_indexed_tags = Set.new(['spectrum_query'])
    @_indexed_tag_keys = {'spectrum_query' => 'spectrum'}
    @_default_id_attr = 'spectrum'
    @_structures_to_flatten = Set.new(['search_score_summary', 'modification_info'])
    # attributes which contain unconverted values
    @_convert_items = {
      'float' => Set.new(['calc_neutral_pep_mass', 'massdiff', 'probability', 'variable', 'static']),
      'int' => Set.new(['start_scan', 'end_scan', 'index', 'num_matched_peptides']),
      'bool' => Set.new(['is_rejected']),
      'floatarray' => Set.new(['all_ntt_prob'])
    }

    @mpx = MultiProcessingXML.new
    @isx = IndexSavingXML.new
  end

  # # from class IndexSavingXML
  # def _read_byte_offsets
  #   File.open(@_byte_offset_filename, 'r') { |f|
  #     index = @_index_class.load(f)
  #     if index.schema_version.nil?
  #       raise TypeError("Legacy Offset Index!")
  #     end
  #     @_offset_index = index
  #   }
  # end
  
  # # from class MultiProcessingXML
  # def _task_map_iterator
  #   @_offset_index[self._default_iter_tag].to_enum
  # end

  def _get_info_smart(element, **kwargs)
    if kwargs.include?('ename')
      name = kwargs.delete('ename')
    else
      name = _local_name(element)
    end
    rec = kwargs.delete('recursive') || nil
    if name == @_root_element
      info = _get_info(element, 'ename' => name, 'recursive' => rec.nil?.! ? rec : false, **kwargs)
    else
      info = _get_info(element, 'ename' => name, 'recursive' => rec.nil?.! ? rec : true, **kwargs)
    end

    def safe_float(s)
      if s.upcase == 'NAN'
        Float::NAN
      elsif s.upcase == '-NAN'
        -Float::NAN
      elsif s.upcase == 'INF'
        Float::INFINITY
      elsif s.upcase == '-INF'
        -Float::INFINITY
      else
        s.to_f
      end
    end

    converters = {
      'float' => lambda { |x| safe_float(x) },
      'int' => lambda { |x| x.to_i },
      'bool' => lambda { |x| ['1', 'true'].include?(x.downcase) },
      'floatarray' => lambda { |x| x[1...-1].split(',').map{ _1.to_f } }
    }
    info.to_h.each do |k, v|
      @_convert_items.each do |t, s|
        if s.include?(k)
          info.delete(k)
          inf[k] = converters[t].call(v)
        end
      end
    end
    ['search_score', 'parameter'].each do |k|
      if info.include?(k) && info[k].instance_of?(Array) && 
        info[k].map{ |x| x.instance_of?(Hash) && x.size == 1 }.all?
        scores = {}
        info[k].each do |score|
          name, value = score.to_a
          scores[name] = value.to_f || value
        end
        info[k] = scores
      end
    end
    if info.include?('search_result') && info['search_result'].size == 1
      info.merge(info['search_result'][0])
      info.delete('search_result')
    end
    if info.include?('protein') && info.include?('peptide')
      info['proteins'] = [{'protein' => info.delete('protein'),
        'protein_descr' => info.delete('protein_descr')}]
      ['peptide_prev_aa', 'peptide_next_aa', 'protein_mw'].each do |add_key|
        if info.include?(add_key)
          info['proteins'][0][add_key] = info.delete(add_key)
        end
      end
      info['proteins'][0]['num_tol_term'] = info.delete('num_tol_term') || 0
      if info.include?('alternative_protein')
        info['proteins'].concat(info['alternative_protein'])
        info.delete('alternative_protein')
      end
    end
    if info.include?('peptide') && info.include?('modified_peptide').!
      info['modified_peptide'] = info['peptide'].dup
    end
    if info.include?('peptide')
      info['modifications'] = info.delete('mod_aminoacid_mass') || []
      if info.include?('mod_nterm_mass')
        info['modifications'].insert(0, {'position' => 0,
          'mass' => info.delete('mod_nterm_mass').to_f})
      end
      if info.include?('mod_cterm_mass')
        info['modifications'] << ({'position' => 1 + info['peptide'].size,
          'mass' => info.delete('mod_cterm_mass').to_f})
      end
    end
    if info.include?('modified_peptide') && info['modified_peptide'] == info['peptide']
      if ['', 0, nil, false, [], {}].include?(info['modifications'])
        info['modifications'] = []
      else
        mp = info['modified_peptide']
        info['modifications'].sort_by{ |m| -m['position'] }.each do |mod|
          if [0, info['peptide'].size + 1].include?(mod['position']).!
            ps = mod['position']
            mp = mp[0...ps] + "[#{mod['mass'].to_i}]" + mp[ps..]
          end
        end
        info['modified_peptide'] = mp
      end
    end
    if info.include?('search_hit')
      info['search_hit'] = info['search_hit'].sort_by{ |x| x['hit_rank'] }
    end
    info
  end

  def method_missing(name, *args)
    [@mpx, @isx].each do |x|
      if x.respond_to?(name)
        raise NoMethodError.new(name)
      end
    end  
  end
end

def read(source, **kwargs)
  read_schema = kwargs['read_schema'] || false
  iterative = kwargs['iterative'] ||  true
  PepXML.new(source, read_schema: read_schema, iterative: iterative)
end

def iterfind(source, path, **kwargs)
  PepXML.new(source, **kwargs).iterfind(path, **kwargs)
end

# Version_info = _make_version_info(PepXML)

def roc_curve(source)
  # parser = etree.XMLParser(remove_comments: true, ns_clean: true)
  tree = REXML::Document.new(source)

  roc_curve = []
  tree.xpath(
    "/*[local-name()='msms_pipeline_analysis'] \
    //*[local-name()='analysis_summary' and @analysis='peptideprophet'] \
    //*[local-name()='peptideprophet_summary'] \
    //*[local-name()='roc_error_data']"
  ).each do |roc_error_data|
    roc_error_data.xpath("*[local-name()='roc_data_point' or local-name()='error_point']").each do |element|
      data_point = element.attrib.to_a
      data_point.each do |key|
        data_point[key] = data_point[key].to_f
      end
      data_point["charge"] = roc_error_data.attrib["charge"]
      data_point["tag"] = etree.QName(element).localname
      roc_curve << data_point
    end
  end
  roc_curve
end

Chain = ChainBase._make_chain('read')

_is_decoy_prefix = lambda do |psm, prefix: 'DECOY_'|
  psm['search_hit'][0]['proteins'].map{ |protein| protein['protein'].start_with(prefix) }.all?
end

_is_decoy_suffix = lambda do |psm, suffix: '_DECOY'|
  psm['search_hit'][0]['proteins'].map{ |protein| protein['protein'].end_with(suffix) }.all?
end

IS_decoy = _is_decoy_prefix
FDR = _make_fdr(_is_decoy_prefix, _is_decoy_suffix)
_key = lambda { |x| x['search_hit'].min{ |sh| sh['search_score']['expect'] } }
Qvalues = _make_qvalues(Chain, _is_decoy_prefix, _is_decoy_suffix, _key)
Filter = _make_filter(Chain, _is_decoy_prefix, _is_decoy_suffix, _key, Qvalues)
Filter.chain = _make_chain(Filter, 'filter', true)

def DataFrame(*args, **kwargs)
  import Pandas
  kwargs = kwargs.dup
  sep = kwargs.pop('sep', None)
  pd_kwargs = kwargs.pop('pd_kwargs', {})

end