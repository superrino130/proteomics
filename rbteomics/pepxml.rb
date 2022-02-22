# from lxml import etree
require 'rexml/document'
# from . import xml, auxiliary as aux, _schema_defaults
require_relative 'xml'
require_relative 'auxiliary/file_helpers'
require_relative 'auxiliary/target_decoy'
require_relative '_schema_defaults'
require 'set'

module Pepxml
  module_function
  class PepXML
    extend File_helpers
    extend File_helpers::IndexSavingMixin
    include Xml::IndexSavingXML
    include Xml::MultiProcessingXML
  
    @@file_format = 'pepXML'
    @@_root_element = 'msms_pipeline_analysis'
    @@_default_schema = Pepxml_schema_defaults
    @@_default_version = '1.15'
    @@_default_iter_tag = 'spectrum_query'
    @@_indexed_tags = Set.new(['spectrum_query'])
    @@_indexed_tag_keys = {'spectrum_query' => 'spectrum'}
    @@_default_id_attr = 'spectrum'
    @@_structures_to_flatten = Set.new(['search_score_summary', 'modification_info'])
    @@_convert_items = {
      'float' => Set.new(['calc_neutral_pep_mass', 'massdiff', 'probability', 'variable', 'static']),
      'int' => Set.new(['start_scan', 'end_scan', 'index', 'num_matched_peptides']),
      'bool' => Set.new(['is_rejected']),
      'floatarray' => Set.new(['all_ntt_prob'])
    }

    def initialize(...)
      __init__(...)
    end
        
    def _get_info_smart(element, **kwargs)
      if kwargs.include?('ename')
        name = kwargs.delete('ename')
      else
        name = Xml._local_name(element)
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
        elsif s.start_with('+-0')
          0
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
          info[k].all?{ |x| x.instance_of?(Hash) && x.size == 1 }
          scores = {}
          info[k].each do |score|
            name, value = score.to_a
            scores[name] = value.to_f || value
          end
          info[k] = scores
        end
      end
      if info.include?('search_result') && info['search_result'].size == 1
        info.merge!(info['search_result'][0])
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
          info['modifications'] << {'position' => 1 + info['peptide'].size,
            'mass' => info.delete('mod_cterm_mass').to_f}
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
  end
    
  Read = lambda do |source, **kwargs|
    read_schema = kwargs['read_schema'] || false
    iterative = kwargs['iterative'] || true
    PepXML.new(source, 'read_schema' => read_schema, 'iterative' => iterative)
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
        data_point = element.attrib.tally
        data_point.each do |key, _|
          data_point[key] = data_point[key].to_f
        end
        data_point["charge"] = roc_error_data.attrib["charge"]
        data_point["tag"] = etree.QName(element).localname
        roc_curve << data_point
      end
    end
    roc_curve
  end
  
  Chain = File_helpers::ChainBase._make_chain('Read', Read)
  
  @_is_decoy_prefix = lambda do |psm, prefix: 'DECOY_'|
    psm['search_hit'][0]['proteins'].all?{ |protein| protein['protein'].start_with(prefix) }
  end
  
  @_is_decoy_suffix = lambda do |psm, suffix: '_DECOY'|
    psm['search_hit'][0]['proteins'].all?{ |protein| protein['protein'].end_with(suffix) }
  end
  
  is_decoy = @_is_decoy_prefix
  fdr = Target_decoy::Make_fdr.call(@_is_decoy_prefix, @_is_decoy_suffix)
  _key = lambda { |x| x['search_hit'].min{ |sh| sh['search_score']['expect'] } }
  @qvalues = Target_decoy::Make_qvalues.call(Chain, @_is_decoy_prefix, @_is_decoy_suffix, _key)
  filter = lambda do |x = nil|
    if x.nil?
      Target_decoy::Make_filter(Chain, @_is_decoy_prefix, @_is_decoy_suffix, _key, @qvalues)
    elsif x == 'chain' || x == :chain
      y = Target_decoy::Make_filter(Chain, @_is_decoy_prefix, @_is_decoy_suffix, _key, @qvalues)
      Target_decoy::Make_chain(y, 'filter', true)
    end
  end
  # Filter = _make_filter(Chain, _is_decoy_prefix, _is_decoy_suffix, _key, Qvalues)
  # Filter.chain = _make_chain(Filter, 'filter', true)
  
  def dataframe(*args, **kwargs)
    require 'pandas'
    kwargs = kwargs.dup
    sep = kwargs.delete('sep') || nil
    pd_kwargs = kwargs.delete('pd_kwargs') || {}
    def gen_item
      f = Chain.new(*args, **kwargs)
      f.each do |item|
        info = {}
        item.each do |k, v|
          if [String, Integer, Float].include?(v.class)
            info[k] = v
          end
        end
        if item.include?('search_hit')
          sh = item['search_hit'][0]
          proteins = sh.delete('proteins')
          prot_dict = {}
          proteins.each do |p|
            p.each do |k, _|
              prot_dict[k] = []
            end
          end
          proteins.each do |p|
            prot_dict.each do |k, v|
              v << p[k]
            end
          end
          if seq.nil?
            info.merge!(prot_dict)
          else
            prot_dict.each do |k, v|
              info[k] = v.map{ |val| val.nil?.! ? val.to_s : '' }.join(sep)
            end
          end
          info.merge!(sh.delete('search_score'))
          mods = sh.delete('modifications') || []
          formatted_mods = mods.map{ |x| "#{sprintf("%.03f", x[mass])}@#{x[position]}" }
          if sep.nil?.!
            info['modifications'] = formatted_mods.join(sep)
          else
            info['modifications'] = formatted_mods
          end
          sh.each do |k, v|
            if [String, Integer, Float].include?(v.class)
              info[k] = v
            end
          end
          if sh.include?('analysis_result')
            sh['analysis_result'].each do |ar|
              if ar['analysis'] == 'peptideprophet'
                if ar['peptideprophet_result'].include?('parameter')
                  info.merge!(ar['peptideprophet_result']['parameter'])
                end
                info['peptideprophet_probability'] = ar['peptideprophet_result']['probability']
                info['peptideprophet_ntt_prob'] = ar['peptideprophet_result']['all_ntt_prob']
              elsif ar['analysis'] == 'interprophet'
                info.merge!(ar['interprophet_result']['parameter'])
                info['interprophet_probability'] = ar['interprophet_result']['probability']
                info['interprophet_ntt_prob'] = ar['interprophet_result']['all_ntt_prob']
              end
            end
          end
        end
        yield info
      end
    end
    Pandas.DataFrame(gen_items, **pd_kwargs)
  end

  def filter_df(*args, **kwargs)
    require 'pandas'
    sep = kwargs['sep']
    kwargs['key'] = 'expect' if kwargs.include?('key').!
    if args.all?{ |arg| arg.instance_of?(Pandas.DataFrame) }
      if args.size > 1
        df = Pandas.concat(args)
      else
        df = args[0]
      end
    else
      read_kw = {}
      kwargs.each do |k, _|
        if ['iterative', 'read_schema', 'sep', 'pd_kwargs'].include(k)
          read_kw[k] = kwargs.delete(k)
        end
      end
      df = dataframe(*args, **read_kw)
    end
    if kwargs.include?('is_decoy').!
      if sep.nil?.!
        if kwargs.include?('decoy_suffix')
          kwargs['is_decoy'] = df['protein'].to_s.split(';').apply(
            lambda { |s| s.all?{ |x| x.end_with(kwargs['decoy_suffix']) } })
        else
          kwargs['is_decoy'] = df['protein'].apply(
            lambda { |s| s.all?{ |x| x.start_with(kwargs['decoy_preffix'] || 'DECOY_') } })
        end
      else
      end
    end
    Target_decoy::filter.call(df, **kwargs)
  end
end
