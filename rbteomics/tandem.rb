# import operator
# from . import xml, auxiliary as aux, _schema_defaults
require 'pandas'
require_relative 'xml'
require_relative '_schema_defaults'
require_relative '../rbteomics/auxiliary/utils'
require_relative '../rbteomics/auxiliary/target_decoy'

module Tandem
  module_function
  class TandemXML
    include Xml::XML
    @@file_format = "TandemXML"
    @@_root_element = "bioml"
    @@_default_schema = Tandem_schema_defaults
    @@_default_iter_path = 'group[@type="model"]'
    @@_structures_to_flatten = Set.new(['domain'])

    def initialize(...)
      __init__(...)
    end
  
    def __init__(*args, **kwargs)
      if kwargs.include?('recursive').!
        super(*args, 'recursive' => true, **kwargs)
      else
        super
      end
    end
  
    # __init__.__doc__ = xml.XML.__init__.__doc__
  
    def _get_info_smart(element, **kw)
      info = _get_info(element, **kw)
      if info['note'].is_a?(Array) && info['note'].size == 1 && info['note'][0].to_set == ['label', 'note'].to_set
        info['note'] = info['note'][0]['note']
      end
      if info.include?('protein') && info.include?('label')
        info.delete('label')
      end
      if info.include?('group')
        info['group'].each do |g|
          label = g.delete('label')
          type_ = g.delete('type')
          info[type_] = {} if info.include?(type_).!
          info[type_][label] = g
        end
        info.delete('group')
      end
      if info.include?('trace')
        info['trace'].each do |t|
          info[t.delete('type')] = t
        end
        info.delete('trace')
      end
      if info['values'].is_a?(Hash)
        info['values'] = info['values']['values']
      end
      if info['attribute'].is_a?(Array)
        info.delete('attribute').each do |a|
          info[a['type']] = a['attribute'].to_f
        end
      end
      if info.include?('support')
        (info['support']['supporting data'] || {}).each do |_, d|
          ['Xdata', 'Ydata'].each do |label|
            d[label]['values'] = d[label]['values'].to_i
            d[label].delete('label')
          end
        end
        if info['support'].include?('fragment ion mass spectrum')
          fims = info['support']['fragment ion mass spectrum']
          fims.merge!(fims.delete('tandem mass spectrum'))
          ['Xdata', 'Ydata'].each do |label|
            info['support']['fragment ion mass spectrum'][label].delete['label']
          end
        end
      end
      if info.include?('charge')
        info['charge'] = info['charge'].to_i
      end
      if info['rt'] == ''
        info['rt'] = nil
      end
      info
    end
  
    def _get_schema_info(read_schema)
      @@_default_schema
    end
  
    def __next__
      n = super
      n.delete('type')
      n
    end
  
    def next
      __next__
    end
  end
  
  def read(source, **kwargs)
    iterative = kwargs['iterative'] || true
  
    TandemXML.new(source, 'read_schema' => false, 'recursive' => true, 'iterative' => iterative)
  end
  
  def iterfind(source, path, **kwargs)
    TandemXML.new(source, **kwargs).iterfing(path, **kwargs)
  end
  
  Chain = ChainBase._make_chain(TandemXML)
  
  @_is_decoy_prefix = lambda do |psm, prefix='DECOY_'|
    psm['protein'].all?{ |prot| prot['label'].start_with?(prefix) }
  end
  
  @_is_decoy_suffix = lambda do |psm, suffix='_DECOY'|
    psm['protein'].all?{ |prot| prot['label'].end_with?(suffix) }
  end
  
  is_decoy = @_is_decoy_prefix
  @qvalues = Target_decoy::Make_qvalues.call(Chain, @_is_decoy_prefix, @_is_decoy_suffix, Itemgetter.new('expect'))
  fdr = Target_decoy::Make_fdr.call(@_is_decoy_prefix, @_is_decoy_suffix)
  filter = lambda do |x = nil|
    if x.nil?
      _make_filter(Chain, @_is_decoy_prefix, @_is_decoy_suffix, Itemgetter.new('expect'), @qvalues)
    elsif x == 'chain' || x == :chain
      y = _make_filter(Chain, @_is_decoy_prefix, @_is_decoy_suffix, Itemgetter.new('expect'), @qvalues)
      _make_chain(t, 'filter', true)
    end
  end
  # filter = _make_filter(chain, _is_decoy_prefix, _is_decoy_suffix, Itemgetter.new('expect'), qvalues)
  # filter.chain = _make_chain(filter, 'filter', true)
  
  def dataframe(*args, **kwargs)
    data = []
    prot_keys = ['id', 'uid', 'label', 'expect']
    pep_keys = ['id', 'pre', 'post', 'start', 'end']
    sep = kwargs.delete('sep') || nil
    pd_kwargs = kwargs.delete('pd_kwargs') || {}
    Chain.new(*args, **kwargs) do |f|
      f.each do |item|
        info = {}
        item.each do |k, v|
          if [String, Integer, Float].include?(v.class)
            info[k] = v
          end
        end
        protein = item['protein'][0]
    
        prot_keys.each do |key|
          vals = item['protein'].map{ |prot| prot[key] }
          if sep.nil?.!
            vals = vals.map{ |val| val.nil?.! ? val.to_s : '' }.join(sep)
          end
          info['protein_' + key] = vals
        end
        pep_keys.each do |key|
          vals = item['protein'].map{ |prot| prot['peptide'][key] }
          if sep.nil?.!
            vals = vals.map{ |val| val.nil?.! ? val.to_s : '' }.join(sep)
          end
          info['peptide_' + key] = vals
        end
        aa = protein['peptide'].delete('aa') || []
        info['modifications'] = aa.map{ |x| "{#{x}[modified]:.3f}@{#{x}[type]}" }.join(',')

        prot_keys.each do |k|
          protein.delete(k)
        end
        pep_keys.each do |k|
          protein['peptide'].delete(k)
        end
        info.merge!(protein['peptide'])
        info['scan'] = item['support']['fragment ion mass spectrum']['note']
        data << info
      end
    end
    Pandas.DataFrame(data, **pd_kwargs)
  end
  
  def filter_df(*args, **kwargs)
    sep = kwargs['sep']
    kwargs['key'] = 'expect' if kwargs.include?('key').!
    if args.all?{ |arg| arg.instance_of?(Pandas.DataFrame) }
      if args.size > 1
        df = Pandas.concat(args)
      else
        df = args[0]
      end
    else
      read_kw = kwargs.select{ |k, _| [0, '', nil, false, [], {}].include?(k).! }.select{ |k| ['iterative', 'read_schema', 'sep', 'pd_kwargs'].include?(k) }.map{ |k| [k, kwargs.delete(k)] }.to_h
      df = dataframe(*args, **read_kw)
    end

    if kwargs.include?('is_decoy').!
      if sep.nil?.!
        if kwargs.include?('decoy_suffix')
          kwargs['is_decoy'] = df['protein_label'].to_s.split(sep).apply(
            lambda { |s| s.all{ |x| x.end_with?(kwargs['decoy_suffix']) } })
        else
          kwargs['is_decoy'] = df['protein_label'].to_s.split(sep).apply(
            lambda { |s| s.all?{ |x| x.start_with(kwargs['decoy_prefix'] || 'DECOY_') } })
        end
      else
        if kwargs.include?('decoy_suffix')
          kwargs['is_decoy'] = df['protein_label'].apply(
            lambda { |s| s.all?{ |x| x.endswith(kwargs['decoy_suffix']) } })
        else
          kwargs['is_decoy'] = df['protein_label'].apply(
            lambda { |s| s.all?{ |x| x.start_with(kwargs['decoy_prefix'] || 'DECOY_') } })
        end
      end
    end

    Target_decoy::Filter(df, **kwargs)
  end  
end
