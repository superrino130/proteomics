# import operator
# from . import xml, auxiliary as aux, _schema_defaults
require 'pandas'
require_relative 'xml'
require_relative '_schema_defaults'
require_relative '../rbteomics/auxiliary/utils'
require_relative '../rbteomics/auxiliary/target_decoy'

class TandemXML < XML
  def initialize(...)
    @file_format = "TandemXML"
    @_root_element = "bioml"
    @_default_schema = Tandem_schema_defaults
    @_default_iter_path = 'group[@type="model"]'
    @_structures_to_flatten = Set.new(['domain'])
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
          # d[label]['values'] = d[label]['values'].astype(int)
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
    @_default_schema
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

def tandem_read(source, **kwargs)
  iterative = kwargs['iterative'] || true

  TandemXML.new(source, 'read_schema' => false, 'recursive' => true, 'iterative' => iterative)
end

def iterfind(source, path, **kwargs)
  # Not started
end

TD_chain = ChainBase._make_chain(TandemXML)

@_is_decoy_prefix = lambda do |psm, prefix='DECOY_'|
  # Not started
end

@_is_decoy_suffix = lambda do |psm, suffix='_DECOY'|
  # Not started
end

is_decoy = @_is_decoy_prefix
@qvalues = _make_qvalues(TD_chain, @_is_decoy_prefix, @_is_decoy_suffix, Itemgetter.new('expect'))
fdr = _make_fdr(@_is_decoy_prefix, @_is_decoy_suffix)
filter = lambda do |x = nil|
  if x.nil?
    _make_filter(TD_chain, @_is_decoy_prefix, @_is_decoy_suffix, Itemgetter.new('expect'), @qvalues)
  elsif x == 'chain' || x == :chain
    y = _make_filter(TD_chain, @_is_decoy_prefix, @_is_decoy_suffix, Itemgetter.new('expect'), @qvalues)
    _make_chain(t, 'filter', true)
  end
end
# filter = _make_filter(TD_chain, _is_decoy_prefix, _is_decoy_suffix, Itemgetter.new('expect'), qvalues)
# filter.chain = _make_chain(filter, 'filter', true)

def dataframe(*args, **kwargs)
  data = []
  prot_keys = ['id', 'uid', 'label', 'expect']
  pep_keys = ['id', 'pre', 'post', 'start', 'end']
  sep = kwargs.delete('sep') || nil
  pd_kwargs = kwargs.delete('pd_kwargs') || {}
  f = TD_chain(*args, **kwargs)
  f.each do |item|
    info = {}
    item.each do |k, v|
      if [String, Integer, Float].include?(v.class)
        info[k] = v
      end
    end
    protein = info['protein'][0]

    prot_keys.each do |key|
      vals = item['protein'].map{ |prot| prot[key] }
      if sep.nil?.!
        vals = vals.map{ |val| val.nil?.! ? val.to_s : '' }.join(sep)
      end
      info['peptide_' + key] = vals
    end
    aa = protein['peptide'].delete('aa') || []
    info.merge!(protein['peptide'])
    data << info
  end
  Pandas.DataFrame(data, **pd_kwargs)
end

def filter_df(*args, **kwargs)
  # Not started
end
