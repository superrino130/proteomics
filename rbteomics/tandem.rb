# import operator
# from . import xml, auxiliary as aux, _schema_defaults
require 'pandas'
require_relative 'xml'
require_relative '_schema_defaults'
require_relative '../rbteomics/auxiliary/utils'
require_relative '../rbteomics/auxiliary/target_decoy'

class TandemXML < XML
  # Not started
end

def read(source, **kwargs)
  iterative = kwargs['iterative'] || true

  # Not started
end

def iterfind(source, path, **kwargs)
  # Not started
end

TD_chain = ChainBase._make_chain(TandemXML)

_is_decoy_prefix = lambda do |psm, prefix='DECOY_'|
  # Not started
end

_is_decoy_suffix = lambda do |psm, suffix='_DECOY'|
  # Not started
end

is_decoy = _is_decoy_prefix
qvalues = _make_qvalues(TD_chain, _is_decoy_prefix, _is_decoy_suffix, Itemgetter.new('expect'))
filter = _make_filter(TD_chain, _is_decoy_prefix, _is_decoy_suffix, Itemgetter.new('expect'), qvalues)
fdr = _make_fdr(_is_decoy_prefix, _is_decoy_suffix)
filter.chain = _make_chain(filter, 'filter', true)

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

