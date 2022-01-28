# from __future__ import absolute_import
# import re
# import operator as op
# import math

# try:
#     basestring
# except NameError:
#     basestring = (str, bytes)
BaseString ||= String

# try:
#     from collections.abc import Container, Sized
# except ImportError:
#     from collections import Container, Sized
# from bisect import bisect_right
# from contextlib import contextmanager


# from .structures import PyteomicsError
# from .file_helpers import _keepstate, IteratorContextManager, _make_chain, ChainBase, TableJoiner
# from .patch import pd
require_relative 'structures'
require_relative 'file_helpers'
# require_relative 'patch'

require 'numpy'
require 'pandas'
require 'set'

def _fix_docstring(f, **defaults)
  defaults.each do |argname, v|
    if !!v
      # f.__doc__ = re.sub('{} : .*'.format(argname),
      #                          lambda m: m.group() + ', optional', f.__doc__)
      f.__doc__ = ""
    end
  end
end

def _calculate_qvalues(scores, isdecoy, **kwargs)
  peps = kwargs['peps'] || false
  correction = kwargs.delete('correction') || 0
  ratio = kwargs.delete('ratio') || 1
  remove_decoy = kwargs['remove_decoy'] || false
  formula = kwargs.delete('formula') || [2, 1][remove_decoy]
  if [1, 2].include?(formula).!
    raise PyteomicsError.new("'formula' must be either 1 or 2")
  end

  cumsum = isdecoy.cumsum(dtype: Numpy.float64)
  tfalse = cumsum.copy()
  ind = Numpy.arange(1.0, scores.shape[0] + 1.0, dtype: Numpy.float64)

  if peps
    q = cumsum / ind
  else
    if correction.instance_of?(Integer)
      if correction == 1
        tfalse += 1
      elsif
        p = 1.0 / (1.0 + ratio)
        targ = ind - cumsum
        tfalse.size.times do |i|
          tfalse[i] = _expectation(cumsum[i], targ[i], p)
        end
      end
    elsif 0 < correction && correction < 1
      p = 1.0 / (1.0 + ratio)
      targ = ind - cumsum
      tfalse.size.times do |i|
        tfalse[i] = _confidence_value(correction, cumsum[i], targ[i], p)
      end
    elsif correction
      raise PyteomicsError.new("Invalid value for 'correction'.")
    end
    if formula == 1
      q = tfalse / (ind - cumsum) / ratio
    else
      q = (cumsum + tfalse / ratio) / ind
    end
  end

  (scores.size - 1...0).each do |i|
    if scores[i] == scores[i - 1] || q[i - 1] > q[i]
      q[i - 1] = q[i]
    end
  end
  q
end

def _qvalues_df(psms, keyf, isdecoy, **kwargs)
  full = kwargs['full_output'] || false
  remove_decoy = kwargs['remove_decoy'] || false
  peps = kwargs['pep']
  decoy_or_pep_label = _decoy_or_pep_label(**kwargs)
  q_label = kwargs.key?('q_label') ? kwargs['q_label'] : kwargs['q_label'] = 'q'
  score_label = kwargs.key?('score_label') ? kwargs['score_label'] : kwargs['score_label'] = 'score'
  keyf = psms.apply(keyf, axis=1) if defined?(keyf)
  isdecoy = psms.apply(isdecoy, axis=1) if defined?(isdecoy)
  if keyf.instance_of?(String)
    if psms.shape[0]
      psms[score_label] = keyf
    else
      psms[score_label] = []
    end
    keyf = kwargs['score_label']
  end
  if isdecoy.instance_of?(String)
    if psms.shape[0]
      psms[decoy_or_pep_label] = isdecoy
    else
      psms[decoy_or_pep_label] = []
    end
    isdecoy = decoy_or_pep_label
  end
  reverse = kwargs['reverse'] || false

  if full.!
    if peps.nil?
      fields = [[keyf, Numpy.float64], [isdecoy, Numpy.bool_], [q_label, Numpy.float64]]
    else
      fields = [[isdecoy, Numpy.float64], [q_label, Numpy.float64]]
    end
    dtype = Numpy.dtype(fields)
  end
  psms.sort_values([keyf, isdecoy], ascending: [reverse.!, true], inplace: true)

  if psms.shape[0].!
    if full
      psms[q_label] = []
      return psms
    else
      return Numpy.array([], dtype=dtype)
    end
  end

  q = _calculate_qvalues(psms[keyf].values, psms[isdecoy].values, peps.nil?.!, **kwargs)
  if remove_decoy
    q = q[~psms[isdecoy].values]
    psms = psms[~psms[isdecoy]].copy()
  end
  if full.!
    psms_ = Numpy.empty_like(q, dtype=dtype)
    psms_[keyf] = psms[keyf] if peps.nil?
    psms_[isdecoy] = psms[isdecoy]
    psms_[q_label] = q
    psms = psms_
  else
    q_label = kwargs['q_label']
    psms[q_label] = q
  end
  psms
end

def _decoy_or_pep_label(**kwargs)
  peps = kwargs['pep']
  return peps.nil? ? (kwargs['decoy_label'] || 'is decoy') : (kwargs['pep_label'] || peps.instance_of?(BaseString) ? peps : 'PEP')
end

def _construct_dtype(*args, **kwargs)
  full = kwargs.delete('full_output') || false
  peps = kwargs['pep']
  q_label = kwargs.include?('q_label') ? kwargs['q_label'] : kwargs['q_label'] = 'q'
  score_label = kwargs.include?('score_label') ? kwargs['score_label'] : kwargs['score_label'] = 'score'

  fields = [
    [score_label, Numpy.float64],
    [_decoy_or_pep_label(**kwargs), peps.nil? ? Numpy.bool_ : Numpy.float64],
    [q_label, Numpy.float64]
  ]
  
  if full
    dtypes = Set.new(args.map{ getattr(_1, 'dtype', nil) })
    if dtypes.size == 1 || dtypes.all?{ _1 }
      psm_dtype = dtype.to_a[-1]
      dtype.delete(dtype.to_a[-1])
    else
      psm_dtype = Numpy.object_
    end
    dtype = Numpy.dtype(fields + [['psm', psm_dtype]])
  else
    dtype = Numpy.dtype(fields)
  end
  dtype
end

def _make_qvalues(read, is_decoy_prefix, is_decoy_suffix, key)
  @mq_is_decoy_prefix = is_decoy_prefix
  @mq_key = key
  def qvalues(*args, **kwargs)
    #@_keepstate
    def get_scores(*args, **kwargs)
      scores = []
      f = read(*args, **kwargs)
      f.each_with_index do |psm, i|
        row = []
        [keyf, isdecoy].each do |func|
          if defined?(func)
            row << func(psm)
          elsif func.instance_of?(String)
            row << psm[func]
          else
            row << func[i]
          end
        end
        row << nil
        row << psm if full
        scores << [row]
      end
      scores
    end

    peps = kwargs['pep'] || nil
    if peps.nil?.!
      x = Set.new(['is_decoy', 'remove_decoy', 'formula',
        'ratio', 'correction']) & kwargs.keys
      if x.empty?.!
        raise PyteomicsError.new("Can't use these parameters with 'pep': " + x.join(', '))
      end
    end
    keyf = kwargs.delete('key') || @mq_key
    reverse = kwargs['reverse'] || false
    if keyf.nil?
      keyf = peps
      if reverse
        raise PyteomicsError.new('reverse = True when using PEPs for sorting')
      end
    end

    if defined?(keyf).! && [Sized, Container].include?(keyf.class).!
      keyf = Numpy.array(keyf.to_a)
    end

    if peps.nil?
      if kwargs.include?('is_decoy').!
        if kwargs.include?('decoy_suffix')
          isdecoy = lambda { |x| is_decoy_suffix(x, kwargs['decoy_suffix']) }
        elsif kwargs.include?('decoy_prefix')
          isdecoy = lambda { |x| is_decoy_prefix(x, kwargs['decoy_prefix']) }
        else
          isdecoy = @mq_is_decoy_prefix
        end
      else
        isdecoy = kwargs['is_decoy']
      end
    else
      isdecoy = peps
    end
    
    if defined?(isdecoy).! and [Sized, Container].include?(isdecoy.class).!
      isdecoy = Numpy.array(isdecoy.to_a)
    end
    remove_decoy = kwargs['remove_decoy'] || false
    decoy_or_pep_label = _decoy_or_pep_label(**kwargs)
    score_label = kwargs.key?('score_label') ? kwargs['score_label'] : kwargs['score_label'] = 'score'
    q_label = kwargs.key?('q_label') ? kwargs['q_label'] : kwargs['q_label'] = 'q'
    dtype = _construct_dtype(*args, **kwargs)
    full = kwargs['full_output'] || false
    arr_flag = false
    psms = nil

    if args.map{ _1.instance_of?(Pandas.DataFrame) }.all?{ _1 }
      psms = Pandas.concat(args)
      return _qvalues_df(psms, keyf, isdecoy, **kwargs)
    end

    if args.map{ _1.instance_of?(Numpy.ndarray)}.all{ _1 }.!
      keyf = Itemgetter.new(keyf) if keyf.instance_of?(String)
      isdecoy = Itemgetter.new(isdecoy) if isdecoy.instance_of?(String)
      peps = Itemgetter.new(peps) if peps.instance_of?(String)
    end

    if defined?(keyf) || defined?(isdecoy)
      kwargs.delete('full_output')
      scores = Numpy.array(get_scores(*args, **kwargs), dtype=dtype)
    else
      if args.map{ _1.instance_of?(Numpy.ndarray) }.all?{ _ 1 }
        psms = Numpy.concatenate(args)
      end
      if keyf.instance_of?(String).!
        keyf = Numpy.array(keyf)
        arr_flag = true
      end
      if isdecoy.instance_of?(String).!
        isdecoy = Numpy.array(isdecoy)
        arr_flag = true
      end

      if arr_flag
        scores = Numpy.empty(hasattr(keyf, 'size') ? keyf.size : isdecoy.size, dtype=dtype)
        [keyf, isdecoy].zip([score_label, decoy_or_pep_label]).each do |func, label|
          if func.instance_of?(String).!
            scores[label] = func
          else
            scores[label] = psms[func]
          end
        end
      else
        scores = Numpy.empty(psms.shape[0], dtype=dtype)
        scores[score_label] = psms[keyf]
        scores[decoy_or_pep_label] = psms[isdecoy]
      end
    end

    if scores.size != 0
      if full && !!psms
        return psms
      end
      return scores
    end

    if reverse.!
      keys = scores[decoy_or_pep_label], scores[score_label]
    else
      keys = scores[decoy_or_pep_label], -scores[score_label]
    end
    lexsort = Numpy.lexsort(keys)
    scores = scores[lexsort]
    psms = psms[lexsort] if !!psms

    scores[q_label] = _calculate_qvalues(scores[score_label], scores[decoy_or_pep_label], !!peps, **kwargs)
    if remove_decoy
      psms = psms[~scores[decoy_or_pep_label]] if !!psms
      scores = scores[~scores[decoy_or_pep_label]]
    end

    if full && psms.nil?.!
      if psms.instance_of?(Numpy.ndarray)
        # fields = sorted(psms.dtype.fields, key=lambda x: psms.dtype.fields[x][1])
        fields = psms.dtype.fields.sort_by{ |x| psms.dtype.fields[x][1] }
        extra = []
        [keyf, isdecoy].zip(['score', decoy_or_pep_label]).each do |func, label|
          if (func.instance_of?(BaseString) || psms.dtype.fields.include?(label)).!
            extra << label
          elsif psms.dtype.fields.include?(label)
            psms[label] = scores[label]
          end
        end
        newdt = fields.map{ [_1, psms.dtype.fields[_1][0]] } + 
          extra.map{ [_1, Numpy.float64] } + [q_label, Numpy.float64]
        psms_ = psms
        psms = Numpy.empty_like(psms_, dtype: newdt)
        fields.each do |f|
          psms[f] = psms_[f]
        end
        extra.each do |f|
          psms[f] = scores[f]
        end
      else
        [keyf, isdecoy].zip(['score', decoy_or_pep_label]).each do |func, label|
          psms[label] = scores[label] if label.instance_of?(String).!
        end
      end
      psms[q_label] = scores[q_label]
      return psms
    end
    scores
  end

  _fix_docstring(qvalues, 'is_decoy' => is_decoy_prefix, 'key' => key)
  if read == _iter
    qvalues.__doc__ = qvalues.__doc__.replace(
      "positional args : file or str
      Files to read PSMs from. All positional arguments are treated as
      files.",
      "positional args : iterables
      Iterables to read PSMs from. All positional arguments are chained."
      ).replace(
      "\n            .. warning::
      The default function may not work
      with your files, because format flavours are diverse.
      decoy_prefix : str, optional
      If the default `is_decoy` function works for you, this parameter specifies which
      protein name prefix to use to detect decoy matches. If you provide your own
      `is_decoy`, or if you specify `decoy_suffix`, this parameter has no effect.
      Default is `DECOY_`.
      decoy_suffix : str, optional
      If the default `is_decoy` function works for you, this parameter specifies which
      protein name suffix to use to detect decoy matches. If you provide your own
      `is_decoy`, this parameter has no effect. Mutually exclusive with `decoy_prefix`.\n",
      "")
  end
  qvalues
end

def _make_filter(read, is_decoy_prefix, is_decoy_suffix, key, qvalues)
  def filter(*args, **kwargs)
    begin
      fdr = kwargs.delete('fdr')
    rescue => exception
      raise PyteomicsError.new('Keyword argument required: fdr')
    end

    args = args.map{ [Container, Sized].instance_of?(_1.class).! ? arg.to_a : _1 }
    peps = kwargs['pep']
    if peps.nil?
      remove_decoy = kwargs.delete('remove_decoy') || true
      scores = qvalues(*args, remove_decoy=remove_decoy, **kwargs)
    else
      scores = qvalues(*args, **kwargs)
    end
    keyf = kwargs.delete('key') || key
    keyf = peps if keyf.nil?
    reverse = kwargs.delete('reverse') || false
    better = [op.lt, op.gt][reverse == 0 ? 0 : 1]
    if kwargs.include?('is_decoy').!
      if kwargs.include?('decoy_suffix')
        isdecoy = lambda { |x| is_decoy_suffix(x, kwargs['decoy_suffix']) }
      elsif kwargs.include?('decoy_prefix')
        isdecoy = lambda { |x| is_decoy_prefix(x, kwargs['decoy_prefix']) }
      else
        isdecoy = is_decoy_prefix
      end
    else
      isdecoy = kwargs['is_decoy']
    end
    kwargs.delete('formula')
    decoy_or_pep_label = _decoy_or_pep_label(**kwargs)
    score_label = kwargs.include?('score_label') ? kwargs['score_label'] : kwargs['score_label'] = 'score'
    q_label = kwargs['q_label'] || 'q'

    begin
      i = scores[q_label].searchsorted(fdr, side: 'right')
      i = i[0] if i.instance_of?(Sized)
    rescue => exception
      i = scores['q'].bsearch_index{ _1 > fdr } || scores['q'].size
    end
    if kwargs.delete('full_output')
      if scores.instance_of?(Pandas.DataFrame)
        return scores.iloc[0..i]
      elsif defined?(keyf) || defined?(isdecoy)
        return scores['psm'][0...i]
      else
        return scores[0...i]
      end
    elsif scores.size != 0
      return []
    end
    if peps.nil?
      label = score_label
    else
      label = decoy_or_pep_label
    end
    cutoff = i < scores.size ? scores[label][i] : scores[label][-1] + [1, -1][reverse == 0 ? 0 : 1]
    def out()
      f = read(*args, **kwargs)
      f.zip(scores).each do |p1, s|
        if peps.nil?.! || remove_decoy.! || s[decoy_or_pep_label].!
          if better(s[label], cutoff)
            p1
            yield
          end
        end
      end
    end
    return out()
  end

  def _filter(*args, **kwargs)
    if kwargs.key?('full_output')
      if [0, '', nil, []].include?(kwargs['full_output']).!
        kwargs.delete('full_output')
        return filter(*args, full_output=True, **kwargs)
      end
    end
    return IteratorContextManager.new(*args, parser_func=filter, **kwargs)
  end

  _fix_docstring(_filter, is_decoy: is_decoy_prefix, key: key)
  if read == _iter
    _filter.__doc__ = _filter.__doc__.replace(
      "positional args : file or str
      Files to read PSMs from. All positional arguments are treated as
      files.",
      "positional args : iterables
      Iterables to read PSMs from. All positional arguments are chained."
      ).replace(
      "\n            .. warning::
      The default function may not work
      with your files, because format flavours are diverse.
      decoy_prefix : str, optional
      If the default `is_decoy` function works for you, this parameter specifies which
      protein name prefix to use to detect decoy matches. If you provide your own
      `is_decoy`, or if you specify `decoy_suffix`, this parameter has no effect.
      Default is `DECOY_`.
      decoy_suffix : str, optional
      If the default `is_decoy` function works for you, this parameter specifies which
      protein name suffix to use to detect decoy matches. If you provide your own
      `is_decoy`, this parameter has no effect. Mutually exclusive with `decoy_prefix`.\n",
      "")
  end
  _filter
end

#@contextmanager
def _itercontext(x, **kw)
  begin
    yield x.iterrows().map{ |_, row| row }
  rescue => exception
    yield x
  end
end

# _iter = ChainBase._make_chain(_itercontext)
_iter = ''
qvalues = _make_qvalues(_iter, nil, nil, nil)

filter = _make_filter(_iter, nil, nil, nil, qvalues)
filter.chain = _make_chain(filter, 'filter', true)

begin
  _precalc_fact = Numpy.log((0...20).map{ |m| (1..m).inject(1, :*) })

  def log_factorial(x)
    x = Numpy.array(x)
    pf = _precalc_fact
    m = (x >= pf.size)
    out = Numpy.empty(x.shape)
    out[~m] = pf[x[~m].astype(int)]
    x = x[m]
    out[m] = x * Numpy.log(x) - x + 0.5 * Numpy.log(2 * Numpy.pi * x)
    out
  end

  def _expectation(d, t, p: 0.5)
    return d + 1 if t.nil?
    t = Numpy.array(t, dtype: int)
    m = Numpy.arange(t.max + 1, dtype: int)
    pi = Numpy.exp(_log_pi(d, m, p))
    ((m * pi).cumsum() / pi.cumsum())[t]
  end

  def _confidence_value(conf, d, t, p: 0.5)
    if t.nil?.!
      t = Numpy.array(t, dtype: int)
      m = Numpy.arange(t.max + 1, dtype: int)
    else
      m = Numpy.arange(max(50 * d, 10000))
    end
    log_pi = _log_pi(d, m, p)
    pics = Numpy.exp(log_pi).cumsum()
    Numpy.searchsorted(pics, conf * (t.nil?.! ? pics[T] : 1))
  end
rescue => exception
  def log_factorial(n)
    if n > 10
      return n * Math.log(n) - n + 0.5 * Math.log(2 * Math.PI * n)
    else
      return math.log((1..n).inject(1, :*))
    end
  end

  def _expectation(*a, **k)
    raise NotImplementedError('NumPy required')
  end

  def _confidence_value(*a, **k)
    raise NotImplementedError('NumPy required')
  end
end

def _log_pi_r(d, k, p: 0.5)
  k * Math.log(p) + log_factorial(k + d) - log_factorial(k) - log_factorial(d)
end

def _log_pi(d, k, p: 0.5)
  _log_pi_r(d, k, p) + (d + 1) * Math.log(1 - p)
end

def _count_psms(psms, is_decoy, pep, decoy_prefix, decoy_suffix, is_decoy_prefix, is_decoy_suffix)
  total, decoy = 0, 0
  if pep.nil?.!
    is_decoy = pep
  elsif is_decoy.nil?
    if decoy_suffix.nil?.!
      is_decoy = lambda { |x| is_decoy_suffix(x, decoy_suffix) }
    else
      is_decoy = lambda { |x| is_decoy_prefix(x, decoy_prefix) }
    end
  end
  if is_decoy.instance_of?(String)
    decoy = psms[is_decoy].sum()
    total = psms.shape[0]
  elsif defined?(is_decoy)
    psms.each do |psm|
      total += 1
      d = is_decoy(psm)
      decoy += pep.nil?.! ? d : [0, '', nil, []].include?(d) ? 0 : 1
    end
  else
    if [Sized, Container].include?(is_decoy.class).!
      is_decoy = is_decoy.to_a
    end
    if !!pep
      decoy = is_decoy.sum
    else
      decoy = is_decoy.select{ _1 != 0 }.size
    end
    total = is_decoy.size
  end
  [decoy, total]
end

def _make_fdr(is_decoy_prefix, is_decoy_suffix)
  def fdr(psms: nil, formula: 1, is_decoy: nil, ratio: 1, correction: 0, pep: nil, decoy_prefix: 'DECOY_', decoy_suffix: nil)
    if [1, 2].include?(formula).!
      raise PyteomicsError.new("'formula' must be either 1 or 2.")
    end

    decoy, total = _count_psms(psms, is_decoy, pep, decoy_prefix, decoy_suffix, is_decoy_prefix, is_decoy_suffix)
    return decoy.to_f / total if pep.nil?.!
    tfalse = decoy
    if correction == 1 || (correction == 2 && total.to_f / decoy > 10)
      tfalse += 1
    elsif correction == 2
      p = 1 / (1.0 + ratio)
      tfalse = _expectation(decoy, total - decoy, p)
    elsif 0 < correction && correction < 1
      p = 1 / (1.0 + ratio)
      tfalse = _confidence_value(correction, decoy, total - decoy, p)
    end
    tfalse.to_f / (total - decoy) / ratio if formula == 1
    return (decoy.to_f + tfalse / ratio) / total
  end

  _fix_docstring(fdr, is_decoy: is_decoy_prefix)
  if is_decoy_prefix.nil?
    fdr.__doc__ = fdr.__doc__.replace(
      "\n            .. warning::
      The default function may not work
      with your files, because format flavours are diverse.
      decoy_prefix : str, optional
      If the default `is_decoy` function works for you, this parameter specifies which
      protein name prefix to use to detect decoy matches. If you provide your own
      `is_decoy`, or if you specify `decoy_suffix`, this parameter has no effect.
      Default is `DECOY_`.
      decoy_suffix : str, optional
      If the default `is_decoy` function works for you, this parameter specifies which
      protein name suffix to use to detect decoy matches. If you provide your own
      `is_decoy`, this parameter has no effect. Mutually exclusive with `decoy_prefix`.\n",
      "")
  end
  fdr
end

fdr = _make_fdr(nil, nil)

def _sigma_T(decoy, ratio)
  Math.sqrt((decoy + 1) * (ratio + 1) / (ratio * ratio).to_f)
end

def sigma_T(psms, is_decoy, ratio: 1)
  decoy, total = _count_psms(psms, is_decoy, nil, nil, nil, nil, nil)
  _sigma_T(decoy, ratio)
end

def sigma_fdr(psms: nil, formula: 1, is_decoy: nil, ratio: 1)
  if [1, 2].include?(formula).!
    raise PyteomicsError.new("'formula' must be either 1 or 2.")
  end
  decoy, total = _count_psms(psms, is_decoy, nil, nil, nil, nil, nil)
  sigmaT = _sigma_T(decoy, ratio)
  return sigmaT.to_f / (total - decoy) / ratio if formula == 1
  sigmaT.to_f / total / ratio
end