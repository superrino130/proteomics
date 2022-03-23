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

module Target_decoy
  module_function
  def _fix_docstring(f, **defaults)
    # defaults.each do |argname, v|
    #   if v.nil?.!
        # f.__doc__ = re.sub('{} : .*'.format(argname),
        #                          lambda m: m.group() + ', optional', f.__doc__)
    #   end
    # end
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
    tfalse = cumsum.dup
    ind = Numpy.arange(1.0, scores.shape[0] + 1.0, dtype: Numpy.float64)
  
    if peps
      q = cumsum / ind
    else
      if correction.instance_of?(Integer)
        if correction == 1
          tfalse += 1
        elsif correction == 2
          ps = 1.0 / (1.0 + ratio)
          targ = ind - cumsum
          tfalse.size.times do |i|
            tfalse[i] = _expectation(cumsum[i], targ[i], ps)
          end
        end
      elsif 0 < correction && correction < 1
        ps = 1.0 / (1.0 + ratio)
        targ = ind - cumsum
        tfalse.size.times do |i|
          tfalse[i] = _confidence_value(correction, cumsum[i], targ[i], ps)
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
  
    (scores.size - 1).dowto(1) do |i|
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
    q_label = kwargs.include?('q_label') ? kwargs['q_label'] : kwargs['q_label'] = 'q'
    score_label = kwargs.include?('score_label') ? kwargs['score_label'] : kwargs['score_label'] = 'score'
    keyf = psms.apply(keyf, axis: 1).to_a if callable?(keyf)
    isdecoy = psms.apply(isdecoy, axis: 1).to_a if callable?(isdecoy)
    if keyf.instance_of?(BaseString)
      if psms.shape[0]
        psms[score_label] = keyf
      else
        psms[score_label] = []
      end
      keyf = kwargs['score_label']
    end
    if isdecoy.instance_of?(BaseString)
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
        fields = [PyCall::Tuple.([keyf, Numpy.float64]), PyCall::Tuple.([isdecoy, Numpy.bool_]), PyCall::Tuple.([q_label, Numpy.float64])]
      else
        fields = [PyCall::Tuple.([isdecoy, Numpy.float64]), PyCall::Tuple.([q_label, Numpy.float64])]
      end
      dtype = Numpy.dtype(fields)
    end
    psms.sort_values([keyf, isdecoy], ascending: [reverse.!, true], inplace: true)
  
    if psms.shape[0].!
      if full
        psms[q_label] = []
        return psms
      else
        return Numpy.array([], dtype: dtype)
      end
    end
  
    q = _calculate_qvalues(psms[keyf].values, psms[isdecoy].values, peps.nil?.!, **kwargs)
    if remove_decoy
      q = q[~psms[isdecoy].values]
      psms = psms[~psms[isdecoy]].dup
    end
    if full.!
      psms_ = Numpy.empty_like(q, dtype: dtype)
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
      PyCall::Tuple.([score_label, Numpy.float64]),
      PyCall::Tuple.([_decoy_or_pep_label(**kwargs), peps.nil? ? Numpy.bool_ : Numpy.float64]),
      PyCall::Tuple.([q_label, Numpy.float64])
    ]
    
    if full
      # dtypes = Set.new(args.map{ getattr(_1, 'dtype', nil) })
      dtypes = Set.new(args.map{ _1.respond_to?('dtype') ? _1.dtype : nil })
      if dtypes.size == 1 || dtypes.all?{ _1 }
        psm_dtype = dtypes.to_a[-1]
        dtypes.delete(dtypes.to_a[-1])
      else
        psm_dtype = Numpy.object_
      end
      dtype = Numpy.dtype(fields + PyCall::Tuple.(['psm', psm_dtype]))
    else
      dtype = Numpy.dtype(fields)
    end
    dtype
  end
  
  Make_qvalues = lambda do |read, is_decoy_prefix, is_decoy_suffix, key|
    @read = read
    @mq_is_decoy_prefix = is_decoy_prefix
    @mq_key = key
    @mqvalues = lambda do |*args, **kwargs|
      #@_keepstate
      def in_get_scores(*args, **kwargs)
        scores = []
        f = @read.new(*args, **kwargs)
          f.each_with_index do |psm, i|
            row = []
            [keyf, isdecoy].each do |func|
              if func.instance_of?(String)
                row << psm[func]
              elsif callable?(func)
                row << func(psm)
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
      def get_scores(*args, **kwargs)
        File_helpers._keepstate(:in_get_scores)
        in_get_scores(*args, **kwargs)
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
  
      if callable?(keyf).! && sizeable?(keyf.class).! && includeable?(keyf.class).!
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
      
      if callable?(isdecoy).! && sizeable?(isdecoy.class).! && includeable?(isdecoy.class).!
        isdecoy = Numpy.array(isdecoy.to_a)
      end
      remove_decoy = kwargs['remove_decoy'] || false
      decoy_or_pep_label = _decoy_or_pep_label(**kwargs)
      score_label = kwargs.include?('score_label') ? kwargs['score_label'] : kwargs['score_label'] = 'score'
      q_label = kwargs.include?('q_label') ? kwargs['q_label'] : kwargs['q_label'] = 'q'
      dtype = _construct_dtype(*args, **kwargs)
      full = kwargs['full_output'] || false
      arr_flag = false
      psms = nil
  
      if args.empty?.! && args.all?{ _1.is_a?(Pandas::DataFrame) }
        psms = Pandas.concat(args)
        return _qvalues_df(psms, keyf, isdecoy, **kwargs)
      end
  
      if args.empty?.! && args.all?{ _1.is_a?(Numpy.ndarray)}.!
        keyf = Itemgetter.new(keyf) if keyf.instance_of?(String)
        isdecoy = Itemgetter.new(isdecoy) if isdecoy.instance_of?(String)
        peps = Itemgetter.new(peps) if peps.instance_of?(String)
      end
  
      if callable?(keyf) || callable?(isdecoy)
        kwargs.delete('full_output')
        scores = Numpy.array(get_scores(*args, **kwargs), dtype: dtype)
      else
        if args.all?{ _1.instance_of?(Numpy.ndarray) }
          psms = Numpy.concatenate(args)
        end
        if keyf.instance_of?(BaseString).!
          keyf = Numpy.array(keyf)
          arr_flag = true
        end
        if isdecoy.instance_of?(BaseString).!
          isdecoy = Numpy.array(isdecoy)
          arr_flag = true
        end
  
        if arr_flag
          scores = Numpy.empty(hasattr(keyf, 'size') ? keyf.size : isdecoy.size, dtype=dtype)
          [keyf, isdecoy].zip([score_label, decoy_or_pep_label]).each do |func, label|
            if func.instance_of?(BaseString).!
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
          newdt = fields.map{ PyCall::Tuple.([_1, psms.dtype.fields[_1][0]]) } + 
            extra.map{ PyCall::Tuple.([_1, Numpy.float64]) } + PyCall::Tuple.([q_label, Numpy.float64])
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
  
    _fix_docstring(@mqvalues, 'is_decoy' => @is_decoy_prefix, 'key' => key)
    if read == @_iter
      # qvalues.__doc__ = qvalues.__doc__.replace(
      #   "positional args : file or str
      #   Files to read PSMs from. All positional arguments are treated as
      #   files.",
      #   "positional args : iterables
      #   Iterables to read PSMs from. All positional arguments are chained."
      #   ).replace(
      #   "\n            .. warning::
      #   The default function may not work
      #   with your files, because format flavours are diverse.
      #   decoy_prefix : str, optional
      #   If the default `is_decoy` function works for you, this parameter specifies which
      #   protein name prefix to use to detect decoy matches. If you provide your own
      #   `is_decoy`, or if you specify `decoy_suffix`, this parameter has no effect.
      #   Default is `DECOY_`.
      #   decoy_suffix : str, optional
      #   If the default `is_decoy` function works for you, this parameter specifies which
      #   protein name suffix to use to detect decoy matches. If you provide your own
      #   `is_decoy`, this parameter has no effect. Mutually exclusive with `decoy_prefix`.\n",
      #   "")
    end
    @mqvalues
  end
  
  Make_filter = lambda do |read, is_decoy_prefix, is_decoy_suffix, key, qvalues|
    @mfqvalues = qvalues
    @is_decoy_prefix = is_decoy_prefix
    filter = lambda do |*args, **kwargs|
      if kwargs.include?('fdr').!
        raise PyteomicsError.new('Keyword argument required: fdr')
      end
  
      args = args.map{ |arg| sizeable?(arg.class).! && includeable?(arg.class).! ? arg.to_a : arg }
      peps = kwargs['pep']
      if peps.nil?
        remove_decoy = kwargs.delete('remove_decoy') || true
        scores = @mfqvalues.call(*args, 'remove_decoy' => remove_decoy, **kwargs)
      else
        scores = @mfqvalues.call(*args, **kwargs)
      end
      keyf = kwargs.delete('key') || key
      keyf = peps if keyf.nil?
      reverse = kwargs.delete('reverse') || false
      better = [op.lt, op.gt][reverse == false ? 0 : 1]
      if kwargs.include?('is_decoy').!
        if kwargs.include?('decoy_suffix')
          isdecoy = lambda { |x| is_decoy_suffix(x, kwargs['decoy_suffix']) }
        elsif kwargs.include?('decoy_prefix')
          isdecoy = lambda { |x| is_decoy_prefix(x, kwargs['decoy_prefix']) }
        else
          isdecoy = @is_decoy_prefix
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
        i = i[0] if sizeable?(i)
      rescue => exception
        i = scores['q'].bsearch_index{ _1 > fdr } || scores['q'].size
      end
      if kwargs.delete('full_output')
        if scores.instance_of?(Pandas.DataFrame)
          return scores.iloc[0..i]
        elsif callable?(keyf) || callable?(isdecoy)
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
        File.open(*args) do |f|
          f.zip(scores).each do |p1, s|
            if peps.nil?.! || remove_decoy.! || s[decoy_or_pep_label].!
              if better(s[label], cutoff)
                p1
                yield
              end
            end
          end
        end
      end
      return out()
    end
  
    _filter = lambda do |*args, **kwargs|
      if kwargs.include?('full_output').! || [0, '', nil, []].include?(kwargs.delete('full_output')).!
        return filter.call(*args, 'full_output' => true, **kwargs)
      end
      return IteratorContextManager.__init__(*args, 'parser_func' => filter, **kwargs)
    end
  
    _fix_docstring(_filter, 'is_decoy' => is_decoy_prefix, 'key' => key)
    # if read == @_iter
      # _filter.__doc__ = _filter.__doc__.replace(
      #   "positional args : file or str
      #   Files to read PSMs from. All positional arguments are treated as
      #   files.",
      #   "positional args : iterables
      #   Iterables to read PSMs from. All positional arguments are chained."
      #   ).replace(
      #   "\n            .. warning::
      #   The default function may not work
      #   with your files, because format flavours are diverse.
      #   decoy_prefix : str, optional
      #   If the default `is_decoy` function works for you, this parameter specifies which
      #   protein name prefix to use to detect decoy matches. If you provide your own
      #   `is_decoy`, or if you specify `decoy_suffix`, this parameter has no effect.
      #   Default is `DECOY_`.
      #   decoy_suffix : str, optional
      #   If the default `is_decoy` function works for you, this parameter specifies which
      #   protein name suffix to use to detect decoy matches. If you provide your own
      #   `is_decoy`, this parameter has no effect. Mutually exclusive with `decoy_prefix`.\n",
      #   "")
    # end
    _filter
  end
  
  #@contextmanager
  Itercontext = lambda do |x, **kw|
    begin
      yield x.iterrows().map{ |_, row| row }
    rescue => exception
      yield x
    end
  end
  
  # @_iter = ChainBase._make_chain('_itercontext')
  @_iter = File_helpers::ChainBase._make_chain(:Itercontext, Itercontext)
  # def qvalues(read = @_iter, key: nil, is_decoy: nil, remove_decoy: nil)
  #   @qvalues = _make_qvalues(read, key, is_decoy, remove_decoy)
  # end
  Qvalues = Make_qvalues.call(@_iter, nil, nil, nil)
  # Filter = lambda do |x = nil|
  #   if x.nil?
  #     _make_filter.call(@_iter, nil, nil, nil, qvalues)
  #   elsif x == 'chain' || x == :chain
  #     y = _make_filter.call(@_iter, nil, nil, nil, qvalues)
  #     _make_chain(y, 'filter', full_output: true)
  #   end
  # end
  Filter = lambda do |x = nil|
    if x.nil?
      Make_filter.call(@_iter, nil, nil, nil, Qvalues)
    elsif x == 'chain' || x == :chain
      y = Make_filter.call(@_iter, nil, nil, nil, Qvalues)
      Male_chain(y, 'filter', 'full_output' => true)
    end
  end  
  # filter = _make_filter.call(@_iter, nil, nil, nil, qvalues)
  # filter.chain = _make_chain(filter, 'filter', full_output: true)
  
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
    if is_decoy.instance_of?(BaseString)
      decoy = psms[is_decoy].sum()
      total = psms.shape[0]
    elsif callable?(is_decoy)
      psms.each do |psm|
        total += 1
        d = is_decoy.call(psm)
        decoy += pep.nil?.! ? d : [0, '', nil, false, []].include?(d) ? 0 : 1
      end
    else
      if sizeable?(is_decoy.class).! && includeable?(is_decoy.class).!
        is_decoy = is_decoy.to_a
      end
      if pep.nil?.!
        decoy = is_decoy.sum
      else
        decoy = is_decoy.select{ ['', 0, nil, false].include?(_1).! }.size
      end
      if is_decoy.is_a?(Enumerator)
        total = is_decoy.map{}.size
      else
        total = is_decoy.size
      end
    end
    [decoy, total]
  end
  
  Make_fdr = lambda do |is_decoy_prefix, is_decoy_suffix|
    @mf_is_decoy_prefix = is_decoy_prefix
    @mf_is_decoy_suffix = is_decoy_suffix
    fdr = lambda do |psms, **kwargs|
      formula = kwargs['formula'] || 1
      is_decoy = kwargs['is_decoy'] || nil
      ratio = kwargs['ratio'] || 1
      correction = kwargs['correction'] ||  0
      pep = kwargs['pep'] || nil
      decoy_prefix = kwargs['decoy_prefix'] || 'DECOY_'
      decoy_suffix = kwargs['decoy_suffix'] || nil
      if [1, 2].include?(formula).!
        raise PyteomicsError.new("'formula' must be either 1 or 2.")
      end
  
      decoy, total = _count_psms(psms, is_decoy, pep, decoy_prefix, decoy_suffix, @mf_is_decoy_prefix, @mf_is_decoy_suffix)
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
      return tfalse.to_f / (total - decoy) / ratio if formula == 1
      return (decoy.to_f + tfalse / ratio) / total
    end
  
    _fix_docstring(fdr, 'is_decoy' => is_decoy_prefix)
    # if is_decoy_prefix.nil?
      # fdr.__doc__ = fdr.__doc__.replace(
      #   "\n            .. warning::
      #   The default function may not work
      #   with your files, because format flavours are diverse.
      #   decoy_prefix : str, optional
      #   If the default `is_decoy` function works for you, this parameter specifies which
      #   protein name prefix to use to detect decoy matches. If you provide your own
      #   `is_decoy`, or if you specify `decoy_suffix`, this parameter has no effect.
      #   Default is `DECOY_`.
      #   decoy_suffix : str, optional
      #   If the default `is_decoy` function works for you, this parameter specifies which
      #   protein name suffix to use to detect decoy matches. If you provide your own
      #   `is_decoy`, this parameter has no effect. Mutually exclusive with `decoy_prefix`.\n",
      #   "")
    # end
    fdr
  end
  
  Fdr = Make_fdr.call(nil, nil)
  
  def _sigma_T(decoy, ratio)
    Math.sqrt((decoy + 1) * (ratio + 1) / (ratio * ratio).to_f)
  end
  
  def sigma_T(psms, is_decoy: nil, ratio: 1)
    decoy, total = _count_psms(psms, is_decoy, nil, nil, nil, nil, nil)
    _sigma_T(decoy, ratio)
  end
  
  def sigma_fdr(psms, formula: 1, is_decoy: nil, ratio: 1)
    if [1, 2].include?(formula).!
      raise PyteomicsError.new("'formula' must be either 1 or 2.")
    end
    decoy, total = _count_psms(psms, is_decoy, nil, nil, nil, nil, nil)
    sigmaT = _sigma_T(decoy, ratio)
    return sigmaT.to_f / (total - decoy) / ratio if formula == 1
    sigmaT.to_f / total / ratio
  end
end
