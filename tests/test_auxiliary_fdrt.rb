require 'minitest/autorun'
require_relative '../rbteomics/auxiliary/math'
require_relative '../rbteomics/auxiliary/structures'
require_relative '../rbteomics/auxiliary/file_helpers'
require_relative '../rbteomics/auxiliary/target_decoy'

PSMS = []
alph = ('A'..'Z').to_a + ('a'..'z').to_a
52.times do |i|
  PSMS << [i, alph[i], 0.01 + i * 0.001]
end
p PSMS

class TestFDRT < Minitest::Test
  def setup
    @is_decoy = lambda { |x| x[1] == x[1].downcase }
    @pep = lambda { |x| x[2] }
  end

  def _run_check(psms, **kwargs)
    is_decoy = kwargs.delete('is_decoy') || @is_decoy
    pep = kwargs.delete('pep') || @pep
    assert_equal Target_decoy::Fdr.call(psms, 'is_decoy' => is_decoy, 'formula' => 1).round(7), 1.0
    assert_equal Target_decoy::Fdr.call(psms, 'is_decoy' => is_decoy, 'formula' => 2).round(7), 1.0
    assert_equal Target_decoy::Fdr.call(psms, 'pep' => pep).round(7), 0.0355
  end

  def test_fdr
    _run_check(PSMS)
  end

  def test_fdr_iter
    # assert_equal Target_decoy::Fdr.call(iter(PSMS), 'is_decoy' => @is_decoy).round(7), 1.0
    # assert_equal Target_decoy::Fdr.call(iter(PSMS), 'pep' => @pep).round(7), 0.0355
    # isd = PSMS.map{ |s, l, p| @is_decoy.call([s, l, p]) }
    # pep = PSMS.map{ |s, l, p| @pep.call([s, l, p]) }
    # assert_equal Target_decoy::Fdr.call(iter(PSMS), 'is_decoy' => iter(isd)).round(7), 1.0
    # assert_equal Target_decoy::Fdr.call(iter(PSMS), 'pep' => iter(pep)).round(7), 0.0355
    
    assert_equal Target_decoy::Fdr.call(PSMS, 'is_decoy' => @is_decoy).round(7), 1.0
    assert_equal Target_decoy::Fdr.call(PSMS, 'pep' => @pep).round(7), 0.0355
    isd = PSMS.map{ |s, l, p| @is_decoy.call([s, l, p]) }
    pep = PSMS.map{ |s, l, p| @pep.call([s, l, p]) }
    assert_equal Target_decoy::Fdr.call(PSMS, 'is_decoy' => isd).round(7), 1.0
    assert_equal Target_decoy::Fdr.call(PSMS, 'pep' => pep).round(7), 0.0355
  end

  def test_fdr_array_str
    # dtype = [PyCall::Tuple.(['score', Numpy.int8]), PyCall::Tuple.(['label', Numpy.str_, 1]), PyCall::Tuple.(['pep', Numpy.float64]), PyCall::Tuple.(['is decoy', Numpy.bool_])]
    dtype = [PyCall::Tuple.(['score', Numpy.int8]), PyCall::Tuple.(['label', Numpy.unicode_, 1]), PyCall::Tuple.(['pep', Numpy.float64]), PyCall::Tuple.(['is decoy', Numpy.bool_])]
    psms_ = Numpy.array(PSMS.map{ |s, l, p| [s, l, p, @is_decoy.call([s, l, p])] }, dtype: dtype)
    _run_check(psms_, 'is_decoy' => 'is decoy', 'pep' => 'pep')
  end
end