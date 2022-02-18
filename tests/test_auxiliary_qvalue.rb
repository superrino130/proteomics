require 'minitest/autorun'
require_relative '../rbteomics/electrochem'
require_relative '../rbteomics/tandem'
require_relative '../rbteomics/version'
require_relative '../rbteomics/auxiliary/structures'

require 'numpy'
require 'pandas'

PSMS = []
alph = ('A'..'Z').to_a + ('a'..'z').to_a
52.times do |i|
  PSMS << [i, alph[i], 0.01 + i * 0.001]
end

class TestQvalue < Minitest::Test
  def setup
    @key = lambda { |x| x[0] }
    @is_decoy = lambda { |x| x[1] == x[1].downcase }
    @pep = lambda { |x| x[2] }
    @psms = PSMS
    @rpsms = Numpy.random.shuffle(@psms)
  end
  
  def _run_check(q, formula)
    assert Numpy.allclose(q['q'][0...26], 0)
    if formula == 2
      assert Numpy.allclose(q['q'][26..], 2 * Numpy.arange(1.0, 27.0) / (26 + Numpy.arange(1, 27)))
    else
      assert Numpy.allclose(q['q'][26..], Numpy.arange(1.0, 27.0) / 26)
    end
    assert Numpy.allclose(q['is decoy'][0...26], 0)
    assert Numpy.allclose(q['is decoy'][26..], 1)
    assert Numpy.allclose(q['score'], Numpy.arange(52))
    setup
    spsms = @psms.sort_by(&key)
    assert Numpy.allclose(spsms.map(&@is_decoy), q['is decoy'])
    assert Numpy.allclose(spsms.map(&@key), q['score'])
    setup
  end

  def _run_check_pep(q)
    assert Numpy.allclose(q['q'], Numpy.arange(0.01, 0.036, 0.0005))
    setup
  end

  def test_qvalues
    q = Target_decoy::Qvalues.call(@psms, 'key' => @key, 'is_decoy' => @is_decoy, 'remove_decoy' => true)
    assert Numpy.allclose(q['q'], 0)
    assert Numpy.allclose(q['is decoy'], 0)
    assert Numpy.allclose(q['score'], Numpy.arange(26))
  end

  # def test_qvalues_pep
  #   q = Target_decoy::Qvalues.call(@psms, 'pep' => @pep)
  #   _run_check_pep(q)
  #   q = Target_decoy::Qvalues.call(@psms, 'pep' => @pep, 'key' = @key)
  #   _run_check_pep(q)
  # end

  # def test_qvalues_with_decoy
  #   q = Target_decoy::Qvalues.call(@psms, 'key' => @key, 'is_decoy' => @is_decoy, 'remove_decoy' => false)
  #   _run_check(q, 2)
  #   q = Target_decoy::Qvalues.call(@psms, 'key' => @key, 'is_decoy' => @is_decoy, 'remove_decoy' => false, 'formula' => 1)
  #   _run_check(q, 1)
  # end

  # def test_qvalues_full_output
  #   q = Target_decoy::Qvalues.call(@psms, 'key' => @key, 'is_decoy' => @is_decoy, 'remove_decoy' => false, 'full_output' => true)
  #   _run_check(q, 2)
  # end

  # def test_qvalues_pep_full_output
  #   q = Target_decoy::Qvalues.call(@psms, 'pep' => @pep, 'full_output' => true)
  #   _run_check_pep(q)
  #   q = Target_decoy::Qvalues.call(@psms, 'key' => @key, 'pep' => @pep, 'full_output' => true)
  #   _run_check_pep(q)
  # end

  # def test_qvalues_from_numpy
  #   dtype = [PyCall::Tuple.(['score', Numpy.int8]), PyCall::Tuple.(['label', Numpy.str_, 1]), PyCall::Tuple.(['pep', Numpy.float64])]
  #   psms = Numpy.array(@psms.to_a, dtype: dtype)
  #   q = Target_decoy::Qvalues.call(psms, 'key' => @key, 'is_decoy' => @is_decoy, 'remove_decoy' => false, 'formula' => 1)
  #   _run_check(q, 1)
  #   q = Target_decoy::Qvalues.call(psms, 'key' => @key, 'is_decoy' => @is_decoy, 'remove_decoy' => false, 'formula' => 1, 'full_output' => true)
  #   _run_check(q, 1)
  #   assert q['psm'].dtype == dtype
  # end

  # def test_qvalues_pep_from_numpy
  #   dtype = [PyCall::Tuple.(['score', Numpy.int8]), PyCall::Tuple.(['label', Numpy.str_, 1]), PyCall::Tuple.(['pep', Numpy.float64])]
  #   psms = Numpy.array(@psms.to_a, dtype: dtype)
  #   q = Target_decoy::Qvalues.call(psms, 'pep' => @pep)
  #   _run_check_pep(q)
  #   q = Target_decoy::Qvalues.call(psms, 'key' => @key, 'pep' => @pep, 'full_output' => true)
  #   _run_check_pep(q)
  #   assert q['psm'].dtype == dtype
  # end

  # def test_qvalues_from_dataframe
  #   dtype = [PyCall::Tuple.(['score', Numpy.int8]), PyCall::Tuple.(['label', Numpy.str_, 1]), PyCall::Tuple.(['pep', Numpy.float64])]
  #   psms = Pandas.DataFrame.new(Numpy.array(@psms.to_a, dtype: dtype))
  #   q = Target_decoy::Qvalues.call(psms, 'key' => @key, 'is_decoy' => @is_decoy, 'remove_decoy' => false, 'formula' => 1)
  #   _run_check(q, 1)
  #   q = Target_decoy::Qvalues.call(psms, 'key' => @key, 'is_decoy' => @is_decoy, 'remove_decoy' => false, 'formula' => 1, 'full_output' => true)
  #   _run_check(q, 1)
  # end

  # def test_qvalues_empty_dataframe
  #   dtype = [PyCall::Tuple.(['score', Numpy.int8]), PyCall::Tuple.(['label', Numpy.str_, 1]), PyCall::Tuple.(['pep', Numpy.float64])]
  #   psms = Pandas.DataFrame.new(Numpy.array([], dtype: dtype))
  #   q = Target_decoy::Qvalues.call(psms, 'key' => @key, 'is_decoy' => @is_decoy, 'remove_decoy' => false, 'formula' => 1)
  #   assert_equal q.shape[0], 0
  #   q = Target_decoy::Qvalues.call(psms, 'key' => @key, 'is_decoy' => @is_decoy, 'remove_decoy' => false, 'formula' => 1, 'full_output' => true)
  #   assert_equal q.shape[0], 0
  # end

  # def test_qvalues_pep_from_dataframe
  #   dtype = [PyCall::Tuple.(['score', Numpy.int8]), PyCall::Tuple.(['label', Numpy.str_, 1]), PyCall::Tuple.(['pep', Numpy.float64])]
  #   psms = Pandas.DataFrame.new(Numpy.array(@psms.to_a, dtype: dtype))
  #   q = Target_decoy.qvalues(psms, 'pep' => @pep)
  #   _run_check_pep(q)
  #   q = Target_decoy.qvalues(psms, 'pep' => @pep, 'full_output' => true)
  #   _run_check_pep(q)
  # end
end