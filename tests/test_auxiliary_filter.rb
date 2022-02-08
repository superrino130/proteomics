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

class TestFDR < Minitest::Test
  def setup
    @key = lambda { |x| x[0] }
    @is_decoy = lambda { |x| x[1] == x[1].downcase }
    @pep = lambda { |x| x[2] }

    @psms = PSMS
    @rpsms = Numpy.random.shuffle(@psms)
  end

  def _run_check(*args, **kwargs)
    key = kwargs['key'] || @key
    is_decoy = kwargs['is_decoy'] || @is_decoy
    f11 = Target_decoy::Filter.call().call(*args, 'key' => key, 'is_decoy' => is_decoy, 'fdr' => 0.5)
    f12 = Target_decoy::Filter.call().call(*args, 'key' => key, 'is_decoy' => is_decoy, 'fdr' => 0.5, 'formula' => 2)
    f21 = Target_decoy::Filter.call().call(*args, 'key' => key, 'is_decoy' => is_decoy, 'fdr' => 0.5, 'remove_decoy' => false, 'formula' => 1)
    f22 = Target_decoy::Filter.call().call(*args, 'key' => key, 'is_decoy' => is_decoy, 'fdr' => 0.5, 'remove_decoy' => false)

    assert_equal f11.shape[0], 26
    assert_equal f12.shape[0], 26
    assert_equal f21.shape[0], 39
    assert_equal f22.shape[0], 34

    Target_decoy::Filter.call(*args, 'key' => key, 'is_decoy' => is_decoy, 'fdr' => 0.5, 'full_output' => false) do |f|
      assert_equal f.to_a.size, 26
    end
    Target_decoy::Filter.call(*args, 'key' => key, 'is_decoy' => is_decoy, 'fdr' => 0.5, 'formula' => 2, 'full_output' => false) do |f|
      assert_equal f.to_a.size, 26
    end
    Target_decoy::Filter.call(*args, 'key' => key, 'is_decoy' => is_decoy, 'fdr' => 0.5, 'remove_decoy' => false, 'formula' => 1, 'full_output' => false) do |f|
      assert_equal f.to_a.size, 39
    end
    Target_decoy::Filter.call(*args, 'key' => key, 'is_decoy' => is_decoy, 'fdr' => 0.5, 'remove_decoy' => false, 'full_output' => false) do |f|
      assert_equal f.to_a.size, 34
    end
  end

  def _run_check_pep(*args, **kwargs)
    key = kwargs.delete('key') || @key
    f11 = Target_decoy::filter(*args, 'key' => @key, 'fdr' => 0.02, **kwargs)
    f12 = Target_decoy::filter(*args, 'fdr' => 0.02, **kwargs)

    assert_equal f11.shape[0], 21
    assert_equal f12.shape[0], 21

    Target_decoy::filter(*args, key=key, fdr=0.02, full_output=False, **kwargs) do |f|
      assert_equal f.to_a.size, 21
    end
    Target_decoy::filter(*args, fdr=0.02, full_output=False, **kwargs) do |f|
      assert_equal f12.to_a.size, 21
    end
  end

  def test_filter
    _run_check(@psms)
  end


end
