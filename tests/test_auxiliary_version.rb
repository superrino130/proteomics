require 'minitest/autorun'
require_relative '../rbteomics/version'

class TestVersion < Minitest::Test
  def setup; end
  
  def test_short_version
    assert_equal Version.versioninfo('1.2'), ['1', '2', nil, nil, nil]
  end

  def test_longer_version
    assert_equal Version.versioninfo('1.2.3'), ['1', '2', '3', nil, nil]
  end

  def test_short_dev_version
    assert_equal Version.versioninfo('1.2dev3'), ['1', '2', nil, 'dev', '3']
  end

  def test_longer_dev_version
    assert_equal Version.versioninfo('1.2.3dev4'), ['1', '2', '3', 'dev', '4']
  end
end