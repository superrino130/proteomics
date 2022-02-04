require 'minitest/autorun'
require_relative '../rbteomics/auxiliary/file_helpers'

class TestOffsetIndex < Minitest::Test
  def setup
    @sequence = (0...10).map{ [_1.to_s, _1] }
    @index = OffsetIndex.new(@sequence)
  end

  def test_index_sequence
    assert_equal @index.index_sequence, @sequence
  end

  def test_find
    assert_equal @index.find('3'), 3
  end

  def test_from_index
    assert_equal @index.from_index(3), '3'
    assert_equal @index.from_index(4, include_value: true), ['4', 4]
  end

  def test_from_slice
    assert_equal @index.from_slice(1...3), ['1', '2']
    assert_equal @index.from_slice(1...3, include_value: true), @sequence[1...3]
  end

  def test_between
    assert_equal @index.between('1', '3'), ['1', '2', '3']
    assert_equal @index.between('1', '3', include_value: true), [['1', 1], ['2', 2], ['3', 3]]
    assert_equal @index.between('3', '1'), ['1', '2', '3']
    assert_equal @index.between(nil, '3'), ['0', '1', '2', '3']
    assert_equal @index.between('8', nil), ['8', '9']
  end
end