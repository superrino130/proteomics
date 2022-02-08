require 'minitest/autorun'
require 'tempfile'
require_relative '../rbteomics/electrochem'
require_relative '../rbteomics/tandem'
require_relative '../rbteomics/version'
require_relative '../rbteomics/auxiliary/structures'
require_relative '../rbteomics/auxiliary/file_helpers'

require 'numpy'
require 'pandas'

PSMS = []
alph = ('A'..'Z').to_a + ('a'..'z').to_a
52.times do |i|
  PSMS << [i, alph[i], 0.01 + i * 0.001]
end

class TestUseIndex < Minitest::Test
  class MockFile
    def initialize(seekable, mode)
      @mode = mode if mode.nil?.!
      @seekable = seekable if seekable.nil?.!
    end

    def seekable

    end
  end


  def _check_file_object(fo, value)

    result = _check_use_index(fo, nil, nil)
    assert_equal result, value
  end

  def test_str_name
    [false, true].each do |ui|
      [false, true].each do |default|
        assert_equal _check_use_index('test.mgf', ui, default), ui
      end
    end
  end

  def test_textfile
    f = './tests/test.fasta'
    _check_file_object(f, false)
  end

  def test_binfile
    f = './tests/test.mgf'
    # test.mgf is not binary file!
    # _check_file_object(f, true)
    _check_file_object(f, false)
  end

  def test_tmpfile_text; end

  def test_tmpfile_bin; end

  def test_stringio
    _check_use_index('test', nil, nil)
  end

  def test_error_not_seekable
    source = MockFile.new(false, 'rb')
    p [67, _check_use_index(source, nil, nil)]
  end

  def test_warning_not_seekable
    source = MockFile.new(false, 'r')
    p [72, _check_use_index(source, true, nil)]
  end

  def test_warning_wrong_mode
    ['rb', 'r'].each do |m|
      source = MockFile.new(true, m)
      p [78, _check_use_index(source, m.include?('b').!, nil)]
    end
  end

  def test_warning_no_mode
    source = MockFile.new(nil, nil)
    p [84, _check_use_index(source, nil, nil)]
  end

  
end