require 'minitest/autorun'
require_relative '../rbteomics/auxiliary/math'
require_relative '../rbteomics/auxiliary/structures'

class TestRegression < Minitest::Test
  def setup
    @x = [1, 2, 3]
    @y = [3, 5, 7]
    @a = 2
    @b = 1
    @r = 1
    @stderr = 0
  end

  def _test_linreg(result)
    a, b, r, stderr = result
    assert_equal a.round(7), @a
    assert_equal b.round(7), @b
    assert_equal r.round(7), @r
    assert_equal stderr.round(7), @stderr
  end

  def test_linear_regression_simple
    result = linear_regression(@x, y: @y)
    _test_linreg(result)
  end

  def test_linear_regression_simple_vertical
    result = linear_regression_vertical(@x, y: @y)
    _test_linreg(result)
  end

  def test_linear_regression_simple_perpendicular
    result = linear_regression_perpendicular(@x, y: @y)
    _test_linreg(result)
  end

  def test_linear_regression_no_y_list
    x = @x.zip(@y)
    result = linear_regression(x)
    _test_linreg(result)
  end

  def test_linear_regression_no_y_list_vertical
    x = @x.zip(@y)
    result = linear_regression_vertical(x)
    _test_linreg(result)
  end

  def test_linear_regression_no_y_list_perpendicular
    x = @x.zip(@y)
    result = linear_regression_perpendicular(x)
    _test_linreg(result)
  end

  def test_linear_regression_no_y_arr
    x = Numpy.array(@x.zip(@y))
    result = linear_regression(x)
    _test_linreg(result)
  end

  def test_linear_regression_no_y_arr_vertical
    x = Numpy.array(@x.zip(@y))
    result = linear_regression_vertical(x)
    _test_linreg(result)
  end

  def test_linear_regression_no_y_arr_perpendicular
    x = Numpy.array(@x.zip(@y))
    result = linear_regression_perpendicular(x)
    _test_linreg(result)
  end

  def test_linear_regression_shape_exception_vertical
    e = assert_raises PyteomicsError do
      linear_regression_vertical(@x)
    end
    assert e.message.include?("If `y` is not given, x.shape should be (N, 2), given:")
  end

  def test_linear_regression_shape_exception_perpendicular
    e = assert_raises PyteomicsError do
      linear_regression_perpendicular(@x)
    end
    assert e.message.include?("If `y` is not given, x.shape should be (N, 2), given:")
  end
end