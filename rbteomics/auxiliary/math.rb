require_relative './structures'
require 'numpy'

def linear_regression_vertical(x, y: nil, a: nil, b: nil)
  x = Numpy.array(x, copy=nil)
  if y
    y = Numpy.array(y, copy=nil)
  else
    if x.shape.size != 2 || x.shape[-1] != 2
      raise PyteomicsError("f `y` is not given, x.shape should be (N, 2), given: #{x.shape}")
    end
    y = x[0..-1, 1] # 本当は[:,1]
    x = x[0..-1, 0] # 本当は[:,0]
  end
  if !!a && !b
    b = (y - a * x).mean()
  elsif !!a && !!b
    # pass
  else
    m = Numpy.polyfit(x, y, 1)
    a = m[0]
    b = m[1]
  end

  r = Numpy.corrcoef(x, y)[0, 1]
  stderr = (y - a * x - b).std()
  [a, b, r, stderr]
end

def linear_regression(...)
  linear_regression_vertical(...)
end

def linear_regression_perpendicular(x, y: nil)
  x = Numpy.array(x, copy: nil)
  if !!y
    y = Numpy.array(y, copy: nil)
    data = Numpy.hstack([x.reshape([-1, 1]), y.reshape([-1, 1])])
  else
    if x.shape.size != 2 || x.shape[-1] != 2
      raise PyteomicsError("If `y` is not given, x.shape should be (N, 2), given: #{x.shape}")
    end
    data = x
  end
  mu = data.mean(axis=0)
  eigenvectors, eigenvalues, v = Numpy.linalg.svd((data - mu).T, full_matrices=false)
  a = eigenvectors[0][1] / eigenvectors[0][0]
  xm, ym = data.mean(axis=0)
  b = ym - a * xm

  r = Numpy.corrcoef(data[0..-1, 0], data[0..-1, 1])[0, 1]
  stderr = ((data[0..-1, 1] - a * data[0..-1, 0] - b) / Numpy.sqrt(a**2 + 1)).std()

  [a, b, r, stderr]
end