# import pylab
# import numpy as np
# from .auxiliary import linear_regression, PyteomicsError
# from . import parser, mass, mgf

# try:
#     import spectrum_utils.spectrum as sus
#     import spectrum_utils.plot as sup
# except ImportError:
#     sus = sup = None

require 'numpy'
require 'matplotlib/pyplot'
require_relative 'auxiliary/structures'
require_relative 'auxiliary/file_helpers'
require_relative 'auxiliary/utils'

require 'set'

module Pylab_aux
  module_function

  def plot_line(a, b, *args, **kwargs)
    xlim = kwargs['xlim'] || nil

    # Not started
  end

  def scatter_trend(x, **kwargs)
    y = kwargs['y'] || nil

    # Not started
  end

  def plot_function_3d(x, y, function, **kwargs)
    # Not started
  end

  def plot_function_contour(x, y, function, **kwargs)
    # Not started
  end

  def plot_qvalue_curve(qvalues, *args, **kwargs)
    # Not started
  end

  Default_plot_spectrum = lambda do |spectrum, *args, **kwargs|
    # Not started
  end

  Spectrum_utils_plot = lambda do |spectrum, *args, **kwargs|
    # Not started
  end

  Spectrum_utils_iplot = lambda do |spectrum, *args, **kwargs|
    # Not started
  end

  @@_plot_backends = {
    'default' => Default_plot_spectrum,
    'spectrum_utils' => Spectrum_utils_plot,
    'spectrum_utils.iplot' => Spectrum_utils_iplot,
  }

  def plot_spectrum(spectrum, *args, **kwargs)
    # Not started
  end

  Default_annotate_spectrum = lambda do |spectrum, peptide, *args, **kwargs|
    # Not started
  end

  def _get_precursor_charge(spectrum)
    # Not started
  end

  def _spectrum_utils_create_spectrum(spectrum, peptide, *args, **kwargs)
    # Not started
  end

  def _spectrum_utils_annotate_spectrum(spectrum, peptide, *args, **kwargs)
    # Not started
  end

  class SpectrumUtilsColorScheme
    # Not s
  end

  def _spectrum_utils_parse_sequence(sequence, aa_mass: nil)
    # Not s
  end

  Spectrum_utils_annotate_plot = lambda do |spectrum, peptide, *args, **kwargs|
    # Not s
  end

  Spectrum_utils_annotate_iplot = lambda do |spectrum, peptide, *args, **kwargs|
    # N s
  end

  @@_annotation_backends = {
    'default' => Default_annotate_spectrum,
    'spectrum_utils' => Spectrum_utils_annotate_plot,
    'spectrum_utils.iplot' => Spectrum_utils_annotate_iplot,
  }

  def annotate_spectrum(spectrum, peptide, *args, **kwargs)
    bname = kwargs.delete('backend') || 'default'
    backend = @@_annotation_backends[bname]
    if backend.nil?
      raise PyteomicsError.new("Unknown backend name: #{bname}. Should be one of: #{@@_annotation_backends.to_a.join('; ')}.")
    end

    Matplotlib::Pyplot.xlabel(kwargs.delete('xlabel') || 'm/z')
    Matplotlib::Pyplot.ylabel(kwargs.delete('ylabel') || 'intensity')
    Matplotlib::Pyplot.title(kwargs.delete('title') || '')
    backend.call(spectrum, peptide, *args, **kwargs)
  end

  Spectrum_utils_mirror = lambda do |spec_top, spec_bottom, **kwargs|
    spectrum_kws = kwargs['spectrum_kws'] || nil
    ax = kwargs['ax'] || nil
    # n s
  end

  Spectrum_utils_iplot_mirror = lambda do |spec_top, spec_bottom, **kwargs|
    spectrum_kws = kwargs['spectrum_kws'] || nil
    # n s
  end

  @@_mirror_backends = {
    'spectrum_utils' => Spectrum_utils_mirror,
    'spectrum_utils.iplot' => Spectrum_utils_iplot_mirror,
  }

  def mirror(spec_top, spec_bottom, **kwargs)
    peptide = kwargs['peptide'] || nil
    spectrum_kws = kwargs['spectrum_kws'] || nil
    ax = kwargs['ax'] || nil

    # n s
  end
end
