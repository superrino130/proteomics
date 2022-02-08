# import warnings
# import re

# from lxml import etree
# from sqlalchemy.ext.declarative import declarative_base, DeclarativeMeta
# from sqlalchemy.orm import relationship, backref, object_session
# from sqlalchemy.ext.associationproxy import association_proxy
# from sqlalchemy import (Numeric, Unicode,
#                         Column, Integer, ForeignKey,
#                         UnicodeText, Boolean, event)
# from sqlalchemy import exc as sa_exc
# from sqlalchemy import create_engine
# from sqlalchemy.orm import sessionmaker

# from . import mass
require './mass'
require 'set'

module Unimod
  @model_registry = Set.new

  class SubclassRegisteringDeclarativeMeta < DeclarativeMeta
    def __new__(cls, name, parents, attrs)
      new_type = super
      model_registry << new_type
      new_type
    end
  end
  
  Base = declarative_base(metaclass: SubclassRegisteringDeclarativeMeta)

  @_unimod_xml_download_url = 'http://www.unimod.org/xml/unimod_tables.xml'
  
  begin
    basestring
  rescue => exception
    # basestring = (str, bytes)
    basestring = [String, Integer]
  end
  
  CompositionType = mass.Composition
  
  def simple_repr()
    template = '{self.__class__.__name__}({d})'
    # d = {'%s=%r' % (k, v) for k, v in self.__dict__.items() if not k.startswith('_')}
    d = self.__dict__.map{ |k, v| "#{k}=#{v}" if k.start_with?('_').! }
    template.format(d: d.join(', '))
  end
  
  Base.__repr__ = simple_repr()
  
  def remove_namespace(doc:, namespace:)
    ns = namespace
    nsl = ns.size
    doc.getiterator().each do |elem|
      elem.tag = elem.tag[nsl..] if elem.tag.start_with?(ns)
    end
  end
  
  def preprocess_xml(doc_path)
    tree = etree.parse(doc_path)
    root = tree.getroot()
    root.nsmap.values.each do |ns|
      remove_namespace(tree, ns)
    end
    tree
  end
  
  def _formula_parser(formula, session)
    composition = Composition.new
    formula.split.each do |token|
      m = /(?P<isotope>\d+)?(?P<elemet>[^\(]+)(?:\((?P<count>-?\d+)\))?/.match(token)
      if m.nil?.!
        # isotope, element, count = m[1], m[2], m[3]
        isotope, element, count = m[0], m[1], m[2]
        if count && count != 0
          count = count.to_i
        else
          count = 1
        end
        if isotope
          name = mass._make_isotope_string(element, isotope)
        else
          name = element
        end
        # is_brick = session.query(Brick).filter(Brick.brick == name).first()
        is_brick = session.query(Brick).filter(Brick.brick == name).first
        if si_brick.nil?
          composition[name] += count
        else
          composition += is_brick.composition * count
        end
      end
    end
    composition
  end
  
  def _composition_listener(attr)
    @event.listens_for(attr, 'set')
    def _update_composition_from_formula(target:, value:, oldvalue:, initiator:)
      session = object_session(target)
      return if value.nil? || value == ''
      return if session.nil?
      target.composition = _formula_parser(value, session)
    end
  
    @event.listens_for(attr.class_, 'load')
    def _update_composition_on_load(target, context)
      value = getattr(target, attr.prop.key)
      return if value.nil? || value == ''
      session = object_session(target)
      target.composition = _formula_parser(value, session)
    end
  end
  
  def has_composition(attr_name)
    def decorator(model)
      _composition_listener(getattr(model, attr_name))
      model
    end
     decorator
  end
  
  class HasFullNameMixin()
    def __eq__(other)
      begin
        return self.full_name == other.full_name
      rescue => exception
        return false
      end
    end
  
    def __ne__(other)
      self != other
    end
  end
  
  class AlternativeName < Base
    @__tablename__ = 'AlternativeName'
  
    @_tag_nme = 'alt_names_row'
  
    @classmethod
    def from_tag(cls:, tag:)
      attrib = tag.attrib
      inst = cls(
        id = attrib['record_id'].to_i,
        alt_name = attrib['alt_name'],
        modification_id = attrib['mod_key'].to_i
      )
      inst
    end
  
    @id = Column.new(Integer, primary_key: true)
    # alt_name = column(Unicode(256), index=true)
    alt_name = column(Unicode(256), index: true)
    modification_id = Column(Integer, ForeignKey('Modification.id'), index: true)
  end
  
  class AminoAcid < Base
    include HasFullNameMixin
  
    @__tblename__ = 'AminoAcid'
  
    @_tag_name_ = 'amino_acid_row'
  
    @classmethod
    def from_tag(cls, tag)
      attrib = tag.attrib
      inst = cls(
        id = attrib['record_id'].to_i,
        full_name = attrib['full_name'],
        one_letter = attrib['one_letter'],
        three_letter = attrib['three_letter'],
        num_H = attrib['num_H'].to_i,
        num_O = attrib['num_O'].to_i,
        num_C = attrib['num_C'].to_i,
        num_N = attrib['num_N'].to_i,
        num_S = attrib['num_S'].to_i,
      )
    end
  
    @id = Column(Integer, primary_key: true)
    @num_H = Column(Integer)
    @num_O = Column(Integer)
    @num_C = Column(Integer)
    @num_N = Column(Integer)
    @num_S = Column(Integer)
    @full_name = Column(Unicode(25), index: true)
    @one_letter = Column(Unicode(10), index: true)
    @three_letter = Column(Unicode(10), index: true)
  end
  
  class Classification < Base
    @__tablename__ = 'Classification'
    @_tag_name = 'classifications_row'
  
    @classmethod
    def from_tag(cls, tag)
      attrib = tag.attrib
      inst = cls(
        id = attrib['record_id'].to_i,
        classification = attrib['classification']
      )
    end
  
    @id = Column(Integer, primary_key: true)
    @classification = Column(Unicode(30), index: true)
  end
  
  class Position < Base
    @__tablename__ = 'Position'
    @_tag_name = 'positions_row'
  
    @classmethod
    def from_tag(cls, tag)
      attrib = tag.attrib
      inst = cls(
        @id = attrib['record_id'].to_i,
        @position = attrib['position']
      )
    end
  
    id = Column(Integer, primary_key: true)
    position = Column(Unicode(20), index: true)
  end
  
  class Brick < Base
    include HasFullNameMixin
    @__tablename__ = 'Brick'
    @_tag_name = 'bricks_row'
  
    @classmethod
    def from_tag(cls, tag)
      attrib = tag.attrib
      inst = cls(
        id = attrib['record_id'].to_i,
        brick = attrib['brick'],
        full_name = attrib['full_name']
      )
    end
  
    id = Column(Integer, primary_key: true)
    brick = Column(Unicode(64), index: true)
    full_name = Column(Unicode(128), index: true)
  
    elements = relationship('BrickToElement')
  
    @property
    def composition()
      composition = CompositionType.new()
      self.elements.each do |element_relation|
        symbol = element_relation.element
        isotope, element = /(?P<isotope>\d+)?(?P<element>\S+)/.match(symbol).captures
        if isotope
          isotope = isotope.to_i
          iso_str = mass._make_isotope_string(element, isotope)
        else
          iso_str = element
        end
        count = element_relation.count()
        composition[iso_str] = count
      end
      composition
    end
  end
  
  class Fragment < Base
    @__tablename__ = 'Fragment'
    @_tag_name = 'fragments_row'
  
    @classmethod
    def from_tag(cls, tag)
      attrib = tag.attrib
      inst = cls(
        id = attrib['record_id'].to_i,
        modification_id = attrib['mod_key'].to_i
      )
    end
  
    @id = Column(Integer, primary_key: true)
    @modification_id = Column(Integer, ForeignKey('Modification.id'), index: true)
  
    @_fragment_composition = relationship('FragmentComposition')
  
    #@property
    def composition
      composition = CompositionType()
      session = object_session()
      self._fragment_composition.each do |fragment_composition_relation|
        symbol = fragment_composition_relation.brick_string
        m = /(?P<isotope>\d+)?(?P<element>\S+)/.match(symbol)
        isotope, element = m[1], m[2]
        count = fragment_composition_relation.count
        if count
          count = count.to_i
        else
          count = 1
        end
        if isotope
          name = mass._make_isotope_string(element, isotope)
        else
          name = element
        end
        # is_brick = session.query(Brick).filter(Brick.brick == name).first()
        is_brick = session.query(Brick).select{ Brick.brick == name }.first
        if is_brick.nil?
          composition[name] += count
        else
          composition += is_brick.composition * count
        end
      end
      composition
    end
  end
  
  class FragmentComposition < Base
    @__tablename__ = 'FragmentComposition'
    @_tag_name = 'fragment_comp_row'
  
    @classmethod
    def from_tag(cls, tag)
      attrib = tag.attrib
      inst = cls(
        id = attrib['record_id'].to_i,
        brick_string = attrib['brick'],
        fragment_id = attrib['fragments_key'].to_i,
        count = attrib['num_brick'].to_i
      )
    end
  
    @id = Column(Integer, primary_key: true)
    @brick_string = Column(Unicode(64), ForeignKey(Brick.brick), index: true)
    @fragment_id = Column(Integer, ForeignKey('Fragment.id'), index: true)
    @count = Column(Integer)
  end
  
  class ModificationToBrick < Base
    @__tablename__ = 'ModificationToBrick'
    @_tag_name = 'mod2brick_row'
  
    @classmethod
    def from_tag(cls, tag)
      attrib = tag.attrib
      inst = cls(
        id = attrib['record_id'].to_i,
        brick_string = attrib['brick'],
        modification_id = attrib['mod_key'].to_i,
        count = attrib['num_brick'].to_i
        )
    end
  
    id = Column(Integer, primary_key: true)
    brick_string = Column(Unicode(64), ForeignKey(Brick.brick), index: true)
    modification_id = Column(Integer, ForeignKey('Modification.id'), index: true)
    count = Column(Integer)
  end
  
  class BrickToElement < Base
    @__tablename__ = 'BrickToElement'
    @_tag_name = 'brick2element_row'
  
    @classmethod
    def from_tag(cls, tag)
      attrib = tag.attrib
      inst = cls(
        id = attrib['record_id'].to_i,
        brick_id = attrib['brick_key'].to_i,
        count = attrib['num_element'].to_i,
        element = attrib['element']
        )
    end
  
    id = Column(Integer, primary_key: true)
    brick_id = Column(Integer, ForeignKey(Brick.id), index: true)
    element = Column(Unicode(16), ForeignKey('Element.element'), index: true)
    element_obj = relationship('Element', uselist: false)
    count = Column(Integer)
  end
  
  class Element < Base
    include HasFullNameMixin
    @__tablename__ = 'Element'
    @_tag_name = 'elements_row'
  
    @classmethod
    def from_tag(cls, tag)
      attrib = tag.attrib
      inst = cls(
        id = attrib['record_id'].to_i,
        average_mass = attrib['avge_mass'].to_f,
        monoisotopic_mass = attrib['mono_mass'].to_f,
        full_name = attrib['full_name'],
        element = attrib['element']
        )
    end
  
    id = Column(Integer, primary_key: true)
    average_mass = Column(Numeric(12, 6, asdecimal: false))
    monoisotopic_mass = Column(Numeric(12, 6, asdecimal: false))
    full_name = Column(Unicode(64), index: true)
    element = Column(Unicode(16), index: true)
  end
  
  #@has_composition('_composition')
  class Modification < Base
    include HasFullNameMixin
    @__tablename__ = 'Modification'
  
    @_tag_name = 'modifications_row'
  
    @id = Column(Integer, primary_key: true)
    @username_of_poster = Column(Unicode(128))
    @average_mass = Column(Numeric(12, 6, asdecimal: false), index: true)
    @ex_code_name = Column(Unicode(64), index: true)
    @monoisotopic_mass = Column(Numeric(12, 6, asdecimal: false), index: true)
    @full_name = Column(Unicode(128), index: true)
    @code_name = Column(Unicode(128), index: true)
    @_composition = Column(Unicode(128), index: true)
    @approved = Column(Boolean, index: true)
  
    @notes = relationship('MiscNotesModifications')
    @specificities = relationship('Specificity')
    @bricks = relationship(ModificationToBrick)
    @_fragments = relationship(Fragment)
  
    @_alt_names = relationship(AlternativeName, backref: backref('modification'))
    @alternative_names = association_proxy('_alt_names', 'alt_name')
    @fragments = association_proxy('_fragments', 'composition')
  
    @classmethod
    def from_tag(cls, tag)
      attrib = tag.attrib
      inst = cls(
        id = attrib['record_id'].to_i,
        username_of_poster = attrib['username_of_poster'],
        average_mass = attrib['avge_mass'].to_f,
        monoisotopic_mass = attrib['mono_mass'].to_f,
        ex_code_name = attrib['ex_code_name'],
        code_name = attrib['code_name'],
        full_name = attrib['full_name'],
        approved = attrib['approved'].to_i == 1 ? true : false,
        _composition = attrib['composition']
      )
      tag.each do |note|
        if note.tag == MiscNotesModifications._tag_name
          model_note = MiscNotesModifications._from_tag(note, inst.id)
          inst.notes << model_note if model_note
        end
      end
      inst
    end
  end
  
  class MiscNotesModifications < Base
    @__tablename__ = 'MiscNotesModifications'
    @_tag_name = 'misc_notes'
  
    @id = Column(Integer, primary_key: true)
    @modification_id = Column(Integer, ForeignKey(Modification.id), index: true)
    @text = Column(UnicodeText)
  
    @classmethod
    def _from_tag(cls, tag, modification_id)
      return if if tag.text.nil?
      cls(text: tag.text, modification_id: modification_id)
    end
  end
  
  class Specificity < Base
    @__tablename__ = 'Specificity'
    @_tag_name = 'specificity_row'
  
    @id = Column(Integer, primary_key: true)
    @position_id = Column(Integer, ForeignKey(Position.id), index: true)
    @classification_id = Column(Integer, ForeignKey(Classification.id), index: true)
    @classification = relationship('Classification', uselist: false)
    @amino_acid = Column(Unicode(10), ForeignKey(AminoAcid.one_letter), index: true)
    @modification_id = Column(Integer, ForeignKey(Modification.id), index: true)
    @hidden = Column(Boolean, index: true)
    @group = Column(Integer, index: true)
    @neutral_losses = relationship('SpecificityToNeutralLoss')
  
    @classmethod
    def from_tag(cls, tag)
      attrib = tag.attrib
      inst = cls(
        id = attrib['record_id'].to_i,
        position_id = attrib['position_key'].to_i,
        classification_id = attrib['classifications_key'].to_i,
        hidden = attrib['hidden'].to_i == 1 ? true : false,
        amino_acid = attrib['one_letter'],
        modification_id = attrib['mod_key'].to_i,
      )
    end
  end
  
  class NeutralLoss < Base
    @__tablename__ = 'NeutralLoss'
    @_tag_name = 'neutral_losses_row'
  
    @classmethod
    def from_tag(cls, tag)
      attrib = tag.attrib
      inst = cls(
        id = attrib['record_id'].to_i,
        brick_string = attrib['brick'],
        count = attrib['num_brick'].to_i,
        specificity_id = attrib['spec_key'].to_i
      )
    end
  
    id = Column(Integer, primary_key: true)
    brick_string = Column(Unicode(64), index: true)
    specificity_id = Column(Integer, ForeignKey(Specificity.id), index: true)
    count = Column(Integer)
  end
  
  #@has_composition('_composition')
  class SpecificityToNeutralLoss < Base
    @__tablename__ = 'SpecificityToNeutralLoss'
    @_tag_name = 'spec2nl_row'
  
    @classmethod
    def from_tag(cls, tag)
      attrib = tag.attrib
      inst = cls(
        id = attrib['record_id'].to_i,
        specificity_id = attrib['spec_key'].to_i,
        monoisotopic_mass = attrib['nl_mono_mass'].to_f,
        average_mass = attrib['nl_avge_mass'].to_f,
        is_required_peptide_neutral_loss = attrib['is_req_pep_nl'].to_i == 1 ? true : false,
        is_peptide_neutral_loss = attrib['is_pep_nl'].to_i == 1 ? true : false,
        is_slave = attrib['is_slave_nl'].to_i == 1 ? true : false,
        _composition = attrib['nl_composition']
      )
    end
  
    id = Column(Integer, primary_key: true)
    specificity_id = Column(Integer, ForeignKey(Specificity.id), index: true)
    specificity = relationship(Specificity, uselist: false)
    monoisotopic_mass = Column(Numeric(12, 6, asdecimal: false), index: true)
    average_mass = Column(Numeric(12, 6, asdecimal: false), index: true)
    _composition = Column(Unicode(128))
    is_slave = Column(Boolean, index: true)
    is_peptide_neutral_loss = Column(Boolean, index: true)
    is_required_peptide_neutral_loss = Column(Boolean, index: true)
  end
  
  class CrossreferenceSource < Base
    @__tablename__ = 'CrossreferenceSource'
    @_tag_name = 'xref_sources_row'
  
    @id = Column(Integer, primary_key: true)
    @source = Column(Unicode(64), index: true)
  
    @classmethod
    def from_tag(cls, tag)
      attrib = tag.attrib
      inst = cls(
        id = attrib['record_id'].to_i,
        source = attrib['xref_source']  
      )
    end
  end
  
  class Crossreference < Base
    @__tablename__ = 'Crossreference'
    @_tag_name = 'xrefs_row'
  
    @id = Column(Integer, primary_key: true)
    @source_id = Column(Integer, ForeignKey(CrossreferenceSource.id), index: true)
    @source = relationship(CrossreferenceSource, uselist: false)
    @url = Column(Unicode(128))
    @modification_id = Column(Integer, ForeignKey(Modification.id), index: true)
    @text = Column(UnicodeText)
  
    @classmethod
    def from_tag(cls, tag)
      attrib = tag.attrib
      inst = cls(
        inst.id = attrib['record_id'].to_i,
        inst.url = attrib['xref_url'],
        inst.source_id = attrib['xref_source_key'].to_i,
        modification_id = attrib['mod_key'].to_i  
      )
      text = []
      tag.getchildren().each do |node|
        if node.tag == 'xref_text'
          text << node.text if node.text
        end
      end
      inst.text = text.join("\n")
      inst
    end
  end
  
  def load(doc_path, output_path: 'sqlite://')
    tree = preprocess_xml(doc_path)
    engine = create_engine(output_path)
    Base.metadata.create_all(engine)
    session = sessionmaker(bind: engine, autoflush: false)()
    model_registry.each do |model|
      if hasattr(model, '_tag_name') && hasattr(model, 'from_tag')
        tree.iterfind('.//' + model._tag_name).each do |tag|
          session << model.from_tag(tag)
        end
      end
    end
    session
  end
  
  def session(path='sqlite:///unimod.db')
    engine = create_engine(path)
    Base.metadata.create_all(engine)
    session = sessionmaker(bind: engine, autoflush: false)()
  end
  
  class Unimod
    def initialize(...)
      __init__(...)
    end
  
    def __init__(path)
      if path.nil?
        @path = nil
        @session = load(_unimod_xml_download_url)
      else
        @path = path
        begin
          @session = session(path)
          raise Exception() if self.session.query(Modification).first.nil?
        rescue => exception
          @session = load(_unimod_xml_download_url, path)
          @session.query(Modification).first
        end
      end
    end
  
    def get(identifier, strict: true)
      if identifier.instance_of?(Integer)
        mod = @session.query(Modification).get(identifier)
        raise KeyError(identifier) if mod.nil?
        return mod
      elsif identifier.instance_of?(basestring)
        if strict
          mod = @session.query(Modification).select{
            (Modification.full_name == identifier) ||
            (Modification.code_name == identifier) ||
            (Modification.ex_code_name == identifier) }.first
          if mod.nil?
            alt_name = @session.query(AlternativeName).select{
              AlternativeName.alt_name == identifier }.first
            raise KeyError(identifier) if alt_name.nil?
            mod = alt_name.modification
          end
          return mod
        else
          qname = "%%#{identifier}%%"
          mod = @session.query(Modification).select{
            (Modification.full_name.like(qname)) ||
            (Modification.code_name.like(qname)) ||
            (Modification.ex_code_name.like(qname)) }.first
          if mod.nil?
            alt_name = @session.query(AlternativeName).select{
              AlternativeName.alt_name.like(qname) }.first
            raise KeyError(identifier) if alt_name.nil?
            mod = alt_name.modification
          end
          return mod
        end
      end
    end
  
    @by_title = by_name = get
    @__getitem__ = get
  
    #@property
    def mods()
      @session.query(Modification).all?
    end
  
    def __iter__()
      # self.session.query(Modification).yield_per(1000) sqlデータを逐次処理
      @session.query(Modification).yield_per(1000)
    end
  end
end
