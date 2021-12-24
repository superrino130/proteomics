# proteomics
[Pyteomics](https://pyteomics.readthedocs.io/en/latest/) の `Ruby` 版を目指します。
（仮リポジトリ）
# gem
次のgemが必要です。
```ruby
gem "bio"
```
`fasta` 形式は `bioruby` を利用しています。
その他のファイル形式についての実装を計画中です。
# 進捗状況

[EXAMPLE 1: UNRAVELLING THE PEPTIDOME - Pyteomics](https://pyteomics.readthedocs.io/en/latest/examples/example_fasta.html) の
```python
print('Cleaving the proteins with trypsin...')
unique_peptides = set()
with gzip.open('yeast.fasta.gz', mode='rt') as gzfile:
    for description, sequence in fasta.FASTA(gzfile):
        new_peptides = parser.cleave(sequence, 'trypsin')
        unique_peptides.update(new_peptides)
print('Done, {0} sequences obtained!'.format(len(unique_peptides)))
```
の部分まで進んでいます。
```ruby
Done, 188701 sequences obtained!
```
と表示されます。

その後
```python
peptides = [{'sequence': i} for i in unique_peptides]

print('Parsing peptide sequences...')
for peptide in peptides:
    peptide['parsed_sequence'] = parser.parse(
        peptide['sequence'],
        show_unmodified_termini=True)
    peptide['length'] = parser.length(peptide['parsed_sequence'])
print('Done!')

peptides = [peptide for peptide in peptides if peptide['length'] <= 100]
```
を経て、
```ruby
peptides.size = 188548
```
更に、
```
print('Calculating the mass, charge and m/z...')
for peptide in peptides:
    peptide['charge'] = int(round(
        electrochem.charge(peptide['parsed_sequence'], pH=2.0)))
    peptide['mass'] = mass.calculate_mass(peptide['parsed_sequence'])
    peptide['m/z'] = mass.calculate_mass(peptide['parsed_sequence'],
        charge=peptide['charge'])
print('Done!')
```
まで動作します。