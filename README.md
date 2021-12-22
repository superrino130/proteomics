# proteomics
[Pyteomics](https://pyteomics.readthedocs.io/en/latest/) の `Ruby` 版を目指します。
（仮リポジトリ）
# gem
次のgemが必要です。
```ruby
gem "bio"
```
`fasta` 形式は `biruby` を利用しています。
その他のファイル形式についての実装を計画中です。
# 進捗状況

[EXAMPLE 1: UNRAVELLING THE PEPTIDOME](https://pyteomics.readthedocs.io/en/latest/examples/example_fasta.html) の
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