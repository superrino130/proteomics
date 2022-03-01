require 'net/http'
require 'json'
require 'matplotlib/pyplot'
plt = Matplotlib::Pyplot

url = 'https://repository.jpostdb.org/proxi/spectra?usi=mzspec:PXD005159:150211tk04-whole_2m8h-3.wizd:scan:2&resultType=full'

uri = URI.parse(url)
res = Net::HTTP.get_response(uri)
raise "Not Connected by Code #{res.code.to_i}." if res.code.to_i != 200

json_data = JSON.parse(res.body)

mzs = []
its = []
json_data[0]['mzs'].zip(json_data[0]['intensities']).each do |m, i|
  case m.to_i
  when 0...1500
    mzs << m.to_i
    its << i.to_i
  end
end

plt.figure()
plt.title('PXD005159:150211tk04-whole_2m8h-3')
plt.xlabel('m/z, Th')
plt.ylabel('Intensity, rel. units')

plt.bar(mzs, its, width: 0.1, linewidth: 2,
  edgecolor: 'red')

plt.show()