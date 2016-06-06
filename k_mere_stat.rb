class KMereStat
	
	attr_accessor :data_path
	attr_accessor :sequences
	attr_accessor :append

	def initialize(path)
		self.data_path = path
		self.sequences = parse_fasta(path)
	end

	def parse_fasta(path)
		fastas = {}
		current_sequence = ""
		current_sequence_name = ""
		File.foreach(path) do |line|
			if line.start_with?(">")
				unless current_sequence.empty?
					fastas[current_sequence_name] = current_sequence
					current_sequence = ""
				end
				current_sequence_name = line.split("|")[1]
			else
				current_sequence << line.strip
			end
		end
		return fastas
	end

	def count_stat(k)
		self.append = append
		stat = Hash.new(0)
		self.sequences.each do |seq_name, seq|
			k_meres = get_all_k_meres(seq, k)
			stat[seq_name] = k_meres
		end
		save_statistics(stat,k)
	end

	def get_uniq_k_meres(stat)
		uniq_k_meres = []
		stat.each do |seq, stats|
			uniq_k_meres = uniq_k_meres | stats.keys
 		end
 		return uniq_k_meres
	end

	def get_all_k_meres(sample, k_mere_length)
		k_meres = Hash.new(0)
		(0..sample.length-k_mere_length).each do |sub_sample|
			k_meres[sample[sub_sample,k_mere_length]] += 1
		end
		return k_meres
	end

	def save_statistics(stat, filename)
		uniq_k_meres = get_uniq_k_meres(stat)
		puts "#{uniq_k_meres.length}"
		header = "k_mer," + self.sequences.keys.map { |x| "0-#{x}" }.join(",")
		header += "," + self.sequences.keys.map { |x| "1-#{x}" }.join(",")
		header += "," + self.sequences.keys.map { |x| "2-#{x}" }.join(",")
		File.open("birds_sorted_#{filename}.csv", "w") do |file|
			file.write("#{header}\n") unless self.append
			uniq_k_meres.sort_by {|k,v| k }.each do |k_mere|
				result = []
				stat.each do |seq, seq_stat|
					result << seq_stat[k_mere]
				end
				k_meres_on_distance_d = get_on_distance_d(k_mere,1)
				stat.each do |seq, seq_stat|
					sum = 0
					k_meres_on_distance_d.each do |k_mere_on_d|
					 	sum += seq_stat[k_mere_on_d]
					end
					result << sum
				end
				k_meres_on_distance_d = get_on_distance_d(k_mere,2)
				stat.each do |seq, seq_stat|
					sum = 0
					k_meres_on_distance_d.each do |k_mere_on_d|
					 	sum += seq_stat[k_mere_on_d]
					end
					result << sum
				end
				file.write("#{k_mere},"+result.join(",") + "\n")
			end
		end
	end

	def get_on_distance_d(sample,d)
		samples = {}
		nucleotides = ["A","T","G","C"]
		places = (0..(sample.length-1)).to_a.permutation(d).to_a
		nucleo_samples = nucleotides.repeated_combination(d).to_a
		places.each do |current_places|
			nucleo_samples.each do |nucleo_sample|
				sample_copy = sample.dup
				current_places.length.times do |place|
					sample_copy[current_places[place]] = nucleo_sample[place]
				end
				samples[sample_copy] = 0
			end
		end
		return samples.keys
	end

end

path = "birds_7.fasta"
k_mere_lengths = [7,5]
k_mere_stat = KMereStat.new(path)

k_mere_lengths.each do |length|
	k_mere_stat.count_stat(length)	
end