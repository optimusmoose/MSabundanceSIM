# Author: Rob Smith
# Released under GNU GPLv3
# Contact Rob Smith at robert.smith@mso.umt.edu for other licensing options.

#Input:
# FASTA with one or more abundances per protein, included on first line of each entry after # and separated with ,
# OR text file with each line containing a label of your choice with abundances in the same format as above
#Output:
# FASTAs with readme (parameters), postfixed with _case or _control
# OR text file with each line formatted similarly

##########################
## PARAMETERS:
#master_fasta = File.open("9606-PA_plasma.fasta","r") # Master fasta file
master_fasta = File.open("test.fasta","r") # Master fasta file

control_n = 5 						# how many control samples to generate

case_n = 5 								# how many case samples to generate

diff_express_percent = 3 	# what percent of proteins to differentially express 
													# between case and control

control_variance = 1	 		# Variance for control samples (max lambda for Poisson distribution). 
													# The higher the value, the more fold change will occur
											 		# among healthy populations. Used only when multiple 
													# abundances are not provided in master fasta. Modification
													# not recommended.

case_variance = 2 				# Variance increase for case samples (max lambda for Poisson distribution). 													# The higher the value, the more fold change will occur. 

##########################

class Integer
  def factorial
    f = 1; for i in 1..self; f *= i; end; f
  end
end

def poisson(lambda, k)
	return lambda**k * Math.exp(-lambda) / k.factorial.to_f
end

def inverse_transform_sample(lambda)
	max = poisson(lambda, lambda)
	y = rand * (max) 
	ks = (0..lambda*4).to_a.shuffle # attempt to capture the entire distribution without wasting tests
	ks.each do |k|
		p_k = poisson(lambda, k)
		return k if y < p_k
	end
	puts "Problem: Should have found a k in inverse_transform_sample"
end

def get_fold_change_clean(a_i, lambda)
	# max fold change at lowest abundance
	raw_fc = inverse_transform_sample(lambda)
	norm_fc = raw_fc 
	sign = [1,-1].sample
	return norm_fc * sign
end
 
def get_fold_change(a_i, lambda)
	# max fold change at lowest abundance
	raw_fc = inverse_transform_sample(lambda)
	norm_fc = raw_fc * (1.0-(a_i[0].to_f / $abundance_max))
	pert_fc = norm_fc * rand
	sign = [1,-1].sample
	return pert_fc * sign
end

def sample_abundance(abundances, fold_change)
	abundance = 0.0
	if abundances.size > 1 # linear interpolation
		idx = [1..abundances.size-1].sample(1)
		next_idx = idx+1 < abundances.size ? idx+1 : abundances[0] # edge case: 2 abundances
		abundance = abundances[idx] + rand * (abundances[next_idx] - abundances[idx])
	else # sample w/variance
		if fold_change >= 0 # positive fold change
			abundance = abundances[0].to_f * fold_change + abundances[0]
		else # negative fold change
			fold_change *= -1
			abundance = abundances[0] / (2.0**fold_change)
		end
	end
	return abundance + [1,-1].sample * rand * 0.1 * abundance
end

entries = []
abundances = []
proteins = [] # [0] is list of fasta lines, [1] is abundances list
$abundance_max = 0
$abundance_min = 9999999999999999999
while line = master_fasta.gets
	line = line.chop
	if line.index(">") != nil #first line of entry
		
		unless abundances.size == 0 # first time
			abundances.sort!

			# process last fasta entry
			proteins << [entries,abundances]
			entries = []
			abundances = []
		end

		# grab intensit[ies] of this entry
		parts = line.split("#")
		abundance = parts[1].to_f
		abundances << abundance
		entries << parts[0]

		$abundance_max = abundance if abundance > $abundance_max
		$abundance_min = abundance if abundance < $abundance_min

		line = parts[0]
	else
		entries << line
	end
end

# generate which proteins will be differentially expressed
diff_expressed_ids = [0..proteins.size-1].sample((proteins.size * diff_express_percent/100.0).to_i)
diff_expressed_signs = Array.new(diff_expressed_ids.size){[-1,1].sample}

sample_n = case_n + control_n
(0..sample_n).each do |n| # for each sample
	type = "control"
	if n < case_n # make a case sample
		type = "case"
	end
	puts "Creating sample #{n} of #{sample_n}"

	# create output file
	outfile = File.open("#{n}_#{type}","w")
	proteins.each_with_index do |protein, idx|
		# put first line of fasta with simulated abundance
		sign = [1,-1].sample
		if type=='case' and diff_expressed_ids.index(idx) != nil
			sign = diff_expressed_signs[idx]
		else
			type = 'control'
		end

		outfile.puts "#{protein[0][0]} + ##{sample_abundance(protein[1], get_fold_change(protein[1], type=='control' ? control_variance : case_variance))}" 
		protein[0][1..-1].each do |additional_line|
			outfile.puts additional_line
		end
	end
end
