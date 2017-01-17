#!/usr/bin/env ruby
#
# Authors: (primary) Rob Smith and (secondary) John T. Prince
# Released under GNU GPLv3
# Contact Rob Smith at robert.smith@mso.umt.edu for other licensing options.

require 'optparse'

# REMOVE THIS
srand(47288)

class MSAbundanceSim
  VERSION = "0.1.0"
end

# puts that respects $VERBOSE
def putsv(*args)
  puts(*args) if $VERBOSE
end

ProteinEntry = Struct.new(:entry_line_wo_abundance, :abundances, :additional_lines) do
  def initialize(*args)
    super(*args)
    self.additional_lines ||= []
  end
end

class Integer
  def factorial
    # https://www.rosettacode.org/wiki/Factorial#Ruby (fastest and most concise)
    (2..self).reduce(1, :*)
  end
end

class MSAbundanceSim
  DEFAULTS = {
    num_control: 5,
    num_case: 5,
    diff_express_percent: 3.0,
    control_variance: 1,
    case_variance: 2,
  }

  class << self
    # event_rate (lambda)
    # num_occurences (k)
    def poisson(event_rate, num_occurences)
      return event_rate**num_occurences * Math.exp(-event_rate) / num_occurences.factorial.to_f
    end

    # probability (y)
    def inverse_transform_sample(event_rate)
      max = poisson(event_rate, event_rate)
      probability = rand * (max)
      list_of_num_occurrences = (0..event_rate*4).to_a.shuffle # attempt to capture the entire distribution without wasting tests
      list_of_num_occurrences.each do |num_occurences|
        p_num_occurences = poisson(event_rate, num_occurences)
        return num_occurences if probability < p_num_occurences
      end
      puts "Problem: Should have found a num_occurences in inverse_transform_sample"
    end

    def get_fold_change(a_i, event_rate, max_abundance)
      # max fold change at lowest abundance
      raw_fc = inverse_transform_sample(event_rate)
      norm_fc = raw_fc * (1.0-(a_i[0].to_f / max_abundance))
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
  end

  def get_protein_entries(filename)
    entries = []
    abundances = []
    proteins = [] # [0] is list of fasta lines, [1] is abundances list
    IO.foreach(filename) do |line|
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

        line = parts[0]
      else
        entries << line
      end
    end
    proteins
  end
end

if __FILE__ == $0
  defaults = MSAbundanceSim::DEFAULTS
  opts = MSAbundanceSim::DEFAULTS
  parser = OptionParser.new do |op|
    op.banner = "usage: #{File.basename(__FILE__)} <file>.fasta ..."
    op.separator "output: <file>_<n>_<case|control>"
    op.separator ""
    op.separator "The file must have one or more abundances per protein entry (the protein"
    op.separator "sequence following the header line is optional)."
    op.separator "The abundance is placed at the end of the line following a ' #',"
    op.separator "with multiple abundances separated with a ','.  Here are two examples:"
    op.separator ""
    op.separator "> SWISSAB|23B The anchor protein #23.2"
    op.separator "> SWISSSPECIAL|24B A green protein #23.2,29.4"
    op.separator ""
    op.separator "notes: Protein sequences are optional and files need not end in '.fasta'"

    op.on(
      "--num-control <#{defaults[:num_control]}>",
      Integer,
      "how many control samples to generate"
    ) {|v| opts[:num_control] = v }

    op.on(
      "--num-case <#{defaults[:num_case]}>",
      Integer,
      "how many case samples to generate"
    ) {|v| opts[:num_case] = v }

    op.on(
      "--diff_express_percent <#{defaults[:diff_express_percent]}>",
      Float,
      "percent of proteins to differentially express between case and control"
    ) {|v| opts[:diff_express_percent] = v }

    op.on(
      "--control-variance <#{defaults[:control_variance]}>",
      Integer,
      "Variance for control samples (max lambda for Poisson distribution). ",
      "The higher the value, the more fold change will occur among healthy populations. ",
      "Used only when multiple abundances are not provided in master fasta. ",
      "Not recommended to modify this parameter."
    ) {|v| opts[:control_variance] = v }

    op.on(
      "--case-variance <#{defaults[:case_variance]}>",
      Integer,
      "Variance increase for case samples (max lambda for Poisson distribution). ",
      "The higher the value, the more fold change will occur."
    ) {|v| opts[:case_variance] = v }
  end


  parser.parse!
  opts[:filenames] = ARGV.to_a

  if opts[:filenames].size == 0
    puts parser
    exit
  end

  #########################################
  diff_express_percent, num_case, num_control, control_variance, case_variance = opts.values_at(:diff_express_percent, :num_case, :num_control, :control_variance, :case_variance)

  simulator = MSAbundanceSim.new

  opts[:filenames].each do |filename|

    basename = filename.chomp(File.extname(filename))

    proteins = simulator.get_protein_entries(filename)
    max_abundance = proteins.max_by {|protein| protein.last.last }.last.last

    # generate which proteins will be differentially expressed
    diff_expressed_ids = [0..proteins.size-1].sample((proteins.size * diff_express_percent/100.0).to_i)
    diff_expressed_signs = Array.new(diff_expressed_ids.size){[-1,1].sample}

    sample_n = num_case + num_control
    (0..sample_n).each do |sample_number| # for each sample
      type = "control"
      if sample_number < num_case # make a case sample
        type = "case"
      end
      puts "Creating sample #{sample_number} of #{sample_n}"

      # create output file
      outfilename = "#{basename}_#{sample_number}_#{type}"
      File.open(outfilename, 'w') do |outfile|
        proteins.each_with_index do |protein, idx|
          # put first line of fasta with simulated abundance
          sign = [1,-1].sample
          if type=='case' and diff_expressed_ids.index(idx) != nil
            sign = diff_expressed_signs[idx]
          else
            type = 'control'
          end

          outfile.puts "#{protein[0][0]} + ##{MSAbundanceSim.sample_abundance(protein[1], MSAbundanceSim.get_fold_change(protein[1], type=='control' ? control_variance : case_variance, max_abundance))}"
          protein[0][1..-1].each do |additional_line|
            outfile.puts additional_line
          end
        end
      end
    end
  end
end
