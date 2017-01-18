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
    # returns a Hash keyed by filenames and pointing to the output of
    # process_file (a list of case and control filenames, keyed by 'case' and
    # 'control')
    def process_files(filenames, opts)
      outputs = filenames.map do |filename|
        MSAbundanceSim.new.process_file(filename, opts)
      end
      filenames.zip(outputs).to_h
    end

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

  def process_file(filename, opts)
    diff_express_percent, num_case, num_control, control_variance, case_variance = opts.values_at(:diff_express_percent, :num_case, :num_control, :control_variance, :case_variance)

    simulator = MSAbundanceSim.new

    basename = filename.chomp(File.extname(filename))

    protein_entries = simulator.get_protein_entries(filename)
    max_abundance = protein_entries.max_by {|entry| entry.abundances.last }.abundances.last

    # generate which proteins will be differentially expressed
    diff_expressed_ids = [0..protein_entries.size-1].sample((protein_entries.size * diff_express_percent/100.0).to_i)
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
        protein_entries.each_with_index do |protein_entry, idx|
          # put first line of fasta with simulated abundance
          sign = [1,-1].sample
          if type=='case' and diff_expressed_ids.index(idx) != nil
            sign = diff_expressed_signs[idx]
          else
            type = 'control'
          end

          outfile.puts "#{protein_entry.entry_line_wo_abundance} + ##{MSAbundanceSim.sample_abundance(protein_entry.abundances, MSAbundanceSim.get_fold_change(protein_entry.abundances, type=='control' ? control_variance : case_variance, max_abundance))}"
          protein_entry.additional_lines.each do |additional_line|
            outfile.puts additional_line
          end
        end
      end
    end
  end

  # returns an array with the beginning part of the entry line (without the
  # abundance (which begins with a '#') and an array of abundances (all the
  # numbers, returned as Floats, after the final '#')
  def parse_entry_line(line)
    octothorpe_index = line.rindex("#")
    entry_line_wo_abundance = line[0...octothorpe_index]
    abundance_str = line[(octothorpe_index+1)..-1]
    abundances = abundance_str.split(",").map(&:to_f).sort
    [entry_line_wo_abundance, abundances]
  end

  def get_protein_entries(filename)
    protein_entries = []
    IO.foreach(filename) do |line|
      line.chomp!
      if line[0] == ">"
        protein_entries << ProteinEntry.new(*parse_entry_line(line))
      else
        protein_entries.last.additional_lines << line
      end
    end
    protein_entries[0...-1]  # incorrect behavior (dropping last protein)
  end

  class Commandline
    class << self
      def run(argv)
        parser, opts = create_parser
        parser.parse!(argv)

        filenames = argv.to_a

        if filenames.size == 0
          puts parser
        else
          MSAbundanceSim.process_files(filenames, opts)
        end
      end

      def create_parser
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
          op.on("--verbose", "talk about it") {|v| $VERBOSE = 3 }
        end
        [parser, opts]
      end
    end
  end
end

if __FILE__ == $0
  MSAbundanceSim::Commandline.run(ARGV)
end


