#!/usr/bin/env ruby
#
# Authors: (primary) Rob Smith and (secondary) John T. Prince
# Released under GNU GPLv3
# Contact Rob Smith at robert.smith@mso.umt.edu for other licensing options.

require 'optparse'
require 'set'
require 'yaml'

class MSAbundanceSim
  VERSION = "0.2.0"
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
    pct_diff_express: 3.0,
    control_variance: 1,
    case_variance: 2,
    output_abundance_separator: " #",
    downshift_min_threshold: 0,
    downshift_probability_threshold: 0.75,
    downshift_amount: 1,
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

    # takes 'value' and subtracts 'amount' from it if 'value' is greater than
    # min_threshold AND a random number is less than the probability_threshold
    # (i.e., a probability of 1 means downshifting will happen to all values).
    def downshift(value, min_threshold, probability_threshold, amount, random=rand())
      if (value > min_threshold) && (random < probability_threshold)
        value - amount
      else
        value
      end
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
    opts = DEFAULTS.merge(opts)

    basename = filename.chomp(File.extname(filename))

    protein_entries = get_protein_entries(filename)
    max_abundance = protein_entries.max_by {|entry| entry.abundances.last }.abundances.last

    # generate which proteins will be differentially expressed
    num_proteins_to_diff_express = (protein_entries.size * opts[:pct_diff_express]/100.0).to_i
    diff_expressed_ids = (0...protein_entries.size).to_a.sample(num_proteins_to_diff_express).to_set

    output = Hash.new {|hash, key| hash[key] = [] }
    case_and_controls_needed = ([:case] * opts[:num_case]) + ([:control] * opts[:num_control])
    case_and_controls_needed.each_with_index do |sample_type, sample_number| # for each sample
      outfilename = "#{basename}_#{sample_number}_#{sample_type}"
      output[sample_type] << outfilename
      File.open(outfilename, 'w') do |outfile|
        protein_entries.each_with_index do |protein_entry, idx|
          # put first line of fasta with simulated abundance
          variance = opts[
            if sample_type == :case && diff_expressed_ids.include?(idx)
              :case_variance
            else
              :control_variance
            end
          ]

          fold_change = MSAbundanceSim.get_fold_change(protein_entry.abundances, variance, max_abundance)

          downshift_params = %w(min_threshold probability_threshold amount)
            .map {|key| 'downshift_' + key }
            .map(&:to_sym)
            .map {|key| opts[key] }

          downshifted_fold_change = MSAbundanceSim.downshift(fold_change, *downshift_params)

          sample_abundance = MSAbundanceSim.sample_abundance(protein_entry.abundances, downshifted_fold_change)
          outfile.puts [protein_entry.entry_line_wo_abundance, sample_abundance].join(opts[:output_abundance_separator])
          outfile.puts protein_entry.additional_lines.join("\n")
        end
      end
    end
    output
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
    protein_entries
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
          output = MSAbundanceSim.process_files(filenames, opts)
          format_and_emit_output(output)
        end
      end

      def format_and_emit_output(output)
        # the type keys are symbols--make them strings
        output_string_types = output.map {|key, value| [key, value.map {|k, v| [k.to_s, v] }.to_h] }.to_h
        # emit yaml without leading "---"
        puts output_string_types.to_yaml.gsub(/\A---\n/, '')
      end

      def create_parser
        defaults = MSAbundanceSim::DEFAULTS
        opts = MSAbundanceSim::DEFAULTS
        parser = OptionParser.new do |op|
          op.banner = "usage: #{File.basename(__FILE__)} [options] <file>.fasta ..."
          op.separator "generates files, each of this form: <file>_<n>_<case|control>"
          op.separator "  where the specified percentage of proteins have altered abundance (for cases)."
          op.separator ""
          op.separator "output: (yaml indicating all the related derivative files)"
          op.separator "  for 2 cases and 1 control the output for file.fasta would be: "
          op.separator "    file.fasta:"
          op.separator "      case:"
          op.separator "      - file_0_case"
          op.separator "      - file_1_case"
          op.separator "      control:"
          op.separator "      - file_2_control"
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
          op.separator ""
          op.separator "options (<default>):"

          op.on(
            "--num-control <#{defaults[:num_control]}>",
            Integer,
            "Number of *control* samples to generate."
          ) {|v| opts[:num_control] = v }

          op.on(
            "--num-case <#{defaults[:num_case]}>",
            Integer,
            "Number of *case* samples to generate."
          ) {|v| opts[:num_case] = v }

          op.on(
            "--pct-diff-express <#{defaults[:pct_diff_express]}>",
            Float,
            "Percentage of proteins to differentially ",
            "  express between case and control."
          ) {|v| opts[:pct_diff_express] = v }

          op.on(
            "--control-variance <#{defaults[:control_variance]}>",
            Integer,
            "Variance for control samples (max lambda for",
            "  Poisson distribution). The higher the value, ",
            "  the more fold change will occur among healthy ",
            "  populations. Used only when multiple ",
            "  abundances are not provided in master fasta. ",
            "  Not recommended to modify this parameter."
          ) {|v| opts[:control_variance] = v }

          op.on(
            "--case-variance <#{defaults[:case_variance]}>",
            Integer,
            "Variance increase for case samples (max lambda ",
            "  for Poisson distribution). The higher the ",
            "  value the more fold change will occur."
          ) {|v| opts[:case_variance] = v }

          op.on(
            "--downshift-min-threshold<#{defaults[:downshift_min_threshold]}>",
            Float,
            "Min threshold for downshifting some fold changes."
          ) {|v| opts[:downshift_min_threshold] = v }

          op.on(
            "--downshift-probability-threshold<#{defaults[:downshift_probability_threshold]}>",
            Float,
            "Min probability threshold for downshifting",
            "some fold changes."
          ) {|v| opts[:downshift_probability_threshold] = v }

          op.on(
            "--downshift-amount<#{defaults[:downshift_amount]}>",
            Float,
            "Amount the fold change will be altered",
            "if meeting the threshold and probability-",
            "threshold."
          ) {|v| opts[:downshift_amount] = v }

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


