require 'test_helper'


describe MSAbundanceSim do
  TESTFILE = TESTFILES + "/test.fasta"

  def setup
    srand(67809)
    @delete_files = []
  end

  def teardown
    @delete_files.each do |file|
      if File.exist?(file)
        File.unlink(file)
      end
    end
  end

  it "has a version number" do
    refute_nil ::MSAbundanceSim::VERSION
  end

  it "calculates a fold change" do
    reply = MSAbundanceSim.get_fold_change([1.3], 3, 10000)
    reply.must_equal 0.885067564263667
  end

  it "samples the inverse transform" do
    reply = MSAbundanceSim.inverse_transform_sample(3)
    reply.must_equal 3
  end

  it "samples a list of abundances" do
    reply = MSAbundanceSim.sample_abundance([0.187], -1.24)
    reply.must_equal 0.08406721711401088
  end

  it "writes cases and controls for a real file" do
    reply = MSAbundanceSim.process_files([TESTFILE])
    reply.keys.size.must_equal 1
    reply.keys.must_equal [TESTFILE]

    file_specific_output = reply.values.first

    # assumes 5 case, 5 control
    # 11 is INCORRECT (should be 10) !!!! The behavior of the original program was off
    roughly_expected = {case: (0...5), control: (5...11)}.map do |key, range|
      [key, range.map {|sample_number| "test_#{sample_number}_#{key}" }]
    end.to_h

    filenames = file_specific_output.values.flatten(1)
    file_specific_output.keys.must_equal roughly_expected.keys
    filenames.zip(roughly_expected.values.flatten(1)) do |actual, roughly_expected|
      assert actual.end_with?(roughly_expected)
    end
    @delete_files = filenames
  end

  it "outputs files with the same number of lines as the original file" do
    num_lines = IO.readlines(TESTFILE).size
    # the behavior is currently incorrect!  The last protein is not represented
    # THIS is an INCORRECT and the line needs to be removed:
    num_lines = num_lines - 4
    reply = MSAbundanceSim.process_file(TESTFILE)

    filenames = reply.values.flatten(1)
    filenames.each do |outfile|
      IO.readlines(outfile).size.must_equal num_lines
    end
    @delete_files = filenames
  end
end
