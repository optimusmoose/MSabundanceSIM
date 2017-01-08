require 'test_helper'

# minitest will respect the SEED env variable and always give consistent
# results even though random values are being used.

describe MSAbundanceSim do
  def setup
    srand(67809)
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
    # working here
  end

end
