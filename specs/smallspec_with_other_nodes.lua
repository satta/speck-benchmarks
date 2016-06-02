describe.region(function(region)
  it("starts at coordinate 1", function()
    expect(region:get_range():get_start()).should_be(1)
  end)
end)

describe.comment(function(comment)
  it("is not longer than 80 characters", function()
    expect(string.len(comment:get_comment())).should_be_smaller_than(80)
  end)
end)

describe.feature("gene", function(gene)
  it("contains a transcript", function()
    expect(gene:has_child_of_supertype("transcript")).should_be(true)
  end)

  it("contains all child features within its coordinates", function()
    for child in gene:children() do
      expect(gene:get_range():overlap(child:get_range())).should_be(true)
    end
  end)

  it("has consistent strands across all children", function()
    for child in gene:children() do
      expect(gene:get_strand()).should_be(child:get_strand())
    end
  end)

  it("is not suspiciously short (>30nt)", function()
    expect(gene:get_range():length()).should_be_larger_than(30)
  end)
end)

describe.feature("mRNA", function(mrna)
  local dnaseq = mrna:extract_sequence("CDS", true, region_mapping):lower()
  local protseq = mrna:extract_and_translate_sequence("CDS", true,
                                                      region_mapping)

  it("consists of less than 50% Ns", function()
    expect(dnaseq:char_count("n")/dnaseq:len()).should_be_smaller_than(0.5)
  end)

  it("has at least one CDS child", function()
    expect(mrna:has_child_of_type("CDS")).should_be(true)
  end)

  it("has only CDS children", function()
    expect(count(mrna:children())-1).should_be(count(mrna:children_of_type("CDS")))
  end)

  it("has a coding sequence >= 3bp", function()
    expect(dnaseq:len()).should_be_larger_than(2)
  end)

  it("has CDS with no internal stop codons", function()
    expect(protseq:sub(1, -2)).should_not_match("[*+#]")
  end)
end)

describe.feature("CDS", function(cds)
  it("should not have children", function()
    expect(#(collect(cds:direct_children()))).should_be(0)
  end)
end)
