describe.feature("gene", function(gene)
  it("contains a transcript", function()
    expect(gene:has_child_of_supertype("transcript")).should_be(true)
  end)

  it("contains all child features within its coordinates", function()
    for child in gene:children() do
      expect(gene:get_range():contains(child:get_range())).should_be(true)
    end
  end)

  it("adheres to a ID naming scheme", function()
    expect(gene:get_attribute("ID")).should_match("^PF3D7_%d%d%d%d")
  end)

  it("has consistent strands across all children", function()
    for child in gene:children() do
      expect(gene:get_strand()).should_be(child:get_strand())
    end
  end)
end)

describe.feature("mRNA", function(mrna)
  local dna = mrna:extract_sequence("CDS", true, region_mapping):lower()
  local prot = mrna:extract_and_translate_sequence("CDS", true, region_mapping)

  it("consists of less than 50% Ns", function()
    expect(dna:char_count("n")/dna:len()).should_be_smaller_than(0.5)
  end)

  it("has CDS with no internal stop codons", function()
    expect(prot:sub(1, -2)).should_not_match("[*+#]")
  end)
end)

derives_from = {}
describe.feature("polypeptide", function(pp)
  it("should derive from a unique mRNA", function()
    local dfrom = pp:get_attribute("Derives_from")
    expect(dfrom).should_not_be(nil)
    expect(derives_from).should_not_have_key(dfrom)
    derives_from[dfrom] = true
  end)

  local go_attrib = pp:get_attribute("full_GO")
  it("must have with/from for ISS/ISA/ISO/ISM evidence", function()
    if go_attrib then
      for _,t in ipairs(gff3_extract_structure(go_attrib)) do
        if t.evidence and t.evidence:match("Sequence") then
          expect(t.with).should_not_be(nil)
          if t.with then
            expect(t.with:len()).should_be_larger_than(0)
          end
        end
      end
    end
  end)
end)
