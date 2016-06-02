derives_from = {}

-- template checks
check_parent = function (n, parent_type)
  it("appears as part of a " .. parent_type, function()
    expect(n:appears_as_child_of_type(parent_type)).should_be(true)
  end)
end
is_a_lone_feature = function (n)
  it("appears as a root node", function()
    expect(n:appears_as_root_node()).should_be(true)
  end)

  it("should not have children", function()
    expect(count(n:direct_children())).should_be(0)
  end)
end

describe.feature("gene", function(gene)
  it("contains a transcript", function()
    expect(gene:has_child_of_supertype("transcript")).should_be(true)
  end)

  it("appears as a root node", function()
    expect(gene:appears_as_root_node()).should_be(true)
  end)

  it("contains all child features within its coordinates", function()
    for child in gene:children() do
      expect(gene:get_range():contains(child:get_range())).should_be(true)
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

describe.feature("pseudogene", function(pseudogene)
  it("contains a pseudogenic_transcript", function()
    expect(pseudogene:has_child_of_type("pseudogenic_transcript")).should_be(true)
    expect(pseudogene:has_child_of_type("mRNA")).should_be(false)
  end)

  it("appears as a root node", function()
    expect(pseudogene:appears_as_root_node()).should_be(true)
  end)

  it("contains all child features within its coordinates", function()
    for child in pseudogene:children() do
      expect(pseudogene:get_range():overlap(child:get_range())).should_be(true)
    end
  end)
end)

describe.feature("pseudogenic_transcript", function(ptranscript)
  check_parent(ptranscript, "pseudogene")

  it("contains at least one pseudogenic_exon", function()
    expect(ptranscript:has_child_of_type("pseudogenic_exon")).should_be(true)
  end)
end)

describe.feature("pseudogenic_exon", function(pexon)
  check_parent(pexon, "pseudogenic_transcript")

  it("should not have children", function()
    expect(count(pexon:direct_children())).should_be(0)
  end)
end)

describe.feature("mRNA", function(mrna)
  local dnaseq = mrna:extract_sequence("CDS", true, region_mapping):lower()
  local protseq = mrna:extract_and_translate_sequence("CDS", true,
                                                      region_mapping)

  check_parent(mrna, "gene")

  it("consists of less than 50% Ns", function()
    expect(dnaseq:char_count("n")/dnaseq:len()).should_be_smaller_than(0.5)
  end)

  it("has at least one CDS child", function()
    expect(mrna:has_child_of_type("CDS")).should_be(true)
  end)

  it("has only CDS, UTR or exon children", function()
    local num_feats = count(mrna:children_of_type("CDS"))
                        + count(mrna:children_of_type("exon"))
                        + count(mrna:children_of_supertype("UTR"))
    expect(count(mrna:children())-1).should_be(num_feats)
  end)

  it("has a coding sequence >= 3bp", function()
    expect(dnaseq:len()).should_be_larger_than(2)
  end)

  it("has CDS with no internal stop codons", function()
    expect(protseq:sub(1, -2)).should_not_match("[*+#]")
  end)

  it("has non-partial CDS ending on a stop codon", function()
    if not (mrna:get_attribute("End_range")) then
      expect(protseq:sub(-1)).should_match("[*+#]")
    end
  end)

  it("agrees exactly with CDS/UTR coordinates of its children", function()
    local rng = nil
    -- collect and join CDS ranges
    for c in mrna:children() do
      if c:get_type() == "CDS" or string.match(c:get_type(), "UTR") then
        if not rng then
          rng = c:get_range()
        else
          rng = rng:join(c:get_range())
        end
      end
    end
    -- should overlap with at least one feature
    expect(rng).should_be_truthy()
    -- check if coordinates match
    if rng then
      expect(rng:get_start() == mrna:get_range():get_start() and
             rng:get_end() == mrna:get_range():get_end()).should_be_truthy()
    end
  end)
end)

describe.feature("CDS", function(cds)
 it("appears as child of an mRNA/ncRNA", function()
    expect(cds:appears_as_child_of_supertype("transcript")).should_be(true)
  end)

 it("should not have children", function()
    expect(#(collect(cds:direct_children()))).should_be(0)
  end)
end)

describe.feature("polypeptide", function(pp)
  it("should derive from a unique mRNA", function()
    local dfrom = pp:get_attribute("Derives_from")
    expect(dfrom).should_not_be(nil)
    expect(derives_from).should_not_have_key(dfrom)
    derives_from[dfrom] = true
  end)

  it("appears as a root node", function()
    expect(pp:appears_as_root_node()).should_be(true)
  end)

  it("has a product attribute", function()
    expect(pp:get_attribute("product")).should_not_be(nil)
  end)

  local go_attrib = pp:get_attribute("full_GO")

  it("has correctly formatted GO ID", function()
    if go_attrib then
      for _,t in ipairs(gff3_extract_structure(go_attrib)) do
        expect(t.GOid).should_match("GO:%d+")
      end
    end
  end)

  it("does not have stray whitespace in with/from and Dbxref", function()
    if go_attrib then
      for _,t in ipairs(gff3_extract_structure(go_attrib)) do
        if t.with then
          expect(t.with).should_not_match("%s")
        end
        if t.dbxref then
          expect(t.dbxref).should_not_match("%s")
        end
      end
    end
  end)

  it("must have GO:005515 appearing with IPI", function()
    if go_attrib then
      for _,t in ipairs(gff3_extract_structure(go_attrib)) do
        if t.GOid == "GO:005515" then
          expect(t.evidence).should_match("Physical Interaction")
        end
      end
    end
  end)


  it("must not have IEP evidence for C and F aspects", function()
    if go_attrib then
      for _,t in ipairs(gff3_extract_structure(go_attrib)) do
        if t.evidence and t.evidence:match("Expression Pattern") then
          expect(t.aspect).should_not_match("[CF]")
        end
      end
    end
  end)

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

  it("must not have with/from for IDA/NAS/ND/TAS/EXP evidence", function()
    if go_attrib then
      for _,t in ipairs(gff3_extract_structure(go_attrib)) do
        if t.evidence:match("Direct Assay")
           or t.evidence:match("Statement")
           or t.evidence:match("Experiment")
           or t.evidence:match("No Biological") then
          expect(t.with == nil or t.with:len() == 0).should_be(true)
        end
      end
    end
  end)

  local overlapping = feature_index:get_features_for_range(pp:get_seqid(),
                                                             pp:get_range())

  it("overlaps at least one transcript", function()
    local num_transcripts = 0
    expect(#overlapping).should_be_larger_than(0)
    if #overlapping > 0 then
      for _,ovl_feat in ipairs(overlapping) do
        if ovl_feat:has_child_of_type("mRNA")
            or ovl_feat:has_child_of_type("pseudogenic_transcript") then
          num_transcripts = num_transcripts + 1
        end
      end
      expect(num_transcripts).should_be_larger_than(0)
    end
  end)

  it("agrees exactly with CDS of at >=1 overlapping transcript", function()
    local nof_possible = 0
    local nof_correct = 0
    -- check every feature in the range
    for _,ovl_feat in ipairs(overlapping) do
      for n in ovl_feat:children_of_type("mRNA") do
        -- locate transcript (no pseudogene etc)
        nof_possible =  nof_possible + 1
        local rng = nil
        -- collect and join CDS ranges
        for c in n:children_of_type("CDS") do
          if not rng then
            rng = c:get_range()
          else
            rng = rng:join(c:get_range())
          end
        end
        -- should overlap with at least one feature
        expect(rng).should_not_be(nil)
        -- check if coordinates match
        if rng then
          if rng:get_start() == pp:get_range():get_start() and
             rng:get_end() == pp:get_range():get_end() then
            nof_correct = nof_correct + 1
          end
        end
      end
    end
    if nof_possible > 0 then
      expect(nof_correct).should_be_larger_than(0)
    end
  end)
end)

describe.feature("ncRNA", function(node)
  check_parent(node, "gene")
end)

describe.feature("tRNA", function(node)
  check_parent(node, "gene")
end)

describe.feature("rRNA", function(node)
  check_parent(node, "gene")
end)

describe.feature("snRNA", function(node)
  check_parent(node, "gene")
end)

describe.feature("snoRNA", function(node)
  check_parent(node, "gene")
end)

describe.feature("gap", function(gap)
  is_a_lone_feature(gap)
end)

describe.feature("contig", function(contig)
  is_a_lone_feature(contig)
end)
