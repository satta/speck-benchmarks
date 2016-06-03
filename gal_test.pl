#!/usr/bin/env perl

use GAL::Annotation;
use Data::Dumper;

my $annot = GAL::Annotation->new($ARGV[0], $ARGV[1]);
my $features = $annot->features;

sub has_child_of_type {
  my $root = shift(@_);
  my $type = shift(@_);
  my @children;
  $root->get_recursive_children(\@children);
  foreach (@children) {
    if ($_->type == $type) {
      return 1;
    }
  }
  return 0;
}

my %fails = ();

sub gene_has_transcript {
  my $gene = shift(@_);
  if ($gene->transcripts->count == 0) {
    $fails{"gene"}{"has_transcript"}++;
  }
}

sub gene_children_inside_coords {
  my $gene = shift(@_);
  my @children;
  $gene->get_recursive_children(\@children);
  foreach (@children) {
    if ($_->start < $gene->start || $_->end > $gene->end) {
      $fails{"gene"}{"children_inside_coords"}++;
    }
  }
}

sub gene_children_consistent_strands {
  my $gene = shift(@_);
  my @children;
  $gene->get_recursive_children(\@children);
  foreach (@children) {
    if ($_->strand != $gene->strand) {
      $fails{"gene"}{"children_consistent_strands"}++;
    }
  }
}

sub gene_not_suspiciously_short {
  my $gene = shift(@_);
  if ($gene->length <= 30) {
    $fails{"gene"}{"not_suspiciously_short"}++;
  }
}

sub mrna_n_content {
  my $mrna = shift(@_);
  if ((($mrna->seq =~ tr/Nn//)/($mrna->genomic_length)) >= 0.5 ) {
    $fails{"mRNA"}{"n_content"}++;
  }
}

sub mrna_has_CDS {
  my $mrna = shift(@_);
  if ($mrna->CDSs->count == 0) {
    $fails{"mRNA"}{"has_CDS"}++;
  }
}

sub mrna_has_only_CDS_children {
  my $mrna = shift(@_);
  my @children;
  $mrna->get_recursive_children(\@children);
  if ($mrna->CDSs->count < @children) {
    $fails{"mRNA"}{"has_only_CDS_children"}++;
  }
}

sub mrna_minimum_length {
  my $mrna = shift(@_);
  if ($mrna->CDS_length < 3 ) {
    $fails{"mRNA"}{"minimum_length"}++;
  }
}

sub mrna_internal_stop_codon {
  my $mrna = shift(@_);
  if ($mrna->protein_seq =~ m/[*#+]/) {
    $fails{"mRNA"}{"internal_stop_codon"}++;
  }
}

sub cds_has_no_children {
  my $cds = shift(@_);
  my @children;
  $cds->get_recursive_children(\@children);
  if (@children > 0) {
    $fails{"CDS"}{"has_no_children"}++;
  }
}

my %checks = ();
$checks{"gene"} = [\&gene_has_transcript,
                   \&gene_children_inside_coords,
                   \&gene_children_consistent_strands,
                   \&gene_not_suspiciously_short];
$checks{"mRNA"} = [\&mrna_n_content,
                   \&mrna_has_CDS,
                   \&mrna_has_only_CDS_children,
                   \&mrna_minimum_length,
                   \&mrna_internal_stop_codon];
$checks{"CDS"} = [\&cds_has_no_children];

# stream features and evaluate checks
while (my $feat = $features->next) {
  if (@{$checks{$feat->type}}) {
    foreach (@{$checks{$feat->type}}) {
       &{$_}($feat);
    }
  }
}

print Dumper(\%fails);

