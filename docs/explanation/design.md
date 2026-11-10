# Scientific design

Each sample is treated as a metagenome. The input contract accepts basecalled ONT
DNA or PacBio HiFi reads. RNA reads are linked evidence for gene calling and do
not enter DNA assembly.

Myloasm produces a sample assembly. Read mapping measures contig depth. QuickBin
uses coverage and composition to form bins; unpaired long-read alignments do not
supply paired-read edges. The workflow retains unbinned contigs. QuickClade uses
an explicit local reference, with remote-server fallback disabled.

Tiara classifies every assembly contig of 1,000 bp or longer as eukaryote,
bacterium, archaeon, prokaryote, organelle, or unknown. It runs on the whole
assembly in every full-mode run and keeps its native table in
`15_contig_classification/tiara.tsv`. Routing sums these calls over each bin's
base pairs, and a contig below Tiara's minimum length counts as unclassified.

CheckEUK screens every bin independently of QuickClade. Bacterial and archaeal
quality estimates, eukaryotic marker recovery, viral-contig estimates, and
giant-virus bin evidence retain their own methods and units. Missing estimates
are null. Conflicting evidence remains visible for review before gene calling.
CheckV contig assessments do not establish the purity of a multi-contig viral bin.

CheckM2 can finish gene calling and DIAMOND with no reference annotations.
The workflow accepts this exit only when the native log, empty DIAMOND files
and predicted proteins for every input bin confirm the condition. It records
`skipped/no_diamond_annotations` and retains the diagnostic files, command and
exit code. CheckM2 quality estimates remain unavailable. Other CheckM2 failures
stop the workflow. A successful prediction can omit a bin with no called proteins
only when its native empty-file warning and regular zero-byte protein file agree
with the missing report entry. Protein files must cover every input bin, and
reported and proven omitted bins together must cover those same bins exactly.
The completed stage records `<bin_id>: no_called_proteins` in its notes and
retains the native logs, exit code, and empty protein files. The bin remains in
the catalog without a CheckM2 QC row; reported estimates remain unchanged.

GVClass uses the database path from the registry. Its native run receipt must
match that resolved path and the declared version before results are accepted.

QuickClade labels remain comparison evidence and cannot establish or veto a
route. Routing compares independent lines of evidence for each bin. A share is a
fraction of the bin's base pairs.

- Eukaryote lines: CheckEUK status `ok` with a resolved lineage; GVClass `d_EUK`;
  a physically linked eukaryotic 18S locus; a Tiara eukaryote plus organelle share of
  at least 0.5, because organelle sequence is eukaryote-derived.
- Prokaryote lines: GVClass `d_BAC` or `d_ARC`; a Tiara prokaryote share (bacteria,
  archaea or prokarya) of at least 0.8; CheckM2 completeness of at least 10 percent.
  Tiara prokaryote shares of bins with eukaryote labels are bimodal, below 0.2 for
  real eukaryotes and at least 0.8 for bacteria, so a partial share does not count.
- Viral support: GVClass `d_NCLDV`, `d_PPV`, `d_MIRUS` or `d_PHAGE`, or whole-contig
  geNomad calls that cover at least 0.5 of the bin.

An 18S locus is linked when its model is RF01960 and it sits on a bin contig whose
`taxonomy_domain` is `Eukaryota`, or whose `reference_taxonomy` or
`blast_candidate_taxonomy` starts with `Eukaryota`. SSUextract 1.2.2 writes
`Unclassified` in `taxonomy_domain` for every 18S locus, so the lineage columns carry
the eukaryotic evidence. CheckEUK completeness alone does not establish taxonomy.
CheckEUK reports a lineage only when its host vote has at least five informative
markers, so small bins often have none. CheckM2 evidence is absent for every bin when
CheckM2 was skipped, and for a bin that has no row in its table.

The first rule that matches decides the route:

1. GVClass `d_EUK-pEVE` leaves the bin unresolved with class `conflicting`.
2. GVClass `d_PLASTID` or `d_MITO` leaves the bin unresolved, because the workflow has
   no organelle gene caller. A Tiara organelle share of at least 0.5 does the same
   when the bin is at most 2 Mb and GVClass is not prokaryotic. Larger bins are
   excluded because short AT-rich ciliate nanochromosomes look like organelle DNA;
   the organelle bins seen in real samples are below 0.6 Mb. The class is
   `conflicting` when eukaryote support or viral support exists.
3. The Tiara prokaryote line and the CheckM2 completeness line both hold, and no 18S
   locus is linked: the bin takes the prokaryotic route, so GTDB-Tk also sees partial
   bacterial genomes. Viral support makes it `conflicting` instead. Any CheckEUK
   lineage or GVClass `d_EUK` is listed in `conflicting_evidence`. The reason is the
   GVClass prokaryotic reason when GVClass is prokaryotic, and otherwise "Tiara and
   CheckM2 support a prokaryotic route". This rule routes single-contig bacterial
   genomes that CheckEUK and GVClass label eukaryotic.
4. Eukaryote support is two or more eukaryote lines. The bin takes the eukaryotic
   route unless viral support exists, GVClass is prokaryotic, or the Tiara prokaryote
   share is at least 0.8; those cases are `conflicting`. The reason names the first two
   lines, for example "CheckEUK and GVClass support a eukaryotic route".
   When CheckEUK status is `ok`, its lineage is blank, its host vote has fewer than
   five markers (`marker_tax_host_lineage_selected_markers` in
   `checkeuk_report_detail.tsv`, a blank value counting as zero) and GVClass assigns
   `d_EUK`, the sparse-CheckEUK fallback also gives eukaryote support, with the reason
   "GVClass supports a eukaryotic route; CheckEUK has too few markers to name a
   lineage".
5. Viral support makes the route viral, unless GVClass is prokaryotic, rule 3 holds, or
   a CheckEUK lineage, GVClass `d_EUK` or a linked 18S locus holds. Those cases are
   `conflicting`.
   The Tiara eukaryote line never conflicts with viral support, because Tiara has no
   virus class and labels giant-virus contigs eukaryotic. It still counts toward the
   two lines of rule 4.
6. GVClass `d_BAC` or `d_ARC` gives the prokaryotic route.
7. Otherwise the bin is unresolved. The reason says whether a single eukaryote line
   did not meet the combined support rule.

`conflicting` means unresolved with candidate class `conflicting`. Whole-contig geNomad
calls that cover less than half of a bin never veto a route. The ledger lists them as
`geNomad-minority:<share>` in `conflicting_evidence`, so one viral contig no longer
removes a bin that CheckEUK and GVClass support as eukaryotic.

The routing ledger `evidence.tsv` reports one taxonomy string per bin in
`taxonomy`, with its origin in `taxonomy_source`. A bin with a CheckEUK lineage
carries that lineage and the source `checkeuk`. A bin routed through the GVClass
fallback carries the GVClass `taxonomy_majority` and the source `gvclass_fallback`.
All other bins leave both fields blank.

The ledger appends nine columns after `conflicting_evidence`: `bin_bp`, the Tiara
shares `tiara_eukaryote_fraction`, `tiara_prokaryote_fraction`,
`tiara_organelle_fraction` and `tiara_unclassified_fraction`, `viral_contig_fraction`,
`checkm2_completeness`, `checkm2_contamination`, and `duplicated_contig_fraction`.
Shares have four decimals. A blank CheckM2 value means no CheckM2 evidence for the
bin. `duplicated_contig_fraction` is the share of bases in contigs that Myloasm marks
`_duplicated-yes`; it informs review and does not change a route. Evidence items
`Tiara:<share>`, `CheckM2:<completeness>/<contamination>` and
`geNomad-minority:<share>` appear in the supporting and conflicting lists where used.
The Tiara item of a eukaryote line shows the eukaryote plus organelle share, while the
columns keep the two shares apart.

The ledger retains the organelle lineage, confidence, and independent evidence.
GVClass `d_PHAGE` calls use the viral routing and conflict rules. Native confidence
labels, including `low_support`, are retained without an additional routing gate.

Route labels select downstream analyses. They do not establish final taxonomy
or MAG acceptance.

The evidence ledger retains every linked geNomad call. Only whole-contig calls
count toward the viral share of a bin. Integrated provirus intervals describe
elements within a host; they do not supply whole-bin viral support. When a host
bin is routed to a cellular gene caller, its integrated intervals are excluded
from the unbinned viral gene input.

Each routed eukaryotic bin gets a BUSCO run against the offline `eukaryota_odb12`
lineage (genome mode, miniprot) and a `tidk explore` search for telomeric repeats.
Both are quality evidence for review; neither changes a route. The catalog stores the
BUSCO complete percentage as `completeness_percent` and the duplicated percentage as
`contamination_percent`. Their bases name BUSCO and its lineage. Duplicated markers
are not a contamination estimate. Fragmented and missing percentages go to the stage
notes. The tidk summary and per-bin tables are retained as artifacts only. A run with
no eukaryotic bins records both stages as skipped with
`no_routed_eukaryotic_bins`.

Dedicated gene calling and functional annotation are optional and off by default.
`--enable-gene-calling` turns them on, and `--skip-gene-calling` states the default;
the two options exclude each other. Without gene calling, the gene and annotation
stages are recorded as skipped with the reason "Gene calling not requested
(skip_gene_calling=true)." A resumed run keeps its saved choice unless a gene-calling
option is explicit.

Repeat modeling and masking are optional and off by default. With
`--enable-gene-calling --enable-repeat-masking`, RepeatModeler supplies the repeat
library unless a softmasked assembly is provided. Preparation requires the Dfam
database and the RepeatModeler environment only in that combination.
The workflow confirms a RepeatScout no-seed failure by rerunning its native child
command on the sampled sequence, then retries RepeatModeler once with `-skipRS`.
The same retry applies to a verified empty-refinement failure when every retained
RepeatScout family has valid range evidence and fewer than five native-parsed instances.
The logs and masking provenance record this RECON-only discovery route.
A completed run with no discovered families retains the input sequence.
A blank family count requires explicit zero-family round output and full input coverage.
No discovered families does not establish that the assembly lacks repeats.
Other native failures stop the workflow. Stage logs survive task scratch cleanup.

BRAKER3 requires compatible RNA or protein evidence. With repeat masking off,
it uses the routed unmasked bin FASTA and `--softmasking_off`. With masking on,
it uses the softmasked bin FASTA and default softmask handling. The stage
provenance records this input choice. Skipping masking supplies no evidence
about repeat content and can change gene predictions.
RNA alignments must belong to the target bin. Supplied proteins need an explicit
lineage label. RNA alone selects ET mode; proteins alone select EP
mode. Long-read RNA with proteins uses the pinned experimental ETP variant.
Short-read RNA with proteins remains pending because that variant uses
StringTie long-read mode. PolyA selection alone does not establish the RNA
sequencing technology.

A valid zero-gene result retains the native output files and source metadata.
When no called proteins remain, functional annotation is skipped with
`no_called_proteins` and the catalog retains the original contigs.

Each registered application clone pins a `git_commit`. Preparation fails unless the
clone is at that commit with a clean working tree, writes the commit to the run
configuration, and includes the clone's commit in the source digest. The resume
identity records the effective ONT minimum read length, also when the default
applies, and a resume must request the same length. A run prepared before the
identity held this key still resumes when the option is omitted.

A tool in a grouped stage runs in its own process. Its stage manifest can name the
saved `command_file` of that tool, and the catalog records that script as the
stage command.

LinkML defines the catalog record contract. Generated Pydantic models validate
values, followed by cross-record checks for identifiers, relationships,
coordinate bounds, and expected output keys. SQLite publication is atomic. The
report reads the validated database and records its identity.

The tests must distinguish process wiring, parser fixtures, actual tool execution,
and biological recovery. A small read subset tests interfaces and resource use;
absence of a low-abundance genome at that depth is not evidence of a pipeline bug.
