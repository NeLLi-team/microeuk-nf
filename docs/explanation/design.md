# Scientific design

Each sample is treated as a metagenome. The input contract accepts basecalled ONT
DNA or PacBio HiFi reads. RNA reads are linked evidence for gene calling and do
not enter DNA assembly.

Myloasm produces a sample assembly. Read mapping measures contig depth. QuickBin
uses coverage and composition to form bins; unpaired long-read alignments do not
supply paired-read edges. The workflow retains unbinned contigs. QuickClade uses
an explicit local reference, with remote-server fallback disabled.

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
stop the workflow. Completed predictions require one result per input bin.

GVClass uses the database path from the registry. Its native run receipt must
match that resolved path and the declared version before results are accepted.

QuickClade prokaryotic calls can support a prokaryotic route. Its eukaryotic
labels stay in the evidence table but do not establish a eukaryotic route or
veto a viral route. A eukaryotic route requires a resolved CheckEUK lineage
plus GVClass eukaryotic support or a linked eukaryotic SSU locus.

GVClass `d_PLASTID` and `d_MITO` calls retain an unresolved route because the
workflow has no organelle gene caller. Supported eukaryotic, prokaryotic, or viral
evidence changes their candidate class from unresolved to conflicting.
The ledger retains the organelle lineage, confidence, and independent evidence.
GVClass `d_PHAGE` calls use the viral routing and conflict rules. Native confidence
labels, including `low_support`, are retained without an additional routing gate.

Route labels select downstream analyses. They do not establish final taxonomy
or MAG acceptance.

The evidence ledger retains every linked geNomad call. Only whole-contig calls
support a viral route for the bin. Integrated provirus intervals describe
elements within a host; they do not supply whole-bin viral support. When a host
bin is routed to a cellular gene caller, its integrated intervals are excluded
from the unbinned viral gene input.

RepeatModeler supplies the repeat library unless a softmasked assembly is provided.
The workflow confirms a RepeatScout no-seed failure by rerunning its native child
command on the sampled sequence, then retries RepeatModeler once with `-skipRS`.
The same retry applies to a verified empty-refinement failure when every retained
RepeatScout family has valid range evidence and fewer than five native-parsed instances.
The logs and masking provenance record this RECON-only discovery route.
A completed run with no discovered families retains the input sequence.
A blank family count requires explicit zero-family round output and full input coverage.
No discovered families does not establish that the assembly lacks repeats.
Other native failures stop the workflow. Stage logs survive task scratch cleanup.

BRAKER3 requires a softmasked assembly and compatible RNA or protein evidence.
RNA alignments must belong to the target bin. Supplied proteins need an explicit
lineage label. RNA alone selects ET mode; proteins alone select EP
mode. Long-read RNA with proteins uses the pinned experimental ETP variant.
Short-read RNA with proteins remains pending because that variant uses
StringTie long-read mode. PolyA selection alone does not establish the RNA
sequencing technology.

A valid zero-gene result retains the native output files and source metadata.
When no called proteins remain, functional annotation is skipped with
`no_called_proteins` and the catalog retains the original contigs.

LinkML defines the catalog record contract. Generated Pydantic models validate
values, followed by cross-record checks for identifiers, relationships,
coordinate bounds, and expected output keys. SQLite publication is atomic. The
report reads the validated database and records its identity.

The tests must distinguish process wiring, parser fixtures, actual tool execution,
and biological recovery. A small read subset tests interfaces and resource use;
absence of a low-abundance genome at that depth is not evidence of a pipeline bug.
