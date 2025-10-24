# HiC data: statistics and properties

The functions in this extension are meant to be used with HiC data
in order to assess the evidence for two contigs being adjacent
to each other in the genome. In the ideal case the result would
be:

1. An estimate of the most likely length seperating two contigs.
2. Some likelihood estimate of linkage at that distance.

There are many existing HiC-scaffolders and manual curation programs
that work well and that are widely used, so why bother with this as well?

My reasons for doing this work are to try to understand the reasons
when some of the usual approaches don't work. In particular I have
been involved in the assembly of a teleost genome as part of a wider
effort. That assembly produced 23 normal sized 'chromosomes' scaffolds
and one rather tiny one. My own observations of the data suggest that
this small chromosome should be linked to one of the long ones. But,
neither the scaffolder nor the manual curation seem to think so.

I suspect that the reasons are related to the presence of long regions
of low complexity sequence. Such sequences have very few unique sequences
that can be used to anchor them to other scaffolds unambiguously. I
don't know how the usual programs deal with such issues, but my own
observations of the HiC-data do seem to indicate linkage. However,
my attempts to estimate distance between scaffolds has not been very
succesful and I suspect that this is related to variance in mappability
that deviates from the nice random models that one can make.

This work looks at how we can incorporate long-range links to increase
the amount of evidence for linkage when it is otherwise scarce. It
is based on the observation that the distribution of log-transformed
link distances appears to be uniform. This uniformity will be disturbed
if the distance between scaffolds (or scaffold regions) is wrong.
Unfortunately, the uniformity is also disturbed by difference in
mapping rates and mappability at both locus regions. To compensate for
that requires an estimate of the density of reads at the two regions;
unfortunately that density should be in terms of the log(distance) which
makes it more difficult to estimate in a linear manner. Here I'm
considering the feasibility of the more or less exhaustive 
