saber - Sequence Analysis of Behavioural Events in R
version 1.1 - 30th Oct 2016
by Andrey Barsky
andreybarsky@gmail.com

This package is designed to process discrete behavioural event data and identify
structural dependencies between events based on the strength of their associations.

The main function is also called saber, and will perform all of the processing steps
to output an R object containing observed and expected frequencies of event
transitions, frequency residuals, and a summary of transitions that occur notably
more often than expected.

Simply running saber() with no arguments will bring up a file select prompt. If
you select your sequence data (as plain text, csv, or excel format) it should run
fine in many cases with the defaults. However, it is advised you tune parameters
according to your specific data.

Here are the default arguments, and an explanation of each one:

Usage

saber(data=NULL, 
      datatype='filename', 
      terminals=TRUE, 
      lag=1, 
      event_sep=' ', 
      seq_sep = '\n', 
      residual='standard', 
      critical='acceleration', 
      min_obs=5, 
      lfe_collapse='default', 
      cluster_height=0, 
      n_recombine=0, 
      recombine_crit='default', 
      recombine_obs='default', 
      quiet=FALSE)
      
Arguments
      
  data: The input data to work from. This is either a path to a file; a string; or
a character vector. If NULL, saber brings up an interactive file prompt to choose a
file directly.
The expected input format in plaintext looks like this:
a b c d e
a e d
a b c e
et cetera, with events separated by spaces, and sequences separated by linebreaks,
although this can be changed with the event_sep and seq_sep arguments. Events need not
be single characters, and can be any length, but avoid using the '@' and '-' characters
(which are used in event recombination).

  datatype: Tells saber what the 'data' argument is referring to. Can be 'filename', 'string', or 'vector'. This will usually be 'filename' but some internal functions
work directly from pre-processed character vectors.

  terminals: If TRUE, each sequence in the input data is modified to begin with a new
event called START, and end with a new event called END. This can be useful, but is
not always necessary depending on what your data means; in particular, it should be
set to FALSE if your sequence data already has explicit start and end events.

  lag: The lag parameter for lag sequential analysis. The default value of 1 means that
an event pair consists of neighbouring events in a sequence (the sequence 'a b c'
produces event pairs 'a-b' and 'b-c'), but you can set it higher for e.g. dyadic conversational data.

  event_sep: Tells saber how to read plaintext input data. If something other than 
space is used to separate events in your sequences, specify it here.

  seq_sep: Tells saber how to read plaintext input data. If something other than 
a linebreak is used to separate sequences in your data file, specify it here.

  residual: Which residual to use when calculating deviation from expected transitional
frequencies in the contingency table. Can be 'standard' or 'adjusted'.
The default is the 'standardised residual', which is given by:
(cell.obs - cell.exp) / sqrt(cell.exp)
Some authors recommend the alternative 'adjusted residual', given by:
(cell.obs - cell.exp) / cell.exp * (1 + row_proportion) * sqrt(1 - column_proportion)
The adjusted residual has some advantages, but has the property of occasionally
hitting divide by zero errors when absolute dependencies exist, which R ordinarily
interprets as infinities in the residual matrix. saber gets around this by replacing
infinity with the highest real number in the residual matrix, although this might
not be strictly mathematically sound. The standard residual is recommended.

  critical: The cutoff residual value for a transition to be considered notably
high-frequency. This is either a positive number, or a method to calculate the
cutoff automatically, which is one of 'acceleration' or 'projection'.
Since the plot of residual values against transitions is effectively a scree plot
as encountered in factor analyis, we can determine a non-arbitrary point at which
transition residuals become unimportant by employing various non-graphical solutions
to the scree test.
The 'acceleration' method attempts to find the elbow of the curve by locating the
global maximum of the second derivative to the smoothed curve. The curve is smoothed
using LOESS (smoothing factor 0.5), and the sharpest turn is used as the critical
residual value.
The 'projection' method uses vector projection to find the point on the curve that is
furthest away from the straight line that joins the first and last point on the curve.
It should work well most of the time, but is sensitive to the length of the plot's tail.
Both methods will, by default, produce a plot that visualises this process.

  min_obs: The minimum number of times a transition must be observed in order to be
considered high-frequency. It is recommended you set this manually, although there is
no particular rule of thumb for this - it is meant mainly to avoid drawing attention
to transitions that have high residuals but low observed frequencies, since those
are likely to be false positives.

  lfe_collapse: It is sometimes useful to recode all events that occur rarely as one
single event, particularly if there are very many event codes. Here you can set a
cutoff point for the number of times an event must be observed to not be considered
low-frequency - e.g. if this is 3, then all events that occur less than 3 times will
be collapsed into the new event 'LFE'. 
If lfe_collapse is 'default', then it automatically sets the cutoff to one standard
deviation below the mean of observed event frequencies.
Note that this might become confusing if your data already contains events called LFE.

  cluster_height: As an alternative (or supplement) to collapsing events based on their
observed frequencies, saber can collapse events that are similar to one another if
they have similar profiles of transitional dependencies.
The logic of this method is that if two events tend to be preceded by similar events,
and followed by similar events, then those two event codes might just be instantiations
of a single 'latent' event. This is essentially dimension reduction for sequences.
saber's clustering method runs hierarchical cluster analysis on a distance matrix
based on the correlations between event antecedents and sequiturs. The hierarchical
clustering dendrogram is then cut by the parameter specified in cluster_height. The
height of the tree is 1, and the height of each clustering is equal to 1 minus the
Pearson correlation between those two events.
If this is 0, no clustering is performed; but if set between 0 and 1, it will collapse
events into new cluster events called X1, X2 etc., based on the cutoff height specified.
This also forces saber to output a 'clusters' variable that shows which original
events belong to which clusters, and an 'event_cor' variable that shows the correlation
matrix.
By default, saber outputs a plot that shows the dendrogram and the cut height.

  n_recombine: By default, saber models the input data as a first-order Markov chain,
where each event's probability is based only on the preceding event (at a distance
determined by the lag parameter). In principle, higher-order dependencies can be
modelled by adding new dimensions to the matrix, but this will usually produce an
extremely sparse and unwieldy matrix due to combinatiorial explosion.
However, there is a clever method to selectively add higher-order dependencies to our
contingency table - simply look at transitions that occur significantly often and
recode those event pairs as a new single event. That is, if the pair 'a b' occurs
frequently in the data, we collapse those into a single 'a-b' event and run the
analysis again. This way, the events 'a', 'b', and 'a-b' all occur independently in
the 2-dimensional matrix of transitional frequencies.
saber uses a recursive function to do this as many times as is required, specified by
this n_recombine argument. Set it to 1 to perform a single recombination loop; more
than 2 should not usually be necessary. 

  recombine_crit: Critical residual value for event recombination for higher-order
transitions. By default, this is double the critical residual value for notability,
although it can be set to any real number.

  recombine_obs: Minimum number of observations for event recombination as above. By
default this is twice the mean of observed nonzero transitional frequencies, although
it can be set to any real number.

  quiet: If FALSE, saber outputs plots for critical residual calculation and
hierarchical clustering, shows observed and residual transition matrices in the
viewer, and outputs the notably high frequency transitions to console. Set this to
TRUE to suppress all of those.


Values

The output is a list object of class 'saber', with the following attributes:

$input: The original input that was fed to the main saber function.

$seqdata: How the input was parsed into a character vector of distinct sequences.

$alphabet: A table of the unique event codes encountered and their frequencies.

$observed: Observed frequency matrix of event transitions.

$expected: Expected frequency matrix of event transitions.

$residuals: Matrix of event transition residuals.

$critical_residual: The critical residual value for notability, either calculated
by default or supplied by user.

$event_cor: If clustering was performed, these are the correlations between event codes.

$clusters: If clustering was performed, these are cluster memberships.

$hft: A data frame of high frequency transitions, ordered by descending residual values.
