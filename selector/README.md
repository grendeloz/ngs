# selector

A Selector is a string of the form operation:subject:pattern which can
be used to define select/filter operations on collections of records.

For example, selectors that could be used with GFF3 records might look
like:

  keep:seqid:^GL
  delete:type:.*_UTR

In general, the effect of all selectors is to drop records from some
collection of records. Selectors with delete operations drop any record
that matches the pattern while keep operations drop any record that does
not match the pattern - in both cases records are dropped.

As a general rule, a keep operation is a blunter tool than a delete
operation because every record you wish to keep must match the keep
pattern. Because delete only drops matching records, multiple delete 
selectors, applied sequentially, with tight patterns can be used to
selectively prune away records that you don't wish to retain.

In cases where multiple selectors are allowed, they shoudl probably be
applied sequentially in the order in which they are specified on the
command line, in the config file, etc.

The colon character ':' must not be used in the subject, operation or
pattern of a selector - it is strictly reserved as a separator for the
selector. The package has no easy way to enforce this requirement so
it's up to the user to watch for too many elements in a selector which
might suggest improper use of the colon.
