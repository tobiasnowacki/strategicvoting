ok, so the problem we have is that small numbers (e.g. 1e-16) are being rounded to effectively zero.
R does not recognise the difference.

Andy suggests normalisation as a solution.
One other idea would be to specify a number format with greater precision using Rmpfr.
We've used that in the past but found it to be awfully slow and inefficient.

Maybe run this again using Rmpfr and benchmark calculation time?
If we're using the server, maybe it's not too bad.
