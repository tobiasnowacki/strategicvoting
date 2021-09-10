code summary
serverside
* run replication_server.R to get main results + data
    - this should also produce conjecture plots...
* repeat with different values
* run replication_random_start_gen.R (with 10 different command line args) to get random iterations

client-side
* run replication_random_start_bind.R to bind these together + produce plots
* run replication_client.R to get distances and main summary plot
    * this is not fully written out yet - I must have stopped at some point
* run replication_v_vecs.R to get ternary paths
* run replication_summary.R to get main results table


* replication_cses.R is obsolete (?)
* what about summary_new.R (?)

* pivotprobs investigation... what is up with IRV prevalence?
    * re-run server_small locally but with old code (w/o normalising)
    * OK, let's see what this brings.
    * when running MC function ks library is not always loaded
    * very weird noise...
    