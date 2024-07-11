'pigun_generator.py', in the AFSrepo repository, is the file to generate the 'piguns.slcio' file containing the pion guns for the analysis.
Default: flat in E between 5 and 300 GeV, flat in theta between 10° and 170° (can be modified)

-------------------------------------------------------------------------

The other files in the repositories are for generating the lctuple output root file used to evaluate the pion reconstruction efficiencies.

To generate events with HTCondor:

- Paste the 'AFSrepo' repository in your AFS repository
- Paste the 'EOSrepo' in your EOS repository
- In the 'submit.sub' file, replace the 'EOSrepo' occurrences with the actual path of the 'EOSrepo' repository
- if you want the .slcio output of reconstruction available, uncomment the LCIOWriter processor in the 'EOSrepo/mucoll-benchmarks/digirecolctuple.xml'; then go to the 'submit.sub' file and, following the syntax as for the 'output_lctuple.root' and 'piguns.slcio' files, set 'output_reco.lcio' as an output in transfer_output_remaps
- In the 'submit.sub' file, set the number of jobs you want to run (default=220)
- In the 'runfile.sh' file, specify the number of events per job (default=500)
- from the 'AFSrepo' repository, launch 'condor_submit submit.sub'
- wait for the jobs to finish
- go to the 'EOSrepo' repository, launch 'hadd -f lctuple_piguns.root Outputs/output_lctuple*' and wait for it to finish
- from the 'EOSrepo' repository, launch 'root -l -b rewrite.C' and wait for it to finish
- to produce the plots, from the 'EOSrepo' repository launch 'root -l -b lctuplePiAnalysis.C'


 