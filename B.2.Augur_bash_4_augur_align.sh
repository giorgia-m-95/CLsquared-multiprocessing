for i in {1..20}
do
	bash Augur_align.sh -i /beegfs/labtoppo/giorgiam/Analisi_Napoli/Outs_Dengue/Outs_parallel_timing_${i}/CL_squared_step_02/Temp_FASTA -r /beegfs/labtoppo/giorgiam/Analisi_Napoli/Outs_Dengue/ref_sequence.fasta -o /beegfs/labtoppo/giorgiam/Analisi_Napoli/Outs_Dengue/Outs_parallel_timing_${i}/ -t 10
done
