###-----------------------------------------------------------------------------------------------------------------###
# Notes on where variables are defined, if not defined in qc.smk
#
# output_directory - workflow/Snakefile
# config           - workflow/Snakefile
# SAMPLES          - workflow/Snakefile
#
###-----------------------------------------------------------------------------------------------------------------###

if config['region_centered_bin_averages']['run']:
	# intervals associated with a given interval (input, in bed format) is the set of intervals positioned with respect to the input interval.
	rule get_regions_for_centered_bin_averages_from_bed:
		input:
			region_bed = config['region_centered_bin_averages']['region_file'],
		params:
			args_list = config['region_centered_bin_averages']['args_list'],
			out_directory = f'{output_directory}/analysis/region_centered_bin_averages/regions/',
		output:
			region_centered_bins = f'{output_directory}/analysis/region_centered_bin_averages/regions/regionCenteredBins.bed'
		log:
			f'{output_directory}/logs/region_centered_bin_averages/build_region.log',
		resources:
			mem_gb = config['hpcParameters']['smallMemoryGb'],
			walltime = config['walltime']['medium'],
		shell:
			"""
			# 
			# ~ mkdir -p {params.out_directory} && touch {output.region_centered_bins}
			python3 workflow/scripts/wzmetagene.py {params.args_list} {input.region_bed} | sort -k1,1 -k2,2n > {output.region_centered_bins}
			"""

	rule region_centered_bin_averages:
		input:
			merged_sampleBed = f'{output_directory}/analysis/pileup/{{sample}}_mergecg.bed.gz',
			region_bins = f'{output_directory}/analysis/region_centered_bin_averages/regions/regionCenteredBins.bed'
		output:
			bed = f'{output_directory}/analysis/region_centered_bin_averages/{{sample}}.bed',
		log:
			f'{output_directory}/logs/region_centered_bin_averages/{{sample}}.log',
		benchmark:
			f'{output_directory}/benchmarks/region_centered_bin_averages/{{sample}}.txt',
		threads: 1,
		resources:
			mem_gb = config['hpcParameters']['smallMemoryGb'],
			walltime = config['walltime']['medium'],
		conda:
			'../envs/general.yaml'
		envmodules:
			config['envmodules']['python3'],
			config['envmodules']['bedtools'],
		shell:
			"""
			# do the intersecting for this sample
			# column 14 has the methylation
			bedtools intersect \
			-a {input.region_bins} \
			-b {input.merged_sampleBed} \
			-sorted -wo | \
			bedtools groupby \
			-g 1-10 \
			-c 14 \
			-o mean | \
			awk '{{print $7,"\\t",$8,"\\t",$9,"\\t",$4,"\\t",$10,"\\t",$11}}' > {output.bed}
			
			"""
