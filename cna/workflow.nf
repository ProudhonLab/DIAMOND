

params.sampleListpath = "${workflow.projectDir}/listSample.csv" // a csv file with Sample Id and if the sample is a healthy control or not
params.dirOutput= "${workflow.projectDir}/Zscore"
params.dirInput = "${workflow.projectDir}/raw"
params.condaPath= "${workflow.projectDir}/environment.yml"


process filteringMapping {
  conda "${params.condaPath}"
  executor = 'slurm'
  memory '25 GB'
  maxForks 9


  input:
  tuple val(sample_id) , val(ctrlStatus), val(amp_id)

  output:
  tuple val(sample_id) , val(ctrlStatus), val(amp_id),  path('*.bed'), emit: put

  script:
  """
	samtools sort ${params.dirInput}/${sample_id}.L1HS_${amp_id}.bam -o sort.${sample_id}.L1HS_${amp_id}.bam
	samtools view -bq 15 sort.${sample_id}.L1HS_${amp_id}.bam > filtered.${sample_id}.L1HS_${amp_id}.bam
	bamToBed -i filtered.${sample_id}.L1HS_${amp_id}.bam > ${sample_id}.L1HS_${amp_id}.bed

  """
}



process IntervalArm {
  conda "${params.condaPath}"
  executor = 'slurm'
  memory '3 GB'
  maxForks 8

  input:
   tuple val(sample_id) , val(ctrlStatus), val(amp_id),  path(x)

  output:
   tuple  val(sample_id) , val(ctrlStatus), val(amp_id),  path('*.csv'), emit: put

  script:
  """
CalculateDensity.r ${sample_id} ${amp_id} ${x} ${ctrlStatus} TRUE ${workflow.projectDir}/PositionArm.csv


  """
}

process fuseAmp {
  executor = 'slurm'
	input:
       tuple val(sample_id) , val(ctrl), path(samplePath)
    output:
	tuple val(sample_id) , val(ctrl), path('*_amp.txt')

	"""
	cat ${samplePath}   > ${sample_id}_amp.txt

	"""
}


process fuse {
  executor = 'slurm'
	input:
       tuple val(sample_id) , val(ctrl), path(samplePath) , path(ctrlPath)

	output:
       tuple val(sample_id) , val(ctrl), path("*_final.txt")
	"""
	cat ${samplePath}   ${ctrlPath} > ${sample_id}_final.txt

	"""
}

process Zscore {
  conda "${params.condaPath}"
  executor = 'slurm'
  	publishDir "${params.dirOutput}" , mode: 'copy', overwrite: true, pattern: '*.csv'
	input:
       tuple val(sample_id) , val(ctrl), path(samplePath)

	output:
       tuple val(sample_id) , path("*Genome_wideZscoreByChrArm.csv"), path("*ZscoreByChrArm.csv")
	"""
ZscoreCalculating.r ${samplePath} ${sample_id} TRUE 7q

	"""
}


workflow {
  AmpList = Channel.from(2, 4, 6, 8, 7 , 9)

 sampleList=Channel.fromPath("${params.sampleListpath}")
		.splitText() { it.strip() }
		.map {tuple( it.split('\\;')[0], it.split('\\;')[1] )}
		.combine(AmpList)
 filteringMapping(sampleList)
 IntervalArm(filteringMapping.out.put)

 fuseAmp(IntervalArm.out.put
 	.groupTuple(by: [0,1])
 	.map {id, ctrl , amp, files -> [id,ctrl, files]})
 	.branch { name, status,  txt ->
        control: status == "TRUE"
			return txt
        exp: status == "FALSE"
		}
       .set { prepFile }


	fileControl = prepFile.control.collectFile(name: 'sample.txt').view()
	Zscore(fuse(prepFile.exp.combine(fileControl)))


}


workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"

}
