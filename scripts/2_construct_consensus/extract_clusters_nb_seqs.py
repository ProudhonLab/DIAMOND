import argparse
import glob

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputDir", type=str, default=None, help="Directory containing consensus sequences to process")
    parser.add_argument("-o", "--outputDir", type=str, default=None, help="Directory to save output")
    args = parser.parse_args()

    inputDir = args.inputDir
    outputDir = args.outputDir

    consensusFastaList = glob.glob(f"{inputDir}/*.consensus.fasta")

    for file in consensusFastaList:
        outputFile_name = file.split("/")[-1].split(".consensus.")[0]
        with open(file, "r") as inputFile, open(f"{outputDir}/{outputFile_name}.consensus.nb_seqs.csv", "w") as outputFile:
            for line in inputFile:
                if line.startswith(">"):
                    nb_reads = line.split(";seqs=")[-1]
                    outputFile.write(nb_reads)