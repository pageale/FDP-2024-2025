import os
import sys

def read_outputs(input_folder, output_csv):
    params_files = os.listdir(input_folder)
    RES = {}

    for f in params_files:
        if f.endswith(".params.txt"):
            sample = f.replace('.params.txt','')
            '''
            here, i could add a conditional for removing _R in sample name
            '''
            RES[sample] = {}
            with open(f'{input_folder}/{f}', 'r') as results:
                for ln in results:
                    ln = ln.rstrip()
                    if "Tumor Fraction" in ln:
                        ln = ln.split(':')[-1].strip().replace("\t", "")
                        RES[sample]["purity"] = ln
                    if "Ploidy" in ln:
                        ln = ln.split(':')[-1].strip().replace("\t", "")
                        RES[sample]["ploidy"] = ln
                    if "Subclone Fraction" in ln:
                        ln = ln.split(':')[-1].strip().replace("\t", "")
                        RES[sample]["scFrac"] = ln
                    if "Fraction Genome Subclonal" in ln:
                        ln = ln.split(':')[-1].strip().replace("\t", "")
                        RES[sample]["scGenome"] = ln
                    if "Fraction CNA Subclonal" in ln:
                        ln = ln.split(':')[-1].strip().replace("\t", "")
                        RES[sample]["scCNA"] = ln

    with open(output_csv, 'w') as out:
        header = ["SampleName", "TumorFraction", "Ploidy", "SubcloneFraction", "SubcloneGenomeFraction", "FractionCNASubclonal"]
        out.write('\t'.join(header) + '\n')

        for sample in sorted(RES):
            args = [
                sample,
                RES[sample].get("purity", "NA"),
                RES[sample].get("ploidy", "NA"),
                RES[sample].get("scFrac", "NA"),
                RES[sample].get("scGenome", "NA"),
                RES[sample].get("scCNA", "NA")
            ]
            out.write('\t'.join(str(x) for x in args) + '\n')

    print(f"Created file: {output_csv}")
    return

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(1)

    input_folder = sys.argv[1]
    output_csv = sys.argv[2]

    read_outputs(input_folder, output_csv)
