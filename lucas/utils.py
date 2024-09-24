import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.colorbar import ColorbarBase
from math import ceil
from glob import glob

plt.rc('font', family='Arial')

RK2 = "MNRTFDRKAYRQELIDAGFSAEDAETIASRTVMRAPRETFQSVGSIVQQATAKIERDSVQLAPPALPAPSAAVERSRRLEQEAAGLAKSMTIDTRGTMTTKKRKTAGEDLAKQVSEAKQAALLKHTKQQIKEMQLSLFDIAPWPDTMRAMPNDTARSALFTTRNKKIPREALQNKVIFHVNKDVKITYTGVELRADDDELVWQQVLEYAKRTPIGEPITFTFYELCQDLGWSINGRYYTKAEECLSRLQATAMGFTSDRVGHLESVSLLHRFRVLDRGKKTSRCQVLIDEEIVVLFAGDHYTKFIWEKYRKLSPTARRMFDYFSSHREPYPLKLETFRLMCGSDSTRVKKWREQVGEACEELRGSGLVEHAWVNDDLVHCKR"
BBR1 = "MATQSREIGIQAKNKPGHWVQTERKAHEAWAGLIARKPTAAMLLHHLVAQMGHQNAVVVSQKTLSKLIGRSLRTVQYAVKDLVAERWISVVKLNGPGTVSAYVVNDRVAWGQPRDQLRLSVFSAAVVVDHDDQDESLLGHGDLRRIPTLYPGEQQLPTGPGEEPPSQPGIPGMEPDLPALTETEEWERRGQQRLPMPDEPCFLDDGEPLEPPTRVTLPRR"
PSA = "MPKNNKAPGHRINEIIKTSLALEMEDAREAGLVGYMARCLVQATMPHTDPKTSYFERTNGIVTLSIMGKPSIGLPYGSMPRTLLAWICTEAVRTKDPVLNLGRSQSEFLQRLGMHTDGRYTATLRNQAQRLFSSMISLAGEQGNDFGIENVVIAKRAFLFWNPKRPEDRALWDSTLTLTGDFFEEVTRSPVPIRIDYLHALRQSPLAMDIYTWLTYRVFLLRAKGRPFVQIPWVALQAQFGSSYGSRARNSPELDDKARERAERAALASFKYNFKKRLREVLIVYPEASDCIEDDGECLRIKSTRLHVTRAPGKGARIGPPPT"
PVS1 = "MSGRKPSGPVQIGAALGDDLVEKLKAAQAAQRQRIEAEARPGESWQAAADRIRKESRQPPAAGAPSIRKPPKGDEQPDFFVPMLYDVGTRDSRSIMDVAVFRLSKRDRRAGEVIRYELPDGHVEVSAGPAGMASVWDYDLVLMAVSHLTESMNRYREGKGDKPGRVFRPHVADVLKFCRRADGGKQKDDLVETCIRLNTTHVAMQRTKKAKNGRLVTVSEGEALISRYKIVKSETGRPEYIEIELADWMYREITEGKNPDVLTVHPDYFLIDPGIGRFLYRLARRAAGKAEARWLFKTIYERSGSAGEFKKFCFTVRKLIGSNDLPEYDLKEEAGQAGPILVMRYRNLIEGEASAGS"
THREADS = 4
SAVE_FIGURES = False

MAGENTA_CMAP = LinearSegmentedColormap.from_list("CustomCmap", ["white", "magenta"] , N=100)

def run_fastqc(in_dir, out_dir):
    os.subprocess(['fastqc',in_dir,'-o','../02-OUTPUT/01-QC/'])
    return None

def run_multiqc(in_dir, out_dir):
    os.subprocess(['multiqc',in_dir,'-o','../02-OUTPUT/01-QC/'])
    return None   

def run_trimmomatic(in_dir,out_dir,params):
    for seq in glob(f"{in_dir}/*R1*.fastq.gz"):
        handle = seq.split('/')[-1].replace("_R1_001.fastq.gz","")
        r1 = handle + "_R1_001"
        r2 = handle + "_R2_001"
        os.system(f"trimmomatic PE -threads 4 -phred33 -trimlog {out_dir}/02-trimmomatic-trim_log_{handle}.txt\
        {in_dir}/{r1}.fastq.gz {in_dir}/{r2}.fastq.gz \
        {out_dir}/{r1}.trimmed.fastq.gz {out_dir}/{r1}un.trimmed.fastq.gz \
        {out_dir}/{r2}.trimmed.fastq.gz {out_dir}/{r2}un.trimmed.fastq.gz \
        {params}")
    return None     

def build_bowtie2_index(in_file, out_file):
    os.subprocess(f"bowtie2-build {in_file} {out_file}")
    return None

def run_bowtie2_alignment(in_dir, out_dir):
    for seq in glob(f"{in_dir}/*R1*.fastq.gz"):
        handle = seq.split('/')[-1].replace("_R1_001.fastq.gz","")
        r1 = handle + "_R1_001"
        r2 = handle + "_R2_001"
        os.system(f"bowtie2 -q --phred33 --no-unal \
        -x {out_dir}/Gent_RK2 \
        -1 {in_dir}/{r1}.trimmed.fastq.gz \
        -2 {in_dir}/{r2}.trimmed.fastq.gz \
        -S {out_dir}/{handle}_mapping.sam")
    return None

def run_pysamstats(in_dir,out_dir):
    for sam in glob(f"{in_dir}/*.sam"):
        handle = sam.replace('.sam','').split('/')[-1]
        os.system(f"samtools view -bS {in_dir}/{handle}.sam > {out_dir}/{handle}.bam")
        os.system(f"samtools sort {out_dir}/{handle}.bam > {out_dir}/{handle}.sorted.bam")
        os.system(f"samtools index {out_dir}/{handle}.sorted.bam")
    for bam in glob("../02-OUTPUT/04-ANALYZE/*.sorted.bam"):
        os.system(f"pysamstats \
        --fasta ../01-INPUT/03-REFERENCE/Gent_RK2.fasta \
        --type variation \
        --max-depth 10000000 \
        {bam} > {bam}.var.txt")

def flatten(input):
    new_list = []
    for i in input:
        for j in i:
            new_list.append(j)
    return new_list

def plot_heatmap(df,col,plot_ax,title,min_cbar,max_cbar,grouped):
    g = sns.heatmap(data=[df[col]],
                    xticklabels=True,
                    yticklabels=False,
                    cmap=MAGENTA_CMAP,
                    ax=plot_ax,
                    vmin=min_cbar,
                    vmax=max_cbar,
                    cbar=None,
                    alpha=1.0)
    for i in range(df.shape[0]):
        rect = Rectangle((i, 0), 1, 1, fill=False, edgecolor='black', lw=0.5)
        g.add_patch(rect)
    if not grouped:
        g.set_xticklabels(labels=list(title),rotation=0)
    else:
        g.set_xticklabels(labels=[],rotation=0)
    g.tick_params(axis='both',which='both',length=0)
    g.set_ylabel(g.get_ylabel(), fontdict={'weight': 'bold','size':12})

def plot_residues(df, col, num_cols, seq, id, grouped=False, group_size=3, out_dir='.'):
    if grouped:
        df['number'] = df.index // (1 * group_size)
        df = df.groupby(['number']).max().reset_index()
        figsize = (12,ceil(len(df)/num_cols))
    else:
        figsize = (10, ceil(len(df)/num_cols)*0.5)

    fig = plt.figure(figsize=figsize)
    gs = GridSpec(ceil(len(df)/num_cols),2,width_ratios=[80, 2],hspace=1,wspace=0.025)

    min_cbar, max_cbar = df[col].min(), df[col].max()
    range_list = list(range(0, len(df), num_cols)) + [len(df)]
    range_tuples = [(range_list[i], range_list[i + 1]) for i in range(len(range_list) - 1)]
    for i,(min_res,max_res) in enumerate(range_tuples):
        subseq = seq[min_res:max_res]
        plot_heatmap(df.iloc[min_res:max_res,:],
                     col=col,
                     plot_ax=fig.add_subplot(gs[i,0]),
                     title=subseq,
                     min_cbar=min_cbar,
                     max_cbar=max_cbar,
                     grouped=grouped)
    cbar = ColorbarBase(fig.add_subplot(gs[:,1]), cmap=MAGENTA_CMAP,orientation='vertical',norm=Normalize(min_cbar, max_cbar))
    cbar.outline.set_linewidth(0.5)

    if SAVE_FIGURES:
        fig.savefig(f'{out_dir}/{id}_cols{num_cols}_grouped{grouped}.png',transparent=False,bbox_inches='tight',dpi=400,facecolor='white')