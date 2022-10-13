# gene family analysis and phylogeny tree construction

## 提前配置

### 软件

|    软件     |     版本     |                             引用                             |
| :---------: | :----------: | :----------------------------------------------------------: |
| orthofinder |    2.3.12    | [Emms, D.M., Kelly, S. OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome Biol 20, 238 (2019)](https://doi.org/10.1186/s13059-019-1832-y)<br />[Emms, D.M., Kelly, S. OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy. Genome Biol 16, 157 (2015)](https://doi.org/10.1186/s13059-015-0721-2) |
| jmodeltest2 |    2.1.10    | Darriba D, Taboada GL, Doallo R, Posada D. 2012. jModelTest 2: more models, new heuristics and parallel computing. Nature Methods 9(8), 772.<br />Guindon S and Gascuel O (2003). A simple, fast and accurate method to estimate large phylogenies by maximum-likelihood". Systematic Biology 52: 696-704. |
|   pal2nal   |     v14      | Suyama, M., Torrents, D., & Bork, P. (2006). PAL2NAL: robust conversion of protein sequence alignments into the corresponding codon alignments. *Nucleic acids research*, *34*(Web Server issue), W609–W612. https://doi.org/10.1093/nar/gkl315 |
|    phyml    | 3.3.20190909 | Guindon, S., & Gascuel, O. (2003). A simple, fast, and accurate algorithm to estimate large phylogenies by maximum likelihood. *Systematic biology*, *52*(5), 696–704. https://doi.org/10.1080/10635150390235520 |
|  raxml-ng   |              | Alexey M. Kozlov, Diego Darriba, Tomáš Flouri, Benoit Morel, and Alexandros Stamatakis (2019) **RAxML-NG: A fast, scalable, and user-friendly tool for maximum likelihood phylogenetic inference.** *Bioinformatics, 35 (21), 4453-4455* doi:[10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305) |
|    paml     |     4.9d     | Yang, Z. 2007. PAML 4: a program package for phylogenetic analysis by maximum likelihood. Molecular Biology and Evolution 24: 1586-1591 |



### Perl 模块安装

需要预先安装Bio::SeqIO模块。可以参考[网上教程](https://www.jianshu.com/p/9e90b3524fe2)，这个还是相对简单的。

---

# 2. 文件格式说明

## 2.1 命名

用于分析的文件命名为：[A-Za-z]{3}.fa

![gene%20famil%2097fb2/Untitled.png](gene%20famil%2097fb2/Untitled.png)

---

## 2.2 序列格式

用于分析的reads的名称是>[a-zA-Z]{3}|.*

![gene%20famil%2097fb2/Untitled%201.png](gene%20famil%2097fb2/Untitled%201.png)

---

---

# 3 蛋白质和cds文件预处理

## 3.1 filter protein and cds by length

```bash
perl filter_fasta_by_length_and_sort.pl LENGTH Prefix.fa Prefix.Passed.fa Prefix.unpassed.fa Prefix.length.txt

# LENGTH: 长度小于这个值的序列会被filter掉
# Prefix.length.txt：每条序列的长度，可以用于检查长度分布
```

## 3.2 add species name to sequence header

```bash
perl add_species_name_to_fa_sequences.pl Prefix.fa
```

---

# 4 获取基因家族信息+画树

## 4.1 cluster gene families by *orthofinder*

### 4.1.1 直接运行

```bash
orthofinder -t 10 -a 5 -S diamond -og -f INPUT_DIR

	-t, Number of parallel sequence search threads [Default = 16]
	-a, Number of parallel analysis threads [Default = 1]
	-S, Sequence search program [Default = blast]. Options: blast, blast_gz, diamond
	-og, Stop after inferring orthogroups
	-f, Input dir containing all protein files
```

### 4.1.2 在原有结果的基础上添加物种

### 4.1.3 在原有结果中删除物种

```bash
orthofinder -b previous_orthofinder_directory -f new_fasta_directory

# This will add each species from the 'new_fasta_directory' to existing set of species, reuse all the previous BLAST results, perform only the new BLAST searches required for the new species and recalculate the orthogroups. The 'previous_orthofinder_directory'  is the OrthoFinder 'WorkingDirectory/' containing the file 'SpeciesIDs.txt'.
```

```bash
orthofinder -b previous_orthofinder_directory

#In the 'WorkingDirectory/' from a previous analysis there is a file called 'SpeciesIDs.txt'. Comment out any species to be removed from the analysis using a '#' character and then run OrthoFinder using the command above. The 'previous_orthofinder_directory'  is the OrthoFinder 'WorkingDirectory/' containing the file 'SpeciesIDs.txt'.
```

### 4.1.4 同时进行添加和删除操作

```bash
orthofinder -b previous_orthofinder_directory -f new_fasta_directory
```

---

## 4.2 construct phylogeny tree

### 4.2.0 预准备

将所有的蛋白质文件cat到一个命名为allpepfiles的文件中。并放到工作目录下

将所有的cds文件放在一个命名为all_cds_dir的文件中

![gene%20famil%2097fb2/Untitled%202.png](gene%20famil%2097fb2/Untitled%202.png)

### 4.2.1 提取单拷贝蛋白质序列

**orthofinder会自动生成这个文件夹，这一步只是保险项。**

Orthogroups.txt来自于*Orthofinder*的输出，allpepfiles是用到的所有物种的蛋白质序列的集合，记得更改脚本里第47行的物种数目。

```bash
perl scripts/1.extract_ortholog.pep.pl ./Orthogroups/Orthogroups.txt allpepfiles Single_Copy_Orthologue_Sequences
```

---

### 4.2.2 *mafft*

注意更改脚本里面对应的单拷贝序列前缀的正则匹配。

```bash
perl scripts/2.create_sh.mafft.pl Single_Copy_Orthologue_Sequences mafft.sh mafftseqs
bash mafft.sh >mafft.log 2>>mafft.err 
mkdir mafftout
mv mafft* mafftout
```

---

---

### 4.2.3 提取蛋白质文件对应的cds序列

记得修改第52行的物种数目

```bash
perl scripts/3.extract_ortholog.cds.pl ./Orthogroups/Orthogroups.txt all_cds_dir/ singlecopycdsdir length_check.txt
```

---

### 4.2.4 cds结果排序

```bash
perl scripts/4.sort_cds_mafft.pl mafftoutdir/mafftseqs singlecopycdsdir sortedcdsdir
```

---

### 4.2.5 *pal2nal*

注意更改脚本中的*pal2nal*软件的路径

```bash
perl scripts/5.create_pal2aln.pl mafftoutdir/mafftseqs sortedcdsdir pal2nalfas pal2nal.sh
nohup bash pal2nal.sh >pal2nal.log 2>>pal2nal.err &
mkdir pal2naloutdir
mv pal2nal* err* pal2naloutdir

# 如果存在protein和cds不对应的情况，应当删除这些
```

---

### 4.2.6 *Gblocks*筛选保守位点（这一步可以选做）

如果不能直接调用gblocks，可以修改line10中指定的路径

```bash
perl scripts/6.create.gblocks_sh.pl pal2naloutdir/pal2nalfas gblockshtms gblocks.sh
bash gblocks.sh
mkdir gblocksoutdir
mkdir gblocksgbs
mv pal2naloutdir/pal2nalfas/*htm gblockshtms
mv pal2naloutdir/pal2nalfas/*gb gblocksgbs
mv gblocks* gblocksoutdir
```

---

### 4.2.7 merge cds

```bash
# 不做gblocks的情况
perl scripts/7.cat_cds.pl pal2naloutdir/pal2nalfas mergedcdss

# 做gblocks的情况
perl scripts/7.cat_cds.pl gblocksoutdir/gblocksgbs mergedcdss
```

---

### 4.2.8 生成phylip文件

这里有个很有趣的事情，就是有时候`convertFasta2Phylip.sh`的换行符会出现变化，然后就无法使用，所以建议直接用包里面的脚本。

```bash
mkdir finalseqs
cat mergedcdss/*fa >finalseqs/final.fa
cd finalseqs
sh ../scripts/convertFasta2Phylip.sh final.fa >final.phy
```

---

### 4.2.9 用*jmodeltest2*寻找最佳模型

```bash
mkdir tree
cd tree
ln -s ../finalseqs/final.phy .
nohup java -jar -XX:ParallelGCThreads=4 -Xmx4g ~/software/jmodeltest-2.1.10/jModelTest.jar -tr 10 -d final.phy -s 11 -f -i -g 8 -AIC -BIC -AICc -DT -p -a -w -o bestmodel &

	-tr,	number of threads to execute (default is 32)
	-d,		input data file
	-s,		number of substitution schemes (e.g., -s 11) (it has to be 3,5,7,11,203; 				default is 3)
	-f,		Include models with unequals base frecuencies
	-i,		include models with a proportion invariable sites (e.g., -i) (default is 				false)
	-g,		include models with rate variation among sites and number of categories 				(e.g., -g 8) (default is false & 4 categories)
	-AIC,	calculate the Akaike Information Criterion (e.g., -AIC) (default is 					false)
	-AICc,	calculate the corrected Akaike Information Criterion (e.g., -AICc) 						(default is false)
	-BIC,	calculate the Bayesian Information Criterion (e.g., -BIC) (default is 					false)
	-DT,	calculate the decision theory criterion (e.g., -DT) (default is false)
	-p,		calculate parameter importances (e.g., -p) (default is false)
	-a,		estimate model-averaged phylogeny for each active criterion (e.g., -a) 					(default is false)
	-w,		write PAUP block (e.g., -w) (default is false)
	-o,		set output file (e.g., -o jmodeltest.out)
```

---

### 4.2.10 *phyml*建树

其他参数是根据*jmodeltest*来的，然后-b是bootstrap值

```bash
phyml -i final.phy -d nt -n 1 -b 1000 --run_id 012345 -m GTR+I+G -f m -v e -c 4 -a e --no_memory_check -o tlr -s BEST

	-i,		seq_file_name
			seq_file_name is the name of the nucleotide or amino-acid sequence file in 				PHYLIP format.
	-d,		data_type
			data_type is 'nt' for nucleotide (default), 'aa' for amino-acid	sequences, 				or 'generic'
	-n,		nb_data_sets
			nb_data_sets is an integer corresponding to the number of data sets to  			 analyse.
	-b,		int
			int >  0: int is the number of bootstrap replicates.
			int =  0: neither approximate likelihood ratio test nor bootstrap values 				are computed.
			int = -1: approximate likelihood ratio test returning aLRT statistics.
			int = -2: approximate likelihood ratio test returning Chi2-based 						parametric branch supports.
			int = -4: SH-like branch supports alone.
			int = -5: (default) approximate Bayes branch supports.
	-m,		model
	-f,		e, m, or fA,fC,fG,fT (具体参照phyml -help)
	-v,		prop_invar
			prop_invar : proportion of invariable sites. Can be a fixed value in the 				[0,1] range or ‘e’ to get the maximum likelihood estimate.
	-c,		nb_subst_cat
			nb_subst_cat : number of relative substitution rate categories. Default : 				nb_subst_cat=4.	Must be a positive integer.
	-a,		gamma
			gamma : distribution of the gamma distribution shape parameter. Can be a 				fixed positive value or ‘e’ to get the maximum likelihood estimate.
	-o,		params
			This option focuses on specific parameter optimisation.
			t: tree topology are optimised; l: branch length are optimised;
			r: rate parameters are optimised; n: no parameter is optimised.
	-s,		move
			Tree topology search operation option. Can be either NNI (default, fast) 				or SPR (a bit slower than NNI) or BEST (best of NNI and SPR search).
```

---

### 4.2.11 *raxml-ng* 建树

```bash
raxml-ng --all -msa final.phy --model GTR+I+G --prefix final --seed 2 --threads 10 --bs-metric fbp,tbe --bs-trees 1000

		--all,			all-in-one (ML search + bootstrapping)
			-msa,			alignment file
		--model,		model specification OR partition file
		--prefix,		prefix for output files (default: MSA file name)	
		--seed,			seed for pseudo-random number generator (default: current time)
	--threads,		number of parallel threads to use (default: 24)
	--bs-metric,	branch support metric: fbp = Felsenstein bootstrap (default), tbe 						= transfer distance	
	  --bs-trees,	number of bootstraps replicates.
```

---

### 4.2.12 *iqtree*建树

```bash
iqtree -s final.phy -st CODON -m GTR+I+G -b 1000 -nt 10

	-s,		Input alignment in PHYLIP/FASTA/NEXUS/CLUSTAL/MSF format
	-st,	BIN, DNA, AA, NT2AA, CODON, MORPH (default: auto-detect)
	-m,		model specification. if '-m MFP -mtree' used, this software will 						automatically search for the best model followed by tree inference
	-bb,	Ultrafast bootstrap (>=1000) (这个是iqtree特有的快速bootstrap)
	-b,		Bootstrap + ML tree + consensus tree (>=100)
	-nt,	number of parallel threads to use
```

---

### 4.2.13 mcmctree计算分歧时间

- mcmctree.tree
  
    ```
    10 1
    (Nnu,((Tha,Ath),((Clo,Bgy)'L(0.478)',((Mes,Rco),(Spu,(Peu,Ptr))))'L(0.48)'))'B(1.196,1.2863)';
    ```
    
- mcmctree.ctl
  
    ```bash
    seed = -1
    seqfile = mcmctree.phy
    treefile = mcmctree.tree
    outfile = mcmctree.out
    ndata = 3
    usedata = 1 * 0: no data; 1:seq; 2:approximation; 3:out.BV (in.BV)
    clock = 3 * 1: global clock; 2: independent; and 3: correlated rates
    RootAge = '<1.2863' * safe constraint on root age, used if no fossil for root.
    model = 4 * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
    alpha = 0 * alpha for gamma rates at sites
    ncatG = 5 * No. categories in discrete gamma
    cleandata = 0 * remove sites with ambiguity data (1:yes, 0:no)?
    BDparas = 1 1 0 * birth, death, sampling
    kappa_gamma = 6 2 * gamma prior for kappa
    alpha_gamma = 1 1 * gamma prior for alpha
    rgene_gamma = 2 2 * gamma prior for rate for genes
    sigma2_gamma = 1 10 * gamma prior for sigma^2 (for clock=2 or 3)
    finetune = 1: .1 .1 .1 .1 .1 .1 * auto (0 or 1) : times, rates, etc.
    print = 1
    burnin = 200000
    sampfreq = 10
    nsample = 50000
    ```
    

```bash
perl scripts/phy2mcmctreephy.pl final.phy
mcmctree mcmctree.ctl
```

---

# 5 脚本和整合的shell

---

## 5.1 整个流程的shell

注意：

1. 脚本放在**`./scripts`**目录下
2. 记得更改`1.extract_ortholog.pep.pl`和`3.extract_ortholog.cds.pl`中的物种数目（分别在第47行和第52行）
3. 记得检查`5.create_pal2aln.pl`中pal2nal的路径（第11行）
4. 记得检查`6.create.gblocks_sh.pl`中gblocks的路径（第9行）
5. 记得更改最后的jModelTest.jar的路径
6. 由于用Windows或者某些软件打开后，换行符那里会添加一个`^M`,linux识别的时候可能会出现问题，建议先用`cat -A final.sh` 检查换行符，然后用`sed -e 's/^M//g' final.sh >1.bak | mv 1.bak final.sh` 来去除^M，注意打^M的方式是`Ctrl+M`而不是分别打出^和M。

```bash
# 1. 脚本放在**`./scripts`**目录下
# 2. 记得更改`1.extract_ortholog.pep.pl`和`3.extract_ortholog.cds.pl`中的物种数目（分别在第47行和第52行）
# 3. 记得检查`5.create_pal2aln.pl`中pal2nal的路径（第11行）
# 4. 记得检查`6.create.gblocks_sh.pl`中gblocks的路径（第9行）
# 5. 记得更改最后的jModelTest.jar的路径
# 6. 由于用Windows或者某些软件打开后，换行符那里会添加一个`^M`,linux识别的时候可能会出现问题，建议先用`cat -A [final.sh](http://final.sh)` 检查换行符，然后用`sed -e 's/^M//g' final.sh >1.bak | mv 1.bak final.sh` 来去除^M，注意打^M的方式是`Ctrl+M`而不是分别打出^和M。
# 7. 由于下载到的序列可能会有蛋白质和cds序列不匹配的情况，所以pal2nal那一步可能会有报错，导致无法往后进行，建议跑shell跑到这一步，然后手动删除文件夹中的空文件，之后手动跑后面的步骤。

# perl scripts/1.extract_ortholog.pep.pl ./Orthogroups/Orthogroups.txt allpepfiles Single_Copy_Orthologue_Sequences &&
perl scripts/1.extract_ortholog.pep.pl ./Orthogroups/Orthogroups.txt allpepfiles Single_Copy_Orthologue_Sequences 
perl scripts/2.create_sh.mafft.pl Single_Copy_Orthologue_Sequences mafft.sh mafftseqs
bash mafft.sh >mafft.log 2>>mafft.err 
mkdir mafftout
mv mafft* mafftout
perl scripts/3.extract_ortholog.cds.pl ./Orthogroups/Orthogroups.txt all_cds_dir/ singlecopycdsdir length_check.txt 
perl scripts/4.sort_cds_mafft.pl mafftoutdir/mafftseqs singlecopycdsdir sortedcdsdir 
perl scripts/5.create_pal2aln.pl mafftoutdir/mafftseqs sortedcdsdir pal2nalfas pal2nal.sh 
bash pal2nal.sh >pal2nal.log 2>>pal2nal.err 
mkdir pal2naloutdir 
mv pal2nal* err* pal2naloutdir 
perl scripts/6.create.gblocks_sh.pl pal2naloutdir/pal2nalfas gblockshtms gblocks.sh &&
bash gblocks.sh >gblocks.log 2>>gblocks.err
mkdir gblocksoutdir
mkdir gblocksgbs
mv pal2naloutdir/pal2nalfas/*htm gblockshtms
mv pal2naloutdir/pal2nalfas/*gb gblocksgbs
mv gblocks* gblocksoutdir
perl scripts/7.cat_cds.pl gblocksoutdir/gblocksgbs mergedcdss
mkdir finalseqs
cat mergedcdss/*fa >finalseqs/final.fa
cd finalseqs
sh ../scripts/convertFasta2Phylip.sh final.fa >final.phy
mkdir tree
cd tree
ln -s ../finalseqs/final.phy .
nohup java -jar -XX:ParallelGCThreads=4 -Xmx4g ~/software/jmodeltest-2.1.10/jModelTest.jar -tr 10 -d final.phy -s 11 -f -i -g 8 -AIC -BIC -AICc -DT -p -a -w -o bestmodel &
```

---

## 5.2 脚本

### add_species_name_to_fa_sequences.pl

```bash
use autodie;

my $name=$ARGV[0];
$name=~s/\.fa//;
open IN,"<$ARGV[0]";
open OUT,">bak";

while (<IN>) {
        if (/^>/) {
                s/\s+$//;
                $a=(split/\s+/,$_)[0];
                $a=~s/^>//;
                $a=">$name|" . "$a";
                print OUT "$a\n";
        } else {
                print OUT;
        }
}
close IN;
close OUT;
`mv bak $ARGV[0]`;
```

---

### filter_fasta_by_length_and_sort.pl

```bash
#Usage: perl *pl Length_threshlod IN.fa passed.fa unpassed.fa Length_of_scaffold.txt
use Bio::SeqIO;

$length_threshold=$ARGV[0];
open CACHE,">>CACHE";
open PASSED,">>$ARGV[2]";
open FILTERED,">>$ARGV[3]";
open LEN,">>$ARGV[4]";

my $in=Bio::SeqIO->new(-format=>'fasta',
                    -file=>"$ARGV[1]");

while (my $seq=$in->next_seq) {
                $id=$seq->id;
                $length=$seq->length;
                $base=$seq->seq;
#print "$id\t$length\t$base\n";
                if ($length<=$length_threshold) {
                        print FILTERED ">$id\n$base\n";
                } else {
                        print CACHE ">$id\n$base\n";
                }
}
close FILTERED;
close CACHE;

my $in=Bio::SeqIO->new(-format=>"fasta",
                                        -file=>"CACHE");

while (my $seq=$in->next_seq) {
                $id=$seq->id;
                $length=$seq->length;
                $base=$seq->seq;

                $id{$length}=$id;
                $base{$id}=$base;
}
@length=sort {$b<=>$a} keys %id;
foreach $length (@length) {
        print PASSED ">$id{$length}\n$base{$id{$length}}\n";
        print LEN "$id{$length}\t$length\n";
}

close LEN;
close PASSED;
close CACHE;
unlink CACHE;
```

---

### 1.extract_ortholog.pep.pl

```bash
use strict;
use Bio::SeqIO;

my $orthomclfile=shift;
open (I,"<$orthomclfile");
my $allpepfiles=shift;       
my $outdir=shift;
system("mkdir $outdir");

my %pep;
my $fa=Bio::SeqIO->new(-file=>$allpepfiles,-format=>'fasta');
while (my $seq_obj=$fa->next_seq){
    my $pep_name=$seq_obj->display_name;
    my $seq=$seq_obj->seq;
    $pep{$pep_name}=$seq;
}
print "pep.fa has done\n";

my %family;
while (<I>){
    chomp;
    my $line=$_;
    my @inf=split/\s+/,$line;
    my $genefamily=shift(@inf);
    $genefamily=~s/://;
    foreach my $gene (@inf){
        $gene=~/^(\w+)\|(.+)/;
        my $species=$1;
        my $id=$gene;
        $family{species}{$genefamily}{$species}++;
        $family{id}{$genefamily}{$species}=$gene;

    }
}
close I;

my @familyes=keys %{$family{species}};
foreach my $genefamily (@familyes){
    my @species=sort keys %{$family{species}{$genefamily}};
    my $count=0;
    foreach my $species (@species){
        $count+=$family{species}{$genefamily}{$species};
    }
    my $mark=scalar(@species);
    if ($count==$mark && $mark==10){

        my $output="$outdir/$genefamily";
        open (O,">$output");
        foreach my $species (@species){
            my $id=$family{id}{$genefamily}{$species};
            my $ortholog=$pep{$id};
            print O ">$id\n$ortholog\n";
        }
        close O;
    }
}
```

---

### 2.create_sh.mafft.pl

```bash
use warnings;
use strict;
my $homolog_dir=shift;
my @files=<$homolog_dir/*>;
my $output_sh=shift;
open (O,">$output_sh");

my $outdir=shift;
system("mkdir $outdir");
foreach my $file (@files){
    $file=~/OG(\d+)/;
    my $sign=$1;
    my $outfile="$outdir/mafft$1";
    print O "mafft --maxiterate 1000 --localpair $file >$outfile\n";
}
close O;
```

---

### 3.extract_ortholog.cds.pl

```bash
use warnings;
use strict;
use Bio::SeqIO;

my $orthomclfile=shift;
open (I,"<$orthomclfile");
my $allnucledir=shift;
my @nuclefiles=<$allnucledir/*>;
my $outdir=shift;
system("mkdir $outdir");
my $genefamilylength_check=shift;
open (F,">$genefamilylength_check");

my %nucle;
foreach my $nuclefile (@nuclefiles){
    my $fa=Bio::SeqIO->new(-file=>$nuclefile,-format=>'fasta');
    while (my $seq_obj=$fa->next_seq){
    my $nucle_name=$seq_obj->display_name;
    my $seq=$seq_obj->seq;
    $nucle{$nucle_name}=$seq;
    }
    print "$nuclefile\n";
}

my %family;
while (<I>){
    chomp;
    my $line=$_;
    my @inf=split/\s+/,$line;
    my $genefamily=shift(@inf);
    $genefamily=~s/://;
    foreach my $gene (@inf){
    $gene=~/^(\w+)\|(.+)/;
    my $species=$1;
    my $id=$2;
    $family{species}{$genefamily}{$species}++;
    $family{id}{$genefamily}{$species}=$gene;

    }
}
close I;

my %familylength;
my @familyes=keys %{$family{species}};
foreach my $genefamily (@familyes){
    my @species=sort keys %{$family{species}{$genefamily}};
    my $count=0;
    foreach my $species (@species){
    $count+=$family{species}{$genefamily}{$species};
    }
    my $mark=scalar(@species);     
    if ($count==$mark && $mark==10){ 
    my $output="$outdir/$genefamily";
    open (O,">$output");
    foreach my $species (@species){
        my $id=$family{id}{$genefamily}{$species};
        my $ortholog=$nucle{$id};
        my $length=length($ortholog);
        if ($length<150){
        $familylength{$genefamily}++;
        }
        print O ">$id\n$ortholog\n";  
    }
    close O;
    }
}

my @familylength=keys %familylength;
foreach my $familylength (@familylength){
    print F "$familylength\n";
}
close F;
```

---

### 4.sort_cds_mafft.pl

```bash
use warnings;
use strict;

use Bio::SeqIO;
my $mafftdir=shift;
my @pepfiles=<$mafftdir/*>;
my $cdsdir=shift;
my $outdir=shift;
system ("mkdir $outdir");

foreach my $pepfile (@pepfiles){
    $pepfile=~/mafft(\d+)/;
    my $sign=$1;
    my $cdsfile="$cdsdir/OG$sign";

    my %pep;
    my $count=0;
    my $fa_P=Bio::SeqIO->new(-file=>$pepfile,-format=>'fasta');
    while (my $seq_obj=$fa_P->next_seq){
    my $pep_id=$seq_obj->display_name;
    $pep_id=~/^(\w+)\|(.+)/;
    my $species=$1;
    my $pepid=$2;
    $count++;
    $pep{$count}=$species;
    }

    my %nucle;
    my $fa_c=Bio::SeqIO->new(-file=>$cdsfile,-format=>'fasta');
    while (my $seq_obj=$fa_c->next_seq){
    my $nuc_id=$seq_obj->display_name;
    my $seq=$seq_obj->seq;
    $nuc_id=~/^(\w+)\|(.+)/;
    my $species=$1;
    my $nucid=$2;
    $nucle{$species}{$nuc_id}=$seq;
    }

    my $output="$outdir/paltoaln$sign";
    open(O,">$output");
    my @order=sort {$a<=>$b} keys %pep;
    foreach my $order (@order){
    my $species=$pep{$order};
    my @id=keys %{$nucle{$species}};   
    foreach my $id (@id){
        print O ">$id\n$nucle{$species}{$id}\n";
    }
    }
    close O;
}
```

---

### 5.create_pal2aln.pl

```bash
use warnings;
use strict;
my $pepdir=shift;
my @mafftfiles=<$pepdir/*>;
my $cdsdir=shift;
my $outdir=shift;
system ("mkdir $outdir");
my $output_sh=shift;
open(O,">$output_sh");

my $pal2alnpath="~/software/pal2nal.v14/pal2nal.pl";
foreach my $pepfile (@mafftfiles){
    $pepfile=~/mafft(\d+)/;
    my $sign=$1;
    print O "perl $pal2alnpath $pepfile $cdsdir/paltoaln$sign -output fasta >$outdir/pal2aln_$sign\n";
}
close O;
```

---

### 6.create.gblocks_sh.pl

```bash
use strict;
my $pepdir=shift;
my @files=<$pepdir/*>;
my $logdir=shift;
system ("mkdir $logdir");
my $output_sh=shift;
open(O,">$output_sh");

my $gblockspath="~/software/Gblocks_0.91b/Gblocks";
foreach my $pepfile (@files){
    print O "$gblockspath $pepfile -t=c\n";
}
close O;
```

---

### 7.cat_cds.pl

```bash
use warnings;
use strict;
use Bio::SeqIO;
my $trimmeddir=shift;
my @pepfiles=<$trimmeddir/*>;
my $pepdir=shift;
system("mkdir $pepdir");

my %pep;
my $count=0;
foreach my $file (sort @pepfiles){
    my $fa=Bio::SeqIO->new(-file=>$file,-format=>'fasta');
    while (my $seq_obj=$fa->next_seq){
        my $pepid=$seq_obj->display_name;
        my $pepseq=$seq_obj->seq;
    $pepid=~/(\w+)\|(.+)/;
    my $species=$1;
    $pep{$species}.=$pepseq;
    }
    $count++;
}
print "total files:$count\n";

my @species=sort keys %pep;
foreach my $species (@species){
    my $output="$pepdir/$species.fa";
    open (O,">$output");
    print O ">$species\n$pep{$species}\n";
    close O;
}
```

---

### convertFasta2Phylip.sh

```bash
#! /bin/sh

if [ $# != 1 ]; then
    echo "USAGE: ./script <fasta-file>"
    exit
fi

numSpec=$(grep -c  ">" $1)
tmp=$(cat $1 | sed "s/>[ ]*\(\w*\).*/;\1</"  | tr -d "\n" | tr -d ' '  | sed 's/^;//' | tr "<" " " )
length=$(($(echo $tmp | sed 's/[^ ]* \([^;]*\);.*/\1/'   | wc -m ) - 1))

echo "$numSpec $length"
echo  $tmp | tr ";" "\n"
```

---

### phy2mcmctreephy.pl

```perl
open PHY,"<$ARGV[0]";
open OUT,">","mcmctree.phy";

@lines=(<PHY>);
$header=shift @lines;
my ($species,$base)=split/\s+/,$header;
$codonbase=$base/3;

for (1..$species) {
        my $seq=shift @lines;
        my ($spe,$codon)=split/\s+/,$seq;
        my @bases=split//,$codon;
        my $count=1;
        foreach my $base (@bases) {
                $codon1{$spe}.="$base" if $count%3==1;
                $codon2{$spe}.="$base" if $count%3==2;
                $codon3{$spe}.="$base" if $count%3==0;
                $count++;
        }
}

foreach my $num (1..3) {
        print OUT "$species $codonbase\n";
        my $hashname="codon" . "$num";
        foreach my $spe (sort keys %codon1) {
                print OUT "$spe  ${$hashname}{$spe}\n";
        }
}

```

---

### mcmctree.ctl

```bash
seed = -1
seqfile = mcmctree.phy
treefile = mcmctree.tree
outfile = mcmctree.out
ndata = 3
usedata = 1 * 0: no data; 1:seq; 2:approximation; 3:out.BV (in.BV)
clock = 3 * 1: global clock; 2: independent; and 3: correlated rates
RootAge = '<1.2863' * safe constraint on root age, used if no fossil for root.
model = 4 * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
alpha = 0 * alpha for gamma rates at sites
ncatG = 5 * No. categories in discrete gamma
cleandata = 0 * remove sites with ambiguity data (1:yes, 0:no)?
BDparas = 1 1 0 * birth, death, sampling
kappa_gamma = 6 2 * gamma prior for kappa
alpha_gamma = 1 1 * gamma prior for alpha
rgene_gamma = 2 2 * gamma prior for rate for genes
sigma2_gamma = 1 10 * gamma prior for sigma^2 (for clock=2 or 3)
finetune = 1: .1 .1 .1 .1 .1 .1 * auto (0 or 1) : times, rates, etc.
print = 1
burnin = 200000
sampfreq = 10
nsample = 50000
```

---