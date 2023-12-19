# mpboot
MPBoot: Fast phylogenetic maximum parsimony tree inference and bootstrap approximation

## **COMPILING INSTRUCTION SINCE 2020**
* Clone the source code, unzip it, and rename to **source**
* Create folder **build** outside folder **source**
* Change directory to **build**
* Run **cmake** command:

**cmake ../source -DIQTREE_FLAGS=sse4 -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++**
* Replace **sse4** by **avx** in above command if you decide to run MPBoot on AVX architecture
* Run **make**
* You will find the executable named **mpboot** once the **make** command is done.

## **PLACEMENT CORE** 
### **Parameter**
* **-pp_on**: enable placement.
* **-pp_n**: number of available samples on tree.
* **-pp_k**: number of missing samples.
* **-pp_tree**: tree file.
* **-pp_origspr**: run spr without changing origin tree.
* **-pp_test_spr**: enable check correct tree.
* **-pp_origin**: origin tree file.

### **Command**
* Add missing samples to existing tree:
  <br>
  ``./mpboot -s <vcf file> -pptree <tree file> -ppon -ppn <existing samples> -ppk <missing samples>``
  <br>
  *Example:*
  <br>
  ``./mpboot -s data/test5/1/added5.vcf -pptree data/test5/1/origin5.fasta.treefile -ppon -ppk 5 -ppn 5``
* Check if tree unchanged
  <br>
  ``./mpboot -s <vcf file> -pptree <tree file> -pp_origin <origin tree file> -ppon -pp_test_spr -ppk <missing samples>``
  *Example:*
  <br>
  ``./mpboot -s data/test5/1/added5.vcf -pptree addedTree.txt -pp_origin data/test5/1/origin5.fasta.treefile -ppon -pp_test_spr -ppk 5``
<hr>
<br><br><br>


> ## **COMPILING INSTRUCTION PRIOR TO 2020**
> * Clone the source code, unzip it, and rename to **source**
> * Change directory to **source**, run following commands to update the sub-repository
> 
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**git submodule init**
> 
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**git submodule update**
> 
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**cd pllrepo**
> 
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**git checkout subufbootmp**
> 
> * Create folder **build** outside folder **source**
> * Change directory to **build**
> * Run **cmake** command:
> 
> **cmake ../source -DIQTREE_FLAGS=sse4 -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++**
> * Replace **sse4** by **avx** in above command if you decide to run MPBoot on AVX architecture
> * Run **make**
> * You will find the executable named **mpboot** once the **make** command is done.
