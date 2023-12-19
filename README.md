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
* **-ppon**: enable placement.
* **-ppn**: number of available samples on tree.
* **-ppk**: number of missing samples.
* **-pptree**: tree file.
* **-pp_test_spr**: enable check correct tree.
* **-pp_origin**: origin tree file.

### **Command**
* Add missing samples to existing tree:
  <br>
  ``./mpboot -s <vcf file> -pptree <tree file> -ppon -ppn <existing samples> -ppk <missing samples>``
* Check if tree unchanged
  ``./mpboot -s <vcf file> -pptree <tree file> -pp_origin <origin tree file> -ppon -pp_test_spr -ppk <missing samples>``

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
