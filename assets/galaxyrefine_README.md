#GalaxyRefine: A program for protein structure model refinement

##0. Remark
The GalaxyRefine distribution version supports only **Linux 64-bit** OS and binary files compiled with serial option.
**Linux 32-bit** OS or binary files compiled with MPI option are not supported.

##1. Installation
1. Download the GalaxyRefine program
 * Download a copy of GalaxyRefine
    * https://github.com/seoklab/GalaxyRefine
    * https://github.com/seoklab/GalaxyRefine/archive/master.zip
 * Or, clone the GalaxyRefine git repository 
     * git clone https://github.com/seoklab/GalaxyRefine.git

2. Unzip and place the downloaded files
 * unzip GalaxyRefine.zip
 * mv GalaxyRefine $GALAXY_HOME  
    (*example*: GALAXY_HOME=/applic/GalaxyRefine)

3. Check the downloaded files
 * There should exist:
  * bin: directory for executables  
    There should be build_initial_model, local_optimize, generate_model, and **GalaxyRefine**
  * data: directory for data files
  * examples: directory for example files

4. Set environment variable
 * export GALAXY_HOME=$GALAXY_HOME  
    (*example*: export GALAXY_HOME=/applic/GalaxyRefine)

##2. How to use GalaxyRefine
1. Prepare the input protein structure in PDB format
 * You have to provide a protein structure in the PDB file format for structure refinement.
 Gaps (or chain breaks) are highly prohibited in the middle of the initial protein structure.

2. Run GalaxyRefine
 * Usage: $GALAXY_HOME/bin/GalaxyRefine [-h] [-p INPUT PDB File] [-t TITLE] [-s N_sample] [-o N_output]
 * Input arguments and options:     
    *  -p or --pdb      : Input protein structure file to refine (**mendatory**)    
    *  -t or --title    : The refinement job title (default: GalaxyRefine)      
    *  -s or --n_sample : The number of simulation trajectories for sampling (default: 16)   
    *  -o or --n_output : The number of refined output models (default: 5)      

3. Output of GalaxyRefine
 * The GalaxyRefine generates a working directory, which is named by the job title (default: GalaxyRefine)
 * In the working directory, the following three directories will be generated:   
    * init:   A working directory for cleaning up the input PDB file after initial side-chain optimization is performed.   
    * refine: A working directory for refining protein structure model with SC perturbations and MD relaxations.     
    * model:  The output of GalaxyRefine, refined protein model structures in PDB format will be placed.      
 * The final refined model will be **${JOB Title}/model/model.pdb**

##3. Release log
* 01 Feb 2016: The first release of GalaxyRefine

##4. References
* L. Heo, H. Park, and C. Seok, GalaxyRefine: Protein structure refinement driven by side-chain repacking, 
 Nucleic Acids Res. 41 (W1), W384-W388 (2013). [[PUBMED]](http://www.ncbi.nlm.nih.gov/pubmed/23737448)
* G. R. Lee, L. Heo, and C. Seok, Effective protein model structure refinement by loop modeling and overall
 relaxation, Proteins: Structure, Function, and Bioinformatics, in press (2015). [[PUBMED]](http://www.ncbi.nlm.nih.gov/pubmed/26172288)

##5. Contact
* chaok@snu.ac.kr
* compbio.galaxy@gmail.com

