import subprocess
import os

def runMAGMa(mgf, r, f='mgf', i='1', e='smiles', p='15', q='0.001', b='3',
             w='1', slow='False', o=',1200,False', g='False',
	     a='None', max_charge='1', n='1', t='None', l='error', call_back_url='None'):
    """Executes MAGMa through docker

    -h, help message
    -z Description of the job (default: )
    -f MS data input format {mass_tree,form_tree_pos,form_tree_neg,mgf}
       (default: mass_tree)
    -i Ionisation mode {-1,1}
       (default: 1)
    -e Output format for ranked compound list {smiles,sdf}
       (default: smiles)
    -p Maximum relative m/z error (ppm)
       (default: 5)
    -q Maximum absolute m/z error (Da) 
       (default: 0.001)
    -b Maximum number of bond breaks to generate substructures 
       (default: 3)
    -w Maximum number of additional water (OH) and/or ammonia (NH2) losses 
       (default: 1)
    --slow Skip fast calculations of molecules up to 64 atoms
       (default: False)
    -s Retrieve molecules from structure database {pubchem,kegg,hmdb}
       (default: )
    -o Specify structure database option: db_filename,max_mim,max_64atoms,incl_halo,min_refscore(only for
                                          PubChem),ids_file 
       (default: ,1200,False)
    -g Get references to PubChem 
       (default: False)
    -a Specify adduct (as comma separated list) for matching at MS1. Positive mode: [Na,K,NH4] 
                                                                     Negative mode: [Cl]
       (default: None)
    -r Read molecules from filename.sdf, from filename.smiles, or from a smiles string
    --max_charge Maximum charge state 
       (default: 1)
    -n Number of parallel cpus to use for annotation
       (default: 1)
    -t Maximum allowed time in minutes 
       (default: None)
    -l  Set logging level {debug,info,warn,error}
       (default: info)
    --call_back_url Call back url 
       (default: None)
    
    Returns
       space sep file with smiles,score,name,refscore,formula,min 
    """
    try:
        p = subprocess.Popen(["docker", "run", "--rm", "-v", os.getcwd()+":/data", "nlesc/magma", "light", "-f", f,
                              "-i", i, "-e", e, "-p", p, "-q", q, "-b", b, "-w", w,
    			  "-o", o, "--max_charge", max_charge,
    			  "-n", n,  "-l", l, '-r', r, mgf], 
    			  stdout=subprocess.PIPE)
        return p.communicate()
    except Exception as e:
        print(e)
