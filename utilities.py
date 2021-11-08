# Copied from https://github.com/rdkit/rdkit/issues/2320
def disable_rdkit_logging():
    """
    Disables RDKit whiny logging.
    """
    import rdkit.rdBase as rkrb
    import rdkit.RDLogger as rkl
    logger = rkl.logger()
    logger.setLevel(rkl.ERROR)
    rkrb.DisableLog('rdApp.error')
    

def standardize(smiles):
    # Sets up a molecule object from a string molecule representation (SMILE)
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        #  Allo precise substructural spefcification
        mol = rdMolStandardize.FragmentParent(mol)
        # Remove any hydorgens 
        mol = uncharger.uncharge(mol)
    return mol


def data_setup(filename="./data/Beta-secretase_1.balanced.tsv"):
    mol_data = pd.read_csv(filename, header=None,sep='\t')
    mol_data["MOL"] = mol_data[mol_data.columns[0]].apply(standardize)
    mol_data.rename(columns={mol_data.columns[0]: "SMILES", mol_data.columns[1]: "Label"}, inplace=True)
    mol_data.drop_duplicates(['SMILES'],inplace=True)
    mol_data.reset_index(drop=True, inplace=True)
    return mol_data
    
    
def create_fingerprints(morgan_radius, morgan_n_bits, mol_data):
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, morgan_radius, nBits=morgan_n_bits) for m in mol_data["MOL"]]
    X_morgan = np.asarray(fps)
    Y_labels = mol_data["Label"].to_numpy()
    return X_morgan, Y_labels